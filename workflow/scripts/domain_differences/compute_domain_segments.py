#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute Domain Segments
=======================

Stage 4a, step 3: segment every gene's DR profile with TGUH and report the
boundaries, the segment levels and the per-segment reliability columns.

The rule is one threshold per gene — ``lam = z * sigma_gene + sigma_add``, with
``sigma_gene`` the robust point-to-point noise of that gene's own profile — and a
gene is reported only when the merge leaves it with two or more segments. There is
no p-value filter: the effect floor is the only filter, and the permutation p's are
computed by a separate rule so they can be skipped.

Two counting conventions, both reported per segment: ``n_sites`` is the statistic
(the independent insertion sites the segmentation used) and ``n_insertions`` is the
evidence (the raw Forward/Reverse readings behind it). They are never equal — one
site normally carries two readings.

Input
-----
- results/4a_domain_differences/{dataset}/_work/insertions.parquet
- results/4a_domain_differences/{dataset}/domain_hand_regions.tsv — coverage log and
  the ``hand_term`` / ``hand_agrees`` columns
- resources/external/pombase/{version}/Protein_features/peptide_stats.tsv — protein
  lengths, for the coverage columns
- results/1c_annotation/.../gene_annotation_reference.protein.parquet — FYPOviability
  and gRNA DR (the same context the gallery titles read)

Output
------
- domain_segments.tsv — one row per (gene, segment)
- domain_summary.tsv — one row per reported gene
- _work/segments.parquet — the segment index ranges, for the permutation-p rule

The summary's classification columns hang on ONE statistic, ``region_contrasts()``:
each segment's mean minus the weighted mean of the rest of its own gene. An
absolute DR cut is not usable here — DR carries a gene-level offset, and only 5% of
segments reach 0, so "tolerant = DR >= 0" would miss most of the hand-curated C
termini. Against that contrast the curated calls come back at 94% (C) / 88% (N).

  tolerant_region          label of the most tolerant segment — N_term / C_term /
                           internal — or ``none`` when even that segment is within
                           the gene's own lambda of the rest
  region_*                 where that segment is and how far above the rest; filled
                           even for ``none``, so a near-miss stays readable
  N_contrast / C_contrast  the same statistic for the two end segments
  terminal_pattern         which termini clear lambda: N_only / C_only / both / neither
  protein_len, coverage_frac, n_gap_aa, c_gap_aa
                           how much of the protein the scan saw. The gaps matter:
                           the median first/last insertion is 50/41 aa from the real
                           terminus, and only ~30% of genes get within 20 aa — a
                           terminus call is only as good as its gap
  median_H, n_single_site_segments
                           reliability: per-segment H (aggregated) and how many
                           segments rest on a single site
  gene_DR, FYPOviability, gRNA_DR, hand_term, hand_agrees
                           context joins. ``hand_agrees`` is the stage's only
                           external check, promoted from a log line to a column

Usage
-----
    python workflow/scripts/domain_differences/compute_domain_segments.py \\
        --insertions results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet \\
        --hand-regions results/4a_domain_differences/HD_DIT_HAP/domain_hand_regions.tsv \\
        --peptide-stats resources/external/pombase/{version}/Protein_features/peptide_stats.tsv \\
        --annotation results/1c_annotation/{pombase}/{sgd}/gene_annotation_reference.protein.parquet \\
        --min-seg 1 --z 1.5 --sigma-add 0.1 --by-site \\
        --output-segments results/4a_domain_differences/HD_DIT_HAP/domain_segments.tsv \\
        --output-summary results/4a_domain_differences/HD_DIT_HAP/domain_summary.tsv \\
        --output-work results/4a_domain_differences/HD_DIT_HAP/_work/segments.parquet

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from domain_differences.core import (  # noqa: E402
    hand_region_map,
    hand_terms,
    ordered_genes,
    read_gene_context,
    read_protein_lengths,
    segment_counts,
)
from domain_differences.segmentation import (  # noqa: E402
    merge_threshold,
    region_contrasts,
    segment_metrics,
    tguh_tuned,
    universe_floor,
)
from io_table import read_parquet, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SEGMENT_COLUMNS = [
    "Systematic ID", "name", "segment", "start_residue", "end_residue",
    "start_frac", "end_frac", "n_insertions", "n_sites", "mean_DR",
    "sigma_seg", "sigma_gene", "H", "se_seg", "z_boundary", "lambda",
]
SUMMARY_COLUMNS = [
    "Systematic ID", "name",
    "n_sites", "n_insertions", "n_segments", "n_boundaries", "min_seg_size",
    "span", "max_boundary_contrast", "strongest_boundary_residue",
    "sigma_gene", "lambda",
    "tolerant_region", "region_start_residue", "region_end_residue",
    "region_contrast", "region_n_sites",
    "N_contrast", "C_contrast", "terminal_pattern",
    "protein_len", "coverage_frac", "n_gap_aa", "c_gap_aa",
    "median_H", "n_single_site_segments",
    "gene_DR", "FYPOviability", "gRNA_DR", "hand_term", "hand_agrees",
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class SegmentsConfig:
    """Runtime parameters for the segmentation."""

    insertions: Path
    hand_regions: Path
    peptide_stats: Path
    annotation: Path
    min_seg: int
    z: float
    sigma_add: float
    by_site: bool
    output_segments: Path
    output_summary: Path
    output_work: Path


@dataclass(kw_only=True, slots=True, frozen=True)
class References:
    """Per-gene lookups the summary joins in, loaded once per run."""

    lengths: dict[str, int]
    viability: dict[str, str]
    grna_dr: dict[str, float]
    hand: dict[str, str]


# =============================================================================
# CORE LOGIC
# =============================================================================
def block_label(index: int, n_segments: int) -> str:
    """Where a segment sits in its profile: ``N_term``, ``C_term`` or ``internal``."""
    if index == 0:
        return "N_term"
    if index == n_segments - 1:
        return "C_term"
    return "internal"


def terminal_pattern(n_contrast: float, c_contrast: float, threshold: float) -> str:
    """Which termini stand above the rest of the gene by more than ``threshold``."""
    n_tol, c_tol = n_contrast >= threshold, c_contrast >= threshold
    if n_tol and c_tol:
        return "both"
    if n_tol:
        return "N_only"
    if c_tol:
        return "C_only"
    return "neither"


def hand_agrees(term: str | None, n_contrast: float, c_contrast: float, threshold: float) -> str:
    """Does the scan flag the terminus the curator called non-essential? "" if unannotated."""
    if term is None:
        return ""
    asserted = (["N"] if term == "N" else ["C"] if term == "C" else ["N", "C"])
    contrasts = {"N": n_contrast, "C": c_contrast}
    return "yes" if all(contrasts[t] >= threshold for t in asserted) else "no"


def segment_gene(sid: str, frame: pd.DataFrame, config: SegmentsConfig,
                 refs: References) -> dict | None:
    """One gene's segment rows, summary row and index ranges; None if it is not reported."""
    drs = frame["DR"].to_numpy(dtype=float)
    residues = frame["Residue_affected"].astype(int).to_numpy()
    fracs = frame["frac"].to_numpy(dtype=float)

    segments = tguh_tuned(drs, config.min_seg, z=config.z, sigma_add=config.sigma_add)
    if len(segments) < 2:
        return None
    metrics = segment_metrics(drs, segments)
    counts = segment_counts(frame, segments)
    lam = merge_threshold(drs, config.z, config.sigma_add)
    sigma_gene = metrics[0]["sigma_gene"]

    segment_rows = []
    for k, (lo, hi) in enumerate(segments):
        n_insertions, n_sites = counts[k]
        segment_rows.append({
            "Systematic ID": sid, "name": frame["name"].iloc[0], "segment": k + 1,
            "start_residue": int(residues[lo]), "end_residue": int(residues[hi - 1]),
            "start_frac": round(float(fracs[lo]), 6), "end_frac": round(float(fracs[hi - 1]), 6),
            "n_insertions": n_insertions, "n_sites": n_sites, "mean_DR": round(float(drs[lo:hi].mean()), 4),
            "sigma_seg": round(metrics[k]["sigma_seg"], 4),
            "sigma_gene": round(sigma_gene, 4),
            "H": round(metrics[k]["H"], 3),
            "se_seg": round(metrics[k]["se_seg"], 4),
            "z_boundary": round(metrics[k]["z_boundary"], 3),
            "lambda": round(lam, 4),
        })

    means = np.array([float(drs[lo:hi].mean()) for lo, hi in segments])
    jumps = np.abs(np.diff(means))
    # contrasts of each segment against the rest of its own gene, and of the two
    # termini specifically; the label is the argmax's, not a separate computation
    # (see the module docstring for why the reference is the gene itself)
    contrasts = region_contrasts(means, np.array([c[1] for c in counts], dtype=float))
    best = int(contrasts.argmax())
    n_contrast, c_contrast = float(contrasts[0]), float(contrasts[-1])
    protein_len = refs.lengths.get(sid)
    first_residue, last_residue = int(residues[0]), int(residues[-1])
    term = refs.hand.get(sid)
    summary = {
        "Systematic ID": sid, "name": frame["name"].iloc[0],
        "n_sites": sum(c[1] for c in counts), "n_insertions": sum(c[0] for c in counts),
        "n_segments": len(segments), "n_boundaries": len(segments) - 1,
        "min_seg_size": min(hi - lo for lo, hi in segments),
        "span": round(float(drs.max() - drs.min()), 4),
        "max_boundary_contrast": round(float(jumps.max()), 4),
        "strongest_boundary_residue": int(residues[segments[int(jumps.argmax())][1] - 1]),
        "sigma_gene": round(sigma_gene, 4), "lambda": round(lam, 4),
        "tolerant_region": block_label(best, len(segments)) if contrasts[best] >= lam else "none",
        "region_start_residue": int(residues[segments[best][0]]),
        "region_end_residue": int(residues[segments[best][1] - 1]),
        "region_contrast": round(float(contrasts[best]), 4),
        "region_n_sites": int(counts[best][1]),
        "N_contrast": round(n_contrast, 4), "C_contrast": round(c_contrast, 4),
        "terminal_pattern": terminal_pattern(n_contrast, c_contrast, lam),
        "protein_len": protein_len if protein_len is not None else np.nan,
        "coverage_frac": (round((last_residue - first_residue + 1) / protein_len, 4)
                          if protein_len else np.nan),
        "n_gap_aa": first_residue - 1,
        "c_gap_aa": protein_len - last_residue if protein_len is not None else np.nan,
        "median_H": round(float(np.median([m["H"] for m in metrics])), 3),
        "n_single_site_segments": sum(1 for c in counts if c[1] == 1),
        "gene_DR": round(float(frame["gene_DR"].iloc[0]), 4),
        "FYPOviability": refs.viability.get(sid, ""),
        "gRNA_DR": refs.grna_dr.get(sid, np.nan),
        "hand_term": term if term is not None else "",
        "hand_agrees": hand_agrees(term, n_contrast, c_contrast, lam),
    }
    ranges = pd.DataFrame({
        "Systematic ID": sid, "segment": range(1, len(segments) + 1),
        "lo": [lo for lo, _ in segments], "hi": [hi for _, hi in segments],
        "n_points": [hi - lo for lo, hi in segments],
    })
    return {"segments": segment_rows, "summary": summary, "ranges": ranges}


@logger.catch(reraise=True)
def run(config: SegmentsConfig) -> int:
    """Segment every gene in the universe, write the two tables plus the ranges."""
    insertions = read_parquet(config.insertions)
    profiles = ordered_genes(insertions, by_site=config.by_site)
    floor = universe_floor(config.min_seg)
    viability, grna_dr = read_gene_context(config.annotation)
    hand = (pd.read_csv(config.hand_regions, sep="\t") if config.hand_regions.exists()
            else None)
    refs = References(
        lengths=read_protein_lengths(config.peptide_stats),
        viability=viability, grna_dr=grna_dr,
        hand=hand_terms(hand) if hand is not None else {},
    )

    segment_rows, summary_rows, range_frames = [], [], []
    universe = set()
    for sid, frame in profiles.items():
        if len(frame) < floor:
            continue
        universe.add(sid)
        result = segment_gene(sid, frame, config, refs)
        if result is None:
            continue
        segment_rows.extend(result["segments"])
        summary_rows.append(result["summary"])
        range_frames.append(result["ranges"])

    if not segment_rows:
        raise ValueError(
            f"no gene reached two segments (universe {len(universe)}, floor {floor}); "
            f"check --z/--sigma-add/--min-seg"
        )

    segments = pd.DataFrame(segment_rows, columns=SEGMENT_COLUMNS)
    summary = (pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
               .sort_values("max_boundary_contrast", ascending=False).reset_index(drop=True))
    ranges = pd.concat(range_frames, ignore_index=True)

    segments.to_csv(config.output_segments, index=False, sep="\t")
    summary.to_csv(config.output_summary, index=False, sep="\t")
    write_parquet(ranges, config.output_work)

    unit = "sites" if config.by_site else "insertions"
    logger.info(f"universe: {len(universe):,} genes with >= {floor} {unit}; "
                f"reported {len(summary):,} genes with >= 2 segments")
    logger.info(f"segments: {len(segments):,}; "
                f"median {summary['n_segments'].median():.0f} segments per gene, "
                f"min segment size {segments.groupby('Systematic ID').size().min()}")
    logger.info(f"thresholds: lambda median {segments['lambda'].median():.3f} "
                f"(z={config.z}, sigma_add={config.sigma_add})")
    non_first = segments["z_boundary"].dropna()
    logger.info(f"|z_boundary| >= 2 in {100 * (non_first.abs() >= 2).mean():.0f}% "
                f"of the {len(non_first):,} non-first segments; "
                f"H median {segments['H'].median():.2f} "
                f"(degenerate sigma_seg: {100 * segments['sigma_seg'].eq(0).mean():.0f}%)")
    log_summary_columns(summary)
    log_hand_coverage(hand, universe, set(summary["Systematic ID"]), summary)
    logger.info(f"wrote {config.output_segments.name}, {config.output_summary.name}, "
                f"{config.output_work.name}")
    return 0


def log_summary_columns(summary: pd.DataFrame) -> None:
    """Report what the classification columns came out as — the summary's own QC."""
    counts = summary["tolerant_region"].value_counts()
    logger.info("tolerant_region: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    logger.info("terminal_pattern: " + ", ".join(
        f"{k} {v}" for k, v in summary["terminal_pattern"].value_counts().items()))
    for col, name in [("n_gap_aa", "N"), ("c_gap_aa", "C")]:
        gap = summary[col].dropna()
        logger.info(f"{name}-terminus gap: median {gap.median():.0f} aa, "
                    f"{100 * (gap <= 20).mean():.0f}% of genes within 20 aa")


def log_hand_coverage(hand: pd.DataFrame | None, universe: set[str], reported: set[str],
                      summary: pd.DataFrame) -> None:
    """Report how many hand-curated genes the scan found — the only external check."""
    if hand is None:
        logger.warning("no hand-region table; skipping the coverage check")
        return
    regions = hand_region_map(hand)
    in_universe = set(regions) & universe
    hit = set(regions) & reported
    logger.info(f"hand-curated: {len(regions)} genes with a usable region, "
                f"{len(in_universe)} inside the universe, {len(hit)} reported "
                f"({100 * len(hit) / max(len(in_universe), 1):.0f}% of the universe's hand genes)")
    scored = summary[summary["hand_agrees"].ne("")]
    if len(scored):
        logger.info(f"hand_agrees: {100 * scored['hand_agrees'].eq('yes').mean():.0f}% of "
                    f"{len(scored)} reported hand genes flag the curated terminus "
                    f"({(scored['hand_term'] == 'C').sum()} C, "
                    f"{(scored['hand_term'] == 'N').sum()} N, "
                    f"{(scored['hand_term'] == 'N+C').sum()} N+C)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Segment gene DR profiles along the CDS (4a)")
    parser.add_argument("--insertions", type=Path, required=True, help="Input insertions.parquet")
    parser.add_argument("--hand-regions", type=Path, required=True,
                        help="Hand-region TSV (coverage log + hand columns)")
    parser.add_argument("--peptide-stats", type=Path, required=True,
                        help="PomBase peptide_stats.tsv (protein lengths)")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="1c annotation reference parquet (per-gene context)")
    parser.add_argument("--min-seg", type=int, required=True, help="Minimum points per segment")
    parser.add_argument("--z", type=float, required=True, help="Noise multiple in lam = z*sigma + sigma_add")
    parser.add_argument("--sigma-add", type=float, required=True, help="Additive floor in the threshold")
    parser.add_argument("--by-site", action=argparse.BooleanOptionalAction, default=True,
                        help="Collapse each site's Forward/Reverse pair before segmenting")
    parser.add_argument("--output-segments", type=Path, required=True, help="Output segment TSV")
    parser.add_argument("--output-summary", type=Path, required=True, help="Output gene summary TSV")
    parser.add_argument("--output-work", type=Path, required=True, help="Output segment-range parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the segmentation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = SegmentsConfig(
            insertions=args.insertions.resolve(),
            hand_regions=args.hand_regions.resolve(),
            peptide_stats=args.peptide_stats.resolve(),
            annotation=args.annotation.resolve(),
            min_seg=args.min_seg, z=args.z, sigma_add=args.sigma_add, by_site=args.by_site,
            output_segments=args.output_segments.resolve(),
            output_summary=args.output_summary.resolve(),
            output_work=args.output_work.resolve(),
        )
        if config.min_seg < 1:
            raise ValueError(f"--min-seg must be at least 1, got {config.min_seg}")
        config.output_segments.parent.mkdir(parents=True, exist_ok=True)
        config.output_work.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
