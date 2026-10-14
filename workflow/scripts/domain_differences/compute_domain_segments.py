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
- results/4a_domain_differences/{dataset}/domain_hand_regions.tsv — coverage log only

Output
------
- domain_segments.tsv — one row per (gene, segment)
- domain_summary.tsv — one row per reported gene
- _work/segments.parquet — the segment index ranges, for the permutation-p rule

Usage
-----
    python workflow/scripts/domain_differences/compute_domain_segments.py \\
        --insertions results/4a_domain_differences/HD_DIT_HAP/_work/insertions.parquet \\
        --hand-regions results/4a_domain_differences/HD_DIT_HAP/domain_hand_regions.tsv \\
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from domain_differences.core import (  # noqa: E402
    hand_region_map,
    ordered_genes,
    segment_counts,
)
from domain_differences.segmentation import (  # noqa: E402
    merge_threshold,
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
    "Systematic ID", "name", "n_sites", "n_insertions", "n_segments", "min_seg_size",
    "span", "max_boundary_contrast", "sigma_gene", "lambda",
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class SegmentsConfig:
    """Runtime parameters for the segmentation."""

    insertions: Path
    hand_regions: Path
    min_seg: int
    z: float
    sigma_add: float
    by_site: bool
    output_segments: Path
    output_summary: Path
    output_work: Path


# =============================================================================
# CORE LOGIC
# =============================================================================
def segment_gene(sid: str, frame: pd.DataFrame, config: SegmentsConfig) -> dict | None:
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

    means = [float(drs[lo:hi].mean()) for lo, hi in segments]
    summary = {
        "Systematic ID": sid, "name": frame["name"].iloc[0],
        "n_sites": sum(c[1] for c in counts), "n_insertions": sum(c[0] for c in counts),
        "n_segments": len(segments), "min_seg_size": min(hi - lo for lo, hi in segments),
        "span": round(float(drs.max() - drs.min()), 4),
        "max_boundary_contrast": round(
            max(abs(means[i] - means[i + 1]) for i in range(len(means) - 1)), 4),
        "sigma_gene": round(sigma_gene, 4), "lambda": round(lam, 4),
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

    segment_rows, summary_rows, range_frames = [], [], []
    universe = set()
    for sid, frame in profiles.items():
        if len(frame) < floor:
            continue
        universe.add(sid)
        result = segment_gene(sid, frame, config)
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
    log_hand_coverage(config, universe, set(summary["Systematic ID"]))
    logger.info(f"wrote {config.output_segments.name}, {config.output_summary.name}, "
                f"{config.output_work.name}")
    return 0


def log_hand_coverage(config: SegmentsConfig, universe: set[str], reported: set[str]) -> None:
    """Report how many hand-curated genes the scan found — the only external check."""
    if not config.hand_regions.exists():
        logger.warning(f"no hand-region table at {config.hand_regions}; skipping the coverage check")
        return
    regions = hand_region_map(pd.read_csv(config.hand_regions, sep="\t"))
    in_universe = set(regions) & universe
    hit = set(regions) & reported
    logger.info(f"hand-curated: {len(regions)} genes with a usable region, "
                f"{len(in_universe)} inside the universe, {len(hit)} reported "
                f"({100 * len(hit) / max(len(in_universe), 1):.0f}% of the universe's hand genes)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Segment gene DR profiles along the CDS (4a)")
    parser.add_argument("--insertions", type=Path, required=True, help="Input insertions.parquet")
    parser.add_argument("--hand-regions", type=Path, required=True, help="Hand-region TSV (coverage log)")
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
