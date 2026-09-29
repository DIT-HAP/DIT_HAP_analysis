#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Incoherence Attribution (WHY a complex is internally dispersed) — Computation Only
===========================================================================================

Complementary to the coherence z-score: the tightest groups are the expected
obligate machines, but the INCOHERENT ones (z > incoherent_z_threshold = more dispersed than random)
are biologically informative too — a complex's members scatter in DR-DL fitness
space when they play different functional roles. This stage scores the candidate
causes of that dispersion for every group and labels each with its most likely
explanation, so the meaningful incoherent complexes can be found and interpreted.

Diagnostic signals (per group, from workflow/src/coherence/attribution.py)
--------------------------------------------------------------------------
- major/minor split: a 2-component GMM on the members' normalized (DR, DL/10)
  points, called a genuine split when its silhouette is high enough — a tight
  essential "core" + a looser dispensable "minority" (e.g. eIF3 core vs eIF3e).
- shared-subunit fraction: fraction of members that also belong to OTHER groups
  of the same source (cross-complex members drag the centroid apart).
- paralog fraction: fraction of members with a paralog (deletion phenotype may be
  buffered by redundancy, dampening DR and pulling the group toward WT).
The label ladder (attribution.py::attribute_incoherence): conditional_module
(split + shared) > major_minor_split > shared_subunits > paralog_buffered >
data_limited > intrinsic_heterogeneity (real spread, no detected cause — for
manual review; also where annotation/technical artefacts fall, which are NOT
auto-labelled).

This is the computation component (ADR-0001); plot_incoherence_attribution.py
renders the figure from these two outputs.

Input
-----
- --metrics: a coherence metrics table (source, group_id, group_name,
  n_scored_members, scored_member_names, median_pairwise_distance_z, q_value, ...).
  One source's Parquet for a per-source run; the pooled de-duplicated
  representatives TSV for the dataset-level run.
- --annotation: the matching group_annotation_long.tsv table(s) (group -> member
  genes). One path per source; they are concatenated, and `source` is what keeps a
  group_id shared by two sources (173 of them) from resolving to the union of two
  different groups' members.
- --fitting-results: upstream fitting_results.tsv (see coherence/io.py).
- --paralogs: Ensembl paralog export TSV (its "Gene stable ID" column lists genes
  with >=1 paralog).

Output
------
- --output-table: incoherence_attribution.tsv — one row per scored group with
  median_pairwise_distance_z, q_value, gmm_silhouette, core_size, minor_size,
  frac_shared_members, paralog_fraction, attribution_label, is_incoherent, and the
  shared/other-group member detail. Sorted by median_pairwise_distance_z
  descending (most incoherent first). TSV: this
  is a final human-facing table, not a pipeline intermediate.
- --output-points: incoherence_split_points.parquet — one row per member of every
  INCOHERENT group: group_id, Systematic ID, norm_DR, norm_DL, component
  ("core" / "minor" / "single") — plus `source` when the input carried one, which
  is what lets the figure tell two same-group_id panels apart. This is what the
  figure colours by, persisted so the plot never re-fits the GMM.

Usage
-----
    python compute_incoherence_attribution.py \\
        --metrics results/3a_coherence/{dataset}/go_macrocomplex/coherence_metrics.parquet \\
        --annotation results/3a_coherence/{dataset}/go_macrocomplex/group_annotation_long.tsv \\
        --fitting-results .../fitting_results.tsv \\
        --paralogs resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --z-threshold 1.0 \\
        --output-table results/3a_coherence/{dataset}/go_macrocomplex/incoherence_attribution.tsv \\
        --output-points results/3a_coherence/{dataset}/go_macrocomplex/incoherence_split_points.parquet

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  2.0.0
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
from coherence.attribution import (  # noqa: E402
    major_minor_split,
    shared_subunits,
    shared_subunit_fractions,
    paralog_fraction,
    attribute_incoherence,
)
from coherence.io import load_fitting_results, load_long_table  # noqa: E402
from io_table import read_file, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Columns the coherence metrics table must carry for attribution.
_REQUIRED_METRIC_COLUMNS = ["group_id", "group_name", "n_scored_members",
                            "median_pairwise_distance_z"]

# Split-point component labels. A group whose members did not split into two GMM
# components gets "single" for every member.
_COMPONENT_CORE = "core"
_COMPONENT_MINOR = "minor"
_COMPONENT_SINGLE = "single"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AttributionConfig:
    """Inputs, outputs, and parameters for incoherence attribution."""
    metrics: Path
    # One path per source: a per-source run passes one, the pooled dedup run passes
    # every registered source's table. They are concatenated, which is what makes
    # the pooled call correct — `source` is what keeps a group_id shared by two
    # sources (173 of them, e.g. GO:0032040 in both go_cc and go_macrocomplex) from
    # resolving to the union of both member sets.
    annotations: tuple[Path, ...]
    fitting_results: Path
    paralogs: Path
    output_table: Path
    output_points: Path
    z_threshold: float = 1.0
    shared_frac_threshold: float = 0.5
    paralog_frac_threshold: float = 0.5

    def validate(self) -> None:
        """Raise ValueError on missing inputs, then make output dirs."""
        paths = [self.metrics, *self.annotations, self.fitting_results, self.paralogs]
        for path in paths:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_table, self.output_points]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC — data loading
# =============================================================================
def load_paralog_ids(paralogs: Path) -> set[str]:
    """Genes that actually HAVE a paralog, from the Ensembl paralog export."""
    # The export has one row per (gene, paralogue) pair, but genes WITHOUT a
    # paralogue still appear with the paralogue columns left blank — so filtering to
    # a non-empty `...paralogue gene stable ID` is essential (otherwise every gene
    # counts as having a paralog and paralog_fraction is a useless 1.0 everywhere).
    par = pd.read_csv(paralogs, sep="\t")
    gene_col = "Gene stable ID"
    para_col = "Schizosaccharomyces pombe paralogue gene stable ID"
    if gene_col not in par.columns or para_col not in par.columns:
        raise ValueError(f"paralog TSV missing '{gene_col}'/'{para_col}' (have: {list(par.columns)[:5]}...)")
    return set(par.loc[par[para_col].notna(), gene_col].dropna().astype(str))


# =============================================================================
# CORE LOGIC — per-group attribution
# =============================================================================
def group_member_points(
    long_table: pd.DataFrame, group_id: str, points: pd.DataFrame, source: str | None = None
) -> tuple[list[str], np.ndarray]:
    """The group's members that have fitness points, as (ids, (n,2) array)."""
    # Uses the same long-table -> point-cloud join coherence used, so the GMM sees
    # exactly the member set the z-score was computed on (members without a fitted
    # DR/DL are dropped, matching compute_coherence's inner merge).
    #
    # `source` must be passed on a POOLED table: group_id alone is ambiguous there
    # (173 group_ids appear in two sources), so without it the GMM would be fitted
    # to the union of two different groups' members. A single-source table makes it
    # optional — the union over one source is just the group.
    rows = long_table
    if source is not None and "source" in long_table.columns:
        rows = long_table[long_table["source"] == source]
    members = rows.loc[rows["group_id"] == group_id, "Systematic ID"].unique()
    ids = [m for m in members if m in points.index]
    X = points.loc[ids].to_numpy(dtype=float) if ids else np.empty((0, 2))
    return ids, X


def component_labels(split: dict, n_members: int) -> list[str]:
    """Per-member GMM component label: "core", "minor", or "single" when unsplit."""
    labels = split.get("labels")
    core_label = split.get("core_label")
    if labels is None or core_label is None:
        return [_COMPONENT_SINGLE] * n_members
    return [
        _COMPONENT_CORE if int(label) == int(core_label) else _COMPONENT_MINOR
        for label in labels
    ]


def attribute_all(
    metrics: pd.DataFrame,
    long_table: pd.DataFrame,
    points: pd.DataFrame,
    paralog_ids: set[str],
    config: AttributionConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attribution rows for every scored group + the split points the figure draws."""
    # Returns (attribution_table, split_points). Split points cover only the
    # INCOHERENT groups — the figure never draws the rest — and carry each member's
    # normalized DR/DL and its GMM component, so the plot never re-fits the GMM.
    #
    # Everything below is keyed on (source, group_id) when the tables carry a
    # `source` column: group_id is not unique across sources, so on the pooled dedup
    # table a group_id-only lookup returns whichever source happens to come first.
    # `sourced` mirrors `attribution.member_pairs`' own column test, so the two key
    # shapes stay in step.
    sourced = "source" in long_table.columns and "source" in metrics.columns
    shared_fractions = shared_subunit_fractions(long_table)

    rows = []
    point_rows = []
    for _, m in metrics.iterrows():
        group_id = m["group_id"]
        source = str(m["source"]) if sourced else None
        key = (source, group_id) if sourced else group_id
        ids, X = group_member_points(long_table, group_id, points, source=source)
        split = major_minor_split(X)

        shared_frac = shared_fractions.get(key, np.nan)
        par_frac = paralog_fraction(ids, paralog_ids)
        label = attribute_incoherence(
            split, shared_frac, par_frac,
            shared_frac_threshold=config.shared_frac_threshold,
            paralog_frac_threshold=config.paralog_frac_threshold,
        )
        is_incoherent = bool(m["median_pairwise_distance_z"] > config.z_threshold)

        # One shared_subunits call serves both the detail column and its count.
        shared_members = shared_subunits(long_table, group_id, source=source)
        sizes = split.get("component_sizes")
        core_label = split.get("core_label")

        rows.append({
            "source": m.get("source"),
            "group_id": group_id,
            "group_name": m["group_name"],
            "n_scored_members": m["n_scored_members"],
            "median_pairwise_distance_z": m["median_pairwise_distance_z"],
            "q_value": m.get("q_value", np.nan),
            "is_incoherent": is_incoherent,
            "gmm_silhouette": split.get("silhouette", np.nan),
            "gmm_is_split": bool(split.get("is_split")),
            "core_size": sizes[core_label] if sizes and core_label is not None else np.nan,
            "minor_size": sizes[1 - core_label] if sizes and core_label is not None else np.nan,
            "frac_shared_members": shared_frac,
            "paralog_fraction": par_frac,
            "attribution_label": label,
            "n_shared_members": len(shared_members),
            "shared_members": "; ".join(shared_members["Systematic ID"]) if not shared_members.empty else "",
        })

        if is_incoherent and ids:
            components = component_labels(split, len(ids))
            point_rows.extend(
                {"source": source, "group_id": group_id, "Systematic ID": gene,
                 "norm_DR": X[i, 0], "norm_DL": X[i, 1], "component": components[i]}
                for i, gene in enumerate(ids)
            )

    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.sort_values("median_pairwise_distance_z", ascending=False).reset_index(drop=True)
    # `source` is carried only when it keys the rows: an unsourced table has no
    # source to record whatsoever, and a NaN column would just read as a real one
    # (the DataFrame constructor drops the extra key against an explicit column list).
    point_columns = ["group_id", "Systematic ID", "norm_DR", "norm_DL", "component"]
    if sourced:
        point_columns = ["source", *point_columns]
    points_table = pd.DataFrame(point_rows, columns=point_columns)
    return table, points_table


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: AttributionConfig) -> None:
    """Load -> per-group attribution -> TSV + split-points Parquet."""
    config.validate()
    metrics = read_file(config.metrics)
    missing = [col for col in _REQUIRED_METRIC_COLUMNS if col not in metrics.columns]
    if missing:
        raise ValueError(f"metrics table missing required column(s) {missing} (have: {list(metrics.columns)})")
    # One long table per source, concatenated. `source` survives the concat and is
    # what keeps a group_id shared by two sources (173 of them) from resolving to
    # the union of two different groups' members.
    long_table = pd.concat(
        [load_long_table(path) for path in config.annotations], ignore_index=True
    )
    logger.info(
        f"annotation: {len(long_table):,} rows from {len(config.annotations)} source table(s); "
        f"sources={sorted(long_table['source'].unique()) if 'source' in long_table.columns else 'none'}"
    )
    # keep="last" mirrors the dict comprehension this replaced: a repeated gene id
    # resolved to its last row, and .loc[] on a duplicated index would fan out.
    points = (
        load_fitting_results(config.fitting_results)
        .drop_duplicates(subset="Systematic ID", keep="last")
        .set_index("Systematic ID")[["norm_DR", "norm_DL"]]
    )
    paralog_ids = load_paralog_ids(config.paralogs)

    if metrics.empty:
        logger.warning("metrics table is empty; writing empty attribution outputs")
        pd.DataFrame().to_csv(config.output_table, sep="\t", index=False)
        write_parquet(pd.DataFrame(
            columns=["group_id", "Systematic ID", "norm_DR", "norm_DL", "component"]
        ), config.output_points)
        return

    table, points_table = attribute_all(metrics, long_table, points, paralog_ids, config)
    table.to_csv(config.output_table, sep="\t", index=False)
    write_parquet(points_table, config.output_points)

    n_incoherent = int(table["is_incoherent"].sum())
    label_counts = table[table["is_incoherent"]]["attribution_label"].value_counts().to_dict()
    logger.success(
        f"{len(table):,} groups scored, {n_incoherent:,} incoherent (z>{config.z_threshold}); "
        f"labels={label_counts}; {len(points_table):,} split points; wrote {config.output_table}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Attribute the cause of coherence incoherence per group")
    parser.add_argument("--metrics", type=Path, required=True, help="A coherence metrics table (one source's Parquet, or the pooled dedup TSV)")
    parser.add_argument("--annotation", type=Path, nargs="+", required=True,
                        help="The matching group_annotation_long.tsv table(s) — one per source, concatenated")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv (DR/DL)")
    parser.add_argument("--paralogs", type=Path, required=True, help="Ensembl paralog export TSV")
    parser.add_argument("--z-threshold", type=float, default=1.0, help="median_pairwise_distance_z above this = incoherent")
    parser.add_argument("--shared-frac-threshold", type=float, default=0.5, help="frac_shared_members >= this triggers the shared-subunit label")
    parser.add_argument("--paralog-frac-threshold", type=float, default=0.5, help="paralog_fraction >= this triggers the paralog-buffered label")
    parser.add_argument("--output-table", type=Path, required=True, help="Output attribution TSV")
    parser.add_argument("--output-points", type=Path, required=True, help="Output incoherent-group split-points Parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run attribution, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = AttributionConfig(
            metrics=args.metrics,
            annotations=tuple(args.annotation),
            fitting_results=args.fitting_results,
            paralogs=args.paralogs,
            output_table=args.output_table,
            output_points=args.output_points,
            z_threshold=args.z_threshold,
            shared_frac_threshold=args.shared_frac_threshold,
            paralog_frac_threshold=args.paralog_frac_threshold,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
