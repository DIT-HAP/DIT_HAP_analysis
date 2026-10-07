#!/usr/bin/env python3

"""
Gene-Group Coherence Analysis (source-agnostic) — Computation Only
==================================================================

Per-dataset x source: for every group (complex / GO term / ...) whose
DR<threshold members (negative DR = depleted) number between --min-size and
--max-size AND whose total annotated membership is <= --max-term-genes, measures
how tightly its member genes cluster in the 2D DIT-HAP fitness space and tests
that tightness against a genome-wide null via a seeded permutation test.

This is the computation component after the computation-plotting decoupling
refactor (ADR-0001). The visualization is handled by plot_coherence.py.

Fitness "points" are the normalized (DR, DL/10) coordinates of each gene.
Coherence = small median pairwise distance (MPD) among members relative to
random draws of the same number of background genes.

Input
-----
- fitting_results.tsv: the upstream per-gene fitting statistics (see
  coherence/io.py::load_fitting_results for the id/legacy-metric normalization).
- group_annotation_long.tsv: the prepared unified long-table from
  prepare_annotation.py, with the contract columns in
  coherence/sources.py::LONG_TABLE_COLUMNS.
- features.tsv (optional): a gene features table. When given, the per-group
  abundance/conservation uniformity columns are computed here and persisted,
  so the figure is a pure renderer of this table.

Output
------
- coherence_metrics.parquet: one row per surviving group. The columns fall in
  four blocks:
  * identity: source, group_id, group_name, n_annotated_members,
    n_measured_members, n_scored_members, scored_member_names (a list of gene
    names). The three `n_*` are the funnel of this analysis and are nested —
    annotated ⊇ measured (have a fitted DR/DL) ⊇ scored (also DR<threshold) — so
    the two gaps separate a coverage shortfall from non-depleted members.
  * geometry: geom_median_DR, geom_median_DL (the Weiszfeld geometric median,
    NOT a centroid) plus the `{method}` pairwise-distance statistics
    (median_/mean_/std_/min_/max_pairwise_distance).
  * test: for each method in _TABLE_ZSCORE_METHODS, `{method}_z` and
    `{method}_p`; n_permutations; and q_value (BH over the per-source p-values
    of the primary method, median_pairwise_distance).
  * annotations: — when --features is given — abundance_cv / conservation_cv
    each with a `_feature` column naming the feature column actually used, plus
    paralog_fraction (the share of the scored members whose feature-matrix
    `paralog_count` is > 0, so it follows `features.paralog_source` rather than
    naming a paralog source of its own).

Usage
-----
    python compute_coherence.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --source go_macrocomplex \\
        --min-size 3 --max-size 300 --max-term-genes 500 --dr-threshold -0.3 \\
        --n-permutations 1000 --random-state 42 \\
        --output results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        [--features results/1b_features/{pombase_version}/pombe_coding_gene_protein_features.tsv]

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (refactor)
Date:     2026-09-03
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party Imports
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import false_discovery_control

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from coherence.fractions import (  # noqa: E402
    paralog_fraction,
    paralog_ids_from_features,
)
from coherence.io import load_fitting_results, load_long_table  # noqa: E402
from coherence.metrics import coherence_metrics, compute_distance_zscores  # noqa: E402
from io_table import write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# The two pairwise methods the metrics table carries, primary first. Median
# pairwise distance is the primary coherence axis: its z/p drive the `q_value`
# correction, the dedup ranking and every downstream plot. Each method emits
# `{method}_z` / `{method}_p`, so the test columns are named after the statistic
# they test and stay derivable from the method key alone. Both are scored in ONE
# permutation pass, so the null draws are shared rather than redrawn per method.
_TABLE_ZSCORE_METHODS = ("median_pairwise_distance", "mean_pairwise_distance")
_PRIMARY_METHOD = _TABLE_ZSCORE_METHODS[0]

# Feature columns for the optional uniformity panels, tried in order — the first
# one present in the features table wins. Computed here rather than at render
# time so the numbers are persisted in the metrics table and can be checked;
# a plot script that derives them cannot be reviewed or reproduced.
_ABUNDANCE_FEATURE_CANDIDATES = [
    "copies_per_cell_EMM_Proliferating_Cell",
    "copies_per_cell_EMMN_Quiescent_Cell",
    "mean_EMM_Proliferating_Cell_RNA_Abundance",
]
_CONSERVATION_FEATURE_CANDIDATES = ["evolutionary_rate"]

# Feature table key -> metrics-table column name. Each gets a sibling
# `{column}_feature` recording which candidate column it was computed from.
_FEATURE_CV_COLUMNS = {
    "abundance_cv": _ABUNDANCE_FEATURE_CANDIDATES,
    "conservation_cv": _CONSERVATION_FEATURE_CANDIDATES,
}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class CoherenceConfig:
    """Inputs, outputs, and parameters for the gene-group coherence analysis."""
    fitting_results: Path
    annotation: Path
    source: str
    output: Path
    min_size: int = 3
    max_size: int = 300
    max_term_genes: int = 500
    dr_threshold: float = -0.3
    n_permutations: int = 1000
    random_state: int = 42
    features: Path | None = None

    def validate(self) -> None:
        """Raise ValueError if inputs are missing or params invalid, then make output dirs."""
        required = [self.fitting_results, self.annotation]
        if self.features is not None:
            required.append(self.features)
        for path in required:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.min_size < 2:
            raise ValueError(f"min_size must be >= 2 (pairwise distances need 2 points): {self.min_size}")
        if self.max_size < self.min_size:
            raise ValueError(f"max_size ({self.max_size}) must be >= min_size ({self.min_size})")
        if self.max_term_genes < self.min_size:
            raise ValueError(f"max_term_genes ({self.max_term_genes}) must be >= min_size ({self.min_size})")
        if self.n_permutations < 1:
            raise ValueError(f"n_permutations must be >= 1: {self.n_permutations}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC — per-group annotations (biology panels)
# =============================================================================
def first_present_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """The first candidate column present in df, or None if none are."""
    for column in candidates:
        if column in df.columns:
            return column
    return None


def member_feature_cv(
    long_table: pd.DataFrame, features: pd.DataFrame, feature_col: str
) -> dict[str, float]:
    """Per group, coefficient of variation of a numeric feature across its members."""
    # Members missing the feature, or carrying a non-positive value, are dropped;
    # a group left with fewer than 2 usable values gets no entry (NaN after the
    # join). Relative spread is the point — a group whose subunits differ wildly in
    # abundance is not one functional unit in the same sense as a uniform one.
    feature_values = features.set_index("gene_systematic_id")[feature_col].to_dict()
    members = long_table[["group_id", "Systematic ID"]].drop_duplicates().copy()
    members["value"] = members["Systematic ID"].map(feature_values)
    members = members[members["value"].notna() & (members["value"] > 0)]

    grouped = members.groupby("group_id")["value"]
    counts = grouped.size()
    cv = grouped.std(ddof=1) / grouped.mean()
    return cv[counts >= 2].to_dict()


def group_annotations(long_table: pd.DataFrame, features: pd.DataFrame | None) -> pd.DataFrame:
    """Per-group derived annotations, indexed by group_id (empty without features)."""
    # The feature-uniformity CVs only exist when a features table is supplied. Each
    # is computed in one pass over every group (not per surviving group), then
    # left-joined onto the metrics table so a group that failed the size filter
    # simply does not appear.
    annotations = pd.DataFrame()
    if features is None:
        return annotations

    for column, candidates in _FEATURE_CV_COLUMNS.items():
        source_column = first_present_column(features, candidates)
        if source_column is None:
            logger.warning(f"No feature column found for {column} (tried {candidates}); column omitted")
            continue
        annotations[column] = pd.Series(member_feature_cv(long_table, features, source_column))
        # Which candidate the CV was actually taken over. The candidates are tried
        # in order and the first present one wins, so the same `abundance_cv`
        # column means a different quantity depending on which features table was
        # passed; recording the winner is what keeps the number interpretable.
        annotations[f"{column}_feature"] = source_column
    return annotations


# =============================================================================
# CORE LOGIC — coherence per group
# =============================================================================
def measured_member_counts(long_table: pd.DataFrame, measured_ids: set[str]) -> dict[str, int]:
    """Per group_id, how many annotated members have a fitted DR/DL."""
    # The middle rung of the funnel: annotated -> measured -> scored. Kept separate
    # from the DR filter so a shortfall in (
    # n_annotated_members - n_measured_members) — the member has no fitness data at
    # all — can be told apart from one in (
    # n_measured_members - n_scored_members) — it has data but is not depleted. The
    # two have different causes (insertion coverage vs biology) and only the second
    # is a statement about the complex.
    members = long_table[["group_id", "Systematic ID"]].drop_duplicates()
    measured = members[members["Systematic ID"].isin(measured_ids)]
    return measured.groupby("group_id")["Systematic ID"].nunique().to_dict()


def build_groups(
    background: pd.DataFrame,
    long_table: pd.DataFrame,
    min_group_size: int,
    max_group_size: int,
    max_term_genes: int,
) -> dict[str, pd.DataFrame]:
    """Map surviving groups (keyed on group_id) -> their DR<threshold member rows."""
    merged = long_table.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for group_id, grp in merged.groupby("group_id"):
        grp = grp.drop_duplicates(subset="Systematic ID")
        n_total = int(grp["n_annotated_members"].iloc[0])
        if n_total > max_term_genes:
            continue
        if min_group_size <= len(grp) <= max_group_size:
            groups[group_id] = grp
    logger.info(
        f"{merged['group_id'].nunique():,} groups with >=1 background member -> "
        f"{len(groups):,} with {min_group_size} <= size <= {max_group_size} "
        f"and n_annotated_members <= {max_term_genes}"
    )
    return groups


def compute_coherence_table(
    groups: dict[str, pd.DataFrame],
    background_points: np.ndarray,
    background_index: dict[str, int],
    n_permutations: int,
    random_state: int,
    annotations: pd.DataFrame,
    n_measured: dict[str, int],
    paralog_ids: set[str] | None = None,
) -> pd.DataFrame:
    """One coherence row per group: identity + metrics + permutation z-scores."""
    rows = []
    for group_id, grp in groups.items():
        member_ids = grp["Systematic ID"].tolist()
        member_points = background_points[[background_index[gid] for gid in member_ids]]

        zscores = compute_distance_zscores(
            member_points,
            background_points,
            _TABLE_ZSCORE_METHODS,
            n_permutations=n_permutations,
            random_state=random_state,
        )

        row = {
            "source": grp["source"].iloc[0],
            "group_id": group_id,
            "group_name": grp["group_name"].iloc[0],
            # The funnel, widest first: annotated ⊇ measured ⊇ scored. n_scored_members
            # is the permutation test's draw size, not the group's size.
            "n_annotated_members": int(grp["n_annotated_members"].iloc[0]),
            "n_measured_members": int(n_measured.get(group_id, 0)),
            "n_scored_members": len(grp),
            # A real list column (not a joined string): the Parquet round-trips it
            # without a parse step, and only the final human-facing TSVs join it.
            "scored_member_names": sorted(grp["Name"].dropna().astype(str)),
            **coherence_metrics(member_points),
        }
        # Keyed on the systematic ids, not the `Name` display column above: the
        # feature matrix is indexed by gene_systematic_id. Left out entirely when
        # no features table was passed, like the CV columns.
        if paralog_ids is not None:
            row["paralog_fraction"] = paralog_fraction(member_ids, paralog_ids)
        # Test columns are generated off the method key, so `{method}` is the
        # observed statistic and `{method}_z` / `{method}_p` are its test.
        for method, (z, p) in zscores.items():
            row[f"{method}_z"] = z
            row[f"{method}_p"] = p
        row["n_permutations"] = n_permutations
        rows.append(row)

    table = pd.DataFrame(rows)
    if not table.empty:
        # BH over the primary method's p-values, per source (each source is its own
        # hypothesis family; combine_metrics.py does not re-correct).
        table["q_value"] = false_discovery_control(
            table[f"{_PRIMARY_METHOD}_p"].to_numpy(), method="bh"
        )
        table = table.merge(annotations, left_on="group_id", right_index=True, how="left")
        table = table.sort_values(f"{_PRIMARY_METHOD}_z").reset_index(drop=True)
    return table


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: CoherenceConfig) -> None:
    """Load -> filter -> per-group coherence + permutation test -> Parquet."""
    config.validate()

    # One read of fitting_results.tsv serves both layers: `fitting` is every gene
    # with a fitted DR/DL (the n_measured_members layer), `background` the
    # DR-filtered cloud the permutation null is drawn from.
    fitting = load_fitting_results(config.fitting_results)
    long_table = load_long_table(config.annotation)
    background = fitting[fitting["DR"] < config.dr_threshold].reset_index(drop=True)
    logger.info(
        f"fitting_results.tsv: {len(fitting):,} fitted genes -> "
        f"{len(background):,} background genes with DR < {config.dr_threshold}"
    )

    # Genome-wide background point cloud + Systematic ID -> row index map.
    background_points = background[["norm_DR", "norm_DL"]].to_numpy(dtype=float)
    background_index = {gid: i for i, gid in enumerate(background["Systematic ID"])}

    groups = build_groups(
        background, long_table, config.min_size, config.max_size, config.max_term_genes
    )
    features = pd.read_csv(config.features, sep="\t") if config.features is not None else None
    annotations = group_annotations(long_table, features)
    table = compute_coherence_table(
        groups, background_points, background_index,
        config.n_permutations, config.random_state, annotations,
        measured_member_counts(long_table, set(fitting["Systematic ID"])),
        paralog_ids=paralog_ids_from_features(features) if features is not None else None,
    )

    # Write output as Parquet
    write_parquet(table, config.output)

    # Not the `coherent` cohort (config.coherence.coherent_*), which is defined on the
    # z/q pair; this is the raw sign split, kept as a sanity number in the log.
    n_tighter = (
        int((table[f"{_PRIMARY_METHOD}_z"] < 0).sum()) if not table.empty else 0
    )
    logger.success(
        f"[{config.source}] Coherence: {len(table):,} groups scored, {n_tighter:,} tighter than random (z<0); "
        f"wrote {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute gene-group coherence metrics")
    parser.add_argument("--fitting-results", type=Path, required=True,
                        help="Upstream fitting_results.tsv (systematic id as index col 0)")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="Prepared group_annotation_long.tsv (long-table from prepare_annotation.py)")
    parser.add_argument("--source", type=str, required=True,
                        help="Grouping-database source name (fan-out dimension)")
    parser.add_argument("--min-size", type=int, default=3,
                        help="Minimum DR<threshold members per group")
    parser.add_argument("--max-size", type=int, default=300,
                        help="Maximum DR<threshold members per group")
    parser.add_argument("--max-term-genes", type=int, default=500,
                        help="Drop groups whose total annotated membership (n_annotated_members) exceeds this")
    parser.add_argument("--dr-threshold", type=float, default=-0.3,
                        help="Keep genes with DR < this (negative = depleted)")
    parser.add_argument("--n-permutations", type=int, default=1000,
                        help="Permutation null draws")
    parser.add_argument("--random-state", type=int, default=42,
                        help="Permutation RNG seed")
    parser.add_argument("--features", type=Path, default=None,
                        help="Optional gene features TSV; adds the abundance/conservation CV and paralog-fraction columns")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output coherence metrics Parquet")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CoherenceConfig(
            fitting_results=args.fitting_results,
            annotation=args.annotation,
            source=args.source,
            output=args.output,
            min_size=args.min_size,
            max_size=args.max_size,
            max_term_genes=args.max_term_genes,
            dr_threshold=args.dr_threshold,
            n_permutations=args.n_permutations,
            random_state=args.random_state,
            features=args.features,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
