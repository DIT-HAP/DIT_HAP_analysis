#!/usr/bin/env python3

"""
Per-View Fractions (paralog + moonlighting) — Computation Only
==============================================================

For one view of the coherence analysis (one source, `combined`, or `dedup`),
collects the two descriptive fractions the view carries and the gene-breadth
profile behind the second one:

- paralog_fraction: per term, the share of its scored members whose feature-matrix
  `paralog_count` is > 0 (already a column of the metrics table — computed by
  compute_coherence.py through coherence/fractions.py).
- gene breadth: per gene, how many of the view's de-duplication groups it sits in.
  A view's groups are the clusters its terms belong to, so the same gene gets a
  different count in a different view; the cut is that per-view distribution's own
  top quartile (coherence/fractions.py::MOONLIGHTING_QUANTILE), and a gene above it
  is flagged `is_moonlighting`.
- moonlighting_fraction: per term, the share of its own members that are flagged.

These are descriptive per-view quantities, not a diagnosis of WHY a group is
dispersed: every view has them, and each view is read against its own breadth
profile rather than against one global number.

Input
-----
- --metrics: combined/coherence_metrics.parquet. Every view is sliced from it, so
  a per-source run and the combined run read the same numbers. It must carry the
  list column `scored_member_names`; a TSV's joined string would have to be
  re-parsed, which is exactly what the list column exists to avoid.
- --dedup-terms: dedup/coherence_terms_deduplicated.tsv, for each term's
  `redundancy_cluster` and `is_representative`.
- --group-members: dedup/coherence_group_members_long.tsv, the (group, gene) long
  table the breadth counts run over.
- --view: a source name, `combined`, or `dedup`.

Output
------
- --output-terms: one row per term of the view — source, group_id, group_name,
  redundancy_cluster, n_scored_members, paralog_fraction, moonlighting_fraction,
  moonlighting_cut_n_groups (the view's cut, carried so the fractions are readable
  on their own).
- --output-genes: one row per gene in the view's groups — gene, n_groups,
  is_moonlighting.

Usage
-----
    python compute_view_fractions.py \\
        --metrics results/3a_coherence/{dataset}/combined/coherence_metrics.parquet \\
        --dedup-terms results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv \\
        --group-members results/3a_coherence/{dataset}/dedup/coherence_group_members_long.tsv \\
        --view go_macrocomplex \\
        --output-terms results/3a_coherence/{dataset}/go_macrocomplex/view_fractions.tsv \\
        --output-genes results/3a_coherence/{dataset}/go_macrocomplex/view_gene_breadth.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 4.6 (implementation)
Date:     2026-10-05
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
from coherence.fractions import MOONLIGHTING_QUANTILE, moonlighting_fraction, view_breadth  # noqa: E402
from io_table import read_file  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The two views that are not a single source. Everything else `--view` accepts is
# a source name, checked against the metrics table rather than against a list, so
# adding a source needs no edit here.
_COMBINED = "combined"
_DEDUP = "dedup"

_KEY_COLUMNS = ["source", "group_id"]

# Column order of the two outputs. `moonlighting_cut_n_groups` rides on every term
# row (like the metrics table's `*_cv_feature` columns): the fraction is only
# readable against the cut it was taken with.
_TERM_COLUMNS = ["source", "group_id", "group_name", "redundancy_cluster",
                 "n_scored_members", "paralog_fraction", "moonlighting_fraction",
                 "moonlighting_cut_n_groups"]
_GENE_COLUMNS = ["gene", "n_groups", "is_moonlighting"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ViewFractionsConfig:
    """Inputs, the view to compute, and the two outputs."""
    metrics: Path
    dedup_terms: Path
    group_members: Path
    view: str
    output_terms: Path
    output_genes: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then make output dirs."""
        for path in [self.metrics, self.dedup_terms, self.group_members]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_terms, self.output_genes]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def view_terms(metrics: pd.DataFrame, dedup: pd.DataFrame, view: str) -> pd.DataFrame:
    """The view's terms from the metrics table, each with its de-duplication cluster."""
    clusters = dedup[_KEY_COLUMNS + ["redundancy_cluster", "is_representative"]].drop_duplicates(
        subset=_KEY_COLUMNS
    )
    if view == _COMBINED:
        selected = metrics
    elif view == _DEDUP:
        # The representative set: one term per cluster. Selected from the dedup
        # table's own flag, so the dedup view and the cluster list agree by
        # construction and there is no second representative definition here.
        selected = metrics.merge(
            clusters.loc[clusters["is_representative"], _KEY_COLUMNS], on=_KEY_COLUMNS, how="inner"
        )
    else:
        selected = metrics[metrics["source"] == view]

    if selected.empty:
        raise ValueError(f"no terms for view {view!r} in the metrics table")
    merged = selected.merge(clusters[_KEY_COLUMNS + ["redundancy_cluster"]], on=_KEY_COLUMNS, how="left")
    n_unclustered = int(merged["redundancy_cluster"].isna().sum())
    if n_unclustered:
        logger.warning(f"{n_unclustered} term(s) of view {view!r} carry no cluster; they contribute no breadth")
    return merged


def fractions_for_view(terms: pd.DataFrame, group_members: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """The view's per-term fractions and per-gene breadth, over the view's own clusters."""
    # The gene sets stay each cluster's own, and the cut is this view's own: two
    # views differ only in which groups they count, never in how a group is counted.
    clusters = set(terms["redundancy_cluster"].dropna())
    per_gene, cut = view_breadth(group_members, clusters)
    flagged = set(per_gene.loc[per_gene["is_moonlighting"], "gene"])

    table = terms.copy()
    table["moonlighting_fraction"] = [
        moonlighting_fraction(members, flagged) for members in table["scored_member_names"]
    ]
    table["moonlighting_cut_n_groups"] = cut
    return table, per_gene, cut


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: ViewFractionsConfig) -> None:
    """Load -> slice the view -> per-term fractions + per-gene breadth -> two TSVs."""
    config.validate()
    metrics = read_file(config.metrics)
    missing = [column for column in _KEY_COLUMNS + ["group_name", "n_scored_members", "scored_member_names"]
               if column not in metrics.columns]
    if missing:
        raise ValueError(f"metrics table missing required column(s) {missing}")
    if "paralog_fraction" not in metrics.columns:
        logger.warning(
            "metrics table has no paralog_fraction column (run compute_coherence with --features); "
            "the terms table will carry NaN there"
        )
        metrics = metrics.assign(paralog_fraction=pd.NA)

    terms = view_terms(metrics, read_file(config.dedup_terms), config.view)
    group_members = read_file(config.group_members)
    table, per_gene, cut = fractions_for_view(terms, group_members)

    table[_TERM_COLUMNS].to_csv(config.output_terms, sep="\t", index=False)
    logger.success(
        f"[{config.view}] {len(table):,} terms, {len(per_gene):,} genes in "
        f"{table['redundancy_cluster'].nunique():,} groups; cut = {cut:g} groups/gene "
        f"(top {1 - MOONLIGHTING_QUANTILE:.0%}), "
        f"{int(per_gene['is_moonlighting'].sum()):,} genes above it; wrote {config.output_terms}"
    )
    per_gene[_GENE_COLUMNS].to_csv(config.output_genes, sep="\t", index=False)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Per-view paralog / moonlighting fractions and gene breadth")
    parser.add_argument("--metrics", type=Path, required=True, help="combined/coherence_metrics.parquet")
    parser.add_argument("--dedup-terms", type=Path, required=True, help="dedup/coherence_terms_deduplicated.tsv")
    parser.add_argument("--group-members", type=Path, required=True, help="dedup/coherence_group_members_long.tsv")
    parser.add_argument("--view", required=True, help="A source name, 'combined', or 'dedup'")
    parser.add_argument("--output-terms", type=Path, required=True, help="Output per-term fractions TSV")
    parser.add_argument("--output-genes", type=Path, required=True, help="Output per-gene breadth TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, compute the view's fractions, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ViewFractionsConfig(
            metrics=args.metrics,
            dedup_terms=args.dedup_terms,
            group_members=args.group_members,
            view=args.view,
            output_terms=args.output_terms,
            output_genes=args.output_genes,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
