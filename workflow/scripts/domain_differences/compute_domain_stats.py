#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Intra-gene DR Heterogeneity (Domain-Difference) Candidate Statistics
====================================================================

Per-dataset: surfaces genes whose in-gene insertions have heterogeneous
positional distribution as a proxy for functional (sub-gene) domains. For every
gene with a high gene-level DR (> DR_THRESHOLD = 0.15), it positions each in-gene
insertion along the CDS (insertion_fraction, 0 = start codon, 1 = stop codon) and
reports per-gene distribution statistics (count, mean, std of insertion_fraction).

Ported from the deterministic section of
DIT_HAP_pipeline/workflow/notebooks/genes_with_domain_differences.ipynb.

IMPORTANT — notebook vs. this script
------------------------------------
The source notebook is a *visualization* notebook: it (1) loads the release
tables, (2) filters to in-gene insertions with the exact IN_GENE_FILTER quirk,
(3) selects genes with gene-level DR (legacy header `um`) > 0.15, and (4)
scatter-plots each in-gene insertion (Residue_affected vs DR) into per-gene PDF
panels. It computes NO per-gene statistics table and writes no deterministic
data output — the "domain" call is done by human review of the PDFs plus a
hand-curated spreadsheet (resources/non_esential_domain_candidates.xlsx), which
is out of scope for a reproducible rule.

compute_domain_candidate_stats is therefore a NEW deterministic distillation of
that visualization: instead of positioning by Residue_affected (protein residue
index, used only for the notebook's x-axis), it positions each insertion by the
release table's nucleotide distances via
insertion_fraction = Distance_to_start_codon / (Distance_to_start_codon +
Distance_to_stop_codon), clamped to [0, 1], then aggregates per gene. The gene
selection (DR > 0.15) and the in-gene filter are byte-faithful to the notebook.

Input
-----
- Insertion-level fitting_results.tsv (MultiIndex [Chr, Coordinate, Strand,
  Target]) — defines the total insertion universe (unique per coordinate).
- Insertion-level annotations.tsv(.gz) (same MultiIndex, plus Type /
  Distance_to_start_codon / Distance_to_stop_codon / Systematic ID / Name).
  This table can carry duplicate index entries (multiple Features per
  coordinate, e.g. CDS + overlapping intron); duplicates are collapsed
  (any Feature passing IN_GENE_FILTER wins) before reindexing onto
  fitting_results, matching the notebook's merge universe without
  double-counting a single insertion.
- Gene-level fitting_results.tsv (Systematic ID, Name, DR, ...). Legacy releases
  ship the pre-rename um/lam headers; normalized to DR/DL on load.

Output
------
- domain_candidate_stats.tsv: one row per high-DR gene that has >=1 in-gene
  insertion, columns [Systematic ID, Name, n_insertions,
  mean_insertion_fraction, std_insertion_fraction, gene_DR], sorted by
  std_insertion_fraction descending (single-insertion genes have NaN std and
  sort last).

Usage
-----
    python compute_domain_stats.py \\
        --fitting-results .../release/insertion_level/fitting_results.tsv \\
        --annotations .../release/insertion_level/annotations.tsv.gz \\
        --gene-level .../release/gene_level/fitting_results.tsv \\
        --output-stats results/domain_differences/{dataset}/domain_candidate_stats.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-21
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
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Gene-selection cutoff: keep genes with gene-level DR strictly greater than
# this value. Byte-faithful to the notebook's `gene_statistics.query("um > 0.15")`
# (um is the legacy header for DR).
DR_THRESHOLD = 0.15

# Byte-faithful to the source notebook's in_gene_insertions query (identical to
# workflow/scripts/coverage/compute_coverage_stats.py): an insertion counts as
# in-gene only if annotated non-intergenic AND at least 5bp upstream of the stop
# codon (the >4 threshold, not >=5, is the notebook's own quirk — kept verbatim).
IN_GENE_FILTER = "Type != 'Intergenic region' and Distance_to_stop_codon > 4"

# Legacy -> current metric column names (same quirk as
# workflow/scripts/utr/classify_utr_insertions.py and
# workflow/scripts/coverage/compute_coverage_stats.py): some datasets'
# fitting_results.tsv still ship the pre-rename um/lam headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class DomainConfig:
    """Inputs, output, and threshold for the domain-difference candidate analysis."""
    fitting_results: Path
    annotations: Path
    gene_level: Path
    output_stats: Path
    dr_threshold: float = DR_THRESHOLD

    def validate(self) -> None:
        """Raise ValueError if any input is missing or the threshold is invalid, then ensure the output dir exists."""
        for path in [self.fitting_results, self.annotations, self.gene_level]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.dr_threshold < 0:
            raise ValueError(f"dr_threshold must be non-negative, got {self.dr_threshold}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def _normalize_legacy_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Rename legacy um/lam metric columns to DR/DL when the new names are absent."""
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in df.columns and new not in df.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        df = df.rename(columns=rename)
    return df


def load_gene_level(gene_level_path: Path) -> pd.DataFrame:
    """Load gene-level fitting statistics, normalizing legacy um/lam -> DR/DL columns."""
    gene_result = pd.read_csv(gene_level_path, sep="\t")
    return _normalize_legacy_metrics(gene_result)


def resolve_duplicate_annotations(annotations: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate-indexed annotation rows to one row per index value.

    The insertion-level annotations table can carry duplicate index entries
    (multiple Features per coordinate, e.g. an overlapping CDS + intron
    record). Among duplicates sharing an index value, the row that passes
    IN_GENE_FILTER wins if any duplicate does (matching the notebook's merge
    universe: a single insertion is one point, positioned by whichever
    annotation makes it in-gene). Uses an explicit ``kind="stable"`` sort so
    the tie-break is deterministic; when no duplicate passes (or there are no
    duplicates) the first row in original file order is kept. Same pattern as
    workflow/scripts/coverage/compute_coverage_stats.py.
    """
    if not annotations.index.duplicated().any():
        return annotations

    n_dup = int(annotations.index.duplicated().sum())
    logger.info(f"Collapsing {n_dup} duplicate-indexed annotation rows (keep in-gene pass if any)")
    passes = annotations.eval(IN_GENE_FILTER)
    return (
        annotations.assign(_passes=passes)
        .sort_values("_passes", ascending=False, kind="stable")
        .loc[lambda df: ~df.index.duplicated(keep="first")]
        .drop(columns="_passes")
    )


def load_insertion_annotations(fitting_results_path: Path, annotations_path: Path) -> pd.DataFrame:
    """Load insertion-level annotations reindexed onto the fitting_results insertion universe.

    fitting_results defines the set of real insertions (unique per coordinate);
    the annotations carry the in-gene call + CDS distances. Duplicate annotation
    rows are collapsed (see resolve_duplicate_annotations) before reindexing so a
    single insertion contributes exactly one point, byte-faithful to the
    notebook's `pd.merge(insertion_statistics, annotations, how="left")` universe
    without inflating per-gene counts. fitting_results' own metric columns (A/DR/
    DL) are not needed for the positional stats and are intentionally not joined.
    """
    fitting_results = pd.read_csv(fitting_results_path, sep="\t", index_col=[0, 1, 2, 3])
    annotations = pd.read_csv(annotations_path, sep="\t", index_col=[0, 1, 2, 3])
    annotations = resolve_duplicate_annotations(annotations)
    return annotations.reindex(fitting_results.index)


# =============================================================================
# CORE LOGIC — primitives (unit-tested)
# =============================================================================
def filter_high_dr_genes(gene_result: pd.DataFrame, dr_threshold: float = DR_THRESHOLD) -> pd.DataFrame:
    """Keep genes whose gene-level DR is strictly greater than ``dr_threshold``.

    Byte-faithful to the notebook's `gene_statistics.query("um > 0.15")`. NaN DR
    rows are dropped (they fail the strict `>` comparison), so uncovered genes
    never enter the candidate set.
    """
    return gene_result[gene_result["DR"] > dr_threshold]


def compute_insertion_fraction(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``df`` with an ``insertion_fraction`` column in [0, 1].

    insertion_fraction = Distance_to_start_codon / (Distance_to_start_codon +
    Distance_to_stop_codon), i.e. the insertion's fractional position along the
    CDS (0 at the start codon, 1 at the stop codon), clamped to [0, 1]. A
    zero-length span (both distances 0) yields NaN rather than +/-inf. The input
    frame is not mutated.
    """
    df = df.copy()
    denom = df["Distance_to_start_codon"] + df["Distance_to_stop_codon"]
    frac = df["Distance_to_start_codon"] / denom.where(denom != 0, np.nan)
    df["insertion_fraction"] = frac.clip(lower=0, upper=1)
    return df


def compute_domain_candidate_stats(
    in_gene_insertions: pd.DataFrame,
    high_dr_genes: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate in-gene insertions into per-gene positional-spread statistics.

    Restricts ``in_gene_insertions`` to the ``high_dr_genes`` set via an inner
    join on ``Systematic ID`` (so only genes with gene-level DR > threshold are
    reported, and only those with >=1 in-gene insertion appear), then computes
    per gene:
      - n_insertions            (count of in-gene insertions)
      - mean_insertion_fraction (mean position along the CDS)
      - std_insertion_fraction  (sample std, ddof=1; NaN for single-insertion genes)
      - gene_DR                 (the gene-level DR carried from high_dr_genes)

    Rows are sorted by std_insertion_fraction descending (genes whose insertions
    are most positionally spread — the strongest intra-gene heterogeneity /
    domain-difference candidates — first); single-insertion genes have NaN std
    and sort last. ``high_dr_genes`` must carry ``Systematic ID``, ``Name``, and
    ``DR`` columns; ``in_gene_insertions`` must carry ``Systematic ID`` and
    ``insertion_fraction``.
    """
    genes = high_dr_genes[["Systematic ID", "Name", "DR"]].rename(columns={"DR": "gene_DR"})

    restricted = in_gene_insertions[["Systematic ID", "insertion_fraction"]].merge(
        genes[["Systematic ID"]], on="Systematic ID", how="inner"
    )

    per_gene = (
        restricted.groupby("Systematic ID")["insertion_fraction"]
        .agg(
            n_insertions="size",
            mean_insertion_fraction="mean",
            std_insertion_fraction="std",
        )
        .reset_index()
    )

    stats = per_gene.merge(genes, on="Systematic ID", how="left")
    stats = stats[[
        "Systematic ID", "Name", "n_insertions",
        "mean_insertion_fraction", "std_insertion_fraction", "gene_DR",
    ]]
    return stats.sort_values(
        "std_insertion_fraction", ascending=False, na_position="last"
    ).reset_index(drop=True)


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: DomainConfig) -> None:
    """Load -> select high-DR genes -> position in-gene insertions -> per-gene stats -> TSV."""
    config.validate()

    gene_result = load_gene_level(config.gene_level)
    annotations = load_insertion_annotations(config.fitting_results, config.annotations)

    high_dr = filter_high_dr_genes(gene_result, config.dr_threshold)
    logger.info(f"Genes with gene-level DR > {config.dr_threshold}: {len(high_dr):,}")

    in_gene = annotations.query(IN_GENE_FILTER)
    logger.info(f"In-gene insertions (IN_GENE_FILTER): {len(in_gene):,}")
    in_gene = compute_insertion_fraction(in_gene)

    stats = compute_domain_candidate_stats(in_gene, high_dr)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    n_multi = int((stats["n_insertions"] > 1).sum())
    logger.success(
        f"Domain candidates: {len(stats):,} high-DR genes with in-gene insertions "
        f"({n_multi:,} with >1 insertion), {int(stats['n_insertions'].sum()):,} insertions total "
        f"-> {config.output_stats}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute intra-gene DR heterogeneity (domain-difference) candidate statistics")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.tsv(.gz)")
    parser.add_argument("--gene-level", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--dr-threshold", type=float, default=DR_THRESHOLD, help="Gene-level DR selection cutoff (default: 0.15)")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output per-gene domain candidate stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DomainConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_level=args.gene_level,
            output_stats=args.output_stats,
            dr_threshold=args.dr_threshold,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
