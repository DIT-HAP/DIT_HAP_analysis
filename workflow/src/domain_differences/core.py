"""
Intra-gene DR Heterogeneity (Domain-Difference) — Core Logic
==============================================================

Shared constants, loaders, and statistics functions for the domain-differences
stage. Per-dataset: surfaces genes whose in-gene insertions have heterogeneous
positional distribution as a proxy for functional (sub-gene) domains. For every
gene with a high gene-level DR (> DR_THRESHOLD = 0.15), it positions each in-gene
insertion along the CDS (insertion_fraction, 0 = start codon, 1 = stop codon) and
reports per-gene distribution statistics (count, mean, std of insertion_fraction).

Ported from the deterministic section of
DIT_HAP_pipeline/workflow/notebooks/genes_with_domain_differences.ipynb, and
factored out of the original single-script port so the stage can be split into
independent Snakemake rules (prepare -> compute), each re-runnable on its own.

IMPORTANT — notebook vs. this module
-------------------------------------
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

Usage
-----
    from workflow.src.domain_differences.core import (
        load_gene_level, load_insertion_annotations, filter_high_dr_genes,
        compute_insertion_fraction, compute_domain_candidate_stats,
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

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
# HELPERS
# =============================================================================
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

    ``n_insertions`` uses group ``size``, i.e. the count of in-gene insertions
    for the gene; a rare insertion with a zero-length CDS span has a NaN
    insertion_fraction (see compute_insertion_fraction) and is still counted
    here, while ``mean``/``std`` skip it — so on such genes n_insertions can
    exceed the number of positioned values. The real HD_DIT_HAP in-gene set has
    no zero-span rows, so this does not arise in practice.
    """
    # drop_duplicates on the gene key guards both merges below: gene-level
    # Systematic ID is unique in the real release, but a duplicated key would
    # otherwise turn each merge into a cross-product and inflate per-gene counts
    # / emit duplicate output rows.
    genes = (
        high_dr_genes[["Systematic ID", "Name", "DR"]]
        .rename(columns={"DR": "gene_DR"})
        .drop_duplicates("Systematic ID")
    )

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
