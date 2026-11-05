"""
Enrichment Table Schema
=======================

Column names and required-column sets for the two tables the enrichment stage
reads and writes. Kept dependency-free (standard library only) so the plotting
scripts can import it in the cnsplots environment without pulling goatools in
through ``cluster_enrichment``.

Two spellings cross the boundary and both live here: ``final_clusters.tsv``
carries the lower-case ``cluster`` id, while every enrichment result frame
carries the capitalised ``Cluster`` column that the concat helpers add.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from enum import StrEnum

# =============================================================================
# COLUMN CONTRACTS
# =============================================================================


class ClusterCol(StrEnum):
    """Columns read from the curated ``final_clusters.tsv``."""

    GENE_ID = "Systematic ID"
    CLUSTER = "cluster"


class EnrichmentCol(StrEnum):
    """Columns shared by the GO/FYPO/MONDO enrichment result frames."""

    CLUSTER = "Cluster"
    NAMESPACE = "namespace"
    TERM_ID = "term_id"
    TERM = "term"
    P_FDR = "p_fdr"
    GENE_RATIO = "gene_ratio"
    STUDY_COUNT = "study_count"
    STUDY_ITEMS = "study_items"
    POP_COUNT = "pop_count"


# Validated once at each script's I/O boundary.
FINAL_CLUSTER_COLUMNS = [ClusterCol.GENE_ID, ClusterCol.CLUSTER]

ENRICHMENT_COLUMNS = [
    EnrichmentCol.CLUSTER,
    EnrichmentCol.NAMESPACE,
    EnrichmentCol.TERM_ID,
    EnrichmentCol.TERM,
    EnrichmentCol.P_FDR,
    EnrichmentCol.GENE_RATIO,
]

GO_FILTER_COLUMNS = [
    EnrichmentCol.CLUSTER,
    EnrichmentCol.NAMESPACE,
    EnrichmentCol.TERM_ID,
    EnrichmentCol.POP_COUNT,
]
