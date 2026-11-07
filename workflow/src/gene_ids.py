"""
Gene Systematic ID Resolution
===============================

Resolves gene names/synonyms to current PomBase systematic IDs, factored out
of the original `workflow/src/utils.py` (DIT_HAP_pipeline). Depends on
`workflow.src.io.read_file` to load the gene metadata table.

Input
-----
- A list of gene identifiers (names, synonyms, or already-current systematic IDs)
- A PomBase gene metadata table (gene_IDs_names_products.tsv) with columns
  gene_systematic_id, gene_name, synonyms, gene_type

Output
------
- A list of the same length with each entry resolved to its current
  systematic ID where a unique match is found, or left unchanged/NaN
  where the id is unknown or ambiguous (logged via print, matching the
  original notebook workflow's use of these log lines for manual review).

Usage
-----
    from gene_ids import update_sys_ids
    resolved = update_sys_ids(["cdc2", "SPBC11B10.09"], gene_meta_file)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 3. Local Imports
from loguru import logger

# 2. Data Processing Imports
import numpy as np
import pandas as pd

from io_table import read_file


# =============================================================================
# CORE LOGIC
# =============================================================================
def _normalize_gene_id(gene: str) -> str:
    """Strip whitespace and normalize a dotted transcript ID to GENE.transcript case."""
    gene = gene.strip()
    if "." in gene:
        gene = gene.split(".")[0].upper() + "." + gene.split(".")[1].lower()
    return gene


def _apply_lookup(gene: str | float, lookup: pd.DataFrame) -> str | float:
    """Keep a unique match; an ambiguous index hit logs a warning and yields NaN."""
    updated = lookup.loc[gene, "gene_systematic_id"]
    if isinstance(updated, str):
        logger.debug(f"{gene} is updated to {updated}")
        return updated
    logger.warning(f"{gene} has multiple updates: {updated}")
    return np.nan


def _resolve_gene_id(
    gene: str | float,
    sys_ids_now: list[str],
    names2ID: pd.DataFrame,
    synonyms2ID: pd.DataFrame,
) -> str | float:
    """Resolve one identifier: NA -> current ID -> gene name -> synonym -> unchanged."""
    if pd.isna(gene):
        logger.debug(f"{gene} is NA")
        return gene
    if gene in sys_ids_now:
        return gene
    for lookup in (names2ID, synonyms2ID):
        if gene in lookup.index:
            return _apply_lookup(gene, lookup)
    logger.debug(f"{gene} is not found in geneid2symbol or synonyms2ID")
    return gene


def update_sys_ids(
    genes: list[str],
    gene_meta_file: Path,
    gene_filter: str = "gene_type == 'protein coding gene'",
) -> list[str | float]:
    """Resolve each gene name/synonym in `genes` to its current systematic ID."""
    gene_meta = read_file(gene_meta_file)
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])

    filtered_genes = gene_meta.query(gene_filter)
    synonyms2ID = (
        filtered_genes.set_index("gene_systematic_id")["synonyms"]
        .str.split(",")
        .explode()
        .str.strip()
        .dropna()
        .reset_index()
        .set_index("synonyms")
    )
    names2ID = (
        filtered_genes.set_index("gene_name")["gene_systematic_id"]
        .drop_duplicates()
        .reset_index()
        .set_index("gene_name")
    )
    sysIDs_now = filtered_genes["gene_systematic_id"].unique().tolist()

    updated_sysIDs = []
    for gene in genes:
        if isinstance(gene, str):
            gene = _normalize_gene_id(gene)
        updated_sysIDs.append(_resolve_gene_id(gene, sysIDs_now, names2ID, synonyms2ID))
    return updated_sysIDs
