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
    from workflow.src.gene_ids import update_sysIDs
    resolved = update_sysIDs(["cdc2", "SPBC11B10.09"], gene_meta_file)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Local Imports
from workflow.src.io import read_file

# =============================================================================
# CORE LOGIC
# =============================================================================
def update_sysIDs(
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
            gene = gene.strip()
            if "." in gene:
                gene = gene.split(".")[0].upper() + "." + gene.split(".")[1].lower()
        if pd.isna(gene):
            updated_sysIDs.append(gene)
            print(f"{gene} is NA")
        elif gene in sysIDs_now:
            updated_sysIDs.append(gene)
        elif gene in names2ID.index:
            updated = names2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        elif gene in synonyms2ID.index:
            updated = synonyms2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        else:
            updated_sysIDs.append(gene)
            print(f"{gene} is not found in geneid2symbol or synonyms2ID")
    return updated_sysIDs
