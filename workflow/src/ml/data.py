"""
ML Modeling-Data Assembly
===========================

Shared modeling-data assembly for the AutoML pipeline: merges the per-gene
feature matrix with curve-fit targets + cluster labels and applies the notebook's
DR < threshold filter (negative DR = depleted). This is the target- and
mode-INDEPENDENT part of
machine_learning_analysis.ipynb, factored out so the four target x mode AutoML
jobs share one prepared table instead of each re-merging (byte-faithful to
train_automl.py's former load_modeling_data).

Input
-----
- Per-gene feature matrix (results/features/{version}/pombe_coding_gene_protein_features.tsv)
- Curated final_clusters.tsv (Systematic ID, A, DR, DL, cluster)

Output
------
- One merged, DR-filtered DataFrame (Systematic_ID + features + A/DR/DL/DIT_HAP_cluster)

Usage
-----
    from ml.data import load_modeling_data
    data = load_modeling_data(feature_matrix, final_clusters, dr_filter=-0.3)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Upstream flipped the DR sign convention on 2026-09-17: negative DR is now the
# depleted end, so the notebook's "meaningfully depleted" filter reads DR < -0.3
# instead of DR > 0.3 (same gene set, mirrored threshold).
DR_FILTER = -0.3


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_modeling_data(feature_matrix: Path, final_clusters: Path, dr_filter: float = DR_FILTER) -> pd.DataFrame:
    """Merge feature matrix + targets, filter to DR < threshold (notebook behavior)."""
    features = pd.read_csv(feature_matrix, sep="\t")
    targets = pd.read_csv(final_clusters, sep="\t").rename(
        columns={"Systematic ID": "Systematic_ID", "cluster": "DIT_HAP_cluster"}
    )[["Systematic_ID", "A", "DR", "DL", "DIT_HAP_cluster"]]

    data = (
        pd.merge(features, targets, left_on="gene_systematic_id", right_on="Systematic_ID", how="left")
        .drop(columns=["Systematic_ID"])
        .rename(columns={"gene_systematic_id": "Systematic_ID"})
        .query(f"DR < {dr_filter}")
    )
    logger.info(f"Modeling data (DR < {dr_filter}): {data.shape}")
    return data
