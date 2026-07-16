#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ML Feature/Target Preparation
=============================

Merges the per-gene feature matrix with DIT-HAP curve-fit targets + cluster
labels, then applies per-feature transformations (PowerTransformer /
StandardScaler / one-hot / binary), emitting transformed feature/target tables
for the `all`, `DR_gt_p35`, `DR_le_p35`, and `nonWT` splits. Deterministic port
of DIT_HAP_pipeline/workflow/notebooks/machine_learning_data_preparation.ipynb
(the canonical version — not the obsolete Feature_organization copy).

No imputation and no train/test split: missing values are dropped
(`dropna(how='any')`), and splitting is deferred to the AutoML stage (Task 7).

Input
-----
- Per-gene feature matrix (results/features/{version}/pombe_coding_gene_protein_features.tsv)
- Curated final_clusters.tsv (Systematic ID, A, DR, DL, revised_cluster)

Output
------
- all_features_with_target_values.csv, missing_value_analysis.csv
- {split}_transformed_{features,targets,features_and_targets}.csv for
  split in {all, DR_gt_p35, DR_le_p35, nonWT}

Usage
-----
    python prepare_features_targets.py \\
        --feature-matrix results/features/2025-10-01/pombe_coding_gene_protein_features.tsv \\
        --final-clusters resources/curated/final_clusters.tsv \\
        --output-dir results/ml/features_targets/{dataset}

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from sklearn.preprocessing import PowerTransformer, StandardScaler

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
DR_SPLIT_THRESHOLD = 0.35
WT_CLUSTER = 9
TARGET_FEATURES = ["A", "DR", "DL", "DIT_HAP_cluster"]

# Per-feature transformation map (byte-faithful to the canonical notebook).
# Commented entries in the notebook (copies_per_cell_*, t1/2) are excluded here too.
SELECTED_FEATURES_AND_TRANSFORMATIONS = {
    "Systematic_ID": "set_index",
    "Chromosome": "OneHotEncoder",
    "Strand": "binary_encoding",
    "Start": "PowerTransformer",
    "End": "PowerTransformer",
    "Abs_distance_from_telomere": "PowerTransformer",
    "Relative_distance_from_telomere": "PowerTransformer",
    "Abs_distance_from_centromere": "PowerTransformer",
    "Relative_distance_from_centromere": "PowerTransformer",
    "Gene_length": "PowerTransformer",
    "GC_content_of_gene": "StandardScaler",
    "CDS_number": "PowerTransformer",
    "GC_content_of_CDS": "StandardScaler",
    "Fraction_of_CDS": "PowerTransformer",
    "GC3": "StandardScaler",
    "Intron_number": "PowerTransformer",
    "GC_content_of_intron": "PowerTransformer",
    "Total_intron_length": "PowerTransformer",
    "Average_intron_length": "PowerTransformer",
    "Length_of_first_intron": "PowerTransformer",
    "GC_contents_of_first_intron": "PowerTransformer",
    "ENC": "PowerTransformer",
    "Peptide_length": "PowerTransformer",
    "Primary_peptide_length": "PowerTransformer",
    "mean_EMM_Nitrogen_Starved_Cell_RNA_Abundance": "PowerTransformer",
    "mean_EMM_Proliferating_Cell_RNA_Abundance": "PowerTransformer",
    "std_EMM_Nitrogen_Starved_Cell_RNA_Abundance": "PowerTransformer",
    "std_EMM_Proliferating_Cell_RNA_Abundance": "PowerTransformer",
    "cv_EMM_Nitrogen_Starved_Cell_RNA_Abundance": "PowerTransformer",
    "cv_EMM_Proliferating_Cell_RNA_Abundance": "PowerTransformer",
    "tAIg": "StandardScaler",
    "mRNA_half_life_minutes": "PowerTransformer",
    "mRNA_synthesis_rate_per_minute": "PowerTransformer",
    "Mass (kDa)": "PowerTransformer",
    "pI": "StandardScaler",
    "Charge": "StandardScaler",
    "Residues": "PowerTransformer",
    "CAI": "StandardScaler",
    "aromaticity": "StandardScaler",
    "aliphatic_index": "StandardScaler",
    "gravy": "StandardScaler",
    "flexibility": "StandardScaler",
    "instability_index": "StandardScaler",
    "aa_percent_Ala": "StandardScaler",
    "aa_percent_Cys": "StandardScaler",
    "aa_percent_Asp": "StandardScaler",
    "aa_percent_Glu": "StandardScaler",
    "aa_percent_Phe": "StandardScaler",
    "aa_percent_Gly": "StandardScaler",
    "aa_percent_His": "StandardScaler",
    "aa_percent_Ile": "StandardScaler",
    "aa_percent_Lys": "StandardScaler",
    "aa_percent_Leu": "StandardScaler",
    "aa_percent_Met": "StandardScaler",
    "aa_percent_Asn": "StandardScaler",
    "aa_percent_Pro": "StandardScaler",
    "aa_percent_Gln": "StandardScaler",
    "aa_percent_Arg": "StandardScaler",
    "aa_percent_Ser": "StandardScaler",
    "aa_percent_Thr": "StandardScaler",
    "aa_percent_Val": "StandardScaler",
    "aa_percent_Trp": "StandardScaler",
    "aa_percent_Tyr": "StandardScaler",
    "mean_pLDDT": "StandardScaler",
    "std_pLDDT": "StandardScaler",
    "cv_pLDDT": "StandardScaler",
    "PFAM_domain_count": "PowerTransformer",
    "japonicus_ortholog_count": "PowerTransformer",
    "cerevisiae_ortholog_count": "PowerTransformer",
    "human_ortholog_count": "PowerTransformer",
    "paralog_count": "PowerTransformer",
    "evolutionary_rate": "StandardScaler",
    "mean.phylop": "StandardScaler",
    "diversity.S": "PowerTransformer",
    "diversity.Pi": "PowerTransformer",
    "diversity.Theta": "StandardScaler",
    "diversity.Tajima_D": "StandardScaler",
    "GO_term_richness": "PowerTransformer",
    "PPI_degree": "PowerTransformer",
    "GI_degree": "PowerTransformer",
}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PrepConfig:
    """Inputs, output dir, and split/threshold parameters for feature/target prep."""
    feature_matrix: Path
    final_clusters: Path
    output_dir: Path
    dr_split_threshold: float = DR_SPLIT_THRESHOLD
    wt_cluster: int = WT_CLUSTER
    transformations: dict = field(default_factory=lambda: dict(SELECTED_FEATURES_AND_TRANSFORMATIONS))

    def validate(self) -> None:
        """Raise ValueError if a required input is missing."""
        for path in [self.feature_matrix, self.final_clusters]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def merge_features_targets(feature_matrix: Path, final_clusters: Path) -> pd.DataFrame:
    """Left-join the feature matrix with A/DR/DL/DIT_HAP_cluster targets from final_clusters.tsv."""
    gene_features = pd.read_csv(feature_matrix, sep="\t")
    metrics = pd.read_csv(final_clusters, sep="\t").rename(columns={"Systematic ID": "Systematic_ID"})

    # DIT_HAP_cluster comes from the curated revised_cluster (NOT the raw cluster).
    target_values = metrics[["Systematic_ID", "A", "DR", "DL", "revised_cluster"]].rename(
        columns={"revised_cluster": "DIT_HAP_cluster"}
    )

    merged = (
        pd.merge(gene_features, target_values, left_on="gene_systematic_id", right_on="Systematic_ID", how="left")
        .drop(columns=["Systematic_ID"])
        .rename(columns={"gene_systematic_id": "Systematic_ID"})
    )
    logger.info(f"Merged feature+target matrix: {merged.shape}")
    return merged


@logger.catch(reraise=True)
def transform_data(data: pd.DataFrame, transformations: dict, target_features: list[str]) -> pd.DataFrame:
    """Apply the per-feature transformation map column-by-column, then append untransformed targets."""
    transformed_data = pd.DataFrame()
    for feature, transformation in transformations.items():
        if feature not in data.columns and transformation != "set_index":
            continue
        if transformation == "set_index":
            transformed_data[feature] = data.index.tolist()
            transformed_data.set_index(feature, inplace=True)
        elif transformation == "OneHotEncoder":
            values = pd.get_dummies(data[feature], prefix=feature, dummy_na=True).astype(int)
            transformed_data = pd.concat([transformed_data, values], axis=1)
        elif transformation == "binary_encoding":
            transformed_data[feature] = data[feature].map({"+": 1, "-": 0})
        elif transformation == "StandardScaler":
            transformed_data[feature] = StandardScaler().fit_transform(data[feature].values.reshape(-1, 1))
        elif transformation == "PowerTransformer":
            transformed_data[feature] = PowerTransformer().fit_transform(data[feature].values.reshape(-1, 1))
        else:
            transformed_data[feature] = data[feature]

    transformed_data[target_features] = data[target_features]
    return transformed_data


@logger.catch
def run_preparation(config: PrepConfig) -> None:
    """Orchestrate merge -> per-split dropna -> transform -> write feature/target CSVs."""
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)

    merged = merge_features_targets(config.feature_matrix, config.final_clusters)
    merged.to_csv(config.output_dir / "all_features_with_target_values.csv", index=False, float_format="%.5f")

    selected_features = list(config.transformations.keys())
    available = [f for f in selected_features if f in merged.columns]
    missing = [f for f in selected_features if f not in merged.columns and f != "Systematic_ID"]
    if missing:
        logger.warning(f"{len(missing)} configured features absent from matrix: {missing}")

    # Missing-value report over the selected columns.
    merged[available + TARGET_FEATURES].isna().sum().sort_values(ascending=False).to_csv(
        config.output_dir / "missing_value_analysis.csv"
    )

    # Build the three primary splits (no imputation — dropna defines each set).
    working = {
        "all": merged[available + TARGET_FEATURES].copy().dropna(axis=0, how="any").set_index("Systematic_ID"),
        "DR_gt_p35": merged.query(f"DR > {config.dr_split_threshold}")[available + TARGET_FEATURES].copy().dropna(axis=0, how="any").set_index("Systematic_ID"),
        "DR_le_p35": merged.query(f"DR <= {config.dr_split_threshold}")[available + TARGET_FEATURES].copy().dropna(axis=0, how="any").set_index("Systematic_ID"),
    }

    transformed = {des: transform_data(data, config.transformations, TARGET_FEATURES) for des, data in working.items()}

    for des, tdf in transformed.items():
        feature_cols = [c for c in tdf.columns if c not in TARGET_FEATURES]
        tdf.to_csv(config.output_dir / f"{des}_transformed_features_and_targets.csv", index=True, float_format="%.3f")
        tdf[feature_cols].to_csv(config.output_dir / f"{des}_transformed_features.csv", index=True, float_format="%.3f")
        tdf[TARGET_FEATURES].to_csv(config.output_dir / f"{des}_transformed_targets.csv", index=True, float_format="%.3f")
        logger.info(f"  split {des}: {tdf.shape}")

    # nonWT split derived from the 'all' transformed frame (DIT_HAP_cluster != WT).
    non_wt = transformed["all"].query(f"DIT_HAP_cluster != {config.wt_cluster}").copy()
    feature_cols = [c for c in non_wt.columns if c not in TARGET_FEATURES]
    non_wt.to_csv(config.output_dir / "nonWT_transformed_features_and_targets.csv", index=True, float_format="%.3f")
    non_wt[feature_cols].to_csv(config.output_dir / "nonWT_transformed_features.csv", index=True, float_format="%.3f")
    non_wt[TARGET_FEATURES].to_csv(config.output_dir / "nonWT_transformed_targets.csv", index=True, float_format="%.3f")
    logger.info(f"  split nonWT: {non_wt.shape}")

    logger.success(f"Feature/target preparation complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare ML feature/target tables from feature matrix + clusters")
    parser.add_argument("--feature-matrix", type=Path, required=True, help="Per-gene feature matrix tsv")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output dir for features_targets tables")
    parser.add_argument("--dr-split-threshold", type=float, default=DR_SPLIT_THRESHOLD, help="DR split threshold (default 0.35)")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT cluster id for nonWT split (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = PrepConfig(
            feature_matrix=args.feature_matrix,
            final_clusters=args.final_clusters,
            output_dir=args.output_dir,
            dr_split_threshold=args.dr_split_threshold,
            wt_cluster=args.wt_cluster,
        )
        run_preparation(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
