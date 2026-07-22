#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ML AutoML Analysis (mljar-supervised)
=====================================

Trains an mljar-supervised AutoML regressor to predict a growth-fitness target
(DR or DL) from ~70 gene/protein features, in one mljar mode (Explain or
Perform). Deterministic port of
DIT_HAP_pipeline/workflow/notebooks/machine_learning_analysis.ipynb.

The merge + DR-filter is done once by prepare_ml_data.py (the shared spine) and
read here from a parquet — the four target x mode jobs share one prepared table
instead of each re-merging. This script still does its own train-only
PowerTransform (leakage-free, different from the Task 6 transformed tables). One
invocation = one target x mode; the Snakemake rule fans out over the combinations.

Determinism: total_time_limit is set generously and random_state is passed
explicitly so the full algorithm list always completes and results are
reproducible (mljar's default 3600s budget can silently skip algorithms).

Input
-----
- modeling_data.parquet: merged, DR-filtered modeling table (from prepare_ml_data.py)

Output
------
- mljar results tree (leaderboard.csv, per-model dirs, ...)
- metrics.tsv (R2/RMSE/MAE/Pearson on original target scale)
- features_importance.csv (aggregated from per-model permutation importance)
- prediction_and_residuals.pdf
- feature_scaler.joblib, target_scaler.joblib

Usage
-----
    python train_automl.py \\
        --modeling-data results/ml/models/{dataset}/{version}/_work/modeling_data.parquet \\
        --target DR --mode Explain \\
        --output-dir results/ml/models/{dataset}/{version}/DR_Explain

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-16
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import joblib
import matplotlib
from loguru import logger
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PowerTransformer

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
TEST_SIZE = 0.2
RANDOM_STATE = 42
# Generous budget so every algorithm in the mode's list finishes (else mljar
# silently drops algorithms on slow runs -> non-reproducible). 4h ceiling.
TOTAL_TIME_LIMIT = 14400

# Feature list from the notebook. Some names differ from our matrix (t1/2 (min)
# -> protein_half_life_minutes); absent features are dropped at runtime with a warning.
FEATURE_COLUMNS = [
    "Abs_distance_from_telomere", "Relative_distance_from_telomere", "Abs_distance_from_centromere", "Relative_distance_from_centromere",
    "Gene_length", "GC_content_of_gene", "CDS_number", "GC_content_of_CDS", "Fraction_of_CDS", "GC3", "ENC",
    "Peptide_length", "Primary_peptide_length",
    "mean_EMM_Nitrogen_Starved_Cell_RNA_Abundance", "mean_EMM_Proliferating_Cell_RNA_Abundance",
    "std_EMM_Nitrogen_Starved_Cell_RNA_Abundance", "std_EMM_Proliferating_Cell_RNA_Abundance",
    "cv_EMM_Nitrogen_Starved_Cell_RNA_Abundance", "cv_EMM_Proliferating_Cell_RNA_Abundance",
    "tAIg", "mRNA_half_life_minutes", "mRNA_synthesis_rate_per_minute",
    "Mass (kDa)", "pI", "Charge", "Residues", "CAI", "aromaticity", "aliphatic_index", "gravy", "flexibility", "instability_index",
    "aa_percent_Ala", "aa_percent_Cys", "aa_percent_Asp", "aa_percent_Glu", "aa_percent_Phe", "aa_percent_Gly", "aa_percent_His",
    "aa_percent_Ile", "aa_percent_Lys", "aa_percent_Leu", "aa_percent_Met", "aa_percent_Asn", "aa_percent_Pro", "aa_percent_Gln",
    "aa_percent_Arg", "aa_percent_Ser", "aa_percent_Thr", "aa_percent_Val", "aa_percent_Trp", "aa_percent_Tyr",
    "copies_per_cell_EMM_Proliferating_Cell", "copies_per_cell_EMMN_Quiescent_Cell", "protein_half_life_minutes",
    "mean_pLDDT", "std_pLDDT", "cv_pLDDT", "PFAM_domain_count",
    "japonicus_ortholog_count", "cerevisiae_ortholog_count", "human_ortholog_count", "paralog_count",
    "evolutionary_rate", "mean.phylop", "diversity.S", "diversity.Pi", "diversity.Theta", "diversity.Tajima_D",
    "GO_term_richness", "PPI_degree", "GI_degree",
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class AutoMLConfig:
    """Inputs, output dir, target/mode, and training parameters."""
    modeling_data: Path
    target: str
    mode: str
    output_dir: Path
    test_size: float = TEST_SIZE
    random_state: int = RANDOM_STATE
    total_time_limit: int = TOTAL_TIME_LIMIT
    feature_columns: list = field(default_factory=lambda: list(FEATURE_COLUMNS))

    def validate(self) -> None:
        """Raise ValueError on missing input or invalid target/mode."""
        if not self.modeling_data.exists():
            raise ValueError(f"Required input not found: {self.modeling_data}")
        if self.target not in {"DR", "DL"}:
            raise ValueError(f"target must be 'DR' or 'DL', got {self.target!r}")
        if self.mode not in {"Explain", "Perform"}:
            raise ValueError(f"mode must be 'Explain' or 'Perform', got {self.mode!r}")


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


@logger.catch(reraise=True)
def load_modeling_data(config: AutoMLConfig) -> pd.DataFrame:
    """Read the shared modeling table (already merged + DR-filtered by prepare_ml_data.py)."""
    data = read_parquet(config.modeling_data)
    logger.info(f"Modeling data: {data.shape}")
    return data


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def train(config: AutoMLConfig, data: pd.DataFrame):
    """Split, PowerTransform (train-only), fit mljar AutoML, return model + test data + predictions + scalers."""
    from supervised.automl import AutoML

    available = [c for c in config.feature_columns if c in data.columns]
    missing = [c for c in config.feature_columns if c not in data.columns]
    if missing:
        logger.warning(f"{len(missing)} configured features absent from matrix: {missing}")

    X = data[available].copy()
    y = data[config.target].copy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=config.test_size, random_state=config.random_state)

    # Train-only PowerTransform (no leakage), persisted so predictions can be inverse-transformed later.
    feats_scaler = PowerTransformer(method="yeo-johnson")
    X_train = pd.DataFrame(feats_scaler.fit_transform(X_train), index=X_train.index, columns=X_train.columns)
    X_test = pd.DataFrame(feats_scaler.transform(X_test), index=X_test.index, columns=X_test.columns)

    target_scaler = PowerTransformer(method="yeo-johnson")
    y_train_t = target_scaler.fit_transform(y_train.values.reshape(-1, 1)).ravel()

    # Clean the results dir first — mljar errors/resumes on a non-empty results_path.
    if config.output_dir.exists():
        shutil.rmtree(config.output_dir)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    automl = AutoML(
        mode=config.mode,
        ml_task="regression",
        mix_encoding=True,
        results_path=str(config.output_dir),
        explain_level=2,
        random_state=config.random_state,
        total_time_limit=config.total_time_limit,
    )
    automl.fit(X_train, y_train_t)

    y_pred_t = automl.predict(X_test)
    y_pred = target_scaler.inverse_transform(np.asarray(y_pred_t).reshape(-1, 1)).ravel()

    return automl, X_test, y_test, y_pred, feats_scaler, target_scaler


@logger.catch(reraise=True)
def evaluate(config: AutoMLConfig, y_test: pd.Series, y_pred: np.ndarray) -> dict:
    """Compute R2/RMSE/MAE/Pearson on the ORIGINAL target scale and save a scatter+residual PDF."""
    r2 = r2_score(y_test, y_pred)
    rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))  # squared=False deprecated
    mae = mean_absolute_error(y_test, y_pred)
    r, p = pearsonr(y_test, y_pred)
    metrics = {"R2": r2, "RMSE": rmse, "MAE": mae, "Pearson_r": r, "Pearson_p": p}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(y_test, y_pred, s=20, alpha=0.5, color="#962955", edgecolor="none")
    lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
    axes[0].plot(lims, lims, "k--", linewidth=1, alpha=0.6)
    axes[0].set_xlabel(f"True {config.target}")
    axes[0].set_ylabel(f"Predicted {config.target}")
    axes[0].set_title(f"Predicted vs True ({config.mode})")
    axes[0].text(0.05, 0.92, f"R2={r2:.3f}\nr={r:.3f}\nRMSE={rmse:.3f}", transform=axes[0].transAxes,
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    residuals = y_test.values - y_pred
    axes[1].scatter(y_pred, residuals, s=20, alpha=0.5, color="#6479cc", edgecolor="none")
    axes[1].axhline(0, color="k", linestyle="--", linewidth=1, alpha=0.6)
    axes[1].set_xlabel(f"Predicted {config.target}")
    axes[1].set_ylabel("Residual (True - Predicted)")
    axes[1].set_title(f"Residuals ({config.mode})")

    plt.suptitle(f"Target: {config.target} | Mode: {config.mode}")
    plt.tight_layout()
    fig.savefig(config.output_dir / "prediction_and_residuals.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    return metrics


@logger.catch(reraise=True)
def aggregate_feature_importance(output_dir: Path) -> pd.DataFrame:
    """Aggregate per-model/per-fold permutation importance CSVs into a mean-importance table.

    The notebook read a features_importance.csv that mljar never writes (dead code);
    we build it here from the learner_fold_*_importance.csv files mljar does write.
    """
    importance_files = list(output_dir.glob("*/learner_fold_*_importance.csv"))
    if not importance_files:
        logger.warning("No per-model importance files found to aggregate")
        return pd.DataFrame()

    frames = []
    for f in importance_files:
        df = pd.read_csv(f, index_col=0)
        # mljar importance files have a single value column; normalize to 'importance'.
        value_col = df.columns[0]
        frames.append(df[[value_col]].rename(columns={value_col: "importance"}))

    combined = pd.concat(frames)
    agg = combined.groupby(combined.index)["importance"].mean().sort_values(ascending=False)
    result = agg.rename_axis("feature").reset_index()
    result.to_csv(output_dir / "features_importance.csv", index=False)
    return result


@logger.catch(reraise=True)
def run_automl_analysis(config: AutoMLConfig) -> None:
    """Orchestrate load -> train -> evaluate -> persist metrics/scalers/importance."""
    config.validate()
    data = load_modeling_data(config)

    automl, X_test, y_test, y_pred, feats_scaler, target_scaler = train(config, data)

    metrics = evaluate(config, y_test, y_pred)
    pd.DataFrame([metrics]).to_csv(config.output_dir / "metrics.tsv", sep="\t", index=False)
    logger.info("metrics: " + ", ".join(f"{k}={v:.4f}" for k, v in metrics.items()))

    joblib.dump(feats_scaler, config.output_dir / "feature_scaler.joblib")
    joblib.dump(target_scaler, config.output_dir / "target_scaler.joblib")

    aggregate_feature_importance(config.output_dir)

    logger.success(f"AutoML analysis complete ({config.target} {config.mode}) -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Train mljar AutoML regressor for one target x mode")
    parser.add_argument("--modeling-data", type=Path, required=True, help="Shared modeling_data parquet (from prepare_ml_data.py)")
    parser.add_argument("--target", required=True, choices=["DR", "DL"], help="Regression target")
    parser.add_argument("--mode", required=True, choices=["Explain", "Perform"], help="mljar AutoML mode")
    parser.add_argument("--output-dir", type=Path, required=True, help="mljar results_path / output dir")
    parser.add_argument("--test-size", type=float, default=TEST_SIZE, help="Test split fraction (default 0.2)")
    parser.add_argument("--random-state", type=int, default=RANDOM_STATE, help="Random seed (default 42)")
    parser.add_argument("--total-time-limit", type=int, default=TOTAL_TIME_LIMIT, help="mljar total_time_limit seconds")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run AutoML analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = AutoMLConfig(
            modeling_data=args.modeling_data,
            target=args.target,
            mode=args.mode,
            output_dir=args.output_dir,
            test_size=args.test_size,
            random_state=args.random_state,
            total_time_limit=args.total_time_limit,
        )
        run_automl_analysis(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
