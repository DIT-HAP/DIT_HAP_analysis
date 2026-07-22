"""Tests for the split ML pipeline: src.ml.data spine + train_automl driver (no mljar training)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.io import write_parquet, read_parquet

from workflow.src.ml.data import load_modeling_data
from workflow.scripts.ml.prepare_ml_data import PrepareConfig
from workflow.scripts.ml.train_automl import (
    AutoMLConfig,
    FEATURE_COLUMNS,
    aggregate_feature_importance,
    evaluate,
    load_modeling_data as load_modeling_data_pickle,
)


def _cfg(tmp_path, target="DR", mode="Explain", **kw):
    return AutoMLConfig(
        modeling_data=tmp_path / "modeling_data.parquet",
        target=target,
        mode=mode,
        output_dir=tmp_path / "out",
        **kw,
    )


def test_config_validate_rejects_bad_target(tmp_path):
    """validate() rejects a target outside {DR, DL}."""
    (tmp_path / "modeling_data.parquet").write_bytes(b"")
    with pytest.raises(ValueError, match="target must be"):
        _cfg(tmp_path, target="A").validate()


def test_config_validate_rejects_bad_mode(tmp_path):
    """validate() rejects a mode outside {Explain, Perform}."""
    (tmp_path / "modeling_data.parquet").write_bytes(b"")
    with pytest.raises(ValueError, match="mode must be"):
        _cfg(tmp_path, mode="Compete").validate()


def test_config_validate_rejects_missing_input(tmp_path):
    """validate() raises when the modeling-data pickle is absent."""
    with pytest.raises(ValueError, match="Required input not found"):
        _cfg(tmp_path).validate()


def test_feature_list_uses_renamed_half_life_column():
    """Our matrix renamed t1/2 (min) -> protein_half_life_minutes; the feature list must match."""
    assert "protein_half_life_minutes" in FEATURE_COLUMNS
    assert "t1/2 (min)" not in FEATURE_COLUMNS


def test_prepare_config_rejects_missing_input(tmp_path):
    """prepare_ml_data PrepareConfig.validate raises when a required input is absent."""
    cfg = PrepareConfig(
        feature_matrix=tmp_path / "f.tsv",
        final_clusters=tmp_path / "c.tsv",
        output=tmp_path / "work" / "modeling_data.parquet",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_src_load_modeling_data_filters_dr_and_maps_cluster(tmp_path):
    """src.ml.data.load_modeling_data left-joins targets, filters DR > 0.3, maps cluster."""
    feat = pd.DataFrame({"gene_systematic_id": ["SPAC1", "SPAC2", "SPAC3"], "GC3": [0.1, 0.2, 0.3]})
    feat.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    clusters = pd.DataFrame(
        {"Systematic ID": ["SPAC1", "SPAC2", "SPAC3"], "A": [1.0, 2.0, 3.0], "DR": [0.5, 0.1, 0.9], "DL": [3.0, 4.0, 5.0], "cluster": [1, 9, 2]}
    )
    clusters.to_csv(tmp_path / "c.tsv", sep="\t", index=False)

    data = load_modeling_data(tmp_path / "f.tsv", tmp_path / "c.tsv", dr_filter=0.3)
    # DR > 0.3 keeps SPAC1 (0.5) and SPAC3 (0.9), drops SPAC2 (0.1).
    assert set(data["Systematic_ID"]) == {"SPAC1", "SPAC3"}
    assert "DIT_HAP_cluster" in data.columns


def test_train_load_modeling_data_reads_pickle(tmp_path):
    """train_automl.load_modeling_data is now a thin read of the prepared pickle."""
    df = pd.DataFrame({"Systematic_ID": ["SPAC1"], "DR": [0.5], "GC3": [0.1]})
    write_parquet(df, tmp_path / "modeling_data.parquet")
    loaded = load_modeling_data_pickle(_cfg(tmp_path))
    pd.testing.assert_frame_equal(loaded, df)


def test_evaluate_computes_metrics_and_writes_pdf(tmp_path):
    """evaluate returns R2/RMSE/MAE/Pearson and writes the scatter+residual PDF."""
    cfg = _cfg(tmp_path)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    y_test = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y_pred = np.array([1.1, 1.9, 3.2, 3.8, 5.1])
    metrics = evaluate(cfg, y_test, y_pred)
    assert set(metrics) == {"R2", "RMSE", "MAE", "Pearson_r", "Pearson_p"}
    assert metrics["R2"] > 0.9  # near-perfect prediction
    assert (cfg.output_dir / "prediction_and_residuals.pdf").exists()


def test_aggregate_feature_importance_means_across_files(tmp_path):
    """aggregate_feature_importance averages per-model importance CSVs into one table."""
    out = tmp_path / "out"
    (out / "3_Linear").mkdir(parents=True)
    (out / "6_RandomForest").mkdir(parents=True)
    pd.DataFrame({"feature": ["f1", "f2"], "importance": [0.8, 0.2]}).set_index("feature").to_csv(out / "3_Linear" / "learner_fold_0_importance.csv")
    pd.DataFrame({"feature": ["f1", "f2"], "importance": [0.6, 0.4]}).set_index("feature").to_csv(out / "6_RandomForest" / "learner_fold_0_importance.csv")

    result = aggregate_feature_importance(out)
    assert list(result.columns) == ["feature", "importance"]
    f1 = result.set_index("feature").loc["f1", "importance"]
    assert abs(f1 - 0.7) < 1e-9  # mean of 0.8 and 0.6
    assert (out / "features_importance.csv").exists()
