"""Tests for workflow/scripts/clustering/generate_candidate_clusters.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.clustering.generate_candidate_clusters import (
    ClusteringConfig,
    DR_CAP,
    DL_DIVISOR,
    BEST_METHOD,
    scale_features,
    perform_clustering_analysis,
    calculate_clustering_metrics,
)


def test_best_method_is_pinned_to_kmeans():
    """The fragile set()[0] selection is replaced by an explicit pin (quirk #3)."""
    assert BEST_METHOD == "kmeans"


def test_scale_features_caps_DR_and_divides_DL():
    """DR above 1.3 is clamped to 1.3; DL is divided by 10 (byte-faithful quirk #1)."""
    df = pd.DataFrame(
        {"DR": [0.5, 1.3, 2.0, 1.29], "DL": [10.0, 20.0, 5.0, 100.0]},
        index=["g1", "g2", "g3", "g4"],
    )
    scaled = scale_features(df, ["DR", "DL"])
    # DR: 0.5 stays, 1.3 stays (strict < so equal is unchanged), 2.0 -> 1.3, 1.29 stays.
    assert list(scaled["DR"]) == [0.5, 1.3, DR_CAP, 1.29]
    # DL divided by 10.
    assert list(scaled["DL"]) == [1.0, 2.0, 0.5, 10.0]


def test_scale_features_drops_nan_rows():
    """dropna defines the clustered set — genes with NaN DR/DL are excluded."""
    df = pd.DataFrame(
        {"DR": [0.5, np.nan, 0.8], "DL": [10.0, 20.0, np.nan]},
        index=["g1", "g2", "g3"],
    )
    scaled = scale_features(df, ["DR", "DL"])
    assert list(scaled.index) == ["g1"]


def test_perform_clustering_analysis_returns_four_methods_0based():
    """All four algorithms return 0-based labels for the requested cluster count."""
    rng = np.random.default_rng(0)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    results = perform_clustering_analysis(data, n_clusters=4, random_state=42)
    assert set(results) == {"kmeans", "hierarchical_agg", "hierarchical_div", "gmm"}
    for method, labels in results.items():
        assert labels.min() == 0, method
        assert len(np.unique(labels)) == 4, method


def test_calculate_clustering_metrics_scores_each_method():
    """Metrics table has one row per method with the three scores + n_clusters."""
    rng = np.random.default_rng(1)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    results = perform_clustering_analysis(data, n_clusters=3, random_state=42)
    metrics = calculate_clustering_metrics(data, results)
    assert set(metrics["method"]) == {"kmeans", "hierarchical_agg", "hierarchical_div", "gmm"}
    assert {"silhouette_score", "calinski_harabasz_score", "davies_bouldin_score", "n_clusters"} <= set(metrics.columns)


def test_kmeans_is_deterministic_given_seed():
    """Same seed + same data -> identical kmeans labels (reproducibility guarantee)."""
    rng = np.random.default_rng(2)
    data = pd.DataFrame(rng.random((80, 2)), columns=["DR", "DL"])
    r1 = perform_clustering_analysis(data, n_clusters=5, random_state=42)["kmeans"]
    r2 = perform_clustering_analysis(data, n_clusters=5, random_state=42)["kmeans"]
    assert np.array_equal(r1, r2)


def test_config_validate_rejects_missing_input(tmp_path):
    """validate() raises ValueError when a required input file is absent."""
    cfg = ClusteringConfig(
        fitting_results=tmp_path / "nope.tsv",
        essentiality_verification_csv=tmp_path / "also_nope.csv",
        output_clusters=tmp_path / "out.tsv",
        output_metrics=tmp_path / "metrics.tsv",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()
