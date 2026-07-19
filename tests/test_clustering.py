"""Tests for the split gene-level clustering pipeline (shared module + driver configs)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.clustering.candidates import (
    BEST_METHOD,
    DR_CAP,
    DL_DIVISOR,
    METHODS,
    cluster_one_method,
    scale_features,
    score_labels,
)
from workflow.scripts.clustering.prepare_clustering_data import PrepareConfig
from workflow.scripts.clustering.cluster_one_method import MethodConfig
from workflow.scripts.clustering.select_candidate_clusters import SelectConfig, combine_metrics


def test_best_method_is_pinned_to_kmeans():
    """The fragile set()[0] selection is replaced by an explicit pin (quirk #3)."""
    assert BEST_METHOD == "kmeans"
    assert METHODS[0] == "kmeans"


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


def test_scale_features_respects_custom_dr_cap_and_dl_divisor():
    """Non-default dr_cap/dl_divisor (as would come from config/analysis.yaml) are honored."""
    df = pd.DataFrame(
        {"DR": [0.5, 1.0, 2.0], "DL": [10.0, 20.0, 5.0]},
        index=["g1", "g2", "g3"],
    )
    scaled = scale_features(df, ["DR", "DL"], dr_cap=1.0, dl_divisor=5)
    assert list(scaled["DR"]) == [0.5, 1.0, 1.0]
    assert list(scaled["DL"]) == [2.0, 4.0, 1.0]


def test_scale_features_drops_nan_rows():
    """dropna defines the clustered set — genes with NaN DR/DL are excluded."""
    df = pd.DataFrame(
        {"DR": [0.5, np.nan, 0.8], "DL": [10.0, 20.0, np.nan]},
        index=["g1", "g2", "g3"],
    )
    scaled = scale_features(df, ["DR", "DL"])
    assert list(scaled.index) == ["g1"]


def test_cluster_one_method_all_four_return_0based_labels():
    """Each of the four methods returns 0-based labels for the requested cluster count."""
    rng = np.random.default_rng(0)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    for method in METHODS:
        labels = cluster_one_method(method, data, n_clusters=4, random_state=42)
        assert labels.min() == 0, method
        assert len(np.unique(labels)) == 4, method


def test_cluster_one_method_rejects_unknown_method():
    """cluster_one_method raises ValueError on an unrecognized method name."""
    data = pd.DataFrame(np.random.default_rng(0).random((20, 2)), columns=["DR", "DL"])
    with pytest.raises(ValueError, match="Unknown clustering method"):
        cluster_one_method("dbscan", data, n_clusters=3, random_state=42)


def test_score_labels_reports_three_scores_and_count():
    """score_labels returns the three sklearn scores plus the cluster count."""
    rng = np.random.default_rng(1)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    labels = cluster_one_method("kmeans", data, n_clusters=3, random_state=42)
    scores = score_labels(data, labels)
    assert set(scores) == {"silhouette_score", "calinski_harabasz_score", "davies_bouldin_score", "n_clusters"}
    assert scores["n_clusters"] == 3


def test_kmeans_is_deterministic_given_seed():
    """Same seed + same data -> identical kmeans labels (reproducibility guarantee)."""
    rng = np.random.default_rng(2)
    data = pd.DataFrame(rng.random((80, 2)), columns=["DR", "DL"])
    r1 = cluster_one_method("kmeans", data, n_clusters=5, random_state=42)
    r2 = cluster_one_method("kmeans", data, n_clusters=5, random_state=42)
    assert np.array_equal(r1, r2)


def test_combine_metrics_stacks_ksweep_then_rounded_methods():
    """combine_metrics concatenates the unrounded k-sweep and 3dp-rounded method rows."""
    ksweep = pd.DataFrame({"k": [2], "inertia": [1.23456], "silhouette": [0.5]})
    methods = pd.DataFrame(
        {"method": ["kmeans"], "silhouette_score": [0.4531111], "n_clusters": [64]}
    )
    combined = combine_metrics(ksweep, methods)
    assert list(combined["table"]) == ["k_sweep", "method_comparison"]
    # k-sweep row keeps full precision; method row is rounded to 3 dp.
    assert combined.loc[0, "inertia"] == 1.23456
    assert combined.loc[1, "silhouette_score"] == 0.453


def test_prepare_config_rejects_missing_input(tmp_path):
    """PrepareConfig.validate raises ValueError when a required input file is absent."""
    cfg = PrepareConfig(
        fitting_results=tmp_path / "nope.tsv",
        essentiality_verification_csv=tmp_path / "also_nope.csv",
        output_annotated=tmp_path / "w" / "ann.pkl",
        output_scaled=tmp_path / "w" / "scaled.pkl",
        output_ksweep=tmp_path / "w" / "ksweep.pkl",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_method_config_rejects_unknown_method(tmp_path):
    """MethodConfig.validate raises ValueError on an unrecognized method name."""
    scaled = tmp_path / "scaled.pkl"
    scaled.write_bytes(b"")
    cfg = MethodConfig(
        method="spectral",
        scaled_data=scaled,
        output_labels=tmp_path / "w" / "l.pkl",
        output_metrics=tmp_path / "w" / "m.pkl",
    )
    with pytest.raises(ValueError, match="Unknown clustering method"):
        cfg.validate()


def test_select_config_rejects_missing_input(tmp_path):
    """SelectConfig.validate raises ValueError when a required input is absent."""
    cfg = SelectConfig(
        annotated_data=tmp_path / "ann.pkl",
        ksweep=tmp_path / "ksweep.pkl",
        best_labels=tmp_path / "labels.pkl",
        method_metrics=[tmp_path / "m.pkl"],
        output_clusters=tmp_path / "out" / "candidate_clusters.tsv",
        output_metrics=tmp_path / "out" / "clustering_metrics.tsv",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()
