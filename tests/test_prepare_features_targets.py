"""Tests for workflow/scripts/ml/prepare_features_targets.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.ml.prepare_features_targets import (
    PrepConfig,
    TARGET_FEATURES,
    SELECTED_FEATURES_AND_TRANSFORMATIONS,
    merge_features_targets,
    transform_data,
)


def test_transform_dict_has_expected_transformation_types():
    """The transformation map uses only the supported transform names."""
    allowed = {"set_index", "OneHotEncoder", "binary_encoding", "StandardScaler", "PowerTransformer", "np.log1p", "No processing"}
    assert set(SELECTED_FEATURES_AND_TRANSFORMATIONS.values()) <= allowed
    assert SELECTED_FEATURES_AND_TRANSFORMATIONS["Systematic_ID"] == "set_index"
    assert SELECTED_FEATURES_AND_TRANSFORMATIONS["Chromosome"] == "OneHotEncoder"
    assert SELECTED_FEATURES_AND_TRANSFORMATIONS["Strand"] == "binary_encoding"


def test_transform_data_applies_each_transformation():
    """transform_data one-hots Chromosome, binary-encodes Strand, scales numeric, keeps targets."""
    df = pd.DataFrame(
        {
            "Chromosome": ["I", "II", "I", "III"],
            "Strand": ["+", "-", "+", "-"],
            "GC_content_of_gene": [0.4, 0.5, 0.45, 0.6],
            "A": [1.0, 2.0, 3.0, 4.0],
            "DR": [0.1, 0.2, 0.3, 0.4],
            "DL": [5.0, 6.0, 7.0, 8.0],
            "DIT_HAP_cluster": [1, 2, 3, 9],
        },
        index=["g1", "g2", "g3", "g4"],
    )
    transformations = {
        "Chromosome": "OneHotEncoder",
        "Strand": "binary_encoding",
        "GC_content_of_gene": "StandardScaler",
    }
    out = transform_data(df, transformations, TARGET_FEATURES)

    # one-hot with dummy_na=True yields a Chromosome_nan column even with no NaN.
    assert {"Chromosome_I", "Chromosome_II", "Chromosome_III", "Chromosome_nan"} <= set(out.columns)
    assert sorted(out["Strand"].unique()) == [0, 1]
    assert abs(out["GC_content_of_gene"].mean()) < 1e-9  # standardized
    # targets appended untransformed.
    assert list(out["A"]) == [1.0, 2.0, 3.0, 4.0]


def test_transform_data_skips_absent_features():
    """A configured feature missing from the data is skipped (not a KeyError)."""
    df = pd.DataFrame(
        {"GC3": [0.1, 0.2], "A": [1.0, 2.0], "DR": [0.1, 0.2], "DL": [1.0, 2.0], "DIT_HAP_cluster": [1, 2]},
        index=["g1", "g2"],
    )
    # set_index first (mirrors the real dict, which always starts with Systematic_ID).
    transformations = {"Systematic_ID": "set_index", "GC3": "StandardScaler", "NonexistentFeature": "PowerTransformer"}
    out = transform_data(df, transformations, TARGET_FEATURES)
    assert "GC3" in out.columns
    assert "NonexistentFeature" not in out.columns


def test_merge_uses_revised_cluster_as_dit_hap_cluster(tmp_path):
    """merge_features_targets maps revised_cluster -> DIT_HAP_cluster via a left join."""
    feat = pd.DataFrame({"gene_systematic_id": ["SPAC1", "SPAC2", "SPAC3"], "GC3": [0.1, 0.2, 0.3]})
    feat_path = tmp_path / "features.tsv"
    feat.to_csv(feat_path, sep="\t", index=False)

    clusters = pd.DataFrame(
        {"Systematic ID": ["SPAC1", "SPAC2"], "A": [1.0, 2.0], "DR": [0.5, 0.1], "DL": [3.0, 4.0], "revised_cluster": [1, 9]}
    )
    clusters_path = tmp_path / "final_clusters.tsv"
    clusters.to_csv(clusters_path, sep="\t", index=False)

    merged = merge_features_targets(feat_path, clusters_path)
    assert "DIT_HAP_cluster" in merged.columns
    assert "Systematic_ID" in merged.columns
    # SPAC3 has no cluster -> NaN target (left join keeps all feature rows).
    assert len(merged) == 3
    assert merged.set_index("Systematic_ID").loc["SPAC3", "DIT_HAP_cluster"] != merged.set_index("Systematic_ID").loc["SPAC3", "DIT_HAP_cluster"]  # NaN


def test_config_validate_rejects_missing_input(tmp_path):
    """validate() raises ValueError when a required input is absent."""
    cfg = PrepConfig(
        feature_matrix=tmp_path / "nope.tsv",
        final_clusters=tmp_path / "nope2.tsv",
        output_dir=tmp_path / "out",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()
