import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd
import pytest


def _long(rows):
    return pd.DataFrame(rows, columns=["source", "group_id", "group_name",
                                       "Systematic ID", "Name", "n_group_genes"])


def test_shared_subunit_fraction_counts_cross_group_members():
    from compute_coherence import shared_subunit_fraction
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
        ("go_cc", "GO:2", "two", "gA", "gA", 1),
    ])
    frac = shared_subunit_fraction(long)
    assert frac["GO:1"] == 0.5   # gA shared, gB not -> 1/2
    assert frac["GO:2"] == 1.0   # gA shared -> 1/1


def test_shared_subunit_fraction_is_per_source():
    # a gene shared across DIFFERENT sources should NOT count as shared
    from compute_coherence import shared_subunit_fraction
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 1),
        ("go_bp", "GO:9", "nine", "gA", "gA", 1),
    ])
    frac = shared_subunit_fraction(long)
    assert frac["GO:1"] == 0.0  # gA only in one go_cc group
    assert frac["GO:9"] == 0.0


def test_member_feature_cv_computes_per_group_cv():
    from compute_coherence import member_feature_cv
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
    ])
    features = pd.DataFrame({"Systematic ID": ["gA", "gB"], "abundance": [10.0, 30.0]})
    cv = member_feature_cv(long, features, "abundance")
    assert cv["GO:1"] == pytest.approx(0.70711, rel=1e-3)  # sample std/mean of [10,30] = 0.707


def test_member_feature_cv_accepts_gene_systematic_id_column():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 2),
                  ("go_cc", "GO:1", "one", "gB", "gB", 2)])
    features = pd.DataFrame({"gene_systematic_id": ["gA", "gB"], "evolutionary_rate": [1.0, 2.0]})
    cv = member_feature_cv(long, features, "evolutionary_rate")
    assert "GO:1" in cv


def test_member_feature_cv_missing_column_returns_empty():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 1)])
    features = pd.DataFrame({"Systematic ID": ["gA"], "other": [1.0]})
    assert member_feature_cv(long, features, "abundance") == {}


def test_member_feature_cv_drops_nan_feature_members():
    from compute_coherence import member_feature_cv
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),  # gB has NaN feature -> dropped -> <2 left
    ])
    features = pd.DataFrame({"Systematic ID": ["gA", "gB"], "abundance": [10.0, np.nan]})
    assert "GO:1" not in member_feature_cv(long, features, "abundance")


def test_member_feature_cv_skips_zero_mean_group():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 2),
                  ("go_cc", "GO:1", "one", "gB", "gB", 2)])
    features = pd.DataFrame({"Systematic ID": ["gA", "gB"], "abundance": [0.0, 0.0]})
    assert "GO:1" not in member_feature_cv(long, features, "abundance")


def test_plot_coherence_axis_counts():
    import matplotlib.pyplot as plt
    from compute_coherence import plot_coherence
    table = pd.DataFrame({
        "source": ["go_cc", "go_cc"], "group_id": ["GO:1", "GO:2"],
        "group_name": ["a", "b"], "term_size": [3, 4],
        "centroid_x": [0.5, 0.6], "centroid_y": [0.1, 0.2],
        "z_score": [-1.0, 0.5], "p_value": [0.1, 0.5], "n_permutations": [50, 50],
    })
    long = _long([("go_cc", "GO:1", "a", "gA", "gA", 3), ("go_cc", "GO:2", "b", "gB", "gB", 4)])
    features = pd.DataFrame({"Systematic ID": ["gA", "gB"], "evolutionary_rate": [1.0, 2.0],
                             "copies_per_cell_EMM_Proliferating_Cell": [10.0, 30.0]})
    fig_wo = plot_coherence(table, long, features=None)
    fig_w = plot_coherence(table, long, features=features)
    assert len(fig_w.axes) > len(fig_wo.axes)  # B/D add panels
    plt.close(fig_wo)
    plt.close(fig_w)
    empty = plot_coherence(table.iloc[0:0], long, features=None)
    assert len(empty.axes) == 1
    plt.close(empty)
