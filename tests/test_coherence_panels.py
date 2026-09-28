import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd
import pytest


def _long(rows):
    return pd.DataFrame(rows, columns=["source", "group_id", "group_name",
                                       "Systematic ID", "Name", "n_annotated_members"])


def test_shared_subunit_fraction_counts_cross_group_members():
    from coherence.attribution import shared_subunit_fractions
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
        ("go_cc", "GO:2", "two", "gA", "gA", 1),
    ])
    frac = shared_subunit_fractions(long)
    assert frac[("go_cc", "GO:1")] == 0.5   # gA shared, gB not -> 1/2
    assert frac[("go_cc", "GO:2")] == 1.0   # gA shared -> 1/1


def test_shared_subunit_fraction_is_per_source():
    # a gene shared across DIFFERENT sources should NOT count as shared
    from coherence.attribution import shared_subunit_fractions
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 1),
        ("go_bp", "GO:9", "nine", "gA", "gA", 1),
    ])
    frac = shared_subunit_fractions(long)
    assert frac[("go_cc", "GO:1")] == 0.0  # gA only in one go_cc group
    assert frac[("go_bp", "GO:9")] == 0.0


def test_shared_subunit_fraction_keys_on_source_when_group_id_collides():
    """One group_id in two sources must keep two fractions, not overwrite one.

    173 group_ids really do appear in more than one source (GO:0032040 is both a
    go_cc term and a go_macrocomplex complex), so a group_id-keyed dict silently
    lets the last source win and reports the other one's fraction. The repr is
    `{('go_cc','GO:1'): 0.0, ('go_bp','GO:1'): 1.0}` — the same bare id, two
    different answers, and only the pair tells them apart.
    """
    from coherence.attribution import shared_subunit_fractions
    long = _long([
        ("go_cc", "GO:1", "cc view", "gA", "gA", 1),
        ("go_bp", "GO:1", "bp view", "gA", "gA", 2),
        ("go_bp", "GO:1", "bp view", "gB", "gB", 2),
        ("go_bp", "GO:2", "other", "gB", "gB", 1),
    ])
    frac = shared_subunit_fractions(long)
    assert frac[("go_cc", "GO:1")] == 0.0   # gA is alone within go_cc
    assert frac[("go_bp", "GO:1")] == 0.5   # gB is shared within go_bp


def test_member_feature_cv_computes_per_group_cv():
    from compute_coherence import member_feature_cv
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
    ])
    features = pd.DataFrame({"gene_systematic_id": ["gA", "gB"], "abundance": [10.0, 30.0]})
    cv = member_feature_cv(long, features, "abundance")
    assert cv["GO:1"] == pytest.approx(0.70711, rel=1e-3)  # sample std/mean of [10,30] = 0.707


def test_member_feature_cv_accepts_gene_systematic_id_column():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 2),
                  ("go_cc", "GO:1", "one", "gB", "gB", 2)])
    features = pd.DataFrame({"gene_systematic_id": ["gA", "gB"], "evolutionary_rate": [1.0, 2.0]})
    cv = member_feature_cv(long, features, "evolutionary_rate")
    assert "GO:1" in cv


def test_group_annotations_omits_columns_with_no_source_feature():
    """A features table carrying none of the candidate columns yields frac_shared_members only.

    The column triage lives in group_annotations (which tries each candidate in
    order and warns), not in member_feature_cv, so that is what this pins.
    """
    from compute_coherence import group_annotations
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 1)])
    features = pd.DataFrame({"gene_systematic_id": ["gA"], "other": [1.0]})
    annotations = group_annotations(long, features)
    assert list(annotations.columns) == ["frac_shared_members"]


def test_member_feature_cv_drops_nan_feature_members():
    from compute_coherence import member_feature_cv
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),  # gB has NaN feature -> dropped -> <2 left
    ])
    features = pd.DataFrame({"gene_systematic_id": ["gA", "gB"], "abundance": [10.0, np.nan]})
    assert "GO:1" not in member_feature_cv(long, features, "abundance")


def test_member_feature_cv_skips_zero_mean_group():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 2),
                  ("go_cc", "GO:1", "one", "gB", "gB", 2)])
    features = pd.DataFrame({"gene_systematic_id": ["gA", "gB"], "abundance": [0.0, 0.0]})
    assert "GO:1" not in member_feature_cv(long, features, "abundance")


def test_plot_coherence_panel_counts():
    """3 descriptive panels + one per available biology column + one FDR panel per encoding.

    cnsplots only exists in the cnsplots rule env, so this skips rather than fails
    where the figure cannot be drawn at all.
    """
    pytest.importorskip("cnsplots")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plot_coherence import biology_panels, plot_coherence

    base = {
        "group_id": ["GO:1", "GO:2", "GO:3"],
        "group_name": ["a", "b", "c"],
        "n_scored_members": [3, 4, 5], "median_pairwise_distance_z": [-1.0, 0.5, -0.5],
        "q_value": [0.01, 0.2, 0.03],
        "geom_median_DR": [0.5, 0.6, 0.4], "geom_median_DL": [0.1, 0.2, 0.3],
    }
    without = pd.DataFrame(base)
    with_biology = pd.DataFrame(base | {
        "frac_shared_members": [0.0, 0.5, 0.2],
        "abundance_cv": [0.1, 0.2, 0.3],
        "conservation_cv": [0.4, 0.5, 0.6],
    })

    assert biology_panels(without) == []
    assert len(biology_panels(with_biology)) == 3

    def drawn_panels(table):
        """Panel count, excluding the centroid panel's inset colourbar."""
        plot_coherence(table)
        count = len([ax for ax in plt.gcf().axes if not ax.get_label().startswith("<")])
        plt.close("all")
        return count

    # 3 descriptive + 3 biology + one FDR panel per x encoding.
    assert drawn_panels(with_biology) == 8
    assert drawn_panels(without) == 5  # 3 descriptive + both FDR panels
    assert drawn_panels(without.iloc[0:0]) == 1  # placeholder
