"""Tests for the per-view fractions (compute_view_fractions.py + its figure).

The semantics pinned here are the point of the stage: breadth is counted PER VIEW
(a view's groups are the de-duplication clusters its terms belong to) and the mode
is each view's own, so the same gene can be moonlighting in one view and not in
another. The gene sets, by contrast, stay each cluster's own — two views differ
only in which groups they count.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd
import pytest


def _metrics():
    """The combined metrics table, sliced by view_terms: three terms, two sources."""
    return pd.DataFrame({
        "source": ["go_cc", "go_bp", "go_bp"],
        "group_id": ["GO:1", "GO:2", "GO:3"],
        "group_name": ["one", "two", "three"],
        "n_scored_members": [2, 2, 2],
        "scored_member_names": [["a", "b"], ["a", "c"], ["c", "d"]],
        "paralog_fraction": [0.5, 0.0, 1.0],
    })


def _dedup():
    """One cluster spanning both sources (all:1) + one go_bp-only cluster (all:2)."""
    return pd.DataFrame({
        "source": ["go_cc", "go_bp", "go_bp"],
        "group_id": ["GO:1", "GO:2", "GO:3"],
        "redundancy_cluster": ["all:1", "all:1", "all:2"],
        "is_representative": [True, False, True],
    })


def _group_members():
    """The (group, gene) long table, already reduced to cluster + gene."""
    return pd.DataFrame({
        "redundancy_cluster": ["all:1", "all:1", "all:1", "all:2", "all:2"],
        "gene": ["a", "b", "c", "c", "d"],
    })


def test_view_terms_slices_each_view_from_the_combined_table():
    from compute_view_fractions import view_terms

    metrics, dedup = _metrics(), _dedup()

    assert list(view_terms(metrics, dedup, "go_cc")["group_id"]) == ["GO:1"]
    assert list(view_terms(metrics, dedup, "combined")["group_id"]) == ["GO:1", "GO:2", "GO:3"]
    # dedup = the representative set, taken from the dedup table's own flag.
    assert list(view_terms(metrics, dedup, "dedup")["group_id"]) == ["GO:1", "GO:3"]
    with pytest.raises(ValueError, match="no terms for view"):
        view_terms(metrics, dedup, "go_bogus")


def test_breadth_is_counted_within_the_view_and_the_cut_is_the_views_own():
    """go_bp's terms touch both clusters, go_cc's only all:1 — so the counts differ.

    Both views cut at their own 75th percentile, but in go_bp `c` (two groups) is
    above it while in go_cc nothing is; that difference IS the per-view rule.
    """
    from compute_view_fractions import fractions_for_view, view_terms

    metrics, dedup, members = _metrics(), _dedup(), _group_members()

    bp, bp_genes, bp_cut = fractions_for_view(view_terms(metrics, dedup, "go_bp"), members)
    # Q75 of {1,1,1,2} (linear interpolation) = 1.25.
    assert bp_cut == pytest.approx(1.25)
    assert bp_genes.set_index("gene")["n_groups"].to_dict() == {"a": 1, "b": 1, "c": 2, "d": 1}
    assert set(bp_genes.loc[bp_genes["is_moonlighting"], "gene"]) == {"c"}
    # Per term, over its own members: GO:2 {a,c} -> 1/2 flagged; GO:3 {c,d} -> 1/2.
    assert bp.set_index("group_id")["moonlighting_fraction"].to_dict() == {
        "GO:2": 0.5, "GO:3": 0.5,
    }
    # paralog_fraction rides through from the metrics table.
    assert bp.set_index("group_id")["paralog_fraction"].to_dict() == {"GO:2": 0.0, "GO:3": 1.0}

    cc, cc_genes, cc_cut = fractions_for_view(view_terms(metrics, dedup, "go_cc"), members)
    assert cc_cut == pytest.approx(1.0)
    assert cc_genes.set_index("gene")["n_groups"].to_dict() == {"a": 1, "b": 1, "c": 1}
    assert not cc_genes["is_moonlighting"].any()
    assert cc["moonlighting_fraction"].iloc[0] == 0.0


def test_plot_fraction_distributions_draws_the_three_panels(tmp_path):
    """Three panels, and the mode is marked on the breadth one."""
    pytest.importorskip("cnsplots")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from coherence.fractions import MOONLIGHTING_QUANTILE
    from plot_fraction_distributions import plot_distributions

    assert MOONLIGHTING_QUANTILE == 0.75

    terms = pd.DataFrame({
        "paralog_fraction": [0.0, 0.5, 1.0, 0.25],
        "moonlighting_fraction": [0.0, 1.0, 1.0, 0.5],
        "moonlighting_cut_n_groups": [2.0, 2.0, 2.0, 2.0],
    })
    genes = pd.DataFrame({"n_groups": [1, 1, 2, 3, 4], "is_moonlighting": [False] * 2 + [True] * 3})

    plot_distributions(terms, genes)
    fig = plt.gcf()
    titles = [ax.get_title() for ax in fig.axes if ax.get_visible()]
    assert len(fig.axes) == 3
    assert titles == ["Paralog buffering", "Gene breadth", "Moonlighting"]
    # The cut line is labelled with its quantile and value — the cut the fraction uses.
    legends = [legend for ax in fig.axes if (legend := ax.get_legend()) is not None]
    assert "P75 = 2" in [text.get_text() for legend in legends for text in legend.get_texts()]
    plt.close("all")
