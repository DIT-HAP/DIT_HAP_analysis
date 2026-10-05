#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Layout regression tests for the coherence figures.

Every defect this file guards against was found by rendering a figure and
looking at it, not by a failing assertion — a figure can render "successfully"
with a title running into the panel next to it, or a row's axis labels sitting on
top of the row below. Three separate instances turned up while porting these
figures to cnsplots:

- a panel title wider than its own 100 px axes (the named-group grid and the
  coherence overview both had 30+ character titles);
- a bottom row's x and tick labels overlapping the next row's titles, because
  multipanel measures only the left and top decorations after a draw;
- the same overflow colliding with the panel label beside it.

They are cheap to check and invisible to a passing run, so they get assertions.
cnsplots only exists in the cnsplots rule env, so this skips elsewhere.
"""

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cnsplots")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Overlaps shallower than this on both axes are grid-abutting artefacts, not defects.
# Measured on the real figure: neighbouring panels D and E show an 8 px tightbbox
# overlap (D's axes frame ends 8 px into E's decoration reserve) while the visible text
# there is 41 px apart, so a 2 px tolerance reports a collision a reader cannot see.
# The defects worth failing on are the ones that are visible — the row collision this
# file was written for was hundreds of pixels deep.
_TOLERANCE_PX = 12.0


def _overlaps(a, b, tolerance: float = _TOLERANCE_PX) -> bool:
    """Whether two bboxes overlap by more than `tolerance` px on BOTH axes."""
    # Neighbouring panels in a grid routinely abut, and multipanel's measured
    # decorations can land a fraction of a pixel inside the next panel's column.
    # The defects worth failing on are the ones a reader can see, which are tens
    # to hundreds of pixels deep, so anything under a couple of pixels is noise.
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx > tolerance and dy > tolerance


def assert_layout_is_clean(fig, *, check_labels: bool = True) -> None:
    """Assert no panel overflows its own axes and no two panels' decorations collide."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = [ax for ax in fig.axes if ax.get_visible() and not ax.get_label().startswith("<")]

    for ax in axes:
        title = ax.title.get_window_extent(renderer=renderer)
        box = ax.get_window_extent(renderer=renderer)
        assert title.width <= box.width, (
            f"title wider than its axes: {ax.get_title().splitlines()[0]!r} "
            f"({title.width:.0f} > {box.width:.0f} px)"
        )

    for left, right in itertools.combinations(axes, 2):
        assert not _overlaps(
            left.get_tightbbox(renderer), right.get_tightbbox(renderer)
        ), f"panel decorations collide: {left.get_title()!r} vs {right.get_title()!r}"

    if not check_labels:
        return
    for ax in axes:
        title = ax.title.get_window_extent(renderer=renderer)
        for text in ax.texts:
            if not text.get_text().strip():
                continue
            assert not _overlaps(text.get_window_extent(renderer=renderer), title), (
                f"panel label {text.get_text()!r} collides with title "
                f"{ax.get_title().splitlines()[0]!r}"
            )


def _coherence_table(n: int = 40) -> pd.DataFrame:
    return pd.DataFrame({
        "group_id": [f"GO:{i:07d}" for i in range(n)],
        "group_name": [f"a fairly long group name {i}" for i in range(n)],
        "n_scored_members": range(3, 3 + n),
        "median_pairwise_distance_z": [(-1) ** i * (i / 10) for i in range(n)],
        "q_value": [0.001 * (i + 1) for i in range(n)],
        "geom_median_DR": [0.1 * (i % 7) for i in range(n)],
        "geom_median_DL": [0.1 * (i % 5) for i in range(n)],
        "abundance_cv": [(i % 4) / 4 for i in range(n)],
        "conservation_cv": [(i % 6) / 6 for i in range(n)],
    })


def test_coherence_overview_layout_is_clean():
    from plot_coherence import plot_coherence

    plot_coherence(_coherence_table())
    assert_layout_is_clean(plt.gcf())
    plt.close("all")


def test_coherence_overview_without_biology_columns_is_clean():
    from plot_coherence import plot_coherence

    plot_coherence(_coherence_table().drop(columns=["abundance_cv", "conservation_cv"]))
    assert_layout_is_clean(plt.gcf())
    plt.close("all")


def test_group_scatter_layout_is_clean():
    """A long group name and id must not overflow the panel or hit the panel label."""
    from plot_group_scatter import plot_group_scatter_figure

    fitting = pd.DataFrame({
        "Systematic ID": [f"g{i}" for i in range(200)],
        "DR": [0.05 * i for i in range(200)],
        "DL": [0.02 * i for i in range(200)],
    })
    metrics = pd.DataFrame({
        "group_id": ["GO:0005762", "GO:0030684", "GO:0032040"],
        "median_pairwise_distance_z": [-3.99, -2.19, -6.13],
        "median_pairwise_distance_p": [0.01, 0.02, 0.01],
        "q_value": [0.0192, 0.0482, 0.0192],
    })
    resolved = [
        ("GO:0005762", "mitochondrial large ribosomal subunit", [f"g{i}" for i in range(47)]),
        ("GO:0030684", "preribosome", [f"g{i}" for i in range(76)]),
        ("GO:0032040", "small-subunit processome", [f"g{i}" for i in range(47)]),
        ("GO:0005681", "spliceosomal complex", [f"g{i}" for i in range(64)]),
        ("GO:9999999", "a very long group name that must be shortened to fit", [f"g{i}" for i in range(5)]),
    ]
    plot_group_scatter_figure(fitting, resolved, metrics)
    assert_layout_is_clean(plt.gcf())
    plt.close("all")


def test_group_scatter_panels_share_both_axes():
    """Every panel draws the same DR/DL space, so all of them must be linked.

    This asserts the sharing MECHANISM, not the resulting numbers: each panel also
    draws the genome-wide background cloud, so even fully autoscaled panels would
    land on the same limits today. What is being pinned is that the axes are
    genuinely shared, so a future change to what a panel draws (a subsampled
    background, a group reaching past the cloud) cannot silently desynchronise them.
    """
    from plot_group_scatter import plot_group_scatter_figure

    rng = np.random.default_rng(0)
    fitting = pd.DataFrame({
        "Systematic ID": [f"g{i}" for i in range(300)],
        "DR": rng.normal(-0.8, 0.4, 300),
        "DL": rng.normal(4.0, 2.0, 300),
    })
    metrics = pd.DataFrame({
        "group_id": ["GO:1", "GO:2", "GO:3"],
        "median_pairwise_distance_z": [-3.0, -1.0, 0.5],
        "median_pairwise_distance_p": [0.01, 0.02, 0.4], "q_value": [0.03, 0.05, 0.5],
    })
    resolved = [
        ("GO:1", "tight", [f"g{i}" for i in range(20)]),
        ("GO:2", "spread", [f"g{i}" for i in range(100, 200)]),
        ("GO:3", "mid", [f"g{i}" for i in range(200, 240)]),
    ]
    plot_group_scatter_figure(fitting, resolved, metrics)
    axes = [ax for ax in plt.gcf().axes if ax.get_visible()]
    assert len(axes) == 3
    for getter in ("get_shared_x_axes", "get_shared_y_axes"):
        siblings = getattr(axes[0], getter)().get_siblings(axes[0])
        assert set(siblings) == set(axes), f"panels are not linked on {getter}"
    plt.close("all")


def test_coherence_biology_panels_share_the_zscore_axis():
    """The biology panels plot z-score on y, so each must show the full z range.

    The limits are asserted exactly, with no autoscale margin, because that is what
    makes them comparable. x is deliberately NOT shared: the panels plot two
    different CVs on unrelated scales.
    """
    from plot_coherence import plot_coherence

    n = 57
    z_scores = np.linspace(-6.0, 3.0, n)
    table = pd.DataFrame({
        "group_id": [f"GO:{i:07d}" for i in range(n)],
        "group_name": [f"group {i}" for i in range(n)],
        "n_scored_members": range(3, 3 + n),
        "median_pairwise_distance_z": z_scores,
        "q_value": np.linspace(0.001, 0.4, n),
        "geom_median_DR": np.linspace(-1.2, 0.0, n),
        "geom_median_DL": np.linspace(0.0, 0.9, n),
        "abundance_cv": np.linspace(0.0, 3.0, n),
        "conservation_cv": np.linspace(0.0, 1.5, n),
    })
    plot_coherence(table)
    biology = [
        ax for ax in plt.gcf().axes
        if ax.get_visible() and ax.get_title().startswith(
            ("Abundance uniformity", "Conservation uniformity")
        )
    ]
    assert len(biology) == 2
    expected = (float(z_scores.min()), float(z_scores.max()))
    assert {ax.get_ylim() for ax in biology} == {expected}
    assert len({ax.get_xlim() for ax in biology}) == 2, "x is meant to stay per-panel"
    plt.close("all")


def _fdr_table(rows: list[tuple[float, float]]) -> pd.DataFrame:
    """(z-score, q_value) pairs -> the minimum columns labelled_extremes reads."""
    return pd.DataFrame({"median_pairwise_distance_z": [z for z, _ in rows],
                         "q_value": [q for _, q in rows]})


def test_labelled_extremes_caps_each_side():
    from plot_coherence import LabelSettings, labelled_extremes

    # 400 significant groups: 5% would be 20, so the cap is what bites.
    table = _fdr_table([(-i / 100.0, 0.001) for i in range(400)])
    settings = LabelSettings(quantile=0.05, q_max=0.05, max_labels=5)
    coherent = labelled_extremes(table, "coherent", settings)
    assert len(coherent) == 5
    # The most negative z first.
    assert list(coherent["median_pairwise_distance_z"]) == sorted(coherent["median_pairwise_distance_z"])


def test_labelled_extremes_quantile_bites_below_the_cap():
    from plot_coherence import LabelSettings, labelled_extremes

    # 60 significant groups, every one past the default z floor: 5% -> 3, fewer
    # than the cap of 5, so the quantile is what limits the count here.
    table = _fdr_table([(-2.5 - i / 10.0, 0.001) for i in range(60)])
    settings = LabelSettings(quantile=0.05, q_max=0.05, max_labels=5)
    assert len(labelled_extremes(table, "coherent", settings)) == 3


def test_labelled_extremes_ignores_fdr_on_the_incoherent_side():
    """No group more dispersed than random can be FDR-significant.

    The coherence p-value is one-sided for tightness, so an incoherent group's q sits
    near 1 by construction. Selecting that end by significance would return nothing,
    which is why it is selected by z alone.
    """
    from plot_coherence import LabelSettings, labelled_extremes

    table = _fdr_table([(-6.0, 0.001), (1.0, 0.99), (2.0, 1.0), (3.0, 0.98)])
    settings = LabelSettings(quantile=0.05, q_max=0.05, max_labels=5)
    assert labelled_extremes(table, "coherent", settings)["median_pairwise_distance_z"].tolist() == [-6.0]
    # quantile=1.0 so the per-side cap is what limits the count, not the 5% share.
    capped = LabelSettings(quantile=1.0, q_max=0.05, max_labels=2)
    assert labelled_extremes(table, "incoherent", capped)["median_pairwise_distance_z"].tolist() == [3.0, 2.0]


def test_place_label_avoids_a_point_it_could_cover():
    """A label is pushed clear of a drawn point when another candidate is free."""
    from plot_coherence import box_for, place_label

    point = (0.9, 0.2)          # the label's own point, right of centre -> text leftwards
    blocker = [(0.8, 0.2)]      # a point sitting exactly where the nearest candidate lands
    anchor_x, anchor_y, to_the_left = place_label(
        point, width=0.3, height=0.1, band=(0.02, 0.98), occupied=[], points=blocker
    )
    box = box_for(anchor_x, anchor_y, 0.3, 0.1, to_the_left)
    assert not (box[0] <= 0.8 <= box[1] and box[2] <= 0.2 <= box[3])


def test_place_label_keeps_its_box_inside_the_axes():
    from plot_coherence import box_for, place_label

    for point in [(0.02, 0.1), (0.5, 0.5), (0.98, 0.9)]:
        anchor_x, anchor_y, to_the_left = place_label(
            point, width=0.35, height=0.1, band=(0.02, 0.98), occupied=[], points=[]
        )
        box = box_for(anchor_x, anchor_y, 0.35, 0.1, to_the_left)
        assert 0.0 <= box[0] and box[1] <= 1.0, f"box escapes the axes for point {point}: {box}"
        assert 0.0 <= box[2] and box[3] <= 1.0, f"box escapes the axes for point {point}: {box}"


def test_fdr_label_text_does_not_overlap():
    """The leader-line labels must not collide with each other.

    Measured with ``Text.get_window_extent``, NOT ``Annotation.get_window_extent``:
    the latter spans the arrow too, so it reports leader lines as if they were text
    and drowns the real collisions in false positives.

    The fixture is sized so BOTH ends reach the per-side cap of 5 — with few labels
    the gap logic is never stressed and the assertion passes vacuously.
    """
    from matplotlib.text import Text

    from plot_coherence import plot_coherence

    n = 120  # per end, so ceil(5% * 120) = 6 > the cap of 5 on both sides
    table = pd.DataFrame({
        "group_id": [f"GO:{i:07d}" for i in range(2 * n)],
        "group_name": [f"a fairly long group name number {i}" for i in range(2 * n)],
        "n_scored_members": [3 + (i % 40) for i in range(2 * n)],
        "median_pairwise_distance_z": np.concatenate([np.linspace(0.5, 8.0, n), np.linspace(-0.5, -8.0, n)]),
        "q_value": np.linspace(0.001, 0.05, 2 * n),
        "geom_median_DR": np.linspace(-1.2, 0.0, 2 * n),
        "geom_median_DL": np.linspace(0.0, 0.9, 2 * n),
        "abundance_cv": np.linspace(0.0, 3.0, 2 * n),
        "conservation_cv": np.linspace(0.0, 1.5, 2 * n),
    })
    plot_coherence(table)
    fig = plt.gcf()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panels = [ax for ax in fig.axes if ax.get_visible() and ax.get_title() == "Coherence vs significance"]
    assert len(panels) == 2  # one per x encoding

    for ax in panels:
        labels = [text for text in ax.texts if text.get_text().strip()]
        assert len(labels) == 11, f"expected 5 + 5 labels plus the panel letter in {ax.get_xlabel()!r}"
        boxes = [Text.get_window_extent(text, renderer) for text in labels]
        for i, j in itertools.combinations(range(len(boxes)), 2):
            assert not _overlaps(boxes[i], boxes[j], tolerance=0.0), (
                f"label text overlaps in {ax.get_xlabel()!r}: "
                f"{labels[i].get_text().splitlines()[0]!r} vs {labels[j].get_text().splitlines()[0]!r}"
            )
        assert_layout_is_clean(fig)
    plt.close("all")


# --- cross-source comparison mode -------------------------------------------
# `color_by="source"` draws every panel once per source instead of once. The two
# ways it can go wrong are both layout defects: panels A/B switch from the house's
# filled histogram to overlaid step outlines with a legend added, and panels D-H go
# from one cns.scatterplot call to one per source — neither change is visible to a
# test that only checks the figure rendered.
def _sourced_coherence_table(sources: list[str], n_per_source: int = 30) -> pd.DataFrame:
    frames = []
    for position, source in enumerate(sources):
        frame = _coherence_table(n_per_source).copy()
        frame["source"] = source
        # Offset so the sources do not land on identical points, which would make a
        # colour check pass by accident on a single visible series.
        frame["median_pairwise_distance_z"] += position * 0.5
        frame["geom_median_DR"] += position * 0.05
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def test_comparison_layout_is_clean():
    from plot_coherence import plot_coherence

    table = _sourced_coherence_table(["go_macrocomplex", "go_cc", "go_bp", "kegg_brite", "kegg_pathway"])
    plot_coherence(table, color_by="source")
    assert_layout_is_clean(plt.gcf())
    plt.close("all")


def test_comparison_panel_c_has_one_colour_per_source_and_no_z_colourbar():
    from plot_coherence import _source_series, plot_coherence
    from coherence.palette import source_colors

    sources = ["go_macrocomplex", "go_cc", "go_bp", "kegg_brite", "kegg_pathway", "dedup"]
    table = _sourced_coherence_table(sources)
    plot_coherence(table, color_by="source")
    fig = plt.gcf()

    assert _source_series(table) == sources  # config order, not order of appearance
    expected = source_colors(sources)
    assert len(set(expected.values())) == len(sources), "two sources share a colour"

    centroid = next(ax for ax in fig.axes if ax.get_title() == "Group centroid positions")
    drawn = {
        tuple(np.round(collection.get_facecolors()[0][:3], 6))
        for collection in centroid.collections
        if len(collection.get_offsets())
    }
    want = {tuple(np.round(matplotlib.colors.to_rgb(expected[s]), 6)) for s in sources}
    assert drawn >= want, f"missing colours: {want - drawn}"
    # The z-score colourbar is what comparison mode replaces with the source key.
    assert not any(ax.get_label() == "<colorbar>" for ax in fig.axes)
    plt.close("all")


def test_single_series_mode_still_draws_without_a_source_column():
    """`color_by="none"` is the default and must not start requiring `source`."""
    from plot_coherence import plot_coherence

    table = _coherence_table()
    assert "source" not in table.columns
    plot_coherence(table)
    fig = plt.gcf()
    centroid = next(ax for ax in fig.axes if ax.get_title() == "Group centroid positions")
    # One data collection, not one per source. (The four others are the empty
    # size-legend handles, which carry no offsets.)
    data_collections = [c for c in centroid.collections if len(c.get_offsets())]
    assert len(data_collections) == 1
    assert len(data_collections[0].get_offsets()) == len(table)
    assert_layout_is_clean(fig)
    plt.close("all")
