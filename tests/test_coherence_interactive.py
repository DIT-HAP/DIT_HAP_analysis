#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the two interactive Altair pages.

Both scripts build a Vega-Lite spec from hand-assembled frames, so the parts worth
asserting are the joins and the graph reduction that feed it: the scatter page's
member table must not invent coordinates for a gene the fitting table never fitted,
and the network's backbone must be a spanning tree whose edges all carry real
member overlap. The chart builders themselves are checked only for being
serializable, since a spec that parses can still be a broken figure and drawing one
needs a browser.
"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "workflow" / "src"))
sys.path.insert(0, str(_REPO_ROOT / "workflow" / "scripts" / "coherence"))

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

pytest.importorskip("altair")
pytest.importorskip("networkx")


# --- interactive scatter ----------------------------------------------------
def _scatter_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics = pd.DataFrame({
        "source": ["s"] * 3,
        "group_id": ["g1", "g2", "g3"],
        "group_name": ["coherent term", "middling term", "loose term"],
        "n_scored_members": [2, 2, 1],
        "median_pairwise_distance_z": [-3.0, 0.5, 2.0],
        "q_value": [0.01, 0.4, 0.9],
    })
    long_table = pd.DataFrame({
        "source": ["s"] * 5,
        "group_id": ["g1", "g1", "g2", "g2", "g3"],
        "group_name": ["coherent term", "coherent term", "middling term", "middling term", "loose term"],
        "Systematic ID": ["a", "b", "c", "unfitted", "e"],
        "Name": ["a", "b", "c", "unfitted", "e"],
    })
    fitting = pd.DataFrame({
        "Systematic ID": ["a", "b", "c", "e"],
        "DR": [-1.0, -1.1, -0.5, 0.0],
        "DL": [0.1, 0.2, 0.3, 0.4],
        "norm_DR": [-1.0, -1.1, -0.5, 0.0],
        "norm_DL": [0.01, 0.02, 0.03, 0.04],
        "R2": [0.9, 0.8, 0.7, 0.6],
        "FYPOviability": ["viable"] * 4,
    })
    return metrics, long_table, fitting


def test_member_table_drops_members_without_a_fit():
    from plot_interactive_scatter import member_table

    metrics, long_table, fitting = _scatter_frames()
    members = member_table(metrics, long_table, fitting, "s")
    # 'unfitted' is annotated into g2 but has no DR/DL, so drawing it would mean
    # inventing a coordinate for it.
    assert "unfitted" not in set(members["Systematic ID"])
    assert set(members["group_name"]) == {"coherent term", "middling term", "loose term"}


def test_term_options_open_on_the_most_coherent_term():
    from plot_interactive_scatter import term_options

    metrics, _, _ = _scatter_frames()
    assert term_options(metrics)[0] == "coherent term"


def test_scatter_spec_is_serializable_and_carries_a_dropdown():
    from plot_interactive_scatter import build_chart, member_table, term_options

    metrics, long_table, fitting = _scatter_frames()
    options = term_options(metrics)
    chart = build_chart(fitting, member_table(metrics, long_table, fitting, "s"), options)
    params = [p for p in chart.to_dict()["params"] if isinstance(p, dict)]
    pickers = [p for p in params if isinstance(p.get("bind"), dict)
               and p["bind"].get("input") == "select"]
    assert len(pickers) == 1, f"expected exactly one dropdown, got {len(pickers)}"
    assert pickers[0]["bind"]["options"] == options
    assert pickers[0]["value"] == "coherent term", "the page must open on a term, not blank"


# --- redundancy network -----------------------------------------------------
def _combined_table() -> pd.DataFrame:
    """Three terms, two of them sharing most of their members."""
    return pd.DataFrame({
        "source": ["go_cc"] * 3 + ["go_cc"],
        "group_id": ["t1", "t2", "t3", "solo"],
        "group_name": ["alpha", "alpha variant", "beta", "loner"],
        "scored_member_names": [
            ["a", "b", "c", "d"],
            ["a", "b", "c", "e"],       # high Jaccard with t1
            ["a", "b", "x", "y"],       # lower Jaccard with both
            ["z"],
        ],
    })


def _dedup_table() -> pd.DataFrame:
    return pd.DataFrame({
        "source": ["go_cc"] * 4,
        "group_id": ["t1", "t2", "t3", "solo"],
        "redundancy_cluster": ["all:0", "all:0", "all:0", "all:1"],
        "is_representative": [True, False, False, True],
        "representative_source": ["auto"] * 4,
        "representative_name": ["alpha"] * 3 + ["loner"],
        "median_pairwise_distance_z": [-4.0, -3.0, -1.0, 0.0],
        "q_value": [0.01, 0.02, 0.5, 0.9],
        "n_scored_members": [4, 4, 4, 1],
        "n_annotated_members": [10, 12, 40, 2],
    })


def test_network_frames_skip_clusters_of_one():
    from plot_redundancy_network import network_frames

    nodes, edges, _ = network_frames(_combined_table(), _dedup_table())
    # 'solo' was never merged with anything, so it has no decision to show.
    assert set(nodes["group_name"]) == {"alpha", "alpha variant", "beta"}
    assert len(edges) == 2, "a 3-term cluster's backbone is 2 edges"


def test_network_frames_name_the_cluster_after_its_representative():
    from plot_redundancy_network import network_frames

    nodes, _, _ = network_frames(_combined_table(), _dedup_table())
    # A cluster id says nothing about what was merged; the representative does.
    assert nodes["cluster_label"].unique().tolist() == ["alpha (3 terms)"]
    assert nodes.loc[nodes["cluster_label"] == "alpha (3 terms)", "node_id"].nunique() == 3


def test_overview_names_only_representatives(tmp_path):
    import json
    import re

    from plot_redundancy_network import network_frames, write_overview

    nodes, edges, _ = network_frames(_combined_table(), _dedup_table())
    out = tmp_path / "overview.html"
    write_overview(nodes, edges, out)
    page = out.read_text()
    drawn = json.loads(re.search(r"nodes\s*=\s*new vis\.DataSet\((\[.*?\])\);", page, re.S).group(1))

    # pyvis substitutes the node id for a falsy label, so an unfixed page comes back
    # with every absorbed term labelled "go_cc:t2" — 1,944 of them in the real run.
    labels = {node["id"]: node["label"] for node in drawn}
    assert labels == {"go_cc:t1": "alpha", "go_cc:t2": "", "go_cc:t3": ""}


def test_overview_disables_physics_and_dismisses_the_loading_bar(tmp_path):
    import json
    import re

    from plot_redundancy_network import network_frames, write_overview

    nodes, edges, _ = network_frames(_combined_table(), _dedup_table())
    out = tmp_path / "overview.html"
    write_overview(nodes, edges, out)
    page = out.read_text()

    # Every node already carries its layout position, so a physics simulation has
    # nothing to do but spin. It also gates pyvis's loading bar: that is dismissed by
    # vis.js's stabilizationIterationsDone handler, which never fires when physics is
    # off — leaving "0%" over a perfectly drawn graph, indefinitely.
    options = json.loads(re.search(r"var options = (\{.*?\});", page).group(1))
    assert options["physics"]["enabled"] is False
    assert "network.fit()" in page, "the view is never framed"
    assert 'getElementById("loadingBar")' in page, "the loading bar is never dismissed"


def test_cluster_places_puts_related_clusters_together():
    import math

    from plot_redundancy_network import cluster_places

    # A and B share most of their genes; C shares none with either. Relatedness is
    # measured on the REPRESENTATIVE terms, which is what these sets are.
    representatives = [
        {"a", "b", "c", "d"},
        {"a", "b", "c", "e"},
        {"x", "y", "z", "w"},
    ]
    places = cluster_places(representatives, [4, 4, 4])
    ab = math.dist(places[0][:2], places[1][:2])
    ac = math.dist(places[0][:2], places[2][:2])
    assert ab < ac, "a related pair was placed further apart than an unrelated one"


def test_cluster_places_leaves_no_overlapping_footprints():
    import math

    from plot_redundancy_network import _OVERVIEW_GAP, cluster_places

    representatives = [{"a"}, {"b"}, {"c"}, {"a"}, {"a", "b"}]
    sizes = [2, 2, 2, 2, 40]
    places = cluster_places(representatives, sizes)
    for i, (x1, y1, r1) in enumerate(places):
        for j, (x2, y2, r2) in enumerate(places):
            if i >= j:
                continue
            # Two discs that touch would still be drawn on top of each other.
            assert math.dist((x1, y1), (x2, y2)) >= r1 + r2 + _OVERVIEW_GAP - 1e-6


def test_pull_in_lonely_only_moves_the_clusters_with_no_relatives():
    import math

    import networkx as nx

    from plot_redundancy_network import pull_in_lonely

    graph = nx.Graph()
    graph.add_nodes_from(range(4))
    graph.add_edge(0, 1)          # 0 and 1 are related; 2 and 3 are not
    centres = [(0.0, 0.0), (10.0, 0.0), (500.0, 0.0), (9.0, 0.0)]
    pulled = pull_in_lonely(centres, graph)

    assert pulled[0] == centres[0] and pulled[1] == centres[1], "a connected cluster moved"
    assert math.dist(pulled[2], centres[2]) > 0, "the lonely cluster was not pulled in"
    assert math.dist(pulled[3], centres[3]) == 0, "a lonely cluster already inside the rim moved"


def test_overview_size_and_weight_follow_the_data():
    from plot_redundancy_network import overview_edge_width, overview_radius

    assert overview_radius(200) > overview_radius(10)
    assert overview_edge_width(0.9) > overview_edge_width(0.5)
    # Both are clamped, so an outlier cannot blow up the canvas.
    assert overview_radius(100000) == overview_radius(495)
    assert overview_edge_width(5.0) == overview_edge_width(1.0)


def test_backbone_edges_are_a_tree_and_keep_the_strongest_link():
    from plot_redundancy_network import backbone_edges

    members = [{"a", "b", "c", "d"}, {"a", "b", "c", "e"}, {"a", "b", "x", "y"}]
    edges = backbone_edges(members, ["alpha", "alpha variant", "beta"])
    assert len(edges) == 2
    # 3/5 for alpha vs alpha variant against 2/6 for either to beta, so the strong
    # pair must survive the reduction.
    strongest = max(edges, key=lambda edge: edge[2])
    assert {strongest[0], strongest[1]} == {0, 1}
    assert strongest[2] == pytest.approx(3 / 5)


def test_scale_positions_stay_inside_the_canvas():
    from plot_redundancy_network import scale_positions

    scaled = scale_positions({0: (0.0, 0.0), 1: (10.0, 4.0), 2: (3.0, 9.0)})
    for x, y in scaled.values():
        assert 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0


def test_network_spec_is_serializable():
    from plot_redundancy_network import build_chart, network_frames

    nodes, edges, labels = network_frames(_combined_table(), _dedup_table())
    spec = build_chart(nodes, edges, labels).to_dict()
    assert "jaccard" in str(spec), "the edge weights never reach the spec"
