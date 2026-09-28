#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Redundancy Network — an Altair HTML ledger of the de-duplication
==========================================================================

deduplicate_terms.py collapses redundant terms into clusters and keeps one
representative each, and all it leaves behind is a cluster id per row. That is
enough to *use* the reduction but not to audit it: you cannot see which terms
were judged redundant with which, how strongly, or what the term that lost out
looked like next to the one that won.

This renders that decision as an interactive network, one cluster at a time. A
dropdown picks the cluster; nodes are its terms, and the edges are the
maximum-Jaccard spanning tree — the backbone of the merge, with each edge
labelled by the Jaccard index of the two terms' scored-member sets. The
representative is coloured and sized apart from the terms it absorbed, and every
node's tooltip carries the evidence the pick was made on: q, z, member count and
its own source.

The spanning tree, not the full similarity graph: a 169-term cluster has ~1,500
pairs at or above the merge threshold, which is a hairball rather than a figure.
The cluster is a connected component, so a spanning tree shows the same
connectivity with n-1 edges. Edges below the threshold still appear when they are
the strongest link available — that is the case the DAG-lineage rule produces.

Clusters of one term were never merged and are not listed. Clusters above
`_MAX_LABELLED` terms are drawn without node names, which at that size would be
unreadable regardless; hover still identifies a node.

Input
-----
- coherence_metrics_combined.parquet: all sources' metrics, joined on
  (source, group_id) for the member sets the Jaccard is computed from.
- coherence_terms_deduplicated.tsv: the annotated table from deduplicate_terms.py
  (redundancy_cluster, is_representative, representative_source, ...).

Output
------
- redundancy_network.html: one cluster at a time, Altair; the dropdown lists clusters
  by their representative's name. Loads vega-embed from a CDN.
- redundancy_overview.html: every cluster on one pan-and-zoom canvas, pyvis/vis.js.
  Loads vis.js from a CDN. Both pages need a network connection the first time.

Usage
-----
    python plot_redundancy_network.py \\
        --combined results/3a_coherence/{dataset}/coherence_metrics_combined.parquet \\
        --deduplicated results/3a_coherence/{dataset}/coherence_terms_deduplicated.tsv \\
        --output tmp/{dataset}/redundancy_network.html \\
        --output-overview tmp/{dataset}/redundancy_overview.html

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-09-22
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import html
import json
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import altair as alt
import networkx as nx
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from pyvis.network import Network

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
# The merge criterion itself, imported from the stage that applied it rather than
# re-implemented: an audit figure whose Jaccard drifted from deduplicate_terms.py's
# would be describing a merge that never happened. The script is import-safe (its
# work is behind a __main__ guard), which is also how tests/test_coherence_dedup.py
# reaches it.
sys.path.append(str(SCRIPT_DIR))
from deduplicate_terms import jaccard_index, member_set  # noqa: E402
from figures import FURNITURE_COLOR, house_colors  # noqa: E402
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
_REQUIRED_COMBINED_COLUMNS = ["source", "group_id", "group_name", "scored_member_names"]
_REQUIRED_DEDUP_COLUMNS = ["source", "group_id", "redundancy_cluster", "is_representative",
                           "representative_source", "median_pairwise_distance_z",
                           "q_value", "n_scored_members"]

# Columns taken from the de-duplicated table. Its `scored_member_names` is a
# comma-joined STRING (it is a human-facing TSV), while the combined Parquet's is a
# real list column — joining on both and letting name suffixes pick a winner would
# hand `member_set` the string and raise. Only the columns this script needs cross
# the join, and the member sets come from the Parquet alone.
_DEDUP_COLUMNS = ["source", "group_id", "redundancy_cluster", "is_representative",
                  "representative_source", "representative_name",
                  "median_pairwise_distance_z", "q_value", "n_scored_members",
                  "n_annotated_members"]

# Above this many terms a cluster is drawn without node names: the labels would
# overlap into a grey block, and the hover tooltip still identifies a node.
_MAX_LABELLED = 30

# Chart footprint in pixels, and the layout normalizes into it. Read on screen, so
# the canvas is larger than any journal panel.
_CHART_WIDTH = 700
_CHART_HEIGHT = 620
_LABEL_FONT_SIZE = 10

# Node footprint by the term's own gene count (annotated members), so a cluster's
# biggest term reads as its biggest term. The square-root keeps a 500-gene term
# from dwarfing a 3-gene one; the range is calibrated on the real span (3..500
# annotated members) so the small terms are not all pinned to the floor.
_NODE_SIZE_MIN, _NODE_SIZE_MAX = 60.0, 620.0

# Edge weight carries the Jaccard on the edge itself, as width. The pyvis overview
# has 1,600 edges on screen at once and needs the weak ones to recede, which it
# does through the alpha baked into each edge's rgba colour.
_EDGE_WIDTH_RANGE = (1.0, 6.0)

# Spring layout is seeded: the figure is a record of a decision, so it has to come
# back the same on every run.
_LAYOUT_SEED = 42
_LAYOUT_ITERATIONS = 60
# Positions are normalized into this margin so no node sits on the canvas edge.
_LAYOUT_MARGIN = 0.06

# Term colour by source, so a cluster that collapsed terms across databases reads as
# one at a glance. Taken from the house palette positions that stay distinguishable
# in greyscale (teal / amber) plus the house red.
_SOURCE_ORDER = ["go_bp", "go_cc", "go_macrocomplex"]
_SOURCE_COLORS = dict(zip(_SOURCE_ORDER, house_colors((1, 2, 0))))
_HOUSE_RED = house_colors((0,))[0]

# --- pyvis overview ---------------------------------------------------------
# One page with every cluster on it. Positions are in vis.js's arbitrary
# coordinate units, and only their RATIO to the node sizes matters: vis.js fits the
# coordinate extent to the viewport and draws node sizes in pixels at the same zoom.
#
# A cluster's footprint grows as the 0.35 power of its term count, and its own
# normalised layout is scaled into that footprint. The gap is what a reader sees as
# the space between one network and the next.
_OVERVIEW_CLUSTER_UNIT = 110.0   # footprint of a single-term cluster
_OVERVIEW_SIZE_EXPONENT = 0.35
_OVERVIEW_GAP = 160.0            # clearance between two clusters' footprints
_OVERVIEW_SPREAD = 6.0           # initial spring box, in units of the radius mass
_OVERVIEW_SEPARATION_PASSES = 60
_OVERVIEW_RADIUS = (6.0, 26.0)   # vis.js node radius, sqrt-scaled by gene count
_OVERVIEW_EDGE_WIDTH = (0.3, 2.6)
_OVERVIEW_HEIGHT_PX = 950

# Columns the Altair page actually reads, per layer. See build_chart.
_NODE_FIELDS = ["cluster_label", "group_name", "group_id", "source", "role",
                "representative_source", "median_pairwise_distance_z", "q_value",
                "n_scored_members", "n_annotated_members", "x", "y"]
_EDGE_FIELDS = ["cluster_label", "from", "to", "jaccard", "x", "y", "x2", "y2",
                "mid_x", "mid_y"]
_LABEL_FIELDS = ["cluster_label", "group_name", "x", "y"]

_NODE_TOOLTIP = [
    alt.Tooltip("group_name:N", title="Term"),
    alt.Tooltip("group_id:N", title="Term ID"),
    alt.Tooltip("source:N", title="Source"),
    alt.Tooltip("role:N", title="Role"),
    alt.Tooltip("representative_source:N", title="Picked by"),
    alt.Tooltip("median_pairwise_distance_z:Q", title="z-score", format=".2f"),
    alt.Tooltip("q_value:Q", title="q (BH)", format=".3g"),
    alt.Tooltip("n_scored_members:Q", title="Scored members"),
    alt.Tooltip("n_annotated_members:Q", title="Annotated genes"),
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Inputs and outputs for the redundancy-network pages."""
    combined: Path
    deduplicated: Path
    output: Path
    output_overview: Path

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then create the output dirs."""
        for path in [self.combined, self.deduplicated]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output, self.output_overview]:
            out.parent.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True, slots=True)
class _Cluster:
    """One redundancy cluster: what it is called, what it holds, and where it goes."""
    label: str
    rows: pd.DataFrame
    names: list[str]
    edges: list[tuple[int, int, float]]
    positions: dict[int, tuple[float, float]]
    representative_members: set[str]


# =============================================================================
# CORE LOGIC — one cluster's backbone
# =============================================================================
def cluster_terms(table: pd.DataFrame, cluster: str) -> pd.DataFrame:
    """The rows of one cluster, representative first so it draws on top."""
    rows = table[table["redundancy_cluster"] == cluster]
    return rows.sort_values("is_representative", ascending=False)


def backbone_edges(members: list[set[str]], names: list[str]) -> list[tuple[int, int, float]]:
    """The maximum-Jaccard spanning tree of a cluster, as (i, j, jaccard) edges."""
    # A complete graph over the cluster, weighted by member-set Jaccard, reduced to
    # its maximum spanning tree. The cluster is a connected component of the
    # above-threshold graph, so the tree reaches every term while drawing n-1
    # edges instead of the ~1,500 a 169-term cluster actually has.
    graph = nx.Graph()
    graph.add_nodes_from(range(len(names)))
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            weight = jaccard_index(members[i], members[j])
            if weight > 0.0:
                # networkx maximizes weight, and a Jaccard of 0 is no evidence of a
                # merge, so zero-weight pairs stay out of the candidate set entirely.
                graph.add_edge(i, j, weight=weight)

    tree = nx.maximum_spanning_tree(graph) if graph.number_of_edges() else graph
    if tree.number_of_edges() and not nx.is_connected(graph):
        logger.warning(
            f"cluster of {len(names)} terms is not connected by shared members; "
            f"drawing the backbone of its largest component"
        )
    return [(i, j, float(data["weight"])) for i, j, data in tree.edges(data=True)]


def cluster_radii(sizes: list[int]) -> list[float]:
    """The footprint each cluster gets on the overview, from its term count."""
    # Sub-linear in the term count: a cluster of 169 needs more room than one of 2,
    # but not 85x more, or the small ones shrink to specks in the gaps between the
    # big ones. The exponent sits between "every cluster the same size" (0) and
    # "area proportional to the term count" (0.5).
    return [_OVERVIEW_CLUSTER_UNIT * float(size) ** _OVERVIEW_SIZE_EXPONENT for size in sizes]


def separate(
    centres: list[tuple[float, float]], radii: list[float]
) -> list[tuple[float, float]]:
    """Push apart overlapping cluster discs, keeping the arrangement they came in."""
    # The spring layout knows nothing about how much room a cluster's own layout
    # needs, so its discs overlap. Resolving them pairwise — same as the classic
    # de-overlap pass — keeps the neighbourhood structure the layout encoded, which
    # is the whole reason for using it.
    places = [list(centre) for centre in centres]
    for _ in range(_OVERVIEW_SEPARATION_PASSES):
        moved = False
        for i in range(len(places)):
            for j in range(i + 1, len(places)):
                dx = places[j][0] - places[i][0]
                dy = places[j][1] - places[i][1]
                distance = float(np.hypot(dx, dy))
                needed = radii[i] + radii[j] + _OVERVIEW_GAP
                if distance >= needed:
                    continue
                # Coincident centres have no direction to push along; nudge them
                # apart deterministically rather than dividing by zero.
                if distance < 1e-9:
                    dx, dy, distance = 1.0, 0.0, 1.0
                shift = (needed - distance) / 2.0 / distance
                places[i][0] -= dx * shift
                places[i][1] -= dy * shift
                places[j][0] += dx * shift
                places[j][1] += dy * shift
                moved = True
        if not moved:
            break
    return [tuple(place) for place in places]


def cluster_places(
    representatives: list[set[str]], sizes: list[int]
) -> list[tuple[float, float, float]]:
    """Where each cluster sits on the overview and how much room it has."""
    # Two clusters are as related as their REPRESENTATIVE terms are, measured the
    # same way the merge itself was — Jaccard over the scored member sets. Using the
    # merge's own yardstick means "near each other" here says the same thing as
    # "close to being merged" did, and needs no new input: a pair that shared enough
    # genes to clear the threshold would have been one cluster already, so what
    # survives is exactly the related-but-distinct neighbourhoods worth seeing.
    #
    # The weights are rescaled to 0..1 for the layout. Raw Jaccard is bounded by the
    # merge threshold — anything above it is already one cluster — so every surviving
    # edge sits in a narrow band of small values, and a spring whose pull is
    # proportional to them barely pulls at all. Against a measured baseline of 0.0036
    # mean Jaccard over all cluster pairs, the 3 nearest clusters come out at 0.0024
    # mean (i.e. no better than chance) with the raw weights and 0.051 with these.
    graph = nx.Graph()
    graph.add_nodes_from(range(len(representatives)))
    weights = []
    for i in range(len(representatives)):
        for j in range(i + 1, len(representatives)):
            weight = jaccard_index(representatives[i], representatives[j])
            if weight > 0.0:
                weights.append(weight)
                graph.add_edge(i, j, weight=weight)
    if weights:
        strongest = max(weights)
        for _, _, data in graph.edges(data=True):
            data["weight"] /= strongest

    radii = cluster_radii(sizes)
    if graph.number_of_edges():
        raw = nx.spring_layout(
            graph, seed=_LAYOUT_SEED, iterations=_LAYOUT_ITERATIONS,
            weight="weight", scale=1.0,
        )
    else:
        raw = nx.circular_layout(graph)
    # The spring layout works in a unit box; scale that box by the room the discs
    # need. Generous, because the arrangement it encodes is the point of the figure
    # and the separation pass below can only blur it: measured on the real table, the
    # similarity of each cluster's three nearest neighbours rises from 4.9x to 14.2x
    # the chance level between 3x and 6x, while the outline stays as compact.
    spread = _OVERVIEW_SPREAD * float(np.sqrt(sum(radius ** 2 for radius in radii)))
    centres = [(float(raw[i][0] - 0.5) * spread, float(raw[i][1] - 0.5) * spread)
               for i in range(len(representatives))]
    centres = pull_in_lonely(centres, graph)
    return [(x, y, radius) for (x, y), radius in zip(separate(centres, radii), radii)]


def pull_in_lonely(
    centres: list[tuple[float, float]], graph: nx.Graph
) -> list[tuple[float, float]]:
    """Bring clusters with no relatives back onto the rim of the connected mass."""
    # A cluster whose representative shares no gene with any other representative has
    # no position the layout can justify, and a repulsion-only spring answers that by
    # pushing it into the void — which costs the page its framing, because the whole
    # canvas has to be fitted around wherever it lands. On the real table that is 5
    # of 324 clusters, and pulling just those back in cuts the canvas from 40,600 by
    # 29,600 units to 13,700 by 20,500 without moving a single other cluster.
    lonely = [index for index in graph.nodes if graph.degree(index) == 0]
    if not lonely:
        return centres
    connected = [centres[i] for i in graph.nodes if graph.degree(i) > 0]
    if not connected:
        return centres
    centre_x = float(np.mean([position[0] for position in centres]))
    centre_y = float(np.mean([position[1] for position in centres]))
    rim = float(np.percentile(
        [np.hypot(x - centre_x, y - centre_y) for x, y in connected], 90,
    ))

    pulled = list(centres)
    for index in lonely:
        x, y = centres[index]
        distance = float(np.hypot(x - centre_x, y - centre_y))
        if distance > rim:
            scale = rim / distance
            pulled[index] = (centre_x + (x - centre_x) * scale,
                             centre_y + (y - centre_y) * scale)
    return pulled


def layout(nodes: int, edges: list[tuple[int, int, float]]) -> dict[int, tuple[float, float]]:
    """Spring-layout positions for one cluster, normalized into [0, 1]."""
    graph = nx.Graph()
    graph.add_nodes_from(range(nodes))
    graph.add_weighted_edges_from((i, j, w) for i, j, w in edges)
    if graph.number_of_edges():
        # Stronger (more redundant) pairs are pulled closer, so distance on the page
        # reads as dissimilarity about as well as a 2-D projection of it can.
        raw = nx.spring_layout(
            graph, seed=_LAYOUT_SEED, iterations=_LAYOUT_ITERATIONS,
            weight="weight", scale=1.0,
        )
    else:
        raw = nx.circular_layout(graph)
    return {index: tuple(map(float, position)) for index, position in raw.items()}


def scale_positions(positions: dict[int, tuple[float, float]]) -> dict[int, tuple[float, float]]:
    """Rescale one cluster's layout into the canvas, preserving its aspect."""
    xs = [x for x, _ in positions.values()]
    ys = [y for _, y in positions.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    centre_x, centre_y = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    span = span * (1.0 + 2 * _LAYOUT_MARGIN)
    return {
        index: (0.5 + (x - centre_x) / span, 0.5 + (y - centre_y) / span)
        for index, (x, y) in positions.items()
    }


# =============================================================================
# CORE LOGIC — assembling the frames the page draws
# =============================================================================
def cluster_backbone(rows: pd.DataFrame) -> tuple[list[str], list[tuple[int, int, float]], dict]:
    """The member sets, the Jaccard backbone, and a layout for one cluster."""
    names = rows["group_name"].astype(str).tolist()
    members = [member_set(value) for value in rows["scored_member_names"]]
    edges = backbone_edges(members, names)
    return names, edges, scale_positions(layout(len(rows), edges))


def network_frames(combined: pd.DataFrame, table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The (nodes, edges, labels) frames for every cluster worth drawing."""
    merged = table[_DEDUP_COLUMNS].merge(
        combined[["source", "group_id", "group_name", "scored_member_names"]],
        on=["source", "group_id"], how="left",
    )
    # A term missing from the combined table cannot contribute a member set, and a
    # cluster built without it would silently misstate the merge.
    missing = int(merged["scored_member_names"].isna().sum())
    if missing:
        logger.warning(f"{missing} deduplicated term(s) absent from the combined table; excluded")
        merged = merged[merged["scored_member_names"].notna()]

    # Every cluster laid out and labelled first, because the overview places them
    # against each other and so needs all of them before it can put any one down.
    clusters = []
    for cluster, rows in merged.groupby("redundancy_cluster", sort=False):
        if len(rows) < 2:
            continue
        rows = cluster_terms(rows, cluster)
        names, edges, positions = cluster_backbone(rows)
        # The representative's name, not the cluster id: "all:15" says nothing about
        # what was merged, and picking which term to audit is the whole point of the
        # list. The size rides along because it is what makes a big cluster big.
        label = f"{rows['representative_name'].iloc[0]} ({len(rows)} terms)"
        # The representative is drawn first by cluster_terms, so index 0 is it — and
        # its member set is what the overview measures cluster relatedness with.
        representative_members = member_set(rows["scored_member_names"].iloc[0])
        clusters.append(_Cluster(label, rows, names, edges, positions, representative_members))

    # Related clusters next to each other, each with room for its own layout. The
    # packing is in `cluster_places`; this only converts the result into per-node
    # coordinates by dropping each cluster's normalised layout into its footprint.
    places = cluster_places(
        [cluster.representative_members for cluster in clusters],
        [len(cluster.rows) for cluster in clusters],
    )

    node_rows, edge_rows, label_rows = [], [], []
    for cluster, (centre_x, centre_y, radius) in zip(clusters, places):
        label, rows, names, edges, positions = (
            cluster.label, cluster.rows, cluster.names, cluster.edges, cluster.positions
        )

        def overview(index: int) -> tuple[float, float]:
            x, y = positions[index]
            return centre_x + (x - 0.5) * 2 * radius, centre_y + (y - 0.5) * 2 * radius

        for index, (_, row) in enumerate(rows.iterrows()):
            x, y = positions[index]
            ox, oy = overview(index)
            node_rows.append({
                "cluster_label": label,
                # Qualified by source: 173 group_ids appear in more than one source,
                # and the overview addresses nodes by a single id.
                "node_id": f"{row['source']}:{row['group_id']}",
                "group_name": row["group_name"],
                "group_id": row["group_id"],
                "source": row["source"],
                "role": "representative" if row["is_representative"] else "member",
                "representative_source": row.get("representative_source"),
                "median_pairwise_distance_z": row["median_pairwise_distance_z"],
                "q_value": row["q_value"],
                "n_scored_members": row["n_scored_members"],
                "n_annotated_members": row["n_annotated_members"],
                "x": x, "y": y, "ox": ox, "oy": oy,
            })
            if len(rows) <= _MAX_LABELLED:
                label_rows.append({
                    "cluster_label": label, "group_name": row["group_name"], "x": x, "y": y,
                })

        for i, j, weight in edges:
            x1, y1 = positions[i]
            x2, y2 = positions[j]
            ox1, oy1 = overview(i)
            ox2, oy2 = overview(j)
            edge_rows.append({
                "cluster_label": label,
                "from": names[i], "to": names[j], "jaccard": weight,
                "from_id": f"{rows['source'].iloc[i]}:{rows['group_id'].iloc[i]}",
                "to_id": f"{rows['source'].iloc[j]}:{rows['group_id'].iloc[j]}",
                "x": x1, "y": y1, "x2": x2, "y2": y2,
                "mid_x": (x1 + x2) / 2, "mid_y": (y1 + y2) / 2,
                "ox": ox1, "oy": oy1, "ox2": ox2, "oy2": oy2,
            })

    return pd.DataFrame(node_rows), pd.DataFrame(edge_rows), pd.DataFrame(label_rows)


def build_chart(nodes: pd.DataFrame, edges: pd.DataFrame, labels: pd.DataFrame) -> alt.Chart:
    """The layered page: lines, their Jaccard labels, then nodes over the top."""
    # Altair embeds a layer's ENTIRE frame in the page, once per layer it appears in,
    # so the overview-only columns (ox/oy/node_id/from_id/to_id) are dropped here
    # rather than shipped twice over 2,000 rows.
    nodes = nodes[_NODE_FIELDS]
    edges = edges[_EDGE_FIELDS]
    labels = labels[_LABEL_FIELDS]
    # Largest cluster first: the dropdown opens on the merge that most needs
    # auditing, and a reader scrolling it wants the big ones at the top. It opens on
    # the largest *labelled* one, though — a first view of 169 unlabelled dots reads
    # as a failure rather than as a deliberate simplification.
    order = nodes["cluster_label"].value_counts().index.tolist()
    labelled = set(labels["cluster_label"].unique())
    picker = alt.selection_point(
        fields=["cluster_label"],
        bind=alt.binding_select(options=order, name="Redundancy cluster  "),
        value=next((option for option in order if option in labelled), order[0]),
    )
    axes = {"x": alt.X("x:Q", axis=None), "y": alt.Y("y:Q", axis=None)}

    links = (
        alt.Chart(edges)
        .transform_filter(picker)
        .mark_rule(color=FURNITURE_COLOR)
        .encode(
            **axes, x2="x2:Q", y2="y2:Q",
            strokeWidth=alt.StrokeWidth(
                "jaccard:Q", title="Jaccard",
                scale=alt.Scale(domain=[0, 1], range=list(_EDGE_WIDTH_RANGE)), legend=None,
            ),
            tooltip=[alt.Tooltip("from:N", title="Term"), alt.Tooltip("to:N", title="with"),
                     alt.Tooltip("jaccard:Q", title="Jaccard", format=".3f")],
        )
    )
    # The Jaccard is what justifies each edge, so on a readable cluster it is printed
    # on the edge rather than left to a hover: the numbers ARE the argument this
    # figure makes. Above _MAX_LABELLED edges they would overlap into a grey fuzz, so
    # there the hover carries them.
    weights = (
        alt.Chart(edges[edges["cluster_label"].isin(labelled)])
        .transform_filter(picker)
        .mark_text(fontSize=9, color=FURNITURE_COLOR, dy=-4)
        .encode(
            # The edge's midpoint, not its endpoints: anchored on a node the number
            # would be drawn under that node's marker and never seen.
            x=alt.X("mid_x:Q", axis=None), y=alt.Y("mid_y:Q", axis=None),
            text=alt.Text("jaccard:Q", format=".2f"),
        )
    )
    names = (
        alt.Chart(labels)
        .transform_filter(picker)
        .mark_text(fontSize=_LABEL_FONT_SIZE, dy=13)
        .encode(**axes, text="group_name:N")
    )
    points = (
        alt.Chart(nodes)
        .transform_filter(picker)
        .mark_circle(opacity=0.95, stroke="white", strokeWidth=1)
        .encode(
            **axes,
            # Size is the term's own gene count, colour its source: together they say
            # which term in the cluster carries the signal and whether the merge
            # crossed databases. The representative is called out separately, since
            # "which one won" is the decision the figure exists to show.
            #
            # The scale is square-rooted as well as clamped to the table's own range:
            # a cluster's terms usually sit in a narrow band of gene counts, and a
            # linear map over 0..500 renders all of them at the same size.
            #
            # Both legends sit outside the frame: inside, they land on the cluster
            # whatever corner they take, because the spring layout fills the panel.
            size=alt.Size("n_annotated_members:Q", title="Annotated genes",
                          scale=alt.Scale(type="sqrt",
                                          domain=[0, float(nodes["n_annotated_members"].max())],
                                          range=[_NODE_SIZE_MIN, _NODE_SIZE_MAX]),
                          legend=alt.Legend(orient="right", tickCount=5)),
            color=alt.Color("source:N", title="Source",
                            scale=alt.Scale(domain=_SOURCE_ORDER,
                                            range=[_SOURCE_COLORS[s] for s in _SOURCE_ORDER]),
                            legend=alt.Legend(orient="right", columns=1)),
            stroke=alt.condition("datum.role === 'representative'",
                                 alt.value(_HOUSE_RED), alt.value("white")),
            strokeWidth=alt.condition("datum.role === 'representative'",
                                      alt.value(3), alt.value(1)),
            tooltip=_NODE_TOOLTIP,
        )
    )
    # The note appears only when some cluster really is drawn unlabelled, which is
    # knowable from the frames: the label frame skips the oversized clusters.
    unlabelled = nodes["cluster_label"].nunique() > labels["cluster_label"].nunique()
    subtitle = (
        "Edge = maximum-Jaccard spanning tree of the merge; its label is the Jaccard "
        "index of the two terms' scored-member sets."
        + (f" Clusters above {_MAX_LABELLED} terms are drawn unlabelled." if unlabelled else "")
    )
    return (
        (links + weights + names + points)
        .add_params(picker)
        .properties(
            width=_CHART_WIDTH, height=_CHART_HEIGHT,
            title=alt.TitleParams("Coherence redundancy clusters", subtitle=subtitle),
        )
        .configure_view(stroke=None)
    )


# =============================================================================
# CORE LOGIC — the one-page overview (pyvis / vis.js)
# =============================================================================
def overview_radius(n_annotated_members: float) -> float:
    """vis.js node radius for a term, from the term's own gene count."""
    scaled = 3.0 + 1.4 * np.sqrt(max(float(n_annotated_members), 1.0))
    return float(min(max(scaled, _OVERVIEW_RADIUS[0]), _OVERVIEW_RADIUS[1]))


def overview_edge_width(jaccard: float) -> float:
    """vis.js edge width from the Jaccard, so a weak link is a thin line."""
    low, high = _OVERVIEW_EDGE_WIDTH
    return float(low + (high - low) * min(max(jaccard, 0.0), 1.0))


def node_tooltip(row: dict) -> str:
    """The hover card for one term, as the HTML vis.js renders in a title."""
    # Escaped: term names carry apostrophes and ampersands, and vis.js sets the
    # title as innerHTML rather than as text.
    return (
        f"<b>{html.escape(str(row['group_name']))}</b><br>"
        f"{html.escape(str(row['group_id']))} · {html.escape(str(row['source']))}<br>"
        f"{row['role']} (picked by {html.escape(str(row['representative_source']))})<br>"
        f"z = {row['median_pairwise_distance_z']:.2f}, q = {row['q_value']:.3g}<br>"
        f"{int(row['n_annotated_members'])} annotated genes, "
        f"{int(row['n_scored_members'])} scored"
    )


def write_overview(nodes: pd.DataFrame, edges: pd.DataFrame, path: Path) -> None:
    """Every cluster on one pan-and-zoom canvas, representatives named."""
    # This is the landscape view the per-cluster page cannot give: which merges
    # happened, how big they are, and where they sit relative to each other. Names
    # are drawn for representatives only — 1,944 labels is not a figure — so each
    # cluster reads as one annotated blob, and every node still answers on hover.
    net = Network(
        height=f"{_OVERVIEW_HEIGHT_PX}px", width="100%", bgcolor="#ffffff",
        font_color="#222222", directed=False, notebook=False, cdn_resources="remote",
    )
    net.toggle_physics(False)

    representatives = []
    for row in nodes.to_dict("records"):
        representative = row["role"] == "representative"
        if representative:
            representatives.append(row["node_id"])
        net.add_node(
            row["node_id"],
            label=row["group_name"] if representative else "",
            title=node_tooltip(row),
            x=row["ox"], y=row["oy"],
            size=overview_radius(row["n_annotated_members"]),
            shape="dot",
            color={
                "background": _SOURCE_COLORS.get(row["source"], FURNITURE_COLOR),
                "border": _HOUSE_RED if representative else "#ffffff",
                "highlight": {"background": _HOUSE_RED, "border": _HOUSE_RED},
            },
            borderWidth=3 if representative else 1,
            physics=False,
        )
    # pyvis substitutes the node id for any falsy label, so every absorbed term
    # would come back labelled "go_bp:GO:0006235" — 1,944 of them. Blank them after
    # the fact, which is the only place the substitution can be undone.
    for node_id in set(nodes["node_id"]) - set(representatives):
        net.get_node(node_id)["label"] = ""

    for row in edges.to_dict("records"):
        net.add_edge(
            row["from_id"], row["to_id"],
            width=overview_edge_width(row["jaccard"]),
            color=f"rgba(120, 120, 120, {0.25 + 0.6 * row['jaccard']:.2f})",
            title=f"Jaccard {row['jaccard']:.3f}",
        )

    net.set_options(json.dumps({
        # `set_options` REPLACES pyvis's options object rather than merging into it,
        # so this has to be here as well as in toggle_physics() above: calling the
        # latter on its own is silently undone by this call, and the page then runs
        # a physics simulation over 1,944 equally-positioned nodes to no effect.
        "physics": {"enabled": False},
        "interaction": {"hover": True, "tooltipDelay": 120, "navigationButtons": True,
                        "keyboard": {"enabled": False}},
        "nodes": {"font": {"size": 14, "face": "sans-serif"}},
        "edges": {"smooth": False},
    }))
    net.write_html(str(path), notebook=False, open_browser=False)
    close_loading_bar(path)


def close_loading_bar(path: Path) -> None:
    """Dismiss pyvis's loading bar and fit the view, which pyvis leaves to vis.js."""
    # Two things pyvis's template only does inside `stabilizationIterationsDone`, an
    # event that never fires when the physics is off — which is this figure's whole
    # point, since every node carries a fixed position. The bar would sit at "0%"
    # over a perfectly drawn graph forever, and the view would stay wherever vis.js
    # happened to leave it instead of framing the whole layout.
    script = """
<script type="text/javascript">
  window.addEventListener("load", function () {
    var bar = document.getElementById("loadingBar");
    if (bar) { bar.style.display = "none"; }
    if (typeof network !== "undefined" && network) { network.fit(); }
  });
</script>
"""
    page = path.read_text()
    path.write_text(page.replace("</body>", f"{script}</body>"))


# =============================================================================
# MAIN EXECUTION
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Load -> build every cluster's backbone -> save the interactive page."""
    config.validate()
    combined = read_parquet(config.combined)
    missing = [column for column in _REQUIRED_COMBINED_COLUMNS if column not in combined.columns]
    if missing:
        raise ValueError(f"combined metrics missing required column(s) {missing}")
    table = pd.read_csv(config.deduplicated, sep="\t")
    missing = [column for column in _REQUIRED_DEDUP_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"deduplicated table missing required column(s) {missing} (have: {list(table.columns)})")

    nodes, edges, labels = network_frames(combined, table)
    if nodes.empty:
        raise ValueError("no cluster holds more than one term; nothing to draw")
    logger.info(
        f"{table['redundancy_cluster'].nunique():,} clusters, {len(nodes):,} terms over "
        f"{nodes['cluster_label'].nunique():,} merged clusters, {len(edges):,} backbone edges"
    )
    build_chart(nodes, edges, labels).save(config.output)
    logger.success(f"Wrote {config.output}")
    write_overview(nodes, edges, config.output_overview)
    logger.success(f"Wrote {config.output_overview}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Interactive Altair network of the coherence redundancy clusters")
    parser.add_argument("--combined", type=Path, required=True, help="coherence_metrics_combined.parquet")
    parser.add_argument("--deduplicated", type=Path, required=True, help="coherence_terms_deduplicated.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output HTML page: one cluster at a time")
    parser.add_argument("--output-overview", type=Path, required=True,
                        help="Output HTML page: every cluster at once (pyvis)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the pages, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            combined=args.combined,
            deduplicated=args.deduplicated,
            output=args.output,
            output_overview=args.output_overview,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
