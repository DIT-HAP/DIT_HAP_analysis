#!/usr/bin/env python3

"""
Coherence Interactive Scatter — an Altair HTML explorer
=======================================================

The static group_scatter.pdf answers "where do this term's genes sit in DR-DL
space?" for a namelist fixed at config time. Answering it for a term you did not
think of means editing config and re-running the pipeline, and the static figure
cannot say anything about a gene you point at.

This is the same picture as an interactive document, plus the maps that make the
term pickable: three panels side by side.

- **Terms in q–z space** — one marker per group at its coherent-tail q (the x
  axis, pinned 0–1 because q is a probability) against its median pairwise
  distance z, with the marker's area carrying the scored-member count. This is
  the plane the cohort labels are cut on, so the two dashed rules ARE the cuts,
  and it is the interactive counterpart of the static figure's "Coherence vs
  significance" panel.
- **Group centroid positions** — one marker per group, at the mean normalized
  position of its scored members. The selected group is a large solid marker, the
  rest fade back. This is where a reader browses: the axis is the same fitness
  space as the right panel, so an outlying group is visibly outlying.
- **Genes of the selected term in DR–DL space** — the genome-wide cloud with the
  selected term's genes on top. Hovering any point — cloud or highlighted — shows
  that gene's fitness record.

One selection drives all three, and it is set two ways: the dropdown above the
panels, or a click on a marker in the left or middle panel. They are the same
selection, so either route updates every panel — the dropdown's value follows a
click, and both marker highlights follow the dropdown.

A second dropdown next to the first narrows all panels to one cohort end:
`all` (the default, today's view), `coherent` (q <= 0.05 AND z < --coherent-z)
or `incoherent` (z > --incoherent-z) — the same cuts the cohort figures and
workbooks label with. On both term panels the unselected markers carry the
cohort on their colour, so the two ends are findable on sight before any
picking; `other` keeps the original grey.

The term list is every group in the source's metrics table, most coherent first,
so the page opens on the strongest signal. Genes are the term's *scored* members
(the DR<threshold set the coherence z-score was computed on), joined to the
upstream fitting table for the per-gene detail shown on hover.

One script serves all three views. `--source` picks one source's rows; leaving it
out pools everything in `--annotation`, which is how the `combined/` and `dedup/`
folders get the same page over their own metrics table. The pooled run prefixes
each term with its source and joins members on (source, group_id) — a name is only
unique within a source, and 173 group_ids are shared by two of them.

Input
-----
- A metrics table: one source's coherence_metrics.parquet, combined/'s pooled
  table, or dedup/'s coherence_terms_representatives.tsv (source, group_id,
  group_name, n_scored_members, median_pairwise_distance_z, q_value_coherent,
  q_value_incoherent).
- group_annotation_long.tsv: the matching prepared group -> member long table(s);
  one path per source, concatenated.
- fitting_results.tsv: the upstream per-gene fitness table (see coherence/io.py).

Output
------
- interactive_scatter.html: one self-contained page on a fixed 1280x720 stage —
  16:9, the size of a PowerPoint widescreen slide at 96 dpi — so a screenshot of
  the page pastes into a slide without cropping or letterboxing. It loads
  vega-embed from a CDN, so opening it needs a network connection the first time.

Usage
-----
    python plot_interactive_scatter.py \\
        --metrics results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --fitting-results .../gene_level/fitting_results.tsv \\
        --source go_macrocomplex \\
        --output tmp/{source}/interactive_scatter.html

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-09-22 (cohort dropdown + marker colours 2026-10-08;
          q–z term panel 2026-10-22)
Version:  2.3.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import altair as alt
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.io import load_fitting_results, load_long_table
from figures import FURNITURE_COLOR, house_colors
from io_table import read_file
from logging_setup import setup_logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# The genome cloud is context, not data: drawn at the same marker size as the
# highlighted genes — a smaller one made the cloud read as a different kind of
# object rather than as the same measurement in context — but in a light grey and
# translucent enough that a few thousand overlapping markers stay a background.
_BACKGROUND_LIGHT_GREY = "#c9c9c9"
_BACKGROUND_OPACITY = 0.35

# The highlighted genes are the subject of the figure, and a term is 3-300 genes,
# so they can afford to be large. The cloud matches this size.
_MEMBER_SIZE = 70
_MEMBER_OPACITY = 0.85
_BACKGROUND_SIZE = _MEMBER_SIZE

# Floor for the DL axis. Most genes sit at exactly DL = 0, so an axis ending at 0
# draws half of those markers outside the frame.
_DL_FLOOR = -0.05

# BOTH panels are pinned to ONE fixed range on each axis — the data's own global
# extent (norm_DR = DR; norm_DL = DL/10), padded 5% so a marker at the very edge
# is not clipped. Fixed, not autoscaled per panel: autoscale lets the centroid
# panel and the gene panel disagree, and a pan/zoom on one silently rescales only
# that panel, which is exactly the comparison the page is for.
_DR_LIMITS = (-2.14, 0.22)
_DL_LIMITS = (_DL_FLOOR, 1.0)

# The scored/unscored split the coherence computation draws on, shown in the gene
# panel: the genome cloud stays one grey, and a dashed vertical line marks the DR
# threshold — but the highlighted members are recoloured on either side of it,
# since only the term's own genes are asked to read against the cut. Threshold
# arrives via --dr-threshold; this fallback must match config/analysis.yaml's
# `coherence.dr_threshold` and the two static scatter scripts' defaults.
_DR_THRESHOLD = -0.3
# Scored members: the house red (the caller's `highlight`). Not-scored members:
# the house blue-teal, a different hue (not a lightness step) so the two stay
# apart in greyscale print as well.
_UNSCORED_MEMBER_COLOR = house_colors((1,))[0]

# Chart footprint in pixels, per panel. The page itself is a fixed 16:9 stage
# (1280x720 is a PowerPoint widescreen slide at 96 dpi; see _PAGE_CSS), and the
# figure is centred inside it. The two are independent knobs: at 366x522 the three
# panels fill the stage edge to edge, and a smaller panel leaves the figure
# floating on the same slide-sized canvas. Below roughly 250 px wide the axis
# titles and the size legend are wider than the plot, so they start to crowd the
# panel next to them — measure before going smaller.
_CHART_WIDTH = 200
_CHART_HEIGHT = 200

# The centroid panel is a map, not a measurement: one small marker per group and
# the selected group swapped for a larger solid one. The unselected mass stays
# translucent — dimming it is what makes the selected marker findable in a
# 1,400-term source without hiding where the rest of the terms sit. The selected
# size is tuned against the 200x200 panel: at 280 px² it covered a tenth of the
# panel and read as a blob rather than a marker.
_CENTROID_SIZE = 60
_CENTROID_SELECTED_SIZE = 140
_CENTROID_OPACITY = 0.4
_CENTROID_STROKE = "#ffffff"

# Cohort colouring for the term markers on both term panels, so a reader can find
# the coherent and incoherent terms on sight instead of reading z off a tooltip.
# Both ends get a greyscale-separated house colour (positions 2 and 3, the pair
# figures.py curates for that property) — position 0 is the highlight red and is
# therefore already taken. Everything else (the bulk) recedes into the grey the
# panels already use. The same vocabulary the cohort workbooks use
# (export_cohorts.py's Cohort), cut with the same thresholds the figures label
# with (--coherent-z/--incoherent-z, below).
_COHORT_COLORS = dict(
    zip(("coherent", "incoherent"), house_colors((2, 3)))
)

# The coherent cohort's FDR half. The z thresholds arrive through --label-*
# style CLI flags (config/analysis.yaml's `coherence:` block drives them via the
# rules); the q cutoff is the one constant this page shares with the workbooks.
_COHERENT_Q_MAX = 0.05

# The default z cut, used when no caller passes one. The CLI always passes the
# config's `coherence.coherent_z_threshold`; cohort_labels() labels with it and the
# q-z panel draws its dashed rule at it, so the two cannot disagree.
_DEFAULT_COHERENT_Z = -2.0

# --- q–z term panel ---------------------------------------------------------
# q is a probability, so its axis is pinned to the full 0-1 instead of autoscaled
# to the observed range: the significance cut then sits at its true position, and
# pages over different row sets (a source, combined, dedup) stay comparable. This
# is the same convention the static "Coherence vs significance" panel uses.
_Q_LIMITS = (0.0, 1.0)

# Marker AREA — not radius — carries the scored-member count, so a term's ink
# tracks how many genes it was scored on at all (a square-rooted radius would not).
# Log-compressed like the static panel's point_sizes(): counts span 3-207, and a
# linear area mapping lets the biggest terms cover the cloud they are meant to be
# read against, so equal area steps are multiplicative in the count — read the
# legend, not the ratio. The floor keeps a 3-member term clickable and the ceiling
# is tuned to the 200x200 panel: 600 px² put a 24 px disc on a 200 px-wide plot.
_STATS_SIZE_RANGE = (12, 200)

# Coordinates are rounded before they are embedded: the page carries one JSON row
# per (term, member) pair — 58k of them for go_bp — and a full float64 prints
# twice as long as it needs to. 4 dp is below what the tooltips display (3 dp), so
# the rounding is invisible.
_ROUND_DP = 4
_ROUND_COLUMNS = ["DR", "DL", "norm_DR", "norm_DL", "R2"]

_REQUIRED_METRIC_COLUMNS = [
    "source", "group_id", "group_name", "median_pairwise_distance_z",
    "q_value_coherent", "q_value_incoherent",
]

# The cohort dropdown's choices, in order. `all` keeps the page's default view
# unchanged; the two cohort entries narrow the map and the term list to that end.
_COHORT_CHOICES = ["all", "coherent", "incoherent"]

# Per-gene columns carried into the tooltip, in display order. Every one comes
# from the upstream fitting table except the two normalized coordinates.
_GENE_TOOLTIP = [
    alt.Tooltip("Name:N", title="Gene"),
    alt.Tooltip("Systematic ID:N", title="Systematic ID"),
    alt.Tooltip("DR:Q", title="DR", format=".3f"),
    alt.Tooltip("DL:Q", title="DL", format=".3f"),
    alt.Tooltip("norm_DR:Q", title="norm DR", format=".3f"),
    alt.Tooltip("norm_DL:Q", title="norm DL/10", format=".3f"),
    alt.Tooltip("R2:Q", title="R²", format=".3f"),
    alt.Tooltip("FYPOviability:N", title="Viability"),
]

# What a group marker says about the group, not about any one gene.
_TERM_TOOLTIP = [
    alt.Tooltip("term:N", title="Term"),
    alt.Tooltip("cohort:N", title="Cohort"),
    alt.Tooltip("median_pairwise_distance_z:Q", title="z", format=".2f"),
    alt.Tooltip("q_value_coherent:Q", title="q coherent (BH)", format=".3g"),
    alt.Tooltip("q_value_incoherent:Q", title="q incoherent (BH)", format=".3g"),
    alt.Tooltip("n_scored_members:Q", title="scored genes", format="d"),
]

# vega-embed renders the bound <select> with browser defaults at the END of the
# chart container, which reads as a form control tacked onto the bottom of the
# figure. These rules restyle it as a header bar and lift it above the panels:
# `.vega-bindings` is a sibling of the chart inside `.chart-wrapper`, so a flex
# column plus `order: -1` moves it without touching the DOM. The class names are
# vega-embed's, so a future version that renames them degrades to the default
# widget in the default position instead of breaking the page.
_PAGE_CSS = """
<style>
  /* The page is a fixed 16:9 stage the size of a PowerPoint slide (1280x720 at
     96 dpi), so a screenshot — or a browser window at that size — pastes into a
     slide without cropping or letterboxing. The stage is pinned to the TOP of the
     window and the figure to the top of the stage: a panel size that does not
     fill the stage leaves its empty space BELOW the figure, never a band above
     it. Horizontal centring is all that is left to do. */
  html, body { min-height: 100vh; }
  body {
    margin: 0;
    padding: 0;
    background: #ffffff;
    display: flex;
    align-items: flex-start;
  }
  /* `body` outranks vega-embed's own `#vis.vega-embed { width: 100% }`: the stage
     has to stay 1280 px whatever the window is, or the figure would be squeezed
     against the actions rail and clipped on the left. `margin: 0 auto` centres it
     horizontally and, unlike body's justify, collapses to 0 on a narrow window
     instead of clipping the left edge. */
  body #vis.vega-embed {
    width: 1280px;
    height: 720px;
    padding: 0;
    display: block;
    margin: 0 auto;
  }
  #vis > .chart-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
  }
  /* vega-embed's actions <details> is 28 px tall and would otherwise sit after
     the figure, pushing the stage 28 px past the viewport and adding a scrollbar
     to a screenshot that is supposed to be exactly 1280x720. */
  #vis.vega-embed > details { position: absolute; top: 0; right: 0; }
  #vis > .chart-wrapper > .vega-bindings { order: -1; }
  .vega-bindings {
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 10px;
    width: fit-content;
    margin: 0 auto 12px auto;
    padding: 8px 16px;
    background: #f6f7f9;
    border: 1px solid #e3e5e9;
    border-radius: 9px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    font-size: 13px;
  }
  .vega-bind { display: flex; align-items: center; gap: 8px; margin: 0; }
  .vega-bind-name { font-weight: 600; color: #2b2f36; letter-spacing: 0.01em; }
  .vega-bindings select {
    font: inherit;
    color: #22262c;
    background: #ffffff;
    padding: 5px 10px;
    border: 1px solid #c9cdd6;
    border-radius: 6px;
    width: auto;
    max-width: 420px;
    cursor: pointer;
  }
  .vega-bindings select:focus {
    outline: 2px solid #4c78a8;
    outline-offset: 1px;
    border-color: #4c78a8;
  }
</style>
"""


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Inputs and output for the interactive scatter page."""
    metrics: Path
    # One path per source: a per-source run passes one, the pooled `combined/` and
    # `dedup/` runs pass every registered source's table. They are concatenated,
    # which is what makes the pooled call correct — `source` is what keeps a
    # group_id shared by two sources (173 of them) from resolving to the union of
    # both member sets.
    annotations: tuple[Path, ...]
    fitting_results: Path
    # None pools every row of `metrics` (the combined / dedup views).
    source: str | None
    # The cohort cuts, matching config/analysis.yaml's `coherence:` block — the
    # same two numbers the static figures and workbooks label with.
    coherent_z: float
    incoherent_z: float
    # The scored/unscored cut the gene panel's cloud is split and dashed-lined on.
    dr_threshold: float
    output: Path

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then create the output dir."""
        for path in [self.metrics, *self.annotations, self.fitting_results]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def cohort_labels(
    metrics: pd.DataFrame, coherent_z: float, incoherent_z: float
) -> pd.DataFrame:
    """Add the `cohort` column the dropdown filters and the markers colour on.

    Same cuts the figures and workbooks label with (config/analysis.yaml's
    `coherence:` block): coherent = q <= 0.05 AND z < coherent_z; incoherent =
    z > incoherent_z, a descriptive rank cut that carries no FDR (see
    export_cohorts.py's docstring). Everything else is `other`.
    """
    z = metrics["median_pairwise_distance_z"]
    coherent = (metrics["q_value_coherent"] <= _COHERENT_Q_MAX) & (z < coherent_z)
    incoherent = z > incoherent_z
    return metrics.assign(
        cohort=np.where(coherent, "coherent", np.where(incoherent, "incoherent", "other"))
    )


def term_labels(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add the `term` column every panel and the picker key on.

    A name has to identify exactly one group, or the picker selects several at once
    and the centroid panel draws one marker where there are many. Two collisions
    break that, and neither is hypothetical:

    - `group_name` repeats *within* a source — BRITE's generic node labels
      ('Others', 'Large subunit') recur across trees, so 25 names in kegg_brite
      carry 2-7 groups each;
    - once the sources are pooled, the same name appears in two databases.

    So a name shared inside a source is suffixed with its group_id, and a pooled
    table prefixes every term with its source.
    """
    name = metrics["group_name"].astype(str)
    if metrics["source"].nunique() > 1:
        name = metrics["source"].astype(str) + ": " + name
    shared = metrics.groupby(["source", "group_name"])["group_id"].transform("nunique") > 1
    return metrics.assign(
        term=name.where(~shared, name + " (" + metrics["group_id"] + ")")
    )


def gene_detail(fitting: pd.DataFrame) -> pd.DataFrame:
    """The per-gene columns the page shows — the tooltip set, rounded for embedding.

    One frame serves two roles: it is the genome cloud, and it is the lookup table
    the member rows resolve their coordinates in. Passing the same object to both
    is what lets Altair embed it once instead of once per referrer.
    """
    columns = ["Systematic ID", "Name", "FYPOviability", "DR", "DL", "norm_DR", "norm_DL", "R2"]
    detail = fitting[columns].copy()
    detail[_ROUND_COLUMNS] = detail[_ROUND_COLUMNS].round(_ROUND_DP)
    return detail


def member_table(
    metrics: pd.DataFrame, long_table: pd.DataFrame, detail: pd.DataFrame, source: str | None
) -> pd.DataFrame:
    """The term -> gene rows to draw: one row per (group, scored member).

    Deliberately just the two keys. A flat row carrying the gene's own columns
    repeats every field name once per member — 58k times for go_bp — and those
    names, not the digits, are the bottleneck: they are what made the page 24 MB.
    The gene columns are fetched in the browser by `transform_lookup` against the
    cloud's own dataset, which is embedded once.

    `source` is None for the pooled views, where the table holds every source.
    """
    # Inner-join semantics against `Systematic ID`, so a member the fitting table
    # could not fit is dropped rather than drawn at a fabricated coordinate — the
    # same members the coherence z-score was computed on.
    fitted = set(detail["Systematic ID"])
    rows = long_table if source is None else long_table[long_table["source"] == source]
    members = rows[["source", "group_id", "Systematic ID"]]
    members = members[members["Systematic ID"].isin(fitted)].copy()
    # Both columns, not group_id alone: on the pooled views a group_id shared by two
    # sources would otherwise merge two different groups' members. Neither key is
    # kept in the returned frame — the page only ever reads `term` and looks the
    # gene up by `Systematic ID`, and 58k repeated short strings is not nothing.
    # `cohort` rides along because the cohort dropdown filters this layer too —
    # without it an undefined datum would empty the gene panel on every pick.
    merged = members.merge(
        metrics[["source", "group_id", "term", "cohort"]], on=["source", "group_id"], how="inner"
    )
    return merged[["term", "cohort", "Systematic ID"]]


def centroid_table(
    members: pd.DataFrame, detail: pd.DataFrame, metrics: pd.DataFrame
) -> pd.DataFrame:
    """One row per term: where its members sit, on average, in the same space."""
    coords = members.merge(detail[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID")
    means = coords.groupby("term", as_index=False)[["norm_DR", "norm_DL"]].mean()
    stats = metrics[[
        "term", "median_pairwise_distance_z", "q_value_coherent",
        "q_value_incoherent", "n_scored_members", "cohort",
    ]]
    return means.merge(stats, on="term", how="left")


def term_options(centroids: pd.DataFrame) -> list[str]:
    """Dropdown entries: every term, most coherent first."""
    ordered = centroids.sort_values("median_pairwise_distance_z")
    return [str(term) for term in ordered["term"].dropna().unique()]


def padded_extent(values: pd.Series, pad: float = 0.05) -> tuple[float, float]:
    """The column's extent with a proportional pad, so an edge marker is not clipped."""
    lo, hi = float(values.min()), float(values.max())
    margin = pad * (hi - lo) if hi > lo else 1.0
    return lo - margin, hi + margin


def member_extent(centroids: pd.DataFrame) -> tuple[float, float]:
    """The scored-member domain for the size scale; a degenerate page gets a width."""
    lo = float(centroids["n_scored_members"].min())
    hi = float(centroids["n_scored_members"].max())
    return lo, hi if hi > lo else lo + 1.0


def build_chart(
    background: pd.DataFrame, members: pd.DataFrame, centroids: pd.DataFrame,
    options: list[str], cohort_colors: dict[str, str] | None = None,
    dr_threshold: float = _DR_THRESHOLD, coherent_z: float = _DEFAULT_COHERENT_Z,
) -> alt.Chart:
    """The three-panel page: q–z left, the term map middle, the picked genes right."""
    if cohort_colors is None:
        cohort_colors = _COHORT_COLORS
    # ONE selection, two ways to set it. `bind` renders the dropdown; `on="click"`
    # additionally lets a click on a term marker write the same value, so the
    # dropdown and the two term panels can never disagree. The param is declared
    # on those panels' own point layers rather than on the page, which is what
    # scopes the click listener to them — declared a level up, a click on the gene
    # cloud (whose rows carry no `term`) would clear the selection.
    picker = alt.selection_point(
        fields=["term"],
        bind=alt.binding_select(options=options, name="Term  "),
        value=options[0],
        on="click",
        clear=False,
        name="term_pick",
    )
    # A second dropdown narrows every panel to one cohort end, so the coherent /
    # incoherent terms are findable without scrolling a 1,400-entry term list.
    # A plain variable param (`alt.param`), not a point selection: it has no
    # click semantics, and `all` — its initial value — simply disables the
    # filter, which keeps today's default view unchanged.
    cohort_pick = alt.param(
        bind=alt.binding_select(options=_COHORT_CHOICES, name="Cohort  "),
        value=_COHORT_CHOICES[0],
        name="cohort_pick",
    )
    # A plain string, not an `&`-composition: mixing a selection with an expression
    # in one PredicateComposition raises SchemaValidationError at compile time, while
    # separate .transform_filter() calls AND together cleanly.
    cohort_filter = "cohort_pick === 'all' || datum.cohort === cohort_pick"
    highlight = house_colors((0,))[0]
    # One cohort palette for every panel that marks a term: the two flagged ends on
    # their house colours, the bulk on the grey the maps already use.
    cohort_color = alt.Color(
        "cohort:N",
        scale=alt.Scale(
            domain=["coherent", "incoherent", "other"],
            range=[cohort_colors["coherent"], cohort_colors["incoherent"],
                   _BACKGROUND_LIGHT_GREY],
        ),
        legend=None,
    )
    # One fixed domain per axis, shared by the two fitness panels.
    x_scale = alt.Scale(domain=list(_DR_LIMITS))
    y_scale = alt.Scale(domain=list(_DL_LIMITS))
    gene_x = alt.X("norm_DR:Q", title="norm DR", scale=x_scale)
    gene_y = alt.Y("norm_DL:Q", title="norm DL/10", scale=y_scale)
    centroid_x = alt.X("norm_DR:Q", title="centroid norm DR", scale=x_scale)
    centroid_y = alt.Y("norm_DL:Q", title="centroid norm DL/10", scale=y_scale)

    # --- term q-z panel -----------------------------------------------------
    # The two numbers the cohort labels are cut on, on one plane: q on x and z on
    # y. Both arrive on the centroid rows, so this panel costs no extra data — one
    # frame feeds both term panels, and the whole metrics row is also what lets a
    # click here write the same `term` the other panels read.
    z_limits = padded_extent(centroids["median_pairwise_distance_z"])
    stats_x = alt.X("q_value_coherent:Q", title="q coherent (BH)",
                    scale=alt.Scale(domain=list(_Q_LIMITS)))
    stats_y = alt.Y("median_pairwise_distance_z:Q", title="median pairwise distance z",
                    scale=alt.Scale(domain=list(z_limits)))
    size_scale = alt.Scale(
        type="log", domain=list(member_extent(centroids)), range=list(_STATS_SIZE_RANGE)
    )
    # The legend rides on the dimmed layer only: the picked layer needs the same
    # size encoding, or the selected marker would jump to a different scale, but a
    # second copy of the legend would just be a duplicate. It is pinned to the
    # bottom because Vega-Lite hoists a concat's legends to the far RIGHT of the
    # whole page by default, where a size key for the left panel reads as if it
    # belonged to the gene panel.
    stats_size = alt.Size("n_scored_members:Q", scale=size_scale,
                          legend=alt.Legend(title="scored genes", orient="bottom"))
    stats_size_silent = alt.Size("n_scored_members:Q", scale=size_scale, legend=None)
    # The cohort cuts, drawn behind the markers in the same dashed furniture the
    # gene panel uses for the DR threshold: the bottom-left corner — q at or below
    # 0.05, z below the coherent cut — then reads as "coherent" without a legend.
    q_cut = (
        alt.Chart(pd.DataFrame({"q_max": [_COHERENT_Q_MAX]}))
        .mark_rule(strokeDash=[4, 3], color=FURNITURE_COLOR)
        .encode(x=alt.X("q_max:Q", scale=alt.Scale(domain=list(_Q_LIMITS))))
    )
    z_cut = (
        alt.Chart(pd.DataFrame({"coherent_z": [coherent_z]}))
        .mark_rule(strokeDash=[4, 3], color=FURNITURE_COLOR)
        .encode(y=alt.Y("coherent_z:Q", scale=alt.Scale(domain=list(z_limits))))
    )

    def stats_layer(picked: bool, sample: alt.Chart) -> alt.Chart:
        # Same split as the centroid map: the picked marker in a layer of its own,
        # drawn last, so it is never buried under the 1,400-term cloud.
        return sample.mark_circle(
            stroke=_CENTROID_STROKE,
            strokeWidth=0.6,
            opacity=1.0 if picked else _CENTROID_OPACITY,
        ).encode(
            x=stats_x,
            y=stats_y,
            size=stats_size_silent if picked else stats_size,
            tooltip=_TERM_TOOLTIP,
            color=alt.value(highlight) if picked else cohort_color,
        )

    # The params are attached to the dimmed point layer rather than to the layer
    # composition: Altair hands a param added to a LayerChart to its FIRST layer,
    # and here the first two are the cut rules — whose rows carry no `term`, so the
    # click listener would be scoped to a layer nothing can be selected on. Added
    # to the point layer, it lands where the markers are.
    dimmed_stats = stats_layer(
        False,
        alt.Chart(centroids).transform_filter(~picker).transform_filter(cohort_filter),
    ).add_params(picker, cohort_pick)
    picked_stats = stats_layer(
        True,
        alt.Chart(centroids).transform_filter(picker).transform_filter(cohort_filter),
    )
    # No .interactive() here either: a pan gesture would fight the term click.
    stats = (
        (q_cut + z_cut + dimmed_stats + picked_stats)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title=alt.TitleParams(
                text="Terms — q coherent vs z",
                subtitle="size = scored genes; click a term",
            ),
        )
    )

    # Two layers, not one conditional mark: a source with 1,400 terms draws a
    # crowded cloud, and in a single layer the selected marker is drawn wherever its
    # row happens to fall — underneath whichever markers come after it. Splitting on
    # the selection puts the picked marker in a layer of its own, drawn last.
    #
    # Unselected markers are NOT one flat grey any more: the cohort end the term
    # belongs to is carried on the marker colour, so the map answers "where are
    # the coherent terms?" without a click. `other` keeps the old grey and the
    # opacity, so the idle page reads exactly as it did before.
    def centroid_layer(picked: bool, sample: alt.Chart) -> alt.Chart:
        # The picked layer stays `highlight`; the dimmed layer carries the cohort
        # on its colour. `other` keeps the old grey, so the idle page reads as
        # before apart from the two flagged ends.
        color = alt.value(highlight) if picked else cohort_color
        return (
            sample.mark_circle(
                stroke=_CENTROID_STROKE,
                strokeWidth=0.6,
                size=_CENTROID_SELECTED_SIZE if picked else _CENTROID_SIZE,
                opacity=1.0 if picked else _CENTROID_OPACITY,
            )
            .encode(x=centroid_x, y=centroid_y, tooltip=_TERM_TOOLTIP, color=color)
        )

    dimmed = centroid_layer(
        False,
        alt.Chart(centroids).transform_filter(~picker).transform_filter(cohort_filter),
    )
    picked = centroid_layer(
        True,
        alt.Chart(centroids).transform_filter(picker).transform_filter(cohort_filter),
    )
    centroid = (
        (dimmed + picked)
        .add_params(picker, cohort_pick)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title="Group centroids — click a term",
        )
    )
    # No .interactive() here on purpose: a pan gesture on this panel would fight
    # the click that selects a term.
    # The gene panel keeps a single grey cloud with a dashed rule marking the DR
    # threshold; the highlighted members carry the split themselves — scored ones
    # in red, not-scored ones in the blue-teal — via an if/else condition colour.
    cloud = (
        alt.Chart(background)
        .mark_circle(size=_BACKGROUND_SIZE, color=_BACKGROUND_LIGHT_GREY, opacity=_BACKGROUND_OPACITY)
        .encode(x=gene_x, y=gene_y, tooltip=_GENE_TOOLTIP)
    )
    threshold_rule = (
        alt.Chart(pd.DataFrame({"dr_threshold": [dr_threshold]}))
        .mark_rule(strokeDash=[4, 3], color=FURNITURE_COLOR)
        .encode(x=alt.X("dr_threshold:Q", scale=x_scale))
    )
    highlighted = (
        alt.Chart(members)
        .transform_lookup(
            lookup="Systematic ID",
            from_=alt.LookupData(
                data=background,
                key="Systematic ID",
                fields=["Name", "FYPOviability", "DR", "DL", "norm_DR", "norm_DL", "R2"],
            ),
        )
        .transform_calculate(
            unscored=alt.datum.norm_DR >= dr_threshold,
        )
        .transform_filter(picker)
        .transform_filter(cohort_filter)
        .mark_circle(size=_MEMBER_SIZE, opacity=_MEMBER_OPACITY)
        .encode(
            x=gene_x, y=gene_y, tooltip=_GENE_TOOLTIP,
            color=alt.condition(
                alt.datum.unscored,
                alt.value(_UNSCORED_MEMBER_COLOR),
                alt.value(highlight),
            ),
        )
    )
    genes = (
        (cloud + threshold_rule + highlighted)
        .add_params(cohort_pick)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title="Selected term's genes — hover for record",
        )
        .interactive()
    )
    return stats | centroid | genes


def embed_page_css(path: Path) -> None:
    """Splice the binding-bar styling into the page Altair just wrote."""
    html = path.read_text(encoding="utf-8")
    path.write_text(html.replace("</head>", _PAGE_CSS + "</head>", 1), encoding="utf-8")


@logger.catch(reraise=True)
def run(config: PlotConfig) -> None:
    """Load -> build the terms and gene rows -> save the interactive page."""
    config.validate()
    # read_file, not read_parquet: `dedup/` hands this script a TSV while the other
    # two views hand it Parquet, and the extension is what says which.
    metrics = read_file(config.metrics)
    missing = [column for column in _REQUIRED_METRIC_COLUMNS if column not in metrics.columns]
    if missing:
        raise ValueError(f"metrics table missing required column(s) {missing} (have: {list(metrics.columns)})")

    long_table = pd.concat(
        [load_long_table(path) for path in config.annotations], ignore_index=True
    )
    fitting = load_fitting_results(config.fitting_results)

    metrics = term_labels(metrics)
    metrics = cohort_labels(metrics, config.coherent_z, config.incoherent_z)
    detail = gene_detail(fitting)
    members = member_table(metrics, long_table, detail, config.source)
    centroids = centroid_table(members, detail, metrics)
    options = term_options(centroids)
    if not options:
        raise ValueError(f"no groups for source '{config.source}' in {config.metrics}")

    scope = config.source if config.source else f"pooled, {metrics['source'].nunique()} sources"
    logger.info(
        f"{len(metrics):,} groups ({scope}), {len(members):,} member rows, "
        f"{len(fitting):,} background genes; opening on '{options[0]}'"
    )
    chart = build_chart(detail, members, centroids, options, _COHORT_COLORS,
                        config.dr_threshold, config.coherent_z)
    chart.save(config.output)
    embed_page_css(config.output)
    logger.success(f"Wrote {config.output} ({len(members):,} member rows over {len(options):,} terms)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Interactive Altair scatter of a source's coherence terms")
    parser.add_argument("--metrics", type=Path, required=True,
                        help="coherence_metrics.parquet, or dedup/coherence_terms_representatives.tsv")
    parser.add_argument("--annotation", type=Path, nargs="+", required=True,
                        help="The matching group_annotation_long.tsv table(s) — one per source, concatenated")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv")
    parser.add_argument("--source", type=str,
                        help="One source's rows; omit to pool everything in --annotation (combined / dedup)")
    parser.add_argument("--coherent-z", type=float, default=_DEFAULT_COHERENT_Z,
                        help="Coherent cohort: z below this, and the q-z panel's horizontal cut")
    parser.add_argument("--incoherent-z", type=float, default=0.5, help="Incoherent cohort: z above this")
    parser.add_argument("--dr-threshold", type=float, default=_DR_THRESHOLD,
                        help="The coherence scored/unscored cut, shown in the gene panel")
    parser.add_argument("--output", type=Path, required=True, help="Output HTML page")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the page, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotConfig(
            metrics=args.metrics,
            annotations=tuple(args.annotation),
            fitting_results=args.fitting_results,
            source=args.source,
            coherent_z=args.coherent_z,
            incoherent_z=args.incoherent_z,
            dr_threshold=args.dr_threshold,
            output=args.output,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
