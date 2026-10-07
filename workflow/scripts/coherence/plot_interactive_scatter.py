#!/usr/bin/env python3

"""
Coherence Interactive Scatter — an Altair HTML explorer
=======================================================

The static group_scatter.pdf answers "where do this term's genes sit in DR-DL
space?" for a namelist fixed at config time. Answering it for a term you did not
think of means editing config and re-running the pipeline, and the static figure
cannot say anything about a gene you point at.

This is the same picture as an interactive document, plus the map that makes the
term pickable: two panels side by side.

- **Group centroid positions** — one marker per group, at the mean normalized
  position of its scored members. The selected group is a large solid marker, the
  rest fade back. This is where a reader browses: the axis is the same fitness
  space as the right panel, so an outlying group is visibly outlying.
- **Genes of the selected term in DR–DL space** — the genome-wide cloud with the
  selected term's genes on top. Hovering any point — cloud or highlighted — shows
  that gene's fitness record.

One selection drives both, and it is set two ways: the dropdown above the panels,
or a click on a group marker. They are the same selection, so either route
updates both panels — the dropdown's value follows a click, and the marker
highlight follows the dropdown.

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
  group_name, n_scored_members, median_pairwise_distance_z, q_value).
- group_annotation_long.tsv: the matching prepared group -> member long table(s);
  one path per source, concatenated.
- fitting_results.tsv: the upstream per-gene fitness table (see coherence/io.py).

Output
------
- interactive_scatter.html: one self-contained page. It loads vega-embed from a
  CDN, so opening it needs a network connection the first time.

Usage
-----
    python plot_interactive_scatter.py \\
        --metrics results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --fitting-results .../gene_level/fitting_results.tsv \\
        --source go_macrocomplex \\
        --output tmp/{source}/interactive_scatter.html

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-09-22
Version:  2.1.0
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.io import load_fitting_results, load_long_table  # noqa: E402
from figures import house_colors  # noqa: E402
from io_table import read_file  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


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
# draws half of those markers outside the frame; DR has no such pile-up and keeps
# its own limits.
_DL_FLOOR = -0.05

# Chart footprint in pixels, per panel. Two panels share the page, so each is
# smaller than the single-panel version was: side by side they stay inside a
# laptop viewport, which is the whole point of putting them there.
_CHART_WIDTH = 430
_CHART_HEIGHT = 340

# The centroid panel is a map, not a measurement: one small marker per group and
# the selected group swapped for a large solid one. The unselected mass stays
# translucent — dimming it is what makes the selected marker findable in a
# 1,400-term source without hiding where the rest of the terms sit.
_CENTROID_SIZE = 60
_CENTROID_SELECTED_SIZE = 280
_CENTROID_OPACITY = 0.4
_CENTROID_STROKE = "#ffffff"

# Coordinates are rounded before they are embedded: the page carries one JSON row
# per (term, member) pair — 58k of them for go_bp — and a full float64 prints
# twice as long as it needs to. 4 dp is below what the tooltips display (3 dp), so
# the rounding is invisible.
_ROUND_DP = 4
_ROUND_COLUMNS = ["DR", "DL", "norm_DR", "norm_DL", "R2"]

_REQUIRED_METRIC_COLUMNS = ["source", "group_id", "group_name", "median_pairwise_distance_z", "q_value"]

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
    alt.Tooltip("median_pairwise_distance_z:Q", title="z", format=".2f"),
    alt.Tooltip("q_value:Q", title="q (BH)", format=".3g"),
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
  body { margin: 0; padding: 18px; background: #ffffff; }
  #vis > .chart-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
  }
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
def term_labels(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add the `term` column both panels and the picker key on.

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
    merged = members.merge(
        metrics[["source", "group_id", "term"]], on=["source", "group_id"], how="inner"
    )
    return merged[["term", "Systematic ID"]]


def centroid_table(
    members: pd.DataFrame, detail: pd.DataFrame, metrics: pd.DataFrame
) -> pd.DataFrame:
    """One row per term: where its members sit, on average, in the same space."""
    coords = members.merge(detail[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID")
    means = coords.groupby("term", as_index=False)[["norm_DR", "norm_DL"]].mean()
    stats = metrics[[
        "term", "median_pairwise_distance_z", "q_value", "n_scored_members",
    ]]
    return means.merge(stats, on="term", how="left")


def term_options(centroids: pd.DataFrame) -> list[str]:
    """Dropdown entries: every term, most coherent first."""
    ordered = centroids.sort_values("median_pairwise_distance_z")
    return [str(term) for term in ordered["term"].dropna().unique()]


def build_chart(
    background: pd.DataFrame, members: pd.DataFrame, centroids: pd.DataFrame, options: list[str]
) -> alt.Chart:
    """The two-panel page: the term map left, the picked term's genes right."""
    # ONE selection, two ways to set it. `bind` renders the dropdown; `on="click"`
    # additionally lets a click on a centroid marker write the same value, so the
    # dropdown and the map can never disagree. The param is declared on the
    # centroid chart rather than on the page, which is what scopes the click
    # listener to that panel — declared a level up, a click on the gene cloud
    # (whose rows carry no `term`) would clear the selection.
    picker = alt.selection_point(
        fields=["term"],
        bind=alt.binding_select(options=options, name="Term  "),
        value=options[0],
        on="click",
        clear=False,
        name="term_pick",
    )
    highlight = house_colors((0,))[0]
    gene_x = alt.X("norm_DR:Q", title="norm DR")
    gene_y = alt.Y("norm_DL:Q", title="norm DL/10", scale=alt.Scale(domainMin=_DL_FLOOR))
    centroid_x = alt.X("norm_DR:Q", title="centroid norm DR")
    centroid_y = alt.Y("norm_DL:Q", title="centroid norm DL/10", scale=alt.Scale(domainMin=_DL_FLOOR))

    # Two layers, not one conditional mark: a source with 1,400 terms draws a
    # crowded cloud, and in a single layer the selected marker is drawn wherever its
    # row happens to fall — underneath whichever markers come after it. Splitting on
    # the selection puts the picked marker in a layer of its own, drawn last.
    def centroid_layer(picked: bool, sample: alt.Chart) -> alt.Chart:
        return (
            sample.mark_circle(
                stroke=_CENTROID_STROKE,
                strokeWidth=0.6,
                size=_CENTROID_SELECTED_SIZE if picked else _CENTROID_SIZE,
                color=highlight if picked else _BACKGROUND_LIGHT_GREY,
                opacity=1.0 if picked else _CENTROID_OPACITY,
            )
            .encode(x=centroid_x, y=centroid_y, tooltip=_TERM_TOOLTIP)
        )

    dimmed = centroid_layer(False, alt.Chart(centroids).transform_filter(~picker))
    picked = centroid_layer(True, alt.Chart(centroids).transform_filter(picker))
    centroid = (
        (dimmed + picked)
        .add_params(picker)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title="Group centroid positions — click a term",
        )
    )
    # No .interactive() here on purpose: a pan gesture on this panel would fight
    # the click that selects a term.
    cloud = (
        alt.Chart(background)
        .mark_circle(size=_BACKGROUND_SIZE, color=_BACKGROUND_LIGHT_GREY, opacity=_BACKGROUND_OPACITY)
        .encode(x=gene_x, y=gene_y, tooltip=_GENE_TOOLTIP)
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
        .transform_filter(picker)
        .mark_circle(size=_MEMBER_SIZE, color=highlight, opacity=_MEMBER_OPACITY)
        .encode(x=gene_x, y=gene_y, tooltip=_GENE_TOOLTIP)
    )
    genes = (
        (cloud + highlighted)
        .properties(
            width=_CHART_WIDTH,
            height=_CHART_HEIGHT,
            title="Genes of the selected term (hover a point for its gene record)",
        )
        .interactive()
    )
    return centroid | genes


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
    chart = build_chart(detail, members, centroids, options)
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
