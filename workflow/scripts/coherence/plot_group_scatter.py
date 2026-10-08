#!/usr/bin/env python3

"""
Generic Config-Driven Per-Group Scatter Visualization
=====================================================

Per-dataset x source: given a config namelist of term/complex names-or-ids,
draws one feature-space subplot per resolved group — background gene cloud +
that group's members highlighted, annotated with group_name (group_id),
n_members, and the group's coherence median_pairwise_distance_z/_p. Replaces the hardcoded
module-visualization logic from analyze_complex_modules.py; now driven by the
generic long-table + coherence metrics tables produced by prepare_annotation.py
and compute_coherence.py.

The panel used to come from plotting.gene_level.plot_given_genes_on_feature_space.
That helper is still used by the clustering figures, so it stays; this script
draws the same picture with cns.scatterplot instead, which lets the house style
own the styling and drops two defects the helper carried — a gaussian_kde call
whose result was thrown away (a pure "will this raise?" probe) and a bare
`except Exception` that silently recoloured the subset red.

Input
-----
- fitting_results.tsv: the upstream per-gene fitting statistics (see
  coherence/io.py::load_fitting_results). Only the systematic id and DR/DL are read.
- group_annotation_long.tsv: the prepared unified long-table from
  prepare_annotation.py, with the contract columns in
  coherence/sources.py::LONG_TABLE_COLUMNS. Maps groups -> member genes.
- coherence_metrics.parquet: per-group coherence results from compute_coherence.py,
  with at least (source, group_id, group_name, median_pairwise_distance_z,
  median_pairwise_distance_p).

Output
------
- group_scatter.pdf: one feature-space subplot per resolved group from the
  config namelist, annotated with coherence metrics (+ a .review.png sibling via
  save_dual). Empty namelist -> single placeholder "No groups resolved" panel.

Usage
-----
    python plot_group_scatter.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --metrics results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --source go_cc \\
        --groups "['kinetochore', 'GO:0000776']" \\
        --output-figure results/3a_coherence/{dataset}/{source}/group_scatter.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import ast
import math
import sys
from pathlib import Path

# 2. Data Processing Imports
import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.io import load_fitting_results, load_long_table  # noqa: E402
from figures import (  # noqa: E402
    FURNITURE_COLOR,
    PanelShape,
    apply_house_style,
    fit_panels,
    grid_axes,
    house_colors,
    panel_labels,
    save_dual,
)
from io_table import read_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
_MAX_COLUMNS = 3

# Panels are SQUARE; a 100 px box fits about 25 characters at the house title
# size, so the title is split across three short lines rather than one long one.
_TITLE_NAME_WIDTH = 22

# The scored/unscored split the coherence computation draws on, shown in every
# panel: the genome cloud stays a single grey, and a dashed line marks the DR
# threshold — but the highlighted members are recoloured on either side of it,
# since only the term's own genes are asked to read against the cut. A threshold
# arrives via --dr-threshold; this is only the fallback, which must match
# config/analysis.yaml's `coherence.dr_threshold` and compute_coherence.py's.
_DR_THRESHOLD = -0.3
# Scored side (left of the line): the house red. Not-scored side: the house
# blue-teal, a genuinely different hue (not a lightness step) so the two stay
# apart in greyscale print as well.
_UNSCORED_MEMBER_COLOR = house_colors((1,))[0]


# =============================================================================
# CORE LOGIC
# =============================================================================
def parse_groups_arg(raw: str, source: str) -> list[str]:
    """Parse the --groups argument into a namelist for the given source."""
    # Snakemake renders the config list through json.dumps, but accept a Python
    # literal too: ast.literal_eval handles both. A list is used directly as the
    # namelist for `source`; a dict is looked up as parsed.get(source, []).
    # Empty/missing -> empty list (script writes a placeholder figure and exits 0).
    raw = (raw or "").strip()
    if not raw:
        return []
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError) as exc:
        raise ValueError(
            f"Could not parse --groups as a Python list/dict literal: {exc}\n"
            f"Expected: \"['term1', 'term2']\" or \"{{source: ['term1']}}\""
        ) from exc

    if isinstance(parsed, list):
        return [str(x) for x in parsed]
    if isinstance(parsed, dict):
        return [str(x) for x in parsed.get(source, [])]
    raise ValueError(f"--groups must be a list or dict, got {type(parsed).__name__}")


def load_metrics(metrics_path: Path) -> pd.DataFrame:
    """Load the coherence metrics table and check the columns the titles need."""
    metrics = read_parquet(metrics_path)
    for required in ["group_id", "median_pairwise_distance_z", "median_pairwise_distance_p"]:
        if required not in metrics.columns:
            raise ValueError(
                f"metrics table missing required column '{required}' (have: {list(metrics.columns)})"
            )
    return metrics


def resolve_groups(
    long_table: pd.DataFrame, source: str, names: list[str]
) -> list[tuple[str, str, list[str]]]:
    """Resolve a config namelist to (group_id, group_name, sorted member Systematic IDs)."""
    # Each entry matches on `group_name` OR `group_id` within the source; an entry
    # matching several group_ids (same name) contributes each distinct group_id,
    # and an entry matching nothing is skipped rather than an error.
    source_table = long_table[long_table["source"] == source]
    resolved = []
    for name in names:
        matched = source_table[
            (source_table["group_name"] == name) | (source_table["group_id"] == name)
        ]
        if matched.empty:
            logger.warning(f"No group found for '{name}' in source '{source}'; skipping.")
            continue

        # Group by group_id to handle same-name->multiple-ids case
        group_ids_matched = matched.groupby("group_id", sort=False)
        if group_ids_matched.ngroups > 1:
            logger.info(f"Config entry '{name}' matches {group_ids_matched.ngroups} distinct group_ids")
        for group_id, group_df in group_ids_matched:
            group_name = group_df["group_name"].iloc[0]
            members = sorted(set(group_df["Systematic ID"].dropna()))
            resolved.append((str(group_id), str(group_name), members))
            prefix = "  " if group_ids_matched.ngroups > 1 else ""
            logger.info(f"{prefix}Resolved '{name}' -> {group_id} ({group_name}): {len(members)} genes")

    return resolved


def metrics_caption(metrics_row: pd.Series | None) -> str:
    """The metrics half of a panel title, preferring the FDR-corrected q."""
    # Groups missing from the metrics table did not survive the size filter, and
    # say so rather than showing a number that was never computed.
    if metrics_row is None:
        return "not scored"
    q_value = metrics_row.get("q_value", np.nan)
    if pd.notna(q_value):
        return f"z={metrics_row['median_pairwise_distance_z']:.2f}, q={q_value:.3g}"
    return (f"z={metrics_row['median_pairwise_distance_z']:.2f}, "
            f"p={metrics_row['median_pairwise_distance_p']:.3g}")


def panel_title(group_id: str, group_name: str, n_members: int, metrics_row: pd.Series | None) -> str:
    """Panel title: name, id, member count, and the coherence metrics."""
    # Three lines, not one: a one-line "name (id), n=..., z=..." runs to ~37
    # characters and overflows the panel, colliding with the panel label beside it.
    name = group_name if len(group_name) <= _TITLE_NAME_WIDTH else group_name[: _TITLE_NAME_WIDTH - 1] + "…"
    return f"{name}\n({group_id}), n={n_members}\n{metrics_caption(metrics_row)}"


def draw_split_cloud(
    fitting_df: pd.DataFrame, dr_threshold: float, member_color: str
) -> None:
    """The genome cloud on the current axes, in one grey, with the dashed threshold line."""
    cns.scatterplot(fitting_df, "DR", "DL", ax=plt.gca(), color=FURNITURE_COLOR, legend=False)
    plt.axvline(
        dr_threshold, color=FURNITURE_COLOR, linestyle="--",
        linewidth=plt.rcParams["lines.linewidth"] * 0.8, zorder=0.5,
    )


def plot_group_scatter_figure(
    fitting_df: pd.DataFrame,
    resolved_groups: list[tuple[str, str, list[str]]],
    metrics_df: pd.DataFrame,
    dr_threshold: float = _DR_THRESHOLD,
) -> None:
    """One feature-space panel per resolved group, annotated with coherence metrics."""
    # Each panel: the genome-wide background cloud in furniture grey with the
    # group's members over it. An empty namelist draws a single placeholder.
    apply_house_style()

    n_panels = max(len(resolved_groups), 1)
    n_cols = min(_MAX_COLUMNS, n_panels)
    n_rows = math.ceil(n_panels / n_cols)
    # Every panel draws the same DR/DL space, so they must share both ranges:
    # autoscaled panels make a tight group and a sprawled one look identical, which
    # is the one thing these panels exist to distinguish. grid_axes drops the
    # interior tick labels once the range is shared, so the comparison reads cleanly.
    axes = grid_axes(
        n_rows, n_cols, labels=panel_labels(n_panels), shape=PanelShape.SQUARE,
        share_x=True, share_y=True,
    )
    for ax in axes[n_panels:]:
        # grid_axes fills every cell; fit_panels measures only visible axes, so
        # the unfilled ones must be hidden (not just deleted) to be discounted.
        ax.set_visible(False)

    if not resolved_groups:
        axes[0].text(
            0.5, 0.5, "No groups resolved from config namelist",
            ha="center", va="center", transform=axes[0].transAxes,
        )
        axes[0].set_axis_off()
        fit_panels()
        return

    member_color = house_colors((0,))[0]
    scored = metrics_df.set_index("group_id")
    for ax, (group_id, group_name, members) in zip(axes, resolved_groups):
        plt.sca(ax)
        draw_split_cloud(fitting_df, dr_threshold, member_color)
        highlighted = fitting_df[fitting_df["Systematic ID"].isin(members)]
        # The term's own genes are the ones asked to read against the cut, so they
        # carry the threshold split: scored members in red, not-scored ones in the
        # blue-teal — both cns.scatterplot calls with legend=False, no legend needed.
        scored_members = highlighted[highlighted["DR"] < dr_threshold]
        unscored_members = highlighted[highlighted["DR"] >= dr_threshold]
        if not unscored_members.empty:
            cns.scatterplot(unscored_members, "DR", "DL", ax=ax,
                            color=_UNSCORED_MEMBER_COLOR, legend=False)
        cns.scatterplot(scored_members, "DR", "DL", ax=ax, color=member_color, legend=False)

        metrics_row = scored.loc[group_id] if group_id in scored.index else None
        ax.set(
            xlabel="DR",
            ylabel="DL",
            title=panel_title(group_id, group_name, len(members), metrics_row),
        )

    fit_panels()


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(
    fitting_results: Path,
    annotation: Path,
    metrics: Path,
    source: str,
    groups: list[str],
    output_figure: Path,
    dr_threshold: float = _DR_THRESHOLD,
) -> None:
    """Load -> resolve groups -> plot feature-space figure + write PDF."""
    for path in [fitting_results, annotation, metrics]:
        if not path.exists():
            raise ValueError(f"Required input not found: {path}")
    output_figure.parent.mkdir(parents=True, exist_ok=True)

    fitting_df = load_fitting_results(fitting_results)
    long_table = load_long_table(annotation)
    metrics_df = load_metrics(metrics)

    resolved = resolve_groups(long_table, source, groups)
    logger.info(f"Resolved {len(resolved)} groups from namelist of {len(groups)} entries")

    plot_group_scatter_figure(fitting_df, resolved, metrics_df, dr_threshold)
    save_dual(output_figure.with_suffix(""))
    logger.success(f"Wrote {output_figure}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Visualize gene groups in fitness space with coherence annotations"
    )
    parser.add_argument("--fitting-results", type=Path, required=True, help="fitting_results.tsv")
    parser.add_argument("--annotation", type=Path, required=True, help="group_annotation_long.tsv")
    parser.add_argument("--metrics", type=Path, required=True, help="coherence_metrics.parquet")
    parser.add_argument("--source", type=str, required=True, help="Source name (e.g. go_cc)")
    parser.add_argument(
        "--groups", type=str, default="",
        help="List or dict literal of group names/ids (empty -> no groups, placeholder figure)"
    )
    parser.add_argument("--dr-threshold", type=float, default=_DR_THRESHOLD,
                        help="The coherence scored/unscored cut, drawn as a dashed line")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output scatter PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: parse args, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        groups = parse_groups_arg(args.groups, args.source)
        run(
            fitting_results=args.fitting_results,
            annotation=args.annotation,
            metrics=args.metrics,
            source=args.source,
            groups=groups,
            output_figure=args.output_figure,
            dr_threshold=args.dr_threshold,
        )
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
