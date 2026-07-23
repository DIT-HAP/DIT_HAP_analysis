#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generic Config-Driven Per-Group Scatter Visualization
=====================================================

Per-dataset x source: given a config namelist of term/complex names-or-ids,
draws one feature-space subplot per resolved group — background gene cloud +
that group's members highlighted, annotated with group_name (group_id),
n_members, and the group's coherence z_score/p_value. Replaces the hardcoded
module-visualization logic from analyze_complex_modules.py; now driven by the
generic long-table + coherence metrics TSVs produced by prepare_annotation.py
and compute_coherence.py.

Input
-----
- fitting_results.tsv: the upstream per-gene fitting statistics, systematic id
  as the index (column 0), with DR/DL fitness columns. Legacy releases may still
  ship the pre-rename um/lam headers -> normalized to DR/DL. Only the systematic
  id + DR/DL are read here.
- group_annotation_long.tsv: the prepared unified long-table from
  prepare_annotation.py (one row per group-member), with the contract columns
  (source, group_id, group_name, Systematic ID, Name, n_group_genes). Maps
  groups -> member genes.
- coherence_metrics.tsv: per-group coherence results from compute_coherence.py,
  with at least (source, group_id, group_name, z_score, p_value).

Output
------
- group_scatter.pdf: one feature-space subplot per resolved group from the
  config namelist, annotated with coherence metrics. Empty namelist -> single
  placeholder "No groups resolved" panel.

Usage
-----
    python plot_group_scatter.py \\
        --fitting-results .../fitting_results.tsv \\
        --annotation results/coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --metrics results/coherence/{dataset}/{source}/coherence_metrics.tsv \\
        --source go_cc \\
        --groups "['kinetochore', 'GO:0000776']" \\
        --output-figure results/coherence/{dataset}/{source}/group_scatter.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import ast
import sys
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.gene_level import plot_given_genes_on_feature_space  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Legacy -> current metric column names (same quirk as
# workflow/src/clustering/candidates.py's _LEGACY_METRIC_RENAME).
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def parse_groups_arg(raw: str, source: str) -> list[str]:
    """Parse the --groups argument into a namelist for the given source.

    Snakemake renders a config list/dict via str(), so it arrives as a Python
    literal (single quotes) rather than strict JSON — ast.literal_eval handles
    both. Accepts BOTH:
    - a list → use directly as the namelist for source;
    - a dict → take parsed.get(source, []).
    Empty/missing → empty list (script writes a placeholder figure and exits 0).
    """
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
    elif isinstance(parsed, dict):
        return [str(x) for x in parsed.get(source, [])]
    else:
        raise ValueError(f"--groups must be a list or dict, got {type(parsed).__name__}")


def load_fitting_results(fitting_results_path: Path) -> pd.DataFrame:
    """Load upstream fitting_results.tsv, normalizing legacy um/lam -> DR/DL columns.

    The systematic id is the INDEX (column 0), matching how
    workflow/src/clustering/candidates.py reads it; we reset it to a
    `Systematic ID` column (robust to whatever the index was named). Matches
    the loading approach in compute_coherence.py::load_fitting_results.
    """
    fitting = pd.read_csv(fitting_results_path, sep="\t", index_col=0)
    fitting = fitting.reset_index()
    # After reset_index the id column keeps the index's name (e.g. "Systematic ID",
    # "gene_systematic_id", or "index" if the index was unnamed). Normalize it.
    if "Systematic ID" not in fitting.columns:
        first_col = fitting.columns[0]
        logger.info(f"Renaming fitting_results index column '{first_col}' -> 'Systematic ID'")
        fitting = fitting.rename(columns={first_col: "Systematic ID"})

    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in fitting.columns and new not in fitting.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        fitting = fitting.rename(columns=rename)

    for required in ["Systematic ID", "DR", "DL"]:
        if required not in fitting.columns:
            raise ValueError(f"fitting_results.tsv missing required column '{required}' (have: {list(fitting.columns)})")

    # Map +/-inf to NaN before dropna so non-finite DR/DL never reach the
    # feature-space scatter (dropna alone keeps +/-inf).
    return fitting.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"]).copy()


def load_long_table(annotation_path: Path) -> pd.DataFrame:
    """Load the prepared unified long-table (group -> member genes)."""
    annotation = pd.read_csv(annotation_path, sep="\t")
    for required in ["source", "group_id", "group_name", "Systematic ID"]:
        if required not in annotation.columns:
            raise ValueError(
                f"annotation missing required column '{required}' (have: {list(annotation.columns)})"
            )
    return annotation


def load_metrics(metrics_path: Path) -> pd.DataFrame:
    """Load coherence metrics TSV."""
    metrics = pd.read_csv(metrics_path, sep="\t")
    for required in ["group_id", "z_score", "p_value"]:
        if required not in metrics.columns:
            raise ValueError(
                f"metrics TSV missing required column '{required}' (have: {list(metrics.columns)})"
            )
    return metrics


def resolve_groups(
    long_table: pd.DataFrame, source: str, names: list[str]
) -> list[tuple[str, str, list[str]]]:
    """Resolve a config namelist to (group_id, group_name, sorted member Systematic IDs).

    For each entry in `names`, match rows where `group_name == entry` OR
    `group_id == entry` within the given source. If an entry matches multiple
    group_ids (same name), include each distinct group_id. Skip entries with no
    match (don't error). De-dup members within a group.

    Returns:
        List of (group_id, group_name, sorted_member_ids) tuples.
    """
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


# =============================================================================
# PLOTTING
# =============================================================================
def plot_group_scatter_figure(
    fitting_df: pd.DataFrame,
    resolved_groups: list[tuple[str, str, list[str]]],
    metrics_df: pd.DataFrame,
) -> plt.Figure:
    """One feature-space subplot per resolved group, annotated with coherence metrics.

    Each panel: background gene cloud + group's members highlighted, titled with
    group_name (group_id), n_members, and z_score/p_value from metrics_df.
    Groups missing from metrics (didn't survive size filter) show "z=NA, p=NA".
    Empty resolved_groups -> single placeholder "No groups resolved" axis.
    """
    n = max(len(resolved_groups), 1)
    col_num = min(3, n)
    row_num = int(np.ceil(n / col_num))

    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH * col_num, AX_HEIGHT * row_num))
    axes = np.atleast_1d(axes).flatten()

    if not resolved_groups:
        axes[0].text(
            0.5, 0.5, "No groups resolved from config namelist",
            ha="center", va="center", fontsize=14, transform=axes[0].transAxes
        )
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        for j in range(1, len(axes)):
            fig.delaxes(axes[j])
        fig.tight_layout()
        return fig

    # Precompute metrics lookup (O(1) per group)
    metrics_dict = {row["group_id"]: row for _, row in metrics_df.iterrows()}

    for idx, (group_id, group_name, members) in enumerate(resolved_groups):
        # Look up coherence metrics
        metrics_row = metrics_dict.get(group_id)
        if metrics_row is not None:
            z_score = metrics_row["z_score"]
            p_value = metrics_row["p_value"]
            metrics_str = f"z={z_score:.2f}, p={p_value:.3g}"
        else:
            metrics_str = "z=NA, p=NA"

        title = f"{group_name} ({group_id})\nn={len(members)}, {metrics_str}"

        plot_given_genes_on_feature_space(
            ax=axes[idx],
            data_df=fitting_df,
            genes=members,
            gene_column="Systematic ID",
            title=title,
            x_feature="DR",
            y_feature="DL",
            cmap="#9D343C",
            label=group_name,
            title_with_count=False,  # we include count in title ourselves
            s=40,
        )

    for j in range(len(resolved_groups), len(axes)):
        fig.delaxes(axes[j])

    fig.tight_layout()
    return fig


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
) -> None:
    """Load -> resolve groups -> plot feature-space figure + write PDF."""
    # Validate inputs
    for path in [fitting_results, annotation, metrics]:
        if not path.exists():
            raise ValueError(f"Required input not found: {path}")
    output_figure.parent.mkdir(parents=True, exist_ok=True)

    fitting_df = load_fitting_results(fitting_results)
    long_table = load_long_table(annotation)
    metrics_df = load_metrics(metrics)

    resolved = resolve_groups(long_table, source, groups)
    logger.info(f"Resolved {len(resolved)} groups from namelist of {len(groups)} entries")

    fig = plot_group_scatter_figure(fitting_df, resolved, metrics_df)
    with PdfPages(output_figure) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

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
    parser.add_argument("--metrics", type=Path, required=True, help="coherence_metrics.tsv")
    parser.add_argument("--source", type=str, required=True, help="Source name (e.g. go_cc)")
    parser.add_argument(
        "--groups", type=str, default="",
        help="List or dict literal of group names/ids (empty -> no groups, placeholder figure)"
    )
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
        )
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
