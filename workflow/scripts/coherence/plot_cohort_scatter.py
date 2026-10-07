#!/usr/bin/env python3

"""
Cohort Scatter Grids (coherent / incoherent terms) — Visualization
===================================================================

Per view (the five sources + dedup, not combined): one feature-space panel per
term of the view's COHERENT or INCOHERENT cohort, in the same format as
plot_group_scatter.py — the genome-wide DR-DL cloud in furniture grey with the
term's scored members highlighted in the house red, the panel titled with the
term name and its z/q evidence.

The cohort cuts are the figures' and workbooks' own (config/analysis.yaml's
`coherence:` block, read from --q-max/--coherent-z/--incoherent-z):

- coherent   = q_value <= q_max AND median_pairwise_distance_z < coherent_z
- incoherent = median_pairwise_distance_z > incoherent_z (a descriptive rank
  cut; see export_cohorts.py's docstring for why this end carries no FDR)

Panels are ordered most-extreme-first within each cohort, so the grid opens on
the strongest signal. Empty cohort -> single placeholder panel, exit 0.

Input
-----
- --metrics: the view's metrics table — a per-source coherence_metrics.parquet,
  or dedup/coherence_terms_representatives.tsv. Needs median_pairwise_distance_z,
  q_value, group_name, and (for dedup) source.
- --annotation: the matching group_annotation_long.tsv table(s); the member
  genes are resolved through it (Systematic ID), matching how
  compute_coherence.py scored them.
- --fitting-results: the upstream per-gene fitness table, i.e. the plotted
  coordinates (see coherence/io.py::load_fitting_results).
- --cohort: `coherent` or `incoherent`.

Output
------
- --output: cohort_scatter_{cohort}.pdf (+ a .review.png sibling via save_dual).

Usage
-----
    python plot_cohort_scatter.py \\
        --metrics results/3a_coherence/{dataset}/{source}/coherence_metrics.parquet \\
        --annotation results/3a_coherence/{dataset}/{source}/group_annotation_long.tsv \\
        --fitting-results .../gene_level/fitting_results.tsv \\
        --cohort coherent \\
        --q-max 0.05 --coherent-z -2.0 --incoherent-z 0.5 \\
        --output results/3a_coherence/{dataset}/{source}/cohort_scatter_coherent.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 5.5 (implementation)
Date:     2026-10-07
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from enum import StrEnum
from math import ceil
from pathlib import Path

# 2. Data Processing Imports
import cnsplots as cns
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
from io_table import read_file  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
class CohortSide(StrEnum):
    """The two cohort ends the figure can draw."""

    COHERENT = "coherent"
    INCOHERENT = "incoherent"


_Z_COLUMN = "median_pairwise_distance_z"
_Q_COLUMN = "q_value"
_MAX_COLUMNS = 3
# Panel titles truncate name AND id: kegg_brite's fallback group_id is the full
# BRITE path (up to 320 characters), which would run the middle line across the
# panel. Same cap as plot_group_scatter.py's name line.
_TITLE_WIDTH = 22


def _clip(text: str, width: int = _TITLE_WIDTH) -> str:
    """A panel-title line that fits a 100 px panel: truncated with an ellipsis."""
    return text if len(text) <= width else text[: width - 1] + "…"

# Per-gene hover card is not needed here — this is a static grid — but the
# cohort thresholds ARE persisted into the figure path/log only, matching the
# other plot scripts.


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class CohortScatterConfig:
    """Inputs, the cohort to draw, its cuts, and the output figure path."""

    metrics: Path
    annotations: tuple[Path, ...]
    fitting_results: Path
    cohort: CohortSide
    q_max: float
    coherent_z: float
    incoherent_z: float
    output: Path

    def validate(self) -> None:
        """Raise ValueError on missing inputs / a bad cut, then make output dirs."""
        for path in [self.metrics, *self.annotations, self.fitting_results]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.coherent_z >= self.incoherent_z:
            raise ValueError(
                f"coherent_z ({self.coherent_z}) must be below incoherent_z ({self.incoherent_z})"
            )
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def cohort_terms(metrics: pd.DataFrame, config: CohortScatterConfig) -> pd.DataFrame:
    """The view's cohort rows, most extreme first (strongest evidence opens the grid)."""
    z = metrics[_Z_COLUMN]
    if config.cohort == CohortSide.COHERENT:
        eligible = metrics[(metrics[_Q_COLUMN] <= config.q_max) & (z < config.coherent_z)]
        return eligible.sort_values(_Z_COLUMN)
    eligible = metrics[z > config.incoherent_z]
    return eligible.sort_values(_Z_COLUMN, ascending=False)


def member_sets_by_group(long_table: pd.DataFrame) -> dict[str, set[str]]:
    """{group_id: annotated member Systematic IDs} — the sets compute_coherence scored."""
    return long_table.groupby("group_id")["Systematic ID"].apply(set).to_dict()


def draw_grid(
    fitting: pd.DataFrame,
    members_by_group: dict[str, set[str]],
    cohort: pd.DataFrame,
    side: CohortSide,
) -> None:
    """One feature-space panel per cohort term, sharing both axes ranges."""
    apply_house_style()
    n_panels = max(len(cohort), 1)
    n_cols = min(_MAX_COLUMNS, n_panels)
    n_rows = ceil(n_panels / n_cols)
    axes = grid_axes(
        n_rows, n_cols, labels=panel_labels(n_panels), shape=PanelShape.SQUARE,
        share_x=True, share_y=True,
    )
    for ax in axes[n_panels:]:
        # grid_axes fills every cell; fit_panels measures only visible axes.
        ax.set_visible(False)

    if cohort.empty:
        axes[0].text(
            0.5, 0.5, f"No {side.value} terms in this view",
            ha="center", va="center", transform=axes[0].transAxes,
        )
        axes[0].set_axis_off()
        fit_panels()
        return

    member_color = house_colors((0,))[0]
    for ax, (_, row) in zip(axes, cohort.iterrows()):
        members = members_by_group.get(row["group_id"], set())
        cns.scatterplot(fitting, "DR", "DL", ax=ax, color=FURNITURE_COLOR, legend=False)
        highlighted = fitting[fitting["Systematic ID"].isin(members)]
        if not highlighted.empty:
            cns.scatterplot(highlighted, "DR", "DL", ax=ax, color=member_color, legend=False)
        name = _clip(str(row["group_name"]))
        group_id = str(row["group_id"])
        if " > " in group_id:
            # kegg_brite's fallback id is the full tree path (up to 320 chars);
            # the tree segment alone is the short, stable part worth a title line.
            group_id = group_id.split(" > ")[0].split(":")[0]
        ax.set(
            xlabel="DR",
            ylabel="DL",
            title=f"{name}\n({_clip(group_id, _TITLE_WIDTH - 4)}), n={row['n_scored_members']}\n"
                  f"z={row[_Z_COLUMN]:.2f}, q={row[_Q_COLUMN]:.3g}",
        )
    fit_panels()


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: CohortScatterConfig) -> None:
    """Load -> slice the cohort -> resolve members -> draw the grid -> save."""
    config.validate()
    fitting = load_fitting_results(config.fitting_results)
    # read_file dispatches on extension: per-source views hand Parquet, the dedup
    # view hands the representatives TSV.
    metrics = read_file(config.metrics)
    missing = [c for c in (_Z_COLUMN, _Q_COLUMN, "group_id", "group_name", "n_scored_members")
               if c not in metrics.columns]
    if missing:
        raise ValueError(f"metrics table missing required column(s) {missing} (have: {list(metrics.columns)})")

    long_table = pd.concat(
        [load_long_table(path) for path in config.annotations], ignore_index=True
    )
    cohort = cohort_terms(metrics, config)
    logger.info(
        f"[{config.cohort.value}] {len(cohort):,} terms in {config.metrics.name} "
        f"(q<={config.q_max:g} & z<{config.coherent_z:g}"
        if config.cohort == CohortSide.COHERENT else
        f"[{config.cohort.value}] {len(cohort):,} terms in {config.metrics.name} "
        f"(z>{config.incoherent_z:g})"
    )

    draw_grid(fitting, member_sets_by_group(long_table), cohort, config.cohort)
    save_dual(config.output.with_suffix(""))
    logger.success(f"[{config.cohort.value}] {len(cohort):,} panels -> {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="One feature-space scatter panel per coherent/incoherent term of a view",
    )
    parser.add_argument("--metrics", type=Path, required=True,
                        help="The view's metrics table (per-source Parquet or the dedup representatives TSV)")
    parser.add_argument("--annotation", dest="annotations", type=Path, nargs="+", required=True,
                        help="The matching group_annotation_long.tsv table(s), concatenated")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv")
    parser.add_argument("--cohort", required=True, choices=[s.value for s in CohortSide],
                        help="Which cohort end to draw")
    parser.add_argument("--q-max", type=float, default=0.05, help="Coherent cohort: q at or below this")
    parser.add_argument("--coherent-z", type=float, default=-2.0, help="Coherent cohort: z below this")
    parser.add_argument("--incoherent-z", type=float, default=0.5, help="Incoherent cohort: z above this")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output cohort_scatter_{cohort}.pdf")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, draw the cohort grid, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = CohortScatterConfig(
            metrics=args.metrics,
            annotations=tuple(args.annotations),
            fitting_results=args.fitting_results,
            cohort=CohortSide(args.cohort),
            q_max=args.q_max,
            coherent_z=args.coherent_z,
            incoherent_z=args.incoherent_z,
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
