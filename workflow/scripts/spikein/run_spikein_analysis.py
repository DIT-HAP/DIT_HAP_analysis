#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Spike-In Dilution Linearity QC
================================

Standalone (no {dataset} wildcard) analysis of the Spikein calibration
project: extracts read counts at 5 known spike-in insertion sites from the
filtered insertion table, assigns each sample its known dilution ratio by
read-count rank, and fits a log-log linear regression of read ratio vs
dilution ratio to check assay linearity. Ported from
DIT_HAP_pipeline/workflow/notebooks/spike_in.ipynb.

Input
-----
- Filtered raw-reads insertion table (Spikein's pre-release results/13_filtered/
  raw_reads.filtered.tsv — release/ never packages this file, see spikein.smk),
  indexed by [Chr, Coordinate, Strand] (a 4th Target level, if present, is
  dropped), columned by [Sample, Timepoint] (one Timepoint per dilution point).

Output
------
- spike_in_stats.tsv: long-form per-site/per-sample table (Reads, Ratio,
  Relative_Read_Ratio, Relative_Dilution_Ratio).
- spike_in_correlation.pdf: scatter of the 5 sites + combined linear fit.

Usage
-----
    python run_spikein_analysis.py \\
        --raw-reads .../Spikein/results/13_filtered/raw_reads.filtered.tsv \\
        --output-stats results/spikein/spike_in_stats.tsv \\
        --output-figure results/spikein/spike_in_correlation.pdf \\
        --spike-in-sites-json '{"DY215": {"chr": "I", "coord": 3749394, "strand": "-"}, ...}'

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402
from scipy.stats import linregress  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.style import COLORS  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Known dilution series for the 6 Spikein samples (Spikein0..Spikein5), most
# dilute first. Byte-faithful to the source notebook's spike_in_ratio.
SPIKE_IN_RATIO = np.array([1.5, 4, 16, 64, 256, 1024]) / 100000

# The 5 known spike-in insertion coordinates (chr, coord, strand). Overridable
# via --spike-in-sites-json; defaults kept here so the script is runnable
# standalone. Mirrors config/analysis.yaml's spikein.coordinates section.
DEFAULT_SPIKE_IN_SITES = {
    "DY215": {"chr": "I", "coord": 3749394, "strand": "-"},
    "DY217": {"chr": "II", "coord": 3344505, "strand": "-"},
    "DY218": {"chr": "II", "coord": 185161, "strand": "-"},
    "DY339": {"chr": "II", "coord": 1157130, "strand": "-"},
    "DY348": {"chr": "II", "coord": 3065244, "strand": "-"},
}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class SpikeinConfig:
    """Inputs, outputs, and spike-in site coordinates for the linearity QC."""
    raw_reads: Path
    output_stats: Path
    output_figure: Path
    spike_in_sites: dict[str, dict]

    def validate(self) -> None:
        """Raise ValueError if the input is missing, then ensure output dirs exist."""
        if not self.raw_reads.exists():
            raise ValueError(f"Required input not found: {self.raw_reads}")
        for out in [self.output_stats, self.output_figure]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
def assign_ratio_by_order(sub_df: pd.DataFrame, spike_in_ratio: np.ndarray) -> pd.DataFrame:
    """Rank one site's Reads ascending, assign the matching dilution ratio, and
    log2-normalise both reads and ratios relative to their max (in place on a copy).

    Port of the source notebook's assign_ratio_by_order: rank 0 (lowest read
    count) maps to spike_in_ratio[0] (most dilute), etc. The minimum read
    count is subtracted from every sample (background-subtraction), so the
    lowest sample floors at exactly 0.

    Note: the original notebook applied this subtraction via an asymmetric
    `Series.where(Reads == min_val, Reads - min_val)`, which leaves the
    minimum row's value UNCHANGED (a notebook quirk/bug) instead of zeroing
    it. This port applies the subtraction uniformly so the minimum is 0, per
    the intended "subtract min read" behaviour.
    """
    read_rank = (sub_df["Reads"].rank() - 1).astype(int).to_numpy()
    sub_df["Ratio"] = spike_in_ratio[read_rank]
    sub_df["Reads"] = sub_df["Reads"] - sub_df["Reads"].min()
    sub_df["Relative_Read_Ratio"] = np.log2(sub_df["Reads"] / sub_df["Reads"].max())
    sub_df["Relative_Dilution_Ratio"] = np.log2(sub_df["Ratio"] / sub_df["Ratio"].max())
    return sub_df


def build_spike_sites_df(raw_reads: pd.DataFrame, spike_in_sites: dict[str, dict]) -> pd.DataFrame:
    """Extract the known spike-in sites from the filtered reads table.

    `raw_reads` is indexed by [Chr, Coordinate, Strand] (any further index
    levels, e.g. Target, are dropped since they don't disambiguate a site) and
    columned by samples. Returns one row per site with the raw Reads array
    plus Chr/Coordinate/Strand/Strain/Name columns.
    """
    if raw_reads.index.nlevels > 3:
        raw_reads = raw_reads.droplevel(list(range(3, raw_reads.index.nlevels)))

    records = []
    for order, (strain, info) in enumerate(spike_in_sites.items(), start=1):
        key = (info["chr"], info["coord"], info["strand"])
        try:
            reads = raw_reads.loc[key, :].to_numpy()
        except KeyError:
            logger.warning(f"Spike-in site not found: {strain} ({info['chr']}:{info['coord']} {info['strand']})")
            continue
        records.append({
            "Chr": info["chr"],
            "Coordinate": info["coord"],
            "Strand": info["strand"],
            "Strain": strain,
            "Name": f"Spike-in Insertion {order}",
            "Reads": reads,
        })

    if not records:
        raise ValueError("No spike-in sites found in the raw reads table")

    return pd.DataFrame(records)


def build_spike_in_stats(raw_reads: pd.DataFrame, spike_in_sites: dict[str, dict]) -> pd.DataFrame:
    """Build the full long-form spike-in stats table: extract sites, reshape to
    long form (one row per site x sample), then assign dilution ratios per site.
    """
    sites_df = build_spike_sites_df(raw_reads, spike_in_sites)

    n_samples = len(sites_df["Reads"].iloc[0])
    sample_names = [f"Spikein{i}" for i in range(n_samples)]
    wide = pd.DataFrame(sites_df["Reads"].tolist(), columns=sample_names)
    wide[["Chr", "Coordinate", "Strand", "Strain", "Name"]] = sites_df[
        ["Chr", "Coordinate", "Strand", "Strain", "Name"]
    ].values

    long_df = (
        wide.sort_values(["Chr", "Coordinate"])
        .sort_values(["Strain"])
        .reset_index(drop=True)
        .set_index(["Chr", "Coordinate", "Strand", "Name", "Strain"])
        .rename_axis("Sample", axis=1)
        .stack()
        .to_frame("Reads")
    )

    return long_df.groupby("Strain").apply(assign_ratio_by_order, spike_in_ratio=SPIKE_IN_RATIO).droplevel(0, axis=0)


def compute_linear_regression_stats(x: pd.Series, y: pd.Series) -> dict[str, float]:
    """Fit y = slope*x + intercept and return slope/intercept/r_value/p_value/std_err/r2."""
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return {
        "slope": slope,
        "intercept": intercept,
        "r_value": r_value,
        "p_value": p_value,
        "std_err": std_err,
        "r2": r_value ** 2,
    }


def plot_spike_in_correlation(spike_in_stats: pd.DataFrame, stats: dict[str, float], output: Path) -> None:
    """Scatter Relative_Dilution_Ratio vs Relative_Read_Ratio per site + the combined linear fit."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    for idx, (name, sub) in enumerate(spike_in_stats.groupby("Name")):
        ax.scatter(
            sub["Relative_Dilution_Ratio"], sub["Relative_Read_Ratio"],
            label=name, facecolor="none", edgecolor=COLORS[idx % len(COLORS)],
            s=150, lw=1.5, alpha=0.75,
        )

    line_x = np.array([-10, 0])
    line_y = stats["slope"] * line_x + stats["intercept"]
    ax.plot(line_x, line_y, color="black", ls="--", alpha=0.7, lw=2.5)

    ax.set_xlabel("log$_{2}$(relative dilution ratio)")
    ax.set_ylabel("log$_{2}$(relative read ratio)")
    ax.text(
        0.05, 0.95,
        f"Slope={stats['slope']:.2f}\nPCC={stats['r_value']:.2f}\nR$^2$={stats['r2']:.2f}",
        transform=ax.transAxes, ha="left", va="top",
    )
    ax.legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


@logger.catch(reraise=True)
def run(config: SpikeinConfig) -> None:
    """Load raw reads -> extract spike-in sites -> assign ratios -> fit -> save TSV + PDF."""
    config.validate()

    raw_reads = pd.read_csv(config.raw_reads, sep="\t", header=[0, 1], index_col=[0, 1, 2, 3])
    spike_in_stats = build_spike_in_stats(raw_reads, config.spike_in_sites)
    spike_in_stats.to_csv(config.output_stats, sep="\t")

    # The lowest-read sample per site floors at Reads=0 -> Relative_Read_Ratio
    # = log2(0) = -inf (kept in the TSV as the true computed value). Excluded
    # here so the fit/plot aren't skewed or broken by a non-finite point.
    finite = spike_in_stats[np.isfinite(spike_in_stats["Relative_Read_Ratio"])]

    stats = compute_linear_regression_stats(
        finite["Relative_Dilution_Ratio"], finite["Relative_Read_Ratio"]
    )
    plot_spike_in_correlation(finite, stats, config.output_figure)

    logger.success(
        f"Spike-in linearity: slope={stats['slope']:.3f}, R2={stats['r2']:.3f}, "
        f"p={stats['p_value']:.2e} ({len(spike_in_stats)} site x sample rows)"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Spike-in dilution linearity QC")
    parser.add_argument("--raw-reads", type=Path, required=True, help="Filtered raw-reads insertion table (tsv)")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output long-form spike-in stats TSV")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output correlation figure PDF")
    parser.add_argument(
        "--spike-in-sites-json", type=str, default=None,
        help="JSON dict of {strain: {chr, coord, strand}} (default: the 5 hardcoded DY sites)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        spike_in_sites = (
            json.loads(args.spike_in_sites_json) if args.spike_in_sites_json else DEFAULT_SPIKE_IN_SITES
        )
        config = SpikeinConfig(
            raw_reads=args.raw_reads,
            output_stats=args.output_stats,
            output_figure=args.output_figure,
            spike_in_sites=spike_in_sites,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
