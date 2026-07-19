#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PCR / Library-Prep Quality Control Figure
==========================================

Builds the 2x2 library-prep QC figure ported from
DIT_HAP_pipeline/workflow/notebooks/thesis_figures.ipynb ("2. PCR quality
control"). Four panels, all deterministic:
  (a) PBL vs PBR primer reads of one library.
  (b) Technical replicate: one sample processed in two upstream projects.
  (c) Biological replicate: two samples within one project.
  (d) Spike-in dilution linearity (currently a placeholder table; see below).

Panels (a)-(c) use the domain-agnostic create_scatter_correlation_plot; panel
(d) does its own spike-in linregress scatter (specific to this figure).

Input
-----
- Panel (a): one merged reads TSV (results/8_merged/...), columns PBL, PBR, Reads.
- Panels (b)/(c): two merged reads TSVs each, joined on (Chr, Coordinate, Strand).
- Panel (d) spike-in table:
    current: resources/curated/spike_in_results_PLACEHOLDER.tsv  (placeholder copy)
    future : results/spikein/spike_in_results.tsv  (produced by spikein.smk, Phase 3+)

Output
------
- PCR_quality_control.pdf: the assembled 2x2 figure.

Usage
-----
    python plot_pcr_qc.py \\
        --pbl-pbr .../8_merged/LD1328-7_0h_YES.tsv \\
        --tech-rep-1 .../LD_DIT_HAP_generationRAW/.../LD1328-7_0h_YES.tsv \\
        --tech-rep-2 .../Spore2YES6_1328/.../LD1328-7_0h_YES.tsv \\
        --bio-rep-1 .../8_merged/LD1328-4_0h_YES.tsv \\
        --bio-rep-2 .../8_merged/LD1328-8_0h_YES.tsv \\
        --spikein resources/curated/spike_in_results_PLACEHOLDER.tsv \\
        --output results/pcr_qc/PCR_quality_control.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-19
Version:  1.0.0
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
from workflow.src.plotting.generic import create_scatter_correlation_plot  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH, COLORS  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PCRQCConfig:
    """Resolved input/output paths for the 2x2 PCR QC figure."""
    pbl_pbr: Path
    tech_rep_1: Path
    tech_rep_2: Path
    bio_rep_1: Path
    bio_rep_2: Path
    spikein: Path
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any input is missing, then ensure the output dir exists."""
        for path in [self.pbl_pbr, self.tech_rep_1, self.tech_rep_2,
                     self.bio_rep_1, self.bio_rep_2, self.spikein]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def _read_merged(path: Path) -> pd.DataFrame:
    """Read a merged reads TSV indexed by (Chr, Coordinate, Strand) with PBL/PBR/Reads."""
    return pd.read_csv(path, sep="\t", index_col=[0, 1, 2])


def _plot_spikein(ax: plt.Axes, spikein: pd.DataFrame) -> None:
    """Panel (d): spike-in dilution vs read ratio scatter + linear fit (figure-specific)."""
    # Spikein0 is the zero-dilution reference and is excluded from the linearity fit,
    # matching the source notebook.
    spikein = spikein.query("Sample != 'Spikein0'")

    for idx, (name, sub) in enumerate(spikein.groupby("Name")):
        ax.scatter(
            sub["Relative_Dilution_Ratio"], sub["Relative_Read_Ratio"],
            label=name, facecolor="none", edgecolor=COLORS[idx % len(COLORS)],
            s=150, lw=1.5, alpha=0.9,
        )

    slope, intercept, r_value, _p, _se = linregress(
        spikein["Relative_Dilution_Ratio"], spikein["Relative_Read_Ratio"]
    )
    line_x = np.array([-8, 0])
    ax.plot(line_x, slope * line_x + intercept, color="black", ls="--", alpha=0.7, lw=2.5)

    ticks = [-8, -6, -4, -2, 0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel("log$_{2}$(relative dilution ratio)")
    ax.set_ylabel("log$_{2}$(relative read ratio)")
    ax.text(
        0.05, 0.95,
        f"PCC={r_value:.2f}\nR²={r_value**2:.2f}\nSlope={slope:.2f}\nIntercept={intercept:.2f}",
        transform=ax.transAxes, ha="left", va="top",
    )
    ax.legend(loc="lower right", fontsize=16, frameon=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PCRQCConfig) -> None:
    """Assemble the 2x2 QC figure and save it as a PDF."""
    config.validate()

    # Panel (a): PBL vs PBR of a single library.
    pbl_pbr = _read_merged(config.pbl_pbr)

    # Panel (b): technical replicate — same sample, two upstream projects.
    tech = pd.merge(
        _read_merged(config.tech_rep_1), _read_merged(config.tech_rep_2),
        left_index=True, right_index=True, suffixes=("_1", "_2"),
    )

    # Panel (c): biological replicate — two samples, one project.
    bio = pd.merge(
        _read_merged(config.bio_rep_1), _read_merged(config.bio_rep_2),
        left_index=True, right_index=True, suffixes=("_1", "_2"),
    )

    # Panel (d): spike-in linearity.
    spikein = pd.read_csv(config.spikein, sep="\t")

    fig, axes = plt.subplot_mosaic(
        [["(a)", "(b)"], ["(c)", "(d)"]], figsize=(AX_WIDTH * 2, AX_HEIGHT * 2)
    )

    create_scatter_correlation_plot(pbl_pbr["PBL"].values, pbl_pbr["PBR"].values, ax=axes["(a)"], xscale="log", yscale="log")
    axes["(a)"].set_xlabel("PBL Reads")
    axes["(a)"].set_ylabel("PBR Reads")

    create_scatter_correlation_plot(tech["Reads_1"], tech["Reads_2"], ax=axes["(b)"], xscale="log", yscale="log")
    axes["(b)"].set_xlabel("Reads of Technical Replicate 1")
    axes["(b)"].set_ylabel("Reads of Technical Replicate 2")

    create_scatter_correlation_plot(bio["Reads_1"], bio["Reads_2"], ax=axes["(c)"], xscale="log", yscale="log")
    axes["(c)"].set_xlabel("Reads of Biological Replicate 1")
    axes["(c)"].set_ylabel("Reads of Biological Replicate 2")

    _plot_spikein(axes["(d)"], spikein)

    fig.tight_layout(h_pad=2, w_pad=2)
    fig.savefig(config.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Wrote PCR QC figure: {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Build the 2x2 PCR / library-prep QC figure")
    parser.add_argument("--pbl-pbr", type=Path, required=True, help="Panel (a): merged reads TSV (PBL vs PBR)")
    parser.add_argument("--tech-rep-1", type=Path, required=True, help="Panel (b): technical replicate 1 merged reads TSV")
    parser.add_argument("--tech-rep-2", type=Path, required=True, help="Panel (b): technical replicate 2 merged reads TSV")
    parser.add_argument("--bio-rep-1", type=Path, required=True, help="Panel (c): biological replicate 1 merged reads TSV")
    parser.add_argument("--bio-rep-2", type=Path, required=True, help="Panel (c): biological replicate 2 merged reads TSV")
    parser.add_argument("--spikein", type=Path, required=True, help="Panel (d): spike-in results TSV")
    parser.add_argument("--output", type=Path, required=True, help="Output figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the figure, report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PCRQCConfig(
            pbl_pbr=args.pbl_pbr,
            tech_rep_1=args.tech_rep_1,
            tech_rep_2=args.tech_rep_2,
            bio_rep_1=args.bio_rep_1,
            bio_rep_2=args.bio_rep_2,
            spikein=args.spikein,
            output=args.output,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
