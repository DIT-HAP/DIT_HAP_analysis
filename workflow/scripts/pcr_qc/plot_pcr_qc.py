#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PCR / Library-Prep Quality Control Figure
==========================================

Stage 2 of the PCR QC split: read the 4 parquet intermediates produced by
prepare_pcr_qc_data.py and assemble the 2x2 library-prep QC figure ported from
DIT_HAP_pipeline/workflow/notebooks/thesis_figures.ipynb ("2. PCR quality
control"). Four panels, all deterministic:
  (a) PBL vs PBR primer reads of one library.
  (b) Technical replicate: one sample processed in two upstream projects.
  (c) Biological replicate: two samples within one project.
  (d) Spike-in dilution linearity (currently a placeholder table; see below).

Panels (a)-(c) use the domain-agnostic create_scatter_correlation_plot; panel
(d) uses core.plot_spikein_panel (specific to this figure).

Input
-----
- pbl_pbr.parquet: panel (a), columns PBL, PBR, Reads.
- tech.parquet: panel (b), pre-merged technical replicate pair (Reads_1/Reads_2).
- bio.parquet: panel (c), pre-merged biological replicate pair (Reads_1/Reads_2).
- spikein.parquet: panel (d) spike-in dilution table.

Output
------
- PCR_quality_control.pdf: the assembled 2x2 figure.

Usage
-----
    python plot_pcr_qc.py \\
        --pbl-pbr results/pcr_qc/_work/pbl_pbr.parquet \\
        --tech results/pcr_qc/_work/tech.parquet \\
        --bio results/pcr_qc/_work/bio.parquet \\
        --spikein results/pcr_qc/_work/spikein.parquet \\
        --output results/pcr_qc/PCR_quality_control.pdf

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from loguru import logger  # noqa: E402

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import read_parquet  # noqa: E402
from workflow.src.pcr_qc.core import plot_spikein_panel  # noqa: E402
from workflow.src.plotting.generic import create_scatter_correlation_plot  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PlotPCRQCConfig:
    """Resolved parquet input / PDF output paths for the 2x2 PCR QC figure."""
    pbl_pbr: Path
    tech: Path
    bio: Path
    spikein: Path
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any input is missing, then ensure the output dir exists."""
        for path in [self.pbl_pbr, self.tech, self.bio, self.spikein]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PlotPCRQCConfig) -> None:
    """Assemble the 2x2 QC figure and save it as a PDF."""
    config.validate()

    pbl_pbr = read_parquet(config.pbl_pbr)
    tech = read_parquet(config.tech)
    bio = read_parquet(config.bio)
    spikein = read_parquet(config.spikein)

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

    plot_spikein_panel(axes["(d)"], spikein)

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
    parser.add_argument("--pbl-pbr", type=Path, required=True, help="Panel (a): pbl_pbr.parquet")
    parser.add_argument("--tech", type=Path, required=True, help="Panel (b): tech.parquet (merged technical replicate pair)")
    parser.add_argument("--bio", type=Path, required=True, help="Panel (c): bio.parquet (merged biological replicate pair)")
    parser.add_argument("--spikein", type=Path, required=True, help="Panel (d): spikein.parquet")
    parser.add_argument("--output", type=Path, required=True, help="Output figure PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, render the figure, report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PlotPCRQCConfig(
            pbl_pbr=args.pbl_pbr,
            tech=args.tech,
            bio=args.bio,
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
