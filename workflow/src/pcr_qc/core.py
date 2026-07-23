"""
PCR / Library-Prep Quality Control — Core Logic
================================================

Shared loader and figure-panel logic for the PCR QC stage. Ported from
DIT_HAP_pipeline/workflow/notebooks/thesis_figures.ipynb ("2. PCR quality
control") and factored out of the original single-script port so the stage
can be split into independent Snakemake rules (prepare -> plot), each
re-runnable on its own.

Usage
-----
    from workflow.src.pcr_qc.core import read_merged_reads, plot_spikein_panel
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: builders only write PDFs, never display
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import linregress  # noqa: E402

# 4. Local Imports
from workflow.src.plotting.style import COLORS


# =============================================================================
# LOADERS
# =============================================================================
def read_merged_reads(path: Path) -> pd.DataFrame:
    """Read a merged reads TSV indexed by (Chr, Coordinate, Strand) with PBL/PBR/Reads."""
    return pd.read_csv(path, sep="\t", index_col=[0, 1, 2])


# =============================================================================
# FIGURES
# =============================================================================
def plot_spikein_panel(ax: plt.Axes, spikein: pd.DataFrame) -> None:
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
