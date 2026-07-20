"""
Plotting Style
==============

Loads the shared DIT-HAP matplotlib style and exposes the derived layout
constants (axis dimensions, color cycle) used across plotting modules.
Design doc §7: this is the single place that applies `config/DIT_HAP.mplstyle`,
so importing any plotting module gives a consistent look without each caller
re-applying the style.

Usage
-----
    from workflow.src.plotting.style import AX_WIDTH, AX_HEIGHT, COLORS
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Third-party Imports
import matplotlib.pyplot as plt

# =============================================================================
# STYLE APPLICATION
# =============================================================================
# Resolve config/DIT_HAP.mplstyle relative to the repository root (parents[3] =
# workflow/src/plotting/style.py -> repo root).
_STYLE_PATH = Path(__file__).resolve().parents[3] / "config" / "DIT_HAP.mplstyle"
if _STYLE_PATH.exists():
    plt.style.use(str(_STYLE_PATH))

# Derived constants, read once after the style is applied.
AX_WIDTH, AX_HEIGHT = plt.rcParams["figure.figsize"]
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# =============================================================================
# DELETION-LIBRARY PHENOTYPE COLOR MAPS
# =============================================================================
# Byte-faithful to compare_with_deletion_library.ipynb's cell 1 (DIT_HAP_pipeline).
# Two independent palettes for the same phenotype labels, indexed into the same
# COLORS cycle: CATEGORY_COLOR_MAP groups labels into "healthy-looking" (index 2)
# vs "phenotype-affected" (index 0) buckets for boxplot/violin coloring, while
# DONUT_COLOR_MAP gives each category its own distinct wedge color for the donut
# charts. Kept here (not in the verification script) so any plotting module can
# reuse them, matching how COLORS/AX_WIDTH/AX_HEIGHT are already shared.
CATEGORY_COLOR_MAP = {
    "WT": COLORS[2],
    "small colonies": COLORS[2],
    "very small colonies": COLORS[2],
    "E": COLORS[0],
    "E (tiny colonies)": COLORS[0],
    "microcolonies": COLORS[0],
    "germinated": COLORS[0],
    "spores": COLORS[0],
    "Not verified": COLORS[-1],
}

DONUT_COLOR_MAP = {
    "spores": COLORS[7],
    "germinated": COLORS[0],
    "microcolonies": COLORS[-1],
    "E": COLORS[5],
    "E (tiny colonies)": COLORS[5],
    "very small colonies": COLORS[2],
    "small colonies": COLORS[4],
    "WT": COLORS[-4],
}
