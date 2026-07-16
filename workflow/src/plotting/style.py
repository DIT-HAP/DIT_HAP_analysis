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
