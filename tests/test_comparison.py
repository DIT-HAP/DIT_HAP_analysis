"""Tests for large-scale study comparison core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.comparison.compare_large_scale_studies import (
    clip_density_columns,
    compute_pearson_r,
    CLIP_UPPER,
    DENSITY_COLUMNS,
)


def test_clip_upper_constant():
    """clip(upper=200) is the exact value from source notebook."""
    assert CLIP_UPPER == 200


def test_density_columns_include_required_names():
    """Integration density, ipkm, uipkm columns must be clipped."""
    required = {"Integration density, in-vivo (integrations/kb/million inserts)", "ipkm", "uipkm"}
    assert required.issubset(set(DENSITY_COLUMNS))


def test_clip_density_columns_caps_at_200():
    """Values above 200 are clipped to exactly 200."""
    df = pd.DataFrame({
        "Integration density, in-vivo (integrations/kb/million inserts)": [50.0, 250.0, 200.0],
        "ipkm": [100.0, 300.0, 199.0],
        "uipkm": [10.0, 201.0, 5.0],
        "other_col": [1000.0, 2000.0, 3000.0],
    })
    result = clip_density_columns(df)
    assert result["Integration density, in-vivo (integrations/kb/million inserts)"].max() == 200.0
    assert result["ipkm"].max() == 200.0
    assert result["uipkm"].max() == 200.0
    # other_col untouched
    assert result["other_col"].max() == 3000.0


def test_compute_pearson_r_returns_r_and_pvalue():
    """compute_pearson_r returns (r, p_value) for two numeric series."""
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([1.1, 1.9, 3.1, 3.9, 5.1])
    r, p = compute_pearson_r(x, y)
    assert abs(r - 1.0) < 0.05
    assert p < 0.05


def test_compute_pearson_r_ignores_nan_pairs():
    """NaN in either column → drop pair before correlation."""
    x = pd.Series([1.0, 2.0, np.nan, 4.0])
    y = pd.Series([1.0, 2.0, 3.0, np.nan])
    r, p = compute_pearson_r(x, y)
    # Only (1,1) and (2,2) survive — perfect positive correlation
    assert abs(r - 1.0) < 0.01
