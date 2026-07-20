"""Tests for spike-in analysis core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from workflow.scripts.spikein.run_spikein_analysis import (
    assign_ratio_by_order,
    build_spike_sites_df,
    compute_linear_regression_stats,
    SPIKE_IN_RATIO,
)


SPIKE_IN_RATIO_EXPECTED = np.array([1.5, 4, 16, 64, 256, 1024]) / 100000


def test_spike_in_ratio_constant():
    """SPIKE_IN_RATIO constant matches known dilution series."""
    np.testing.assert_allclose(SPIKE_IN_RATIO, SPIKE_IN_RATIO_EXPECTED)


def test_assign_ratio_by_order_basic():
    """Reads rank ascending → lowest read = lowest ratio; relative values are log2-normalised."""
    sub = pd.DataFrame(
        {"Reads": [100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0]},
        index=range(6),
    )
    result = assign_ratio_by_order(sub.copy(), SPIKE_IN_RATIO_EXPECTED)
    # Ratio assigned by rank (0-indexed ascending)
    np.testing.assert_allclose(result["Ratio"].values, SPIKE_IN_RATIO_EXPECTED)
    # Minimum read subtracted: lowest becomes 0
    assert result["Reads"].min() == 0.0
    # Relative_Dilution_Ratio: log2(ratio / max_ratio)
    expected_rel_dil = np.log2(SPIKE_IN_RATIO_EXPECTED / SPIKE_IN_RATIO_EXPECTED.max())
    np.testing.assert_allclose(result["Relative_Dilution_Ratio"].values, expected_rel_dil)


def test_assign_ratio_by_order_monotone_reads():
    """With perfectly ordered reads, the rank assignment is identity."""
    reads = np.array([10.0, 40.0, 160.0, 640.0, 2560.0, 10240.0])
    sub = pd.DataFrame({"Reads": reads}, index=range(6))
    result = assign_ratio_by_order(sub.copy(), SPIKE_IN_RATIO_EXPECTED)
    np.testing.assert_allclose(result["Ratio"].values, SPIKE_IN_RATIO_EXPECTED)


def test_build_spike_sites_df_shape():
    """build_spike_sites_df returns one row per spike-in site with expected columns."""
    mock_index = pd.MultiIndex.from_tuples(
        [("I", 3749394, "-"), ("II", 3344505, "-"), ("II", 185161, "-"),
         ("II", 1157130, "-"), ("II", 3065244, "-")],
        names=["Chr", "Coordinate", "Strand"],
    )
    reads = np.array([10.0, 20.0, 40.0, 80.0, 160.0, 320.0])
    mock_df = pd.DataFrame(
        [reads for _ in range(5)],
        index=mock_index,
        columns=pd.MultiIndex.from_tuples([("S1",), ("S2",), ("S3",), ("S4",), ("S5",), ("S6",)]),
    )
    spike_in_sites = {
        "DY215": {"chr": "I", "coord": 3749394, "strand": "-"},
        "DY217": {"chr": "II", "coord": 3344505, "strand": "-"},
        "DY218": {"chr": "II", "coord": 185161, "strand": "-"},
        "DY339": {"chr": "II", "coord": 1157130, "strand": "-"},
        "DY348": {"chr": "II", "coord": 3065244, "strand": "-"},
    }
    df = build_spike_sites_df(mock_df, spike_in_sites)
    assert len(df) == 5
    assert "Strain" in df.columns


def test_compute_linear_regression_stats():
    """compute_linear_regression_stats returns slope, r_value, p_value, r2."""
    x = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0])
    y = 0.95 * x + 0.1  # near-perfect linear
    stats = compute_linear_regression_stats(pd.Series(x), pd.Series(y))
    assert abs(stats["slope"] - 0.95) < 0.01
    assert stats["r2"] > 0.99
    assert "p_value" in stats
