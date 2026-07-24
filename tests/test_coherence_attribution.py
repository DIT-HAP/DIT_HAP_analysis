"""Tests for the incoherence-attribution diagnostics (workflow/src/coherence/attribution.py).

Pins the behaviour attribute_incoherence.py relies on: the GMM major/minor split
on a synthetic core+minor cloud, the group_id-keyed shared-subunit fraction, the
paralog fraction, and the attribution label priority ladder (including the
CLRC-like split+shared -> conditional_module case).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.coherence.attribution import (
    major_minor_split,
    shared_subunits,
    shared_fraction,
    paralog_fraction,
    attribute_incoherence,
)


# --- GMM major/minor split --------------------------------------------------
def test_major_minor_split_detects_core_plus_minority():
    """A tight core + a far-flung minority is called a genuine split, core = tighter comp."""
    rng = np.random.default_rng(0)
    core = rng.normal(0.0, 0.01, size=(8, 2))          # tight essential core
    minor = rng.normal(1.0, 0.02, size=(3, 2))         # dispersed minority, far away
    X = np.vstack([core, minor])
    res = major_minor_split(X)
    assert res["is_split"] is True
    assert res["silhouette"] > 0.5
    # core = the tighter component -> its size is the 8-point cluster
    assert res["component_sizes"][res["core_label"]] == 8


def test_major_minor_split_too_few_points_is_data_limited():
    """Below MIN_N_FOR_GMM there is no split; reason marks it n<... (data-limited)."""
    X = np.array([[0.0, 0.0], [0.1, 0.1], [0.2, 0.0]])  # n=3 < 6
    res = major_minor_split(X)
    assert res["is_split"] is False
    assert res["reason"].startswith("n<")


def test_major_minor_split_uniform_cloud_not_split():
    """A single diffuse blob has no clean 2-subgroup structure -> not a split."""
    rng = np.random.default_rng(1)
    X = rng.normal(0.5, 0.3, size=(20, 2))
    res = major_minor_split(X)
    assert res["is_split"] is False


# --- shared-subunit fraction (group_id keyed) -------------------------------
def _long(rows):
    return pd.DataFrame(rows, columns=["group_id", "group_name", "Systematic ID"])


def test_shared_subunits_lists_other_groups():
    """A member in >1 group is reported with the OTHER groups it belongs to."""
    long = _long([
        ("C1", "complex one", "g1"), ("C1", "complex one", "g2"), ("C1", "complex one", "g3"),
        ("C2", "complex two", "g1"),   # g1 shared with C2
        ("C3", "complex three", "g1"), # g1 shared with C3 too
        ("C2", "complex two", "g9"),
    ])
    ss = shared_subunits(long, "C1")
    row = ss[ss["Systematic ID"] == "g1"].iloc[0]
    assert row["n_other_groups"] == 2
    assert "complex two" in row["other_groups"] and "complex three" in row["other_groups"]
    # g2, g3 are not shared -> not in the table
    assert set(ss["Systematic ID"]) == {"g1"}


def test_shared_fraction_counts_shared_members():
    """shared_fraction = (#members shared with >=1 other group) / (#members)."""
    long = _long([
        ("C1", "one", "g1"), ("C1", "one", "g2"), ("C1", "one", "g3"), ("C1", "one", "g4"),
        ("C2", "two", "g1"), ("C2", "two", "g2"),  # g1, g2 shared -> 2/4
    ])
    assert shared_fraction(long, "C1") == pytest.approx(0.5)


def test_shared_fraction_empty_group_is_nan():
    long = _long([("C1", "one", "g1")])
    assert np.isnan(shared_fraction(long, "NOPE"))


# --- paralog fraction -------------------------------------------------------
def test_paralog_fraction():
    assert paralog_fraction(["g1", "g2", "g3", "g4"], {"g1", "g3"}) == pytest.approx(0.5)
    assert np.isnan(paralog_fraction([], {"g1"}))
    assert paralog_fraction(["g1", "g2"], set()) == 0.0


# --- attribution label ladder -----------------------------------------------
def test_attribute_split_plus_shared_is_conditional_module():
    """A CLRC-like group (genuine split AND cross-shared) -> conditional_module (top priority)."""
    split = {"is_split": True, "reason": "gmm_2comp"}
    label = attribute_incoherence(split, shared_frac=0.8, paralog_frac=0.1)
    assert label == "conditional_module"


def test_attribute_split_only_is_major_minor():
    split = {"is_split": True, "reason": "gmm_2comp"}
    assert attribute_incoherence(split, shared_frac=0.1, paralog_frac=0.1) == "major_minor_split"


def test_attribute_shared_only_is_shared_subunits():
    split = {"is_split": False, "reason": "gmm_2comp"}
    assert attribute_incoherence(split, shared_frac=0.9, paralog_frac=0.1) == "shared_subunits"


def test_attribute_paralog_only_is_paralog_buffered():
    split = {"is_split": False, "reason": "gmm_2comp"}
    assert attribute_incoherence(split, shared_frac=0.1, paralog_frac=0.8) == "paralog_buffered"


def test_attribute_small_group_is_data_limited():
    split = {"is_split": False, "reason": "n<6"}
    assert attribute_incoherence(split, shared_frac=0.1, paralog_frac=0.1) == "data_limited"


def test_attribute_no_signal_is_intrinsic():
    split = {"is_split": False, "reason": "gmm_2comp"}
    assert attribute_incoherence(split, shared_frac=0.1, paralog_frac=0.1) == "intrinsic_heterogeneity"


def test_attribute_priority_shared_over_paralog():
    """When both shared and paralog fire (no split), shared-subunit wins (higher priority)."""
    split = {"is_split": False, "reason": "gmm_2comp"}
    assert attribute_incoherence(split, shared_frac=0.9, paralog_frac=0.9) == "shared_subunits"
