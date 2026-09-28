"""Tests for the incoherence-attribution diagnostics (workflow/src/coherence/attribution.py).

Pins the behaviour compute_incoherence_attribution.py relies on: the GMM
major/minor split on a synthetic core+minor cloud, the source-scoped
(source, group_id)-keyed shared-subunit fraction, the paralog fraction, and the
attribution label priority ladder (including the CLRC-like split+shared ->
conditional_module case).

The (source, group_id) key is pinned separately below: it is what makes the
pooled de-duplicated attribution correct, because group_id is NOT unique across
sources (173 of them, e.g. GO:0032040, appear in both go_cc and go_macrocomplex).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd
import pytest

from workflow.src.coherence.attribution import (
    major_minor_split,
    shared_subunits,
    shared_subunit_fractions,
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
    """shared fraction = (#members shared with >=1 other group) / (#members)."""
    long = _long([
        ("C1", "one", "g1"), ("C1", "one", "g2"), ("C1", "one", "g3"), ("C1", "one", "g4"),
        ("C2", "two", "g1"), ("C2", "two", "g2"),  # g1, g2 shared -> 2/4
    ])
    fracs = shared_subunit_fractions(long)
    assert fracs["C1"] == pytest.approx(0.5)
    assert fracs["C2"] == pytest.approx(1.0)  # both of C2's members also in C1


def test_shared_fraction_group_with_no_shared_member_is_zero():
    long = _long([
        ("C1", "one", "g1"), ("C1", "one", "g2"),
        ("C2", "two", "g3"),  # disjoint
    ])
    assert shared_subunit_fractions(long)["C1"] == 0.0


def test_shared_subunits_keys_on_group_id_not_name():
    """Two distinct term IDs that share a NAME are still two different groups.

    "Other groups" used to be derived by subtracting the focal group's names, so
    a member of a same-named sibling term read as unshared. GO has real
    same-name/different-ID pairs, and plot_coherence.py's own copy of this metric
    was already group_id-keyed — the two disagreed.
    """
    long = _long([
        ("C1", "one", "g1"), ("C1", "one", "g2"),
        ("C9", "one", "g1"),  # different id, same name
    ])
    ss = shared_subunits(long, "C1")
    assert set(ss["Systematic ID"]) == {"g1"}
    assert ss.iloc[0]["n_other_groups"] == 1
    assert shared_subunit_fractions(long)["C1"] == pytest.approx(0.5)


# --- pooled (all-sources) tables --------------------------------------------
# group_id is not unique across sources, so on the pooled de-duplicated table
# every lookup has to carry the source with it. `_long_pooled` mirrors the real
# long table's column order (source first), which is what the code tests for.
def _long_pooled(rows):
    return pd.DataFrame(rows, columns=["source", "group_id", "group_name", "Systematic ID"])


def test_shared_subunits_scopes_to_the_given_source():
    """The same group_id in two sources is two different groups' member lists."""
    long = _long_pooled([
        ("go_cc", "GO:1", "cc view", "gA"),
        ("go_bp", "GO:1", "bp view", "gB"), ("go_bp", "GO:1", "bp view", "gC"),
        ("go_bp", "GO:2", "other", "gB"),
    ])
    assert set(shared_subunits(long, "GO:1", source="go_cc")["Systematic ID"]) == set()
    assert set(shared_subunits(long, "GO:1", source="go_bp")["Systematic ID"]) == {"gB"}
    # Without the source it falls back to inferring it, which on a collision picks
    # whichever row comes first — the reason the pooled caller must pass it.
    assert set(shared_subunits(long, "GO:1")["Systematic ID"]) == set()


def test_group_member_points_scopes_to_the_given_source():
    """A pooled long table must not union two sources' members under one group_id."""
    from compute_incoherence_attribution import group_member_points

    long = _long_pooled([
        ("go_cc", "GO:1", "cc view", "gA"),
        ("go_bp", "GO:1", "bp view", "gB"), ("go_bp", "GO:1", "bp view", "gC"),
    ])
    points = pd.DataFrame(
        {"norm_DR": [-1.0, -0.5, 0.0], "norm_DL": [0.0, 1.0, 2.0]},
        index=["gA", "gB", "gC"],
    )
    cc_ids, cc_X = group_member_points(long, "GO:1", points, source="go_cc")
    bp_ids, bp_X = group_member_points(long, "GO:1", points, source="go_bp")
    assert cc_ids == ["gA"] and cc_X.shape == (1, 2)
    assert bp_ids == ["gB", "gC"] and bp_X.shape == (2, 2)
    # An unsourced table (every per-source caller) keeps the old behaviour: the
    # group_id is unique there, so omitting the source is not a silent change.
    assert group_member_points(long.drop(columns="source"), "GO:1", points)[0] == ["gA", "gB", "gC"]


def test_attribute_all_pooled_matches_per_source():
    """The pooled run answers each source with its own members, not the union.

    This is the whole reason the dedup attribution is allowed to run once over a
    pooled table instead of five times per source: a group_id's GMM, shared
    fraction and label must come out identical either way.
    """
    from compute_incoherence_attribution import AttributionConfig, attribute_all

    long = _long_pooled([
        ("go_cc", "GO:1", "cc view", "gA"),
        ("go_bp", "GO:1", "bp view", "gB"), ("go_bp", "GO:1", "bp view", "gC"),
        ("go_bp", "GO:2", "other", "gB"),
    ])
    metrics = pd.DataFrame({
        "source": ["go_cc", "go_bp", "go_bp"],
        "group_id": ["GO:1", "GO:1", "GO:2"],
        "group_name": ["cc view", "bp view", "other"],
        "n_scored_members": [1, 2, 1],
        "median_pairwise_distance_z": [2.5, 2.5, -1.0],
    })
    points = pd.DataFrame(
        {"norm_DR": [-1.0, -0.5, 0.0], "norm_DL": [0.0, 1.0, 2.0]},
        index=["gA", "gB", "gC"],
    )
    config = AttributionConfig(
        metrics=Path("unused.parquet"), annotations=(Path("unused.tsv"),),
        fitting_results=Path("unused.tsv"), paralogs=Path("unused.tsv"),
        output_table=Path("unused.tsv"), output_points=Path("unused.parquet"),
    )
    table, split_points = attribute_all(metrics, long, points, set(), config)

    by_key = {(row.source, row.group_id): row for row in table.itertuples()}
    # gA alone in go_cc -> nothing shared; gB shared with GO:2 within go_bp -> 0.5.
    assert by_key[("go_cc", "GO:1")].frac_shared_members == pytest.approx(0.0)
    assert by_key[("go_bp", "GO:1")].frac_shared_members == pytest.approx(0.5)
    # The split points carry the source, which is what lets the figure key on it.
    assert set(split_points.columns) >= {"source", "group_id"}
    assert set(map(tuple, split_points[["source", "group_id"]].drop_duplicates().values)) == {
        ("go_cc", "GO:1"), ("go_bp", "GO:1"),
    }


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
