#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for plot_group_scatter.py
================================

Tests the resolve_groups function that matches group names/ids from the config
namelist against the long-table annotation.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd


def _long():
    return pd.DataFrame({
        "source": ["go_cc"] * 3,
        "group_id": ["GO:1", "GO:1", "GO:2"],
        "group_name": ["alpha", "alpha", "beta"],
        "Systematic ID": ["gA", "gB", "gC"],
        "Name": ["gA", "gB", "gC"],
        "n_group_genes": [2, 2, 1],
    })


def test_resolve_groups_by_name():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["alpha"])
    assert len(out) == 1
    gid, gname, genes = out[0]
    assert gid == "GO:1" and gname == "alpha" and set(genes) == {"gA", "gB"}


def test_resolve_groups_by_id():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["GO:2"])
    assert out[0][0] == "GO:2"


def test_resolve_groups_skips_unknown():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["nonexistent"])
    assert out == []


def test_resolve_groups_same_name_multiple_ids():
    """When a name matches multiple group_ids, all are resolved."""
    import pandas as pd
    long = pd.DataFrame({
        "source": ["go_cc"] * 4,
        "group_id": ["GO:1", "GO:1", "GO:3", "GO:3"],
        "group_name": ["alpha", "alpha", "alpha", "alpha"],
        "Systematic ID": ["gA", "gB", "gC", "gD"],
        "Name": ["gA", "gB", "gC", "gD"],
        "n_group_genes": [2, 2, 2, 2],
    })
    from plot_group_scatter import resolve_groups
    out = resolve_groups(long, "go_cc", ["alpha"])
    assert len(out) == 2  # GO:1 and GO:3
    assert {g[0] for g in out} == {"GO:1", "GO:3"}


def test_resolve_groups_drops_nan_ids():
    """Members with NaN Systematic ID are silently dropped."""
    import numpy as np
    import pandas as pd
    long = pd.DataFrame({
        "source": ["go_cc", "go_cc"],
        "group_id": ["GO:1", "GO:1"],
        "group_name": ["alpha", "alpha"],
        "Systematic ID": ["gA", np.nan],
        "Name": ["gA", "gB"],
        "n_group_genes": [2, 2],
    })
    from plot_group_scatter import resolve_groups
    out = resolve_groups(long, "go_cc", ["alpha"])
    assert out[0][2] == ["gA"]  # only gA, gB (NaN id) dropped
