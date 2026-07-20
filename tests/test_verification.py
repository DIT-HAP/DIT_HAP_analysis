"""Tests for deletion library verification logic."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.verification.compare_deletion_library import (
    merge_deletion_library,
    compute_category_stats,
    CATEGORY_COLOR_MAP,
    DONUT_COLOR_MAP,
)


def _make_gene_results():
    return pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4", "g5"],
        "DR": [0.9, 0.1, 0.5, 0.8, 0.2],
        "DL": [5.0, 1.0, 3.0, 4.0, 2.0],
        "DeletionLibrary_essentiality": ["E", "V", "E", "V", "E"],
    })


def _make_deletion_library():
    return pd.DataFrame({
        "Updated_Systematic_ID": ["g1", "g2", "g3", "g4", "g5"],
        "Category": ["spores", "WT", "germinated", "small colonies", "E"],
    })


def test_category_color_map_has_required_keys():
    """CATEGORY_COLOR_MAP must cover all expected deletion phenotype labels."""
    required = {"WT", "small colonies", "very small colonies", "E",
                "E (tiny colonies)", "microcolonies", "germinated", "spores", "Not verified"}
    assert required.issubset(set(CATEGORY_COLOR_MAP.keys()))


def test_donut_color_map_has_required_keys():
    """DONUT_COLOR_MAP must cover all expected donut chart categories."""
    required = {"spores", "germinated", "microcolonies", "E",
                "E (tiny colonies)", "very small colonies", "small colonies", "WT"}
    assert required.issubset(set(DONUT_COLOR_MAP.keys()))


def test_merge_deletion_library_joins_on_systematic_id():
    """merge_deletion_library joins on Systematic ID / Updated_Systematic_ID."""
    gene = _make_gene_results()
    dl = _make_deletion_library()
    merged = merge_deletion_library(gene, dl)
    assert "Category" in merged.columns
    assert len(merged) == 5


def test_compute_category_stats_returns_counts():
    """compute_category_stats returns count per category in the merged frame."""
    gene = _make_gene_results()
    dl = _make_deletion_library()
    merged = merge_deletion_library(gene, dl)
    stats = compute_category_stats(merged)
    assert "category" in stats.columns
    assert "count" in stats.columns
    assert stats["count"].sum() == 5


def test_category_with_essentiality_flag():
    """small colonies + E essentiality → 'small colonies (E)' label."""
    from workflow.scripts.verification.compare_deletion_library import apply_category_with_essentiality
    row_sc_e = pd.Series({"Category": "small colonies", "DeletionLibrary_essentiality": "E"})
    row_sc_v = pd.Series({"Category": "small colonies", "DeletionLibrary_essentiality": "V"})
    assert apply_category_with_essentiality(row_sc_e) == "small colonies (E)"
    assert apply_category_with_essentiality(row_sc_v) == "small colonies"
