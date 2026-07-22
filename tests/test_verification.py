"""Tests for deletion library verification logic."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.verification.compare_deletion_library import (
    merge_deletion_library,
    compute_category_stats,
    canonicalize_category,
    build_final_merged,
    prepare_verification_data,
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


# =============================================================================
# CRITICAL-GENE ANALYSIS
# =============================================================================
def test_canonicalize_category_folds_schema_drift():
    """WT-like → WT and compound labels fold to their canonical bucket."""
    merged = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4"],
        "Category": ["WT-like", "spores, germinated", "microcolonies, small colonies", "spores"],
    })
    out = canonicalize_category(merged)
    assert out["cat_canon"].tolist() == ["WT", "spores", "small colonies", "spores"]


def _make_critical_fixtures():
    """merged (with cat_canon), full verification, and simplified verification for bucket tests."""
    merged = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4", "g5"],
        "Name": ["a", "b", "c", "d", "e"],
        "DR": [0.9, 0.8, 0.7, 0.5, 0.4],
        "DL": [5.0, 4.0, 3.0, 2.0, 1.0],
        "DeletionLibrary_essentiality": ["V", "V", "V", "V", "V"],
        "Category": ["WT-like", "WT-like", "WT-like", "WT-like", "small colonies"],
    })
    merged = canonicalize_category(merged)
    simplified = pd.DataFrame({
        "Systematic ID": ["g1", "g2"],
        "Verification result": ["E", "small colonies"],
        "Verified essentiality": ["E", "V"],
    })
    verification_full = pd.DataFrame({
        "systematic_id": ["g1", "g2", "g3"],
        "verification_phenotype": ["E", "small colonies", "WT"],
        "verification_essentiality": ["E", "V", "V"],
        "median_area_day3": [0.1, 0.5, 1.0],
    })
    return merged, simplified, verification_full


def test_prepare_verification_data_buckets():
    """Outliers split into verified categories + a Not verified bucket, DR values collected."""
    merged, simplified, verification_full = _make_critical_fixtures()
    final_merged = build_final_merged(merged, verification_full)
    dr_dict, detail = prepare_verification_data(
        merged, final_merged, simplified,
        outlier_filter="cat_canon == 'WT' and DR > 0.35",
    )
    # g1..g4 are WT outliers (DR>0.35); g1 verified E, g2 verified small colonies,
    # g3+g4 unverified. Not-verified DR list is DR-sorted descending.
    assert dr_dict["E"] == [0.9]
    assert dr_dict["small colonies"] == [0.8]
    assert dr_dict["Not verified"] == [0.7, 0.5]
    assert set(detail["Verification result bucket"]) == {"E", "small colonies", "Not verified"}


def test_prepare_verification_data_empty_group():
    """A zero-hit filter returns an empty dict and an empty detail frame without raising."""
    merged, simplified, verification_full = _make_critical_fixtures()
    final_merged = build_final_merged(merged, verification_full)
    dr_dict, detail = prepare_verification_data(
        merged, final_merged, simplified,
        outlier_filter="cat_canon == 'WT' and DR > 100",
    )
    assert dr_dict == {}
    assert len(detail) == 0


def test_build_final_merged_zero_fills_essential_missing_area():
    """An E-verified gene missing day-3 area is zero-filled across area columns."""
    merged = pd.DataFrame({
        "Systematic ID": ["g1"],
        "Name": ["a"],
        "DR": [0.9], "DL": [5.0],
        "DeletionLibrary_essentiality": ["E"],
        "Category": ["spores"],
    })
    verification_full = pd.DataFrame({
        "systematic_id": ["g1"],
        "verification_essentiality": ["E"],
        "median_area_day3": [None],
        "median_area_day6": [None],
    })
    final = build_final_merged(merged, verification_full)
    assert final.loc[0, "median_area_day3"] == 0
    assert final.loc[0, "median_area_day6"] == 0
