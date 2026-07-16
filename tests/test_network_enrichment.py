"""Tests for network enrichment: cache logic + formatting (no live API calls)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.enrichment import pipeline
from workflow.src.enrichment.pipeline import (
    _cache_key,
    _cache_load,
    _cache_store,
    format_string_enrichment_results,
    stringdb_enrichment,
    revigo_analysis,
)


def test_cache_key_is_deterministic_and_order_independent_per_arg():
    """Same parts -> same key; different parts -> different key."""
    assert _cache_key("string", "a,b", "c") == _cache_key("string", "a,b", "c")
    assert _cache_key("string", "a,b", "c") != _cache_key("string", "b,a", "c")


def test_cache_store_and_load_roundtrip(tmp_path):
    """A stored DataFrame is returned verbatim on the next load."""
    df = pd.DataFrame({"term_id": ["GO:1"], "p_fdr": [0.01]})
    _cache_store(tmp_path, "k1", df)
    loaded = _cache_load(tmp_path, "k1")
    pd.testing.assert_frame_equal(loaded, df)


def test_cache_load_miss_returns_none(tmp_path):
    """A missing key or disabled cache returns None."""
    assert _cache_load(tmp_path, "absent") is None
    assert _cache_load(None, "any") is None


def test_stringdb_enrichment_uses_cache_without_network(tmp_path):
    """A pre-seeded cache entry short-circuits the API call entirely."""
    query = ["SPAC1", "SPAC2"]
    bg = ["SPAC1", "SPAC2", "SPBC1"]
    key = _cache_key("string", ",".join(sorted(query)), ",".join(sorted(bg)))
    seeded = pd.DataFrame({"term_id": ["GO:X"], "namespace": ["KEGG Pathways"], "p_fdr": [0.001]})
    _cache_store(tmp_path, key, seeded)

    # If this touches the network, the test environment would fail/slow; cache hit avoids it.
    result = stringdb_enrichment(query, bg, cache_dir=tmp_path)
    pd.testing.assert_frame_equal(result, seeded)


def test_revigo_analysis_uses_cache_without_network(tmp_path):
    """A pre-seeded REVIGO cache entry short-circuits the API call."""
    enrich_df = pd.DataFrame({"term_id": ["GO:0001", "GO:0002"], "p_fdr": [0.01, 0.02]})
    key = _cache_key("revigo", "0.7", "GO:0001:0.01,GO:0002:0.02")
    seeded = pd.DataFrame({"Term ID": ["GO:0001"], "Representative": ["mitosis"]})
    _cache_store(tmp_path, key, seeded)

    result = revigo_analysis(enrich_df, cut_off=0.7, cache_dir=tmp_path)
    pd.testing.assert_frame_equal(result, seeded)


def test_format_string_enrichment_maps_namespaces_and_schema():
    """STRING results are renamed to the shared schema with human-readable namespaces."""
    raw = pd.DataFrame(
        {
            "term": ["KEGG:1"],
            "category": ["KEGG"],
            "description": ["Ribosome"],
            "p_value": [0.001],
            "fdr": [0.01],
            "number_of_genes": [5],
            "number_of_genes_in_background": [50],
            "preferredNames": ["a,b,c,d,e"],
        }
    )
    out = format_string_enrichment_results(raw, ["a", "b"], ["a", "b", "c"])
    assert "term_id" in out.columns and "p_fdr" in out.columns
    assert out["namespace"].iloc[0] == "KEGG Pathways"
    assert out["study_n"].iloc[0] == 2
    assert out["pop_n"].iloc[0] == 3
