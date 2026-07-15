"""Tests for workflow/src/gene_ids.py systematic ID resolution."""

import math
import pandas as pd
from workflow.src.gene_ids import update_sysIDs


def _write_gene_meta(tmp_path):
    """Write a minimal gene_IDs_names_products.tsv fixture."""
    f = tmp_path / "gene_IDs_names_products.tsv"
    df = pd.DataFrame({
        "gene_systematic_id": ["SPBC11B10.09", "SPAC1002.01"],
        "gene_name": ["cdc2", None],
        "synonyms": ["cdk1,cdc28", ""],
        "gene_type": ["protein coding gene", "protein coding gene"],
    })
    df.to_csv(f, sep="\t", index=False)
    return f


def test_already_current_sysid_passes_through(tmp_path):
    """A gene already given as a current systematic ID is returned unchanged."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["SPBC11B10.09"], meta)
    assert result == ["SPBC11B10.09"]


def test_gene_name_resolves_to_sysid(tmp_path):
    """A gene given by its common name resolves to the systematic ID."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["cdc2"], meta)
    assert result == ["SPBC11B10.09"]


def test_synonym_resolves_to_sysid(tmp_path):
    """A gene given by a synonym resolves to the systematic ID."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["cdk1"], meta)
    assert result == ["SPBC11B10.09"]


def test_unknown_gene_passes_through_unchanged(tmp_path):
    """A gene not found anywhere is returned unchanged (for manual review)."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["not_a_real_gene"], meta)
    assert result == ["not_a_real_gene"]


def test_na_input_passes_through(tmp_path):
    """A NaN input is passed through as NaN, not resolved."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs([float("nan")], meta)
    assert len(result) == 1 and math.isnan(result[0])


def test_dotted_transcript_id_is_case_normalized(tmp_path):
    """A lowercase dotted transcript ID is normalized to GENE.transcript case before lookup."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["spac1002.01"], meta)
    assert result == ["SPAC1002.01"]
