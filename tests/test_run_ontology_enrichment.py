"""Tests for workflow/scripts/enrichment/run_ontology_enrichment.py (fast, no goatools run)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.enrichment.run_ontology_enrichment import (
    EnrichmentConfig,
    _concat_by_cluster,
    write_gene_lists,
)


def test_config_validate_rejects_missing_input(tmp_path):
    """validate() raises ValueError when a required input is absent."""
    cfg = EnrichmentConfig(
        final_clusters=tmp_path / "nope.tsv",
        pombase_dir=tmp_path / "nope_dir",
        deletion_library_xlsx=tmp_path / "nope.xlsx",
        output_dir=tmp_path / "out",
        intermediate_dir=tmp_path / "gaf",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_ontology_dir_property(tmp_path):
    """ontology_dir resolves under the PomBase version directory."""
    cfg = EnrichmentConfig(
        final_clusters=tmp_path / "f.tsv",
        pombase_dir=tmp_path / "pombase" / "2025-10-01",
        deletion_library_xlsx=tmp_path / "d.xlsx",
        output_dir=tmp_path / "out",
        intermediate_dir=tmp_path / "gaf",
    )
    assert cfg.ontology_dir == tmp_path / "pombase" / "2025-10-01" / "ontologies_and_associations"


def test_write_gene_lists_emits_all_files(tmp_path):
    """write_gene_lists writes the all-genes file, per-cluster lists, and the matrix."""
    cluster_genes = {1: ["SPAC1", "SPAC2"], 2: ["SPBC1"]}
    bg = ["SPAC1", "SPAC2", "SPBC1"]
    write_gene_lists(cluster_genes, bg, tmp_path)

    assert (tmp_path / "DIT_HAP_all_genes.txt").read_text().split() == bg
    assert (tmp_path / "DIT_HAP_cluster_1_genes.txt").read_text().split() == ["SPAC1", "SPAC2"]
    assert (tmp_path / "DIT_HAP_cluster_2_genes.txt").read_text().split() == ["SPBC1"]
    assert (tmp_path / "DIT_HAP_cluster_genes_matrix.txt").exists()


def test_concat_by_cluster_adds_cluster_column():
    """_concat_by_cluster stacks per-cluster frames and adds a Cluster column."""
    by_cluster = {
        1: pd.DataFrame({"term_id": ["GO:1"], "p_fdr": [0.01]}),
        2: pd.DataFrame({"term_id": ["GO:2"], "p_fdr": [0.02]}),
    }
    out = _concat_by_cluster(by_cluster)
    assert "Cluster" in out.columns
    assert set(out["Cluster"]) == {1, 2}


def test_concat_by_cluster_empty_safe():
    """_concat_by_cluster returns an empty frame when every cluster result is empty."""
    assert _concat_by_cluster({1: pd.DataFrame(), 2: None}).empty
