"""Tests for the GO/FYPO/MONDO enrichment modules and their driver configs."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.enrichment.enrich_one_ontology import OntologyConfig
from workflow.scripts.enrichment.finalize_enrichment import FinalizeConfig
from workflow.scripts.enrichment.prepare_genesets import PrepareConfig
from workflow.src.enrichment.cluster_enrichment import (
    CLUSTER_COLUMN,
    ONTOLOGIES,
    _concat_by_cluster,
    filter_go_full,
    load_cluster_genesets,
    resolve_ontology,
    write_gene_lists,
)


def _final_clusters(rows: list[dict]) -> pd.DataFrame:
    """Build a minimal final_clusters frame with Systematic ID + cluster."""
    return pd.DataFrame(rows)


def test_load_cluster_genesets_splits_background_and_nonwt():
    """Background = all genes, nonWT = clusters below wt; keys are int and sorted."""
    clusters = _final_clusters(
        [
            {"Systematic ID": "SPAC1", "cluster": 1},
            {"Systematic ID": "SPAC2", "cluster": 9},
            {"Systematic ID": "SPBC1", "cluster": 2},
        ]
    )
    gs = load_cluster_genesets(clusters, cluster_column="cluster", wt_cluster=9)
    assert set(gs.bg_genes) == {"SPAC1", "SPAC2", "SPBC1"}
    assert set(gs.nonwt_bg_genes) == {"SPAC1", "SPBC1"}  # cluster 9 excluded
    assert list(gs.cluster_genes) == [1, 2, 9]  # int-keyed, sorted


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


def test_filter_go_full_drops_mf_and_high_popcount():
    """filter_go_full keeps pop_count < max and drops the MF namespace, sorted deterministically."""
    go_full = pd.DataFrame(
        {
            "Cluster": [1, 1, 2],
            "namespace": ["BP", "MF", "CC"],
            "term_id": ["GO:2", "GO:9", "GO:1"],
            "pop_count": [100, 50, 500],  # row 3 exceeds max
        }
    )
    out = filter_go_full(go_full, pop_count_max=400)
    assert list(out["namespace"]) == ["BP"]  # MF dropped, high pop_count dropped
    assert out.iloc[0]["term_id"] == "GO:2"


def test_filter_go_full_empty_safe():
    """filter_go_full returns an empty frame for an empty input."""
    assert filter_go_full(pd.DataFrame(), pop_count_max=400).empty


def test_resolve_ontology_rejects_unknown_name(tmp_path):
    """resolve_ontology names the valid ontologies when handed an unknown one."""
    with pytest.raises(ValueError, match="Unknown ontology"):
        resolve_ontology("KEGG", tmp_path, tmp_path, fdr_threshold=0.05)


def test_prepare_config_rejects_missing_input(tmp_path):
    """PrepareConfig.validate raises ValueError when a required input is absent."""
    cfg = PrepareConfig(
        final_clusters=tmp_path / "nope.tsv",
        pombase_dir=tmp_path / "nope_dir",
        deletion_library_xlsx=tmp_path / "nope.xlsx",
        output_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        wt_cluster=9,
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_ontology_config_dir_property(tmp_path):
    """ontology_dir resolves under the PomBase version directory."""
    cfg = OntologyConfig(
        ontology="GO",
        genesets=tmp_path / "g.pkl",
        id2name=tmp_path / "n.pkl",
        pombase_dir=tmp_path / "pombase" / "2025-10-01",
        output_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        wt_cluster=9,
        fdr_threshold=0.05,
    )
    assert cfg.ontology_dir == tmp_path / "pombase" / "2025-10-01" / "ontologies_and_associations"


def test_finalize_config_rejects_missing_frames(tmp_path):
    """FinalizeConfig.validate raises when a per-ontology frame pickle is absent."""
    work = tmp_path / "work"
    work.mkdir()
    cfg = FinalizeConfig(work_dir=work, output_dir=tmp_path / "out", pop_count_max=400)
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_ontologies_constant_order():
    """ONTOLOGIES lists GO, FYPO, MONDO in the canonical processing order."""
    assert ONTOLOGIES == ["GO", "FYPO", "MONDO"]


def test_cluster_column_default_is_cluster():
    """Pin the final-contract column name (guards the production default)."""
    assert CLUSTER_COLUMN == "cluster"
