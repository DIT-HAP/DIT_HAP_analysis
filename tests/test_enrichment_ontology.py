"""Tests for workflow/src/enrichment/ontology.py — OBO/GAF loading."""

import pytest
from pathlib import Path
from workflow.src.enrichment.ontology import OntologyDataConfig

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")
ONTOLOGY_DIR = POMBASE_DIR / "ontologies_and_associations"


def test_validate_paths_raises_on_missing_obo(tmp_path):
    """validate_paths raises FileNotFoundError naming the missing OBO file."""
    cfg = OntologyDataConfig(
        ontology_obo=tmp_path / "missing.obo",
        ontology_association_gaf=tmp_path / "missing.gaf",
        slim_terms_table=[],
    )
    with pytest.raises(FileNotFoundError, match="Gene ontology file not found"):
        cfg.validate_paths()


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_load_data_concatenates_slim_tables():
    """load_data() concatenates all three GO slim tables into one Term/Description frame."""
    cfg = OntologyDataConfig(
        ontology_obo=ONTOLOGY_DIR / "go-basic.obo",
        ontology_association_gaf=ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            ONTOLOGY_DIR / "bp_go_slim_terms.tsv",
            ONTOLOGY_DIR / "mf_go_slim_terms.tsv",
            ONTOLOGY_DIR / "cc_go_slim_terms.tsv",
        ],
    )
    data = cfg.load_data()
    assert list(data.slim_term_dataframe.columns) == ["Term", "Description"]
    assert len(data.slim_term_dataframe) > 0


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_load_ontology_data_returns_gene2go_for_known_gene():
    """load_ontology_data's gene2go dict has an entry for a known coding gene."""
    from workflow.src.enrichment.ontology import load_ontology_data

    cfg = OntologyDataConfig(
        ontology_obo=ONTOLOGY_DIR / "go-basic.obo",
        ontology_association_gaf=ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[ONTOLOGY_DIR / "bp_go_slim_terms.tsv"],
    )
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(
        cfg.load_data(),
        relationships={"is_a", "part_of"},
        propagate_counts=True,
        load_obsolete=False,
        prt=None,
    )
    assert "SPAC1002.01" in gene2go
    assert len(gene2go["SPAC1002.01"]) > 0
