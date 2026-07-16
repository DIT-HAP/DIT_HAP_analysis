"""Tests for workflow/src/enrichment/pipeline.py and the new ontology.py helpers."""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.enrichment.ontology import (
    GeneMetaConfig,
    assign_term_name,
    format_phaf_file,
)
from workflow.src.enrichment.pipeline import (
    get_slim_ns2assoc,
    ontology_enrichment_pipeline,
)

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")
ONTOLOGY_DIR = POMBASE_DIR / "ontologies_and_associations"
GENE_META = POMBASE_DIR / "Gene_metadata" / "gene_IDs_names_products.tsv"
DELETION_XLSX = Path("resources/curated/deletion_library_categories.xlsx")


class _FakeTerm:
    def __init__(self, name):
        self.name = name


class _FakeDag:
    def __init__(self, mapping):
        self._m = mapping

    def __contains__(self, k):
        return k in self._m

    def __getitem__(self, k):
        return _FakeTerm(self._m[k])


def test_assign_term_name_hits_and_misses():
    """assign_term_name returns the DAG name, or a placeholder for unknown terms."""
    dag = _FakeDag({"GO:0001": "mitosis"})
    assert assign_term_name("GO:0001", dag) == "mitosis"
    assert assign_term_name("GO:9999", dag) == "No record for GO:9999"


def test_gene_meta_config_validate_rejects_missing(tmp_path):
    """GeneMetaConfig.validate_paths raises FileNotFoundError on a missing file."""
    cfg = GeneMetaConfig(
        gene_IDs_names_products=tmp_path / "nope.tsv",
        deletion_library_essentiality=tmp_path / "nope.xlsx",
    )
    with pytest.raises(FileNotFoundError, match="Gene metadata file not found"):
        cfg.validate_paths()


@pytest.mark.skipif(not (GENE_META.exists() and DELETION_XLSX.exists()), reason="requires PomBase metadata + curated xlsx")
def test_gene_meta_loads_and_joins_essentiality():
    """GeneMetaData joins deletion-library essentiality and fills gene_name from id."""
    cfg = GeneMetaConfig(gene_IDs_names_products=GENE_META, deletion_library_essentiality=DELETION_XLSX)
    meta = cfg.load_data()
    assert "Gene dispensability. This study" in meta.gene_info_with_essentiality.columns
    assert len(meta.id2name) > 0
    # gene_name has no NaN (filled with systematic id).
    assert meta.gene_info_with_essentiality["gene_name"].notna().all()


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires PomBase ontology files")
def test_format_phaf_file_uses_fixed_date_no_today_stamp():
    """format_phaf_file writes a fixed date header (reproducible), not date.today()."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "phaf_go_style.tsv"
        result = format_phaf_file(
            ONTOLOGY_DIR / "fypo-simple-pombase.obo",
            ONTOLOGY_DIR / "pombase_phenotype_annotation.phaf.tsv",
            out,
        )
        assert result == out
        header = out.read_text().splitlines()[:5]
        assert any("!date-generated: 2025-09-01" in h for h in header)
        # 17-column GO-style GAF body.
        body = pd.read_csv(out, sep="\t", comment="!", header=None)
        assert body.shape[1] == 17


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires PomBase ontology files")
def test_pipeline_runs_and_returns_formatted_tables():
    """ontology_enrichment_pipeline returns full+slim DataFrames and computes slim assoc without error."""
    from workflow.src.enrichment.ontology import OntologyDataConfig

    go_cfg = OntologyDataConfig(
        ontology_obo=ONTOLOGY_DIR / "go-basic.obo",
        ontology_association_gaf=ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            ONTOLOGY_DIR / "bp_go_slim_terms.tsv",
            ONTOLOGY_DIR / "mf_go_slim_terms.tsv",
            ONTOLOGY_DIR / "cc_go_slim_terms.tsv",
        ],
    )
    data = go_cfg.load_data()
    gaf = pd.read_csv(
        ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        sep="\t", comment="!", header=None, usecols=[1], names=["gene"],
    )
    all_genes = gaf["gene"].dropna().unique().tolist()
    bg = all_genes[:600]
    query = all_genes[:30]

    full_df, slim_df, dag, objanno = ontology_enrichment_pipeline(
        data, query, bg,
        enrichment_kwargs={"alpha": 0.05, "methods": ["fdr_bh"], "propagate_counts": True, "relationships": {"is_a", "part_of"}},
        format_kwargs={"itemid2name": None},
    )
    # Both are DataFrames (possibly empty for a small arbitrary gene set); no exception is the key check.
    assert isinstance(full_df, pd.DataFrame)
    assert isinstance(slim_df, pd.DataFrame)
