"""Contract tests for coherence source adapters (unified long-table schema)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.coherence.sources import (
    LONG_TABLE_COLUMNS,
    load_macrocomplex,
)

LONG_TABLE_COLUMNS_EXPECTED = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]


def _write_macrocomplex(tmp_path: Path) -> Path:
    """A tiny macromolecular_complex_annotation.tsv: 2 complexes, one member with no symbol."""
    df = pd.DataFrame(
        {
            "complex_term_id": ["GO:0001", "GO:0001", "GO:0002", "GO:0002"],
            "GO_term_name": ["alpha complex", "alpha complex", "beta complex", "beta complex"],
            "systematic_id": ["SPAC1", "SPAC2", "SPBC1", "SPBC2"],
            "symbol": ["gene1", None, "gene3", "gene4"],
        }
    )
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    path = tmp_path / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    df.to_csv(path, sep="\t", index=False)
    return tmp_path


def test_long_table_columns_constant():
    assert list(LONG_TABLE_COLUMNS) == LONG_TABLE_COLUMNS_EXPECTED


def test_macrocomplex_returns_contract_columns(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["source"]) == {"go_macrocomplex"}


def test_macrocomplex_fills_missing_name_with_systematic_id(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    row = out[out["Systematic ID"] == "SPAC2"].iloc[0]
    assert row["Name"] == "SPAC2"  # no symbol -> filled with systematic id


def test_macrocomplex_n_group_genes_is_per_term_total(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    alpha = out[out["group_name"] == "alpha complex"]
    assert (alpha["n_group_genes"] == 2).all()


def _write_go_fixture(tmp_path: Path) -> Path:
    """A minimal go-basic.obo + gene_ontology_annotation.gaf.tsv with 1 CC + 1 BP term."""
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True, exist_ok=True)
    obo = d / "go-basic.obo"
    # GO:0000101 is_a GO:0000100 so a child-only annotation propagates up to the parent.
    obo.write_text(
        "format-version: 1.2\n\n"
        "[Term]\nid: GO:0000100\nname: test cc complex\nnamespace: cellular_component\n\n"
        "[Term]\nid: GO:0000101\nname: test cc subcomplex\nnamespace: cellular_component\nis_a: GO:0000100\n\n"
        "[Term]\nid: GO:0000200\nname: test bp process\nnamespace: biological_process\n\n"
    )
    gaf = d / "gene_ontology_annotation.gaf.tsv"
    # GAF 2.1: 17 tab-separated columns; col2=DB_Object_ID(gene), col5=GO_ID, col9=Aspect.
    lines = ["!gaf-version: 2.1"]
    def row(gene, go, aspect):
        cols = ["PomBase", gene, gene, "", go, "PMID:1", "IDA", "", aspect,
                "", "", "gene", "taxon:4896", "20250101", "PomBase", "", ""]
        return "\t".join(cols)
    lines += [row("SPAC1", "GO:0000100", "C"), row("SPAC2", "GO:0000100", "C"),
              row("SPAC3", "GO:0000101", "C"),  # child-only annotation, should propagate to GO:0000100
              row("SPBC1", "GO:0000200", "P"), row("SPBC2", "GO:0000200", "P")]
    gaf.write_text("\n".join(lines) + "\n")
    return tmp_path


def test_gaf_namespace_cc_only_returns_cc_terms(tmp_path):
    pytest.importorskip("goatools")
    from workflow.src.coherence.sources import load_gaf_namespace
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_gaf_namespace(pombase_dir, "CC")
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["group_id"]) == {"GO:0000100", "GO:0000101"}
    assert set(out["source"]) == {"go_cc"}


def test_gaf_namespace_bp_only_returns_bp_terms(tmp_path):
    pytest.importorskip("goatools")
    from workflow.src.coherence.sources import load_gaf_namespace
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_gaf_namespace(pombase_dir, "BP")
    assert set(out["group_id"]) == {"GO:0000200"}
    assert set(out["source"]) == {"go_bp"}


def test_gaf_namespace_propagates_child_gene_to_parent(tmp_path):
    pytest.importorskip("goatools")
    from workflow.src.coherence.sources import load_gaf_namespace
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_gaf_namespace(pombase_dir, "CC")
    parent_members = set(out.loc[out["group_id"] == "GO:0000100", "Systematic ID"])
    assert "SPAC3" in parent_members  # child annotation propagates up via is_a


def test_source_loaders_registry_has_three_sources():
    from workflow.src.coherence.sources import SOURCE_LOADERS
    assert set(SOURCE_LOADERS) == {"go_macrocomplex", "go_cc", "go_bp"}
