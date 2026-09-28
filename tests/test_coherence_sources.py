"""Contract tests for coherence source adapters (unified long-table schema)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.coherence.sources import (
    LONG_TABLE_COLUMNS,
    load_gaf_namespace,
    load_gene_names,
    load_macrocomplex,
    load_source,
)

LONG_TABLE_COLUMNS_EXPECTED = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_annotated_members",
]

# systematic id -> gene name, as PomBase's gene_IDs_names_products.tsv carries it.
# Every gene the GAF fixture annotates is in here, so a go_cc / go_bp Name must come
# from this table (and never fall back to the systematic id); SPBC3 is the one gene
# left out, and it is the macrocomplex fixture's systematic-id fallback case.
GENE_NAMES = {"SPAC1": "gene1", "SPAC2": "gene2", "SPAC3": "gene5",
              "SPBC1": "gene3", "SPBC2": "gene4"}


def _write_gene_metadata(pombase_dir: Path, names: dict[str, str]) -> None:
    """Write a minimal Gene_metadata/gene_IDs_names_products.tsv under `pombase_dir`."""
    directory = pombase_dir / "Gene_metadata"
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "gene_systematic_id": list(names),
        "gene_name": [names[key] for key in names],
    }).to_csv(directory / "gene_IDs_names_products.tsv", sep="\t", index=False)


def _write_macrocomplex(tmp_path: Path) -> Path:
    """A tiny macromolecular_complex_annotation.tsv: 2 complexes, symbols on all but SPAC2/SPBC3."""
    df = pd.DataFrame(
        {
            "complex_term_id": ["GO:0001", "GO:0001", "GO:0002", "GO:0002", "GO:0002"],
            "GO_term_name": ["alpha complex", "alpha complex", "beta complex",
                             "beta complex", "beta complex"],
            "systematic_id": ["SPAC1", "SPAC2", "SPBC1", "SPBC2", "SPBC3"],
            "symbol": ["gene1", None, "gene3", "gene4", None],
        }
    )
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    path = d / "macromolecular_complex_annotation.tsv"
    df.to_csv(path, sep="\t", index=False)
    _write_gene_metadata(tmp_path, GENE_NAMES)
    return tmp_path


def test_long_table_columns_constant():
    assert list(LONG_TABLE_COLUMNS) == LONG_TABLE_COLUMNS_EXPECTED


def test_macrocomplex_returns_contract_columns(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_source("go_macrocomplex", pombase_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["source"]) == {"go_macrocomplex"}


def test_macrocomplex_keeps_its_own_symbol(tmp_path):
    """A symbol in the annotation wins over the gene-name table."""
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_source("go_macrocomplex", pombase_dir)
    assert out[out["Systematic ID"] == "SPAC1"].iloc[0]["Name"] == "gene1"


def test_macrocomplex_fills_a_missing_symbol_from_the_gene_name_table(tmp_path):
    """No symbol in the annotation -> the gene-name table supplies it."""
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_source("go_macrocomplex", pombase_dir)
    assert out[out["Systematic ID"] == "SPAC2"].iloc[0]["Name"] == "gene2"


def test_macrocomplex_falls_back_to_systematic_id_when_unnamed(tmp_path):
    """No symbol AND no gene name -> the systematic id, never NaN."""
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_source("go_macrocomplex", pombase_dir)
    assert out[out["Systematic ID"] == "SPBC3"].iloc[0]["Name"] == "SPBC3"
    assert out["Name"].notna().all()


def test_macrocomplex_n_annotated_members_is_per_term_total(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_source("go_macrocomplex", pombase_dir)
    alpha = out[out["group_name"] == "alpha complex"]
    assert (alpha["n_annotated_members"] == 2).all()


def test_load_gene_names_drops_unnamed_genes(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    assert load_gene_names(pombase_dir) == GENE_NAMES


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
    _write_gene_metadata(tmp_path, GENE_NAMES)
    return tmp_path


def test_gaf_namespace_cc_only_returns_cc_terms(tmp_path):
    pytest.importorskip("goatools")
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_source("go_cc", pombase_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["group_id"]) == {"GO:0000100", "GO:0000101"}
    assert set(out["source"]) == {"go_cc"}


def test_gaf_namespace_bp_only_returns_bp_terms(tmp_path):
    pytest.importorskip("goatools")
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_source("go_bp", pombase_dir)
    assert set(out["group_id"]) == {"GO:0000200"}
    assert set(out["source"]) == {"go_bp"}


def test_gaf_namespace_propagates_child_gene_to_parent(tmp_path):
    pytest.importorskip("goatools")
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_source("go_cc", pombase_dir)
    parent_members = set(out.loc[out["group_id"] == "GO:0000100", "Systematic ID"])
    assert "SPAC3" in parent_members  # child annotation propagates up via is_a


def test_gaf_namespace_resolves_gene_names(tmp_path):
    """The GAF carries no symbol, so every Name must come from the gene-name table.

    This is the regression: go_cc / go_bp used to write Name = pd.NA, which
    _finalize then filled with the systematic id, so the whole column read
    "SPAC1" instead of "gene1".
    """
    pytest.importorskip("goatools")
    pombase_dir = _write_go_fixture(tmp_path)
    # Both namespaces, because the fixture's genes are split across CC (SPAC*) and
    # BP (SPBC*), and each adapter only returns its own namespace's terms.
    out = pd.concat([load_source("go_cc", pombase_dir), load_source("go_bp", pombase_dir)])
    names = dict(zip(out["Systematic ID"], out["Name"]))
    assert names["SPAC1"] == "gene1"
    assert names["SPAC3"] == "gene5"   # present only via the metadata table
    assert names["SPBC2"] == "gene4"
    # The GAF fixture annotates no gene outside GENE_NAMES, so nothing falls back.
    assert not (out["Name"] == out["Systematic ID"]).any()


def test_gaf_namespace_rejects_unknown_namespace(tmp_path):
    pombase_dir = _write_go_fixture(tmp_path)
    with pytest.raises(ValueError, match="namespace must be one of"):
        load_gaf_namespace(pombase_dir, "MF", GENE_NAMES)


def test_source_loaders_registry_has_every_configured_source():
    from workflow.src.coherence.sources import SOURCE_LOADERS
    assert set(SOURCE_LOADERS) == {
        "go_macrocomplex", "go_cc", "go_bp", "kegg_brite", "kegg_pathway",
    }


def test_load_source_rejects_unknown_source(tmp_path):
    with pytest.raises(ValueError, match="unknown source"):
        load_source("not_a_source", tmp_path)


def _write_kegg_dir(tmp_path: Path) -> Path:
    """A tiny kegg_dir + gene-name table: one pathway, two BRITE trees.

    Gene_Symbol carries what KEGG actually puts there: PomBase's gene name when the
    gene has one (gene1), and the bare systematic id when it does not (SPAC2). The
    third row's symbol resolves to nothing and must be dropped.

    The BRITE fixture mirrors the real table's two id regimes: `spo00001` labels its
    nodes with a KEGG id (so `Level_D_ID` decides the group), `spo04131` labels none,
    so kegg_parser repeats the label back as the "id" and the path decides instead.
    Its rows also put two DIFFERENT `Level_B` branches under one `Level_D`, which is
    the case the id has to keep apart.
    """
    _write_gene_metadata(tmp_path, GENE_NAMES)
    kegg = tmp_path / "kegg"
    kegg.mkdir()
    pd.DataFrame({
        "Gene_Symbol": ["gene1", "SPAC2", "nosuchgene"],
        "Pathway_ID": ["spo00010"] * 3,
        "Pathway_Name": ["Glycolysis / Gluconeogenesis"] * 3,
    }).to_csv(kegg / "pathway_gene_mapping.tsv", sep="\t", index=False)
    pd.DataFrame({
        "BRITE_ID": ["spo00001", "spo00001", "spo04131", "spo04131"],
        # Two kingdom branches that both end at the same leaf label.
        "Level_A": ["Metabolism", "Metabolism", "Ribosomal proteins", "Ribosomal proteins"],
        "Level_B": ["Carbohydrate metabolism", "Carbohydrate metabolism",
                    "Eukaryotes", "Archaea"],
        "Level_C": ["Glycolysis / Gluconeogenesis", "Glycolysis / Gluconeogenesis",
                    "Large subunit", "Large subunit"],
        "Level_D": ["Glycolysis / Gluconeogenesis", "Glycolysis / Gluconeogenesis",
                    "Large subunit", "Large subunit"],
        # Non-empty everywhere, but the last two are kegg_parser's label fallback.
        "Level_D_ID": ["spo00010", "spo00010", "Large subunit", "Large subunit"],
        "Gene_Symbol": ["gene1", "gene3", "gene4", "gene4"],
    }).to_csv(kegg / "brite_flat.tsv", sep="\t", index=False)
    return kegg


def test_kegg_pathway_groups_by_pathway_and_drops_unmapped_genes(tmp_path):
    kegg_dir = _write_kegg_dir(tmp_path)
    out = load_source("kegg_pathway", tmp_path, kegg_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["source"]) == {"kegg_pathway"}
    assert set(out["group_id"]) == {"spo00010"}
    assert set(out["Systematic ID"]) == {"SPAC1", "SPAC2"}  # nosuchgene dropped
    assert set(out["Name"]) == {"gene1", "SPAC2"}
    assert (out["n_annotated_members"] == 2).all()


def test_kegg_brite_prefers_the_kegg_id_and_falls_back_to_the_path(tmp_path):
    kegg_dir = _write_kegg_dir(tmp_path)
    out = load_source("kegg_brite", tmp_path, kegg_dir)
    ids = dict(zip(out["Systematic ID"], out["group_id"]))
    # spo00001 labels its node with a KEGG id, so that id IS the group id.
    assert ids["SPAC1"] == "spo00010"
    assert ids["SPBC1"] == "spo00010"
    # spo04131 labels none, so the group id is the tree + the root-to-node path.
    # SPBC2 sits on both branches, hence a set rather than one value.
    assert set(out.loc[out["Systematic ID"] == "SPBC2", "group_id"]) == {
        "spo04131:Ribosomal proteins > Archaea > Large subunit > Large subunit",
        "spo04131:Ribosomal proteins > Eukaryotes > Large subunit > Large subunit",
    }
    # Only the deepest node is emitted — the shallower levels are not groups of their own.
    assert set(out["group_id"]) == {
        "spo00010",
        "spo04131:Ribosomal proteins > Archaea > Large subunit > Large subunit",
        "spo04131:Ribosomal proteins > Eukaryotes > Large subunit > Large subunit",
    }


def test_kegg_brite_keeps_same_named_branches_apart(tmp_path):
    """Two branches that both end at "Large subunit" are two groups, not one.

    The leaf label is not unique within a tree (four kingdom branches share "Large
    subunit" in the real ribosome tree), and `Level_D_ID` does not rescue it: where
    KEGG has no id for a node, kegg_parser repeats the label back as the id. Keying
    on either merged unrelated branches — 68 `Level_D_ID`s cover more than one node
    in the real table, 17 more collide across trees.
    """
    kegg_dir = _write_kegg_dir(tmp_path)
    out = load_source("kegg_brite", tmp_path, kegg_dir)
    brite = out[out["source"] == "kegg_brite"]
    shared_label = brite[brite["group_name"] == "Large subunit"]
    # gene4 is annotated under both the Eukaryotes and the Archaea branch, so the
    # two groups share a member and must NOT have been collapsed into one.
    assert shared_label["group_id"].nunique() == 2
    assert set(shared_label["group_id"]) == {
        "spo04131:Ribosomal proteins > Eukaryotes > Large subunit > Large subunit",
        "spo04131:Ribosomal proteins > Archaea > Large subunit > Large subunit",
    }


def test_kegg_sources_require_a_kegg_dir(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    for source in ("kegg_brite", "kegg_pathway"):
        with pytest.raises(ValueError, match="--kegg-dir"):
            load_source(source, pombase_dir)
