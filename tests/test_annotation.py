"""Tests for workflow/src/annotation/core.py gene annotation assembly."""

import pandas as pd

import pytest

from workflow.src.annotation.core import (
    annotate_table,
    assemble_annotation_reference,
    build_complex_block,
    build_go_slim_block,
    build_grna_block,
    build_hs_ortholog_block,
    build_pombe_block,
    build_sc_essentiality,
    build_sc_ortholog_block,
    parse_ortholog_field,
    read_sgd_phenotype_data,
    summarise_match,
)


def _sgd_phenotype_frame(rows: list[tuple[str, str, str]]) -> pd.DataFrame:
    """Build a minimal SGD phenotype_data frame from (feature, mutant_type, phenotype) rows."""
    return pd.DataFrame(
        [
            {"feature_name": feature, "mutant_type": mutant_type, "phenotype": phenotype}
            for feature, mutant_type, phenotype in rows
        ]
    )


# =============================================================================
# ORTHOLOG STRING PARSING
# =============================================================================
def test_single_ortholog_parses_to_one_id():
    """A plain single ortholog yields one id and no fragment marker."""
    assert parse_ortholog_field("YCL030C") == ["YCL030C"]


def test_pipe_separated_orthologs_split_into_independent_ids():
    """Pipe-separated orthologs are independent genes and split into separate ids."""
    assert parse_ortholog_field("YEL066W|YPR193C") == ["YEL066W", "YPR193C"]


def test_fusion_ortholog_keeps_both_fragments_as_one_group():
    """A '+'-joined fusion is one ortholog relation, so its fragments stay in one group."""
    assert parse_ortholog_field("YBR265W(N)+YDR302W(C)") == ["YBR265W+YDR302W"]


def test_fusion_fragment_markers_are_stripped_from_ids():
    """(N)/(C) fragment markers are stripped from ids so they can be looked up in SGD."""
    parsed = parse_ortholog_field("YGR036C(C)")
    assert parsed == ["YGR036C"]


def test_none_yields_empty_list():
    """PomBase's NONE sentinel means no ortholog and yields an empty list."""
    assert parse_ortholog_field("NONE") == []


def test_missing_value_yields_empty_list():
    """A NaN ortholog field yields an empty list rather than raising."""
    assert parse_ortholog_field(float("nan")) == []


def test_multi_pipe_orthologs_preserve_input_order():
    """Order is preserved so downstream per-ortholog columns stay positionally aligned."""
    assert parse_ortholog_field("YIL123W|YJL116C|YKR042W|YNL066W") == [
        "YIL123W",
        "YJL116C",
        "YKR042W",
        "YNL066W",
    ]


# =============================================================================
# S. CEREVISIAE ESSENTIALITY (null mutant viability)
# =============================================================================
def test_only_inviable_evidence_yields_inviable():
    """A gene whose null mutant is only ever reported inviable is called inviable."""
    pheno = _sgd_phenotype_frame([("YAL001C", "null", "inviable")])
    result = build_sc_essentiality(pheno)
    assert result.loc["YAL001C", "essentiality"] == "inviable"


def test_only_viable_evidence_yields_viable():
    """A gene whose null mutant is only ever reported viable is called viable."""
    pheno = _sgd_phenotype_frame([("YAL002W", "null", "viable")])
    result = build_sc_essentiality(pheno)
    assert result.loc["YAL002W", "essentiality"] == "viable"


def test_both_labels_yield_conflicting():
    """A gene with both inviable and viable null-mutant reports is called conflicting."""
    pheno = _sgd_phenotype_frame(
        [
            ("YAL003W", "null", "inviable"),
            ("YAL003W", "null", "viable"),
        ]
    )
    result = build_sc_essentiality(pheno)
    assert result.loc["YAL003W", "essentiality"] == "conflicting"


def test_evidence_column_counts_each_label():
    """The evidence column reports how many records support each label, inviable first."""
    pheno = _sgd_phenotype_frame(
        [
            ("YAL004W", "null", "inviable"),
            ("YAL004W", "null", "inviable"),
            ("YAL004W", "null", "inviable"),
            ("YAL004W", "null", "viable"),
        ]
    )
    result = build_sc_essentiality(pheno)
    assert result.loc["YAL004W", "essentiality_evidence"] == "inviable:3|viable:1"


def test_non_null_mutant_types_are_ignored():
    """Only null mutants speak to essentiality, so conditional/overexpression rows are dropped."""
    pheno = _sgd_phenotype_frame(
        [
            ("YAL005C", "conditional", "inviable"),
            ("YAL005C", "overexpression", "inviable"),
        ]
    )
    result = build_sc_essentiality(pheno)
    assert "YAL005C" not in result.index


def test_graded_viability_phenotypes_are_ignored():
    """'viability: decreased' is a graded phenotype, not a null-mutant viable/inviable call."""
    pheno = _sgd_phenotype_frame(
        [
            ("YAL006W", "null", "viability: decreased"),
            ("YAL006W", "null", "viability: increased"),
        ]
    )
    result = build_sc_essentiality(pheno)
    assert "YAL006W" not in result.index


def test_gene_with_no_viability_records_is_absent():
    """A gene with only unrelated phenotypes gets no essentiality row (left empty on join)."""
    pheno = _sgd_phenotype_frame([("YAL007C", "null", "heat sensitivity: increased")])
    result = build_sc_essentiality(pheno)
    assert "YAL007C" not in result.index


# =============================================================================
# READING SGD's RAGGED PHENOTYPE FILE
# =============================================================================
def test_reads_headerless_phenotype_file_into_named_columns(tmp_path):
    """SGD's phenotype_data.tab has no header row, so column names are supplied positionally."""
    f = tmp_path / "phenotype_data.tab"
    f.write_text(
        "YAL001C\tORF\tTFC3\tS000000001\tSGD_REF:S1|PMID:1\tclassical genetics\t"
        "null\ttfc3-1\tS288C\tinviable\t\t\t\t\n"
    )
    result = read_sgd_phenotype_data(f)
    assert result.loc[0, "feature_name"] == "YAL001C"
    assert result.loc[0, "mutant_type"] == "null"
    assert result.loc[0, "phenotype"] == "inviable"


def test_reads_rows_carrying_an_extra_trailing_field(tmp_path):
    """39 real rows have a 15th field; the reader must keep them, not abort the whole file."""
    f = tmp_path / "phenotype_data.tab"
    f.write_text(
        "YAL001C\tORF\tTFC3\tS000000001\tSGD_REF:S1|PMID:1\tclassical genetics\t"
        "null\ttfc3-1\tS288C\tinviable\t\t\t\t\t\n"
    )
    result = read_sgd_phenotype_data(f)
    assert len(result) == 1
    assert result.loc[0, "phenotype"] == "inviable"


def test_reads_rows_missing_trailing_fields(tmp_path):
    """19 real rows are truncated to 13 fields; missing trailing columns become NaN."""
    f = tmp_path / "phenotype_data.tab"
    f.write_text(
        "YAL001C\tORF\tTFC3\tS000000001\tSGD_REF:S1|PMID:1\tclassical genetics\t"
        "null\ttfc3-1\tS288C\tinviable\t\t\t\n"
    )
    result = read_sgd_phenotype_data(f)
    assert len(result) == 1
    assert result.loc[0, "phenotype"] == "inviable"


# =============================================================================
# JOINING ANNOTATION ONTO A USER TABLE
# =============================================================================
@pytest.fixture
def annotation_reference() -> pd.DataFrame:
    """A two-gene annotation reference indexed by pombe systematic ID."""
    return pd.DataFrame(
        {"gene_name": ["cdc2", "mrx11"], "Sc_ortholog_name": ["CDC28", "MRX11"]},
        index=pd.Index(["SPBC11B10.09", "SPAC1002.01"], name="gene_systematic_id"),
    )


def test_annotation_columns_are_appended(annotation_reference):
    """Annotation columns land alongside the caller's own columns."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"], "my_score": [1.5]})
    result = annotate_table(user_table, annotation_reference, gene_column="gene")
    assert result.loc[0, "my_score"] == 1.5
    assert result.loc[0, "gene_name"] == "cdc2"
    assert result.loc[0, "Sc_ortholog_name"] == "CDC28"


def test_duplicate_gene_ids_do_not_change_row_count(annotation_reference):
    """A gene repeated across rows must annotate each row without fanning out the table."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09", "SPBC11B10.09"], "sample": ["a", "b"]})
    result = annotate_table(user_table, annotation_reference, gene_column="gene")
    assert len(result) == 2
    assert list(result["gene_name"]) == ["cdc2", "cdc2"]


def test_unmatched_rows_are_kept_with_empty_annotation(annotation_reference):
    """Unmatched genes are kept by default — dropping them would be silent data loss."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09", "not_a_gene"]})
    result = annotate_table(user_table, annotation_reference, gene_column="gene")
    assert len(result) == 2
    assert pd.isna(result.loc[1, "gene_name"])


def test_drop_unmatched_removes_unmatched_rows(annotation_reference):
    """With drop_unmatched, rows without annotation are removed."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09", "not_a_gene"]})
    result = annotate_table(
        user_table, annotation_reference, gene_column="gene", drop_unmatched=True
    )
    assert list(result["gene"]) == ["SPBC11B10.09"]


def test_row_order_is_preserved(annotation_reference):
    """The caller's row order survives the join so results line up with their input."""
    user_table = pd.DataFrame({"gene": ["SPAC1002.01", "SPBC11B10.09"]})
    result = annotate_table(user_table, annotation_reference, gene_column="gene")
    assert list(result["gene_name"]) == ["mrx11", "cdc2"]


def test_selected_columns_only(annotation_reference):
    """Passing columns restricts the appended annotation to those fields."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"]})
    result = annotate_table(
        user_table, annotation_reference, gene_column="gene", columns=["gene_name"]
    )
    assert "gene_name" in result.columns
    assert "Sc_ortholog_name" not in result.columns


def test_missing_gene_column_raises(annotation_reference):
    """A wrong --gene-column must fail loudly rather than produce an unannotated table."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"]})
    with pytest.raises(KeyError, match="wrong_column"):
        annotate_table(user_table, annotation_reference, gene_column="wrong_column")


def test_requesting_absent_annotation_column_raises(annotation_reference):
    """Asking for a column the reference lacks fails rather than silently omitting it."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"]})
    with pytest.raises(KeyError, match="Sp_growth_category"):
        annotate_table(
            user_table,
            annotation_reference,
            gene_column="gene",
            columns=["gene_name", "Sp_growth_category"],
        )


def test_colliding_column_is_suffixed_not_overwritten(annotation_reference):
    """A caller column of the same name is preserved; the annotation gets a suffix."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"], "gene_name": ["my_own_label"]})
    result = annotate_table(user_table, annotation_reference, gene_column="gene")
    assert result.loc[0, "gene_name"] == "my_own_label"
    assert result.loc[0, "gene_name_annotation"] == "cdc2"


# =============================================================================
# MATCH SUMMARY (so unmatched ids are never silent)
# =============================================================================
def test_summary_raises_a_helpful_error_for_a_wrong_gene_column(annotation_reference):
    """A wrong --gene-column must name the available columns, not surface a bare pandas KeyError."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09"], "score": [1]})
    with pytest.raises(KeyError, match="Available columns"):
        summarise_match(user_table, annotation_reference, gene_column="WrongName")


def test_summary_reports_matched_count_and_unmatched_ids(annotation_reference):
    """The caller needs both numbers to notice a wrong column or a non-coding gene slipping in."""
    user_table = pd.DataFrame({"gene": ["SPBC11B10.09", "not_a_gene", "also_missing"]})
    matched, unmatched = summarise_match(
        user_table, annotation_reference, gene_column="gene"
    )
    assert matched == 1
    assert unmatched == ["also_missing", "not_a_gene"]


def test_summary_deduplicates_unmatched_ids(annotation_reference):
    """A repeated unmatched id is reported once, so the log stays readable."""
    user_table = pd.DataFrame({"gene": ["nope", "nope", "nope"]})
    matched, unmatched = summarise_match(
        user_table, annotation_reference, gene_column="gene"
    )
    assert matched == 0
    assert unmatched == ["nope"]


# =============================================================================
# ORTHOLOG ANNOTATION BLOCK
# =============================================================================
@pytest.fixture
def sc_gene_info() -> pd.DataFrame:
    """SGD per-ORF standard name, qualifier and description, indexed by systematic name."""
    return pd.DataFrame(
        {
            "standard_name": ["HIS4", None, "VPS54"],
            "qualifier": ["Verified", "Dubious", "Verified"],
            "description": ["histidinol dehydrogenase", "unknown", "GARP complex subunit"],
        },
        index=pd.Index(["YCL030C", "YBR265W", "YDR027C"], name="systematic_name"),
    )


@pytest.fixture
def sc_essentiality() -> pd.DataFrame:
    """Per-ORF null-mutant essentiality calls, indexed by systematic name."""
    return pd.DataFrame(
        {
            "essentiality": ["inviable", "viable"],
            "essentiality_evidence": ["inviable:2", "viable:5"],
        },
        index=pd.Index(["YCL030C", "YDR027C"], name="feature_name"),
    )


def _ortholog_frame(pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """Build a PomBase curated_orthologs frame from (pombe_id, ortholog_field) pairs."""
    return pd.DataFrame(pairs, columns=["gene_systematic_id", "orthologs"])


def test_ortholog_block_resolves_standard_name(sc_gene_info, sc_essentiality):
    """A single ortholog gets its SGD standard name."""
    orthologs = _ortholog_frame([("SPBC1711.13", "YCL030C")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPBC1711.13", "Sc_ortholog_id"] == "YCL030C"
    assert result.loc["SPBC1711.13", "Sc_ortholog_name"] == "HIS4"


def test_ortholog_without_standard_name_falls_back_to_systematic_id(sc_gene_info, sc_essentiality):
    """~1300 SGD ORFs have no common name, so the systematic id is used instead of a blank."""
    orthologs = _ortholog_frame([("SPAC0001.01", "YBR265W")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPAC0001.01", "Sc_ortholog_name"] == "YBR265W"


def test_ortholog_block_carries_orf_qualifier(sc_gene_info, sc_essentiality):
    """A Dubious ORF makes the ortholog relation weak evidence, so the qualifier is surfaced."""
    orthologs = _ortholog_frame([("SPAC0001.01", "YBR265W")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPAC0001.01", "Sc_ortholog_qualifier"] == "Dubious"


def test_multiple_orthologs_stay_positionally_aligned(sc_gene_info, sc_essentiality):
    """Per-ortholog columns are pipe-joined in the same order so position i refers to one gene."""
    orthologs = _ortholog_frame([("SPAC0002.01", "YCL030C|YDR027C")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    row = result.loc["SPAC0002.01"]
    assert row["Sc_ortholog_id"] == "YCL030C|YDR027C"
    assert row["Sc_ortholog_name"] == "HIS4|VPS54"
    assert row["Sc_ortholog_qualifier"] == "Verified|Verified"
    assert row["Sc_essentiality"] == "inviable|viable"


def test_ortholog_count_counts_independent_orthologs(sc_gene_info, sc_essentiality):
    """Count reflects independent ortholog genes, not fusion fragments."""
    orthologs = _ortholog_frame(
        [
            ("SPAC0002.01", "YCL030C|YDR027C"),
            ("SPCC1450.15", "YBR265W(N)+YCL030C(C)"),
            ("SPAC0003.01", "NONE"),
        ]
    )
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPAC0002.01", "Sc_ortholog_count"] == 2
    assert result.loc["SPCC1450.15", "Sc_ortholog_count"] == 1
    assert result.loc["SPAC0003.01", "Sc_ortholog_count"] == 0


def test_fusion_ortholog_reports_both_fragment_names(sc_gene_info, sc_essentiality):
    """A fusion's fragments are one relation, so their names join with '+' inside one group."""
    orthologs = _ortholog_frame([("SPCC1450.15", "YBR265W(N)+YCL030C(C)")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPCC1450.15", "Sc_ortholog_name"] == "YBR265W+HIS4"


def test_raw_ortholog_field_is_preserved(sc_gene_info, sc_essentiality):
    """The original string keeps the (N)/(C) terminus information that parsing strips."""
    orthologs = _ortholog_frame([("SPCC1450.15", "YBR265W(N)+YCL030C(C)")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    assert result.loc["SPCC1450.15", "Sc_ortholog_raw"] == "YBR265W(N)+YCL030C(C)"


def test_ortholog_lacking_essentiality_record_is_blank_not_dropped(sc_gene_info, sc_essentiality):
    """An ortholog with no SGD viability record leaves that slot empty but keeps alignment."""
    orthologs = _ortholog_frame([("SPAC0004.01", "YCL030C|YBR265W")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    row = result.loc["SPAC0004.01"]
    assert row["Sc_ortholog_id"] == "YCL030C|YBR265W"
    assert row["Sc_essentiality"] == "inviable|"


def test_gene_with_no_ortholog_has_empty_id_and_zero_count(sc_gene_info, sc_essentiality):
    """A NONE row yields empty ortholog fields rather than being dropped from the reference."""
    orthologs = _ortholog_frame([("SPAC0003.01", "NONE")])
    result = build_sc_ortholog_block(orthologs, sc_gene_info, sc_essentiality)
    row = result.loc["SPAC0003.01"]
    assert row["Sc_ortholog_id"] == ""
    assert row["Sc_ortholog_count"] == 0


def test_human_ortholog_block_gives_symbols_and_count():
    """Human orthologs come straight from PomBase symbols, with a count column."""
    orthologs = _ortholog_frame(
        [("SPAC1002.05c", "KDM5A|KDM5B|KDM5C"), ("SPAC0003.01", "NONE")]
    )
    result = build_hs_ortholog_block(orthologs)
    assert result.loc["SPAC1002.05c", "Hs_ortholog_symbol"] == "KDM5A|KDM5B|KDM5C"
    assert result.loc["SPAC1002.05c", "Hs_ortholog_count"] == 3
    assert result.loc["SPAC0003.01", "Hs_ortholog_count"] == 0


# =============================================================================
# POMBE-SIDE BLOCK (identity + two independent essentiality sources)
# =============================================================================
def test_pombe_block_uses_gene_name_and_product():
    """Gene name and product description come from PomBase gene metadata."""
    gene_meta = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1002.01"],
            "gene_name": ["mrx11"],
            "gene_product": ["MIOREX component Mrx11"],
            "synonyms": ["SPAC1610.05"],
        }
    )
    result = build_pombe_block(gene_meta, pd.Series(dtype=str), pd.DataFrame())
    assert result.loc["SPAC1002.01", "gene_name"] == "mrx11"
    assert result.loc["SPAC1002.01", "gene_product"] == "MIOREX component Mrx11"
    assert result.loc["SPAC1002.01", "synonyms"] == "SPAC1610.05"


def test_unnamed_pombe_gene_falls_back_to_systematic_id():
    """An unnamed gene shows its systematic id so the column is never blank."""
    gene_meta = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1002.02"],
            "gene_name": [None],
            "gene_product": ["conserved protein"],
            "synonyms": [None],
        }
    )
    result = build_pombe_block(gene_meta, pd.Series(dtype=str), pd.DataFrame())
    assert result.loc["SPAC1002.02", "gene_name"] == "SPAC1002.02"


def test_pombe_block_carries_fypo_viability():
    """FYPO viability is genome-wide but mostly 'unknown', so it is kept as its own column."""
    gene_meta = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1002.04c"],
            "gene_name": ["abc1"],
            "gene_product": ["p"],
            "synonyms": [None],
        }
    )
    viability = pd.Series({"SPAC1002.04c": "inviable"})
    result = build_pombe_block(gene_meta, viability, pd.DataFrame())
    assert result.loc["SPAC1002.04c", "Sp_FYPO_viability"] == "inviable"


def test_pombe_block_carries_deletion_library_columns():
    """The Hayles deletion library is a separate measured source from FYPO, so both are kept."""
    gene_meta = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1002.04c"],
            "gene_name": ["abc1"],
            "gene_product": ["p"],
            "synonyms": [None],
        }
    )
    deletion = pd.DataFrame(
        {
            "Systematic ID": ["SPAC1002.04c"],
            "Gene dispensability. This study": ["E"],
            "Phenotypic classification used for analysis": ["misshapen essential"],
            "Category": ["microcolonies"],
        }
    )
    result = build_pombe_block(gene_meta, pd.Series(dtype=str), deletion)
    row = result.loc["SPAC1002.04c"]
    assert row["Sp_deletion_essentiality"] == "E"
    assert row["Sp_deletion_phenotype"] == "misshapen essential"
    assert row["Sp_growth_category"] == "microcolonies"


def test_gene_absent_from_deletion_library_keeps_other_columns():
    """The deletion library covers only 4843 genes; the rest keep identity/FYPO columns."""
    gene_meta = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1002.05c"],
            "gene_name": ["xyz1"],
            "gene_product": ["p"],
            "synonyms": [None],
        }
    )
    result = build_pombe_block(
        gene_meta, pd.Series({"SPAC1002.05c": "viable"}), pd.DataFrame()
    )
    assert result.loc["SPAC1002.05c", "gene_name"] == "xyz1"
    assert result.loc["SPAC1002.05c", "Sp_FYPO_viability"] == "viable"
    assert pd.isna(result.loc["SPAC1002.05c", "Sp_deletion_essentiality"])


# =============================================================================
# FUNCTIONAL ANNOTATION BLOCK (complex membership)
# =============================================================================
def test_complex_membership_is_pipe_joined():
    """A gene in several complexes gets all of them, so partial membership is visible."""
    complexes = pd.DataFrame(
        {
            "systematic_id": ["SPBC1815.01", "SPBC1815.01", "SPAC1002.01"],
            "GO_term_name": ["enolase complex", "glycolytic complex", "MIOREX complex"],
        }
    )
    result = build_complex_block(complexes)
    assert result.loc["SPBC1815.01", "complex"] == "enolase complex|glycolytic complex"
    assert result.loc["SPAC1002.01", "complex"] == "MIOREX complex"


def test_duplicate_complex_annotations_are_deduplicated():
    """The same complex annotated from two evidence sources is reported once."""
    complexes = pd.DataFrame(
        {
            "systematic_id": ["SPBC1815.01", "SPBC1815.01"],
            "GO_term_name": ["enolase complex", "enolase complex"],
        }
    )
    result = build_complex_block(complexes)
    assert result.loc["SPBC1815.01", "complex"] == "enolase complex"


# =============================================================================
# GO-SLIM BLOCK
# =============================================================================
@pytest.fixture
def slim_term_names() -> dict[str, str]:
    """GO-slim term id -> readable term name."""
    return {
        "GO:0006281": "DNA repair",
        "GO:0006310": "DNA recombination",
        "GO:0005634": "nucleus",
        "GO:0003677": "DNA binding",
    }


def test_go_slim_block_splits_namespaces_into_columns(slim_term_names):
    """BP/CC/MF land in separate columns so each reads as one coherent functional axis."""
    ns2slim = {
        "BP": {"SPAC1002.01": {"GO:0006281"}},
        "CC": {"SPAC1002.01": {"GO:0005634"}},
        "MF": {"SPAC1002.01": {"GO:0003677"}},
    }
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert result.loc["SPAC1002.01", "GO_slim_BP"] == "DNA repair"
    assert result.loc["SPAC1002.01", "GO_slim_CC"] == "nucleus"
    assert result.loc["SPAC1002.01", "GO_slim_MF"] == "DNA binding"


def test_go_slim_terms_are_names_not_ids(slim_term_names):
    """Readable names are the whole point of using slim terms, so ids are resolved."""
    ns2slim = {"BP": {"SPAC1002.01": {"GO:0006310"}}}
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert result.loc["SPAC1002.01", "GO_slim_BP"] == "DNA recombination"


def test_multiple_slim_terms_are_sorted_for_stable_output(slim_term_names):
    """Terms arrive as an unordered set, so they are sorted to keep the table reproducible."""
    ns2slim = {"BP": {"SPAC1002.01": {"GO:0006310", "GO:0006281"}}}
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert result.loc["SPAC1002.01", "GO_slim_BP"] == "DNA recombination|DNA repair"


def test_gene_with_no_slim_terms_in_a_namespace_is_empty(slim_term_names):
    """A gene annotated only in BP gets empty CC/MF rather than a missing column."""
    ns2slim = {"BP": {"SPAC1002.01": {"GO:0006281"}}, "CC": {"SPAC1002.01": set()}}
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert result.loc["SPAC1002.01", "GO_slim_CC"] == ""


def test_unknown_slim_term_id_is_skipped(slim_term_names):
    """A term absent from the slim name table is dropped rather than emitting a bare GO id."""
    ns2slim = {"BP": {"SPAC1002.01": {"GO:0006281", "GO:9999999"}}}
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert result.loc["SPAC1002.01", "GO_slim_BP"] == "DNA repair"


def test_all_three_slim_columns_exist_even_when_namespace_absent(slim_term_names):
    """Downstream code selects these columns by name, so all three always exist."""
    ns2slim = {"BP": {"SPAC1002.01": {"GO:0006281"}}}
    result = build_go_slim_block(ns2slim, slim_term_names)
    assert list(result.columns) == ["GO_slim_BP", "GO_slim_CC", "GO_slim_MF"]


# =============================================================================
# ASSEMBLING THE FULL REFERENCE FROM BLOCKS
# =============================================================================
def _pombe_block(gene_ids: list[str]) -> pd.DataFrame:
    """A minimal pombe-side block for assembly tests."""
    return pd.DataFrame(
        {"gene_name": [f"name_{g}" for g in gene_ids]},
        index=pd.Index(gene_ids, name="gene_systematic_id"),
    )


def test_assembly_is_keyed_on_the_pombe_gene_set():
    """The pombe block defines the row set; annotation blocks only contribute columns."""
    pombe = _pombe_block(["SPAC0001.01", "SPAC0002.01"])
    sc = pd.DataFrame(
        {"Sc_ortholog_id": ["YAL001C"]},
        index=pd.Index(["SPAC0001.01"], name="gene_systematic_id"),
    )
    result = assemble_annotation_reference(pombe, [sc])
    assert list(result.index) == ["SPAC0001.01", "SPAC0002.01"]


def test_blocks_covering_extra_genes_do_not_add_rows():
    """An ortholog file listing a gene absent from the pombe gene set must not introduce it."""
    pombe = _pombe_block(["SPAC0001.01"])
    sc = pd.DataFrame(
        {"Sc_ortholog_id": ["YAL001C", "YAL002W"]},
        index=pd.Index(["SPAC0001.01", "SPAC_RETIRED.01"], name="gene_systematic_id"),
    )
    result = assemble_annotation_reference(pombe, [sc])
    assert list(result.index) == ["SPAC0001.01"]


def test_genes_missing_from_a_block_get_empty_values():
    """A gene with no ortholog record keeps its row with blanks in that block's columns."""
    pombe = _pombe_block(["SPAC0001.01", "SPAC0002.01"])
    sc = pd.DataFrame(
        {"Sc_ortholog_id": ["YAL001C"]},
        index=pd.Index(["SPAC0001.01"], name="gene_systematic_id"),
    )
    result = assemble_annotation_reference(pombe, [sc])
    assert pd.isna(result.loc["SPAC0002.01", "Sc_ortholog_id"])


def test_all_block_columns_are_present_after_assembly():
    """Every block's columns survive assembly, in block order after the pombe columns."""
    pombe = _pombe_block(["SPAC0001.01"])
    block_a = pd.DataFrame(
        {"Sc_ortholog_id": ["YAL001C"]},
        index=pd.Index(["SPAC0001.01"], name="gene_systematic_id"),
    )
    block_b = pd.DataFrame(
        {"complex": ["enolase complex"]},
        index=pd.Index(["SPAC0001.01"], name="gene_systematic_id"),
    )
    result = assemble_annotation_reference(pombe, [block_a, block_b])
    assert list(result.columns) == ["gene_name", "Sc_ortholog_id", "complex"]


# =============================================================================
# gRNA-LEVEL DEPLETION PARAMETERS
# =============================================================================
def _grna_frame(rows: list[tuple[str, float, float]]) -> pd.DataFrame:
    """Build a curated gRNA fitted-parameters frame from (systematic_id, um, lam) rows."""
    return pd.DataFrame(
        [
            {"Systematic ID": gene, "um": um, "lam": lam, "A": 7.0, "R2": 0.99}
            for gene, um, lam in rows
        ]
    )


def test_grna_block_renames_legacy_um_and_lam():
    """The curated table still uses upstream's legacy um/lam names for DR/DL."""
    grna = _grna_frame([("SPAC3A12.11c", 1.048, 1.103)])
    result = build_grna_block(grna)
    assert result.loc["SPAC3A12.11c", "gRNA_DR"] == 1.048
    assert result.loc["SPAC3A12.11c", "gRNA_DL"] == 1.103


def test_grna_block_emits_only_the_two_depletion_columns():
    """Only DR/DL are wanted, so the other fitted parameters are dropped."""
    grna = _grna_frame([("SPAC3A12.11c", 1.048, 1.103)])
    result = build_grna_block(grna)
    assert list(result.columns) == ["gRNA_DR", "gRNA_DL"]


def test_grna_columns_are_prefixed_to_avoid_colliding_with_gene_level_dr_dl():
    """Gene-level DR/DL are different measurements (DL correlates only ~0.55), so names must differ."""
    grna = _grna_frame([("SPAC3A12.11c", 1.048, 1.103)])
    result = build_grna_block(grna)
    assert "DR" not in result.columns
    assert "DL" not in result.columns


def test_grna_block_accepts_a_table_already_using_dr_dl_names():
    """If upstream ever ships DR/DL directly, the same loader must keep working."""
    grna = pd.DataFrame({"Systematic ID": ["SPAC3A12.11c"], "DR": [1.048], "DL": [1.103]})
    result = build_grna_block(grna)
    assert result.loc["SPAC3A12.11c", "gRNA_DR"] == 1.048
    assert result.loc["SPAC3A12.11c", "gRNA_DL"] == 1.103


def test_duplicate_gene_in_grna_table_raises():
    """The curated table is one row per gene; a duplicate would silently fan out the reference."""
    grna = _grna_frame([("SPAC3A12.11c", 1.0, 1.0), ("SPAC3A12.11c", 2.0, 2.0)])
    with pytest.raises(ValueError, match="duplicate"):
        build_grna_block(grna)


def test_missing_depletion_columns_raise():
    """A table with neither um/lam nor DR/DL is the wrong file and must fail loudly."""
    grna = pd.DataFrame({"Systematic ID": ["SPAC3A12.11c"], "something_else": [1.0]})
    with pytest.raises(KeyError, match="um/lam"):
        build_grna_block(grna)


def test_count_columns_stay_integer_after_joining_partial_blocks():
    """Genes missing from the ortholog block must not turn counts into floats like '1.0'."""
    pombe = _pombe_block(["SPAC0001.01", "SPAC0002.01"])
    sc = pd.DataFrame(
        {"Sc_ortholog_count": [1]},
        index=pd.Index(["SPAC0001.01"], name="gene_systematic_id"),
    )
    result = assemble_annotation_reference(pombe, [sc])
    assert result.loc["SPAC0001.01", "Sc_ortholog_count"] == 1
    assert str(result["Sc_ortholog_count"].dtype) == "Int64"
    assert pd.isna(result.loc["SPAC0002.01", "Sc_ortholog_count"])


def test_duplicate_index_in_a_block_raises():
    """A duplicated gene in a source block would silently fan out rows, so it fails loudly."""
    pombe = _pombe_block(["SPAC0001.01"])
    bad_block = pd.DataFrame(
        {"Sc_ortholog_id": ["YAL001C", "YAL002W"]},
        index=pd.Index(["SPAC0001.01", "SPAC0001.01"], name="gene_systematic_id"),
    )
    with pytest.raises(ValueError, match="duplicate"):
        assemble_annotation_reference(pombe, [bad_block])
