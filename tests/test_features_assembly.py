"""Tests for the split feature-collection layout: assembly helpers + per-level driver configs."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.io_table import write_parquet, read_parquet

from workflow.src.features.assembly import (
    count_paralogs,
    get_ortholog_counts,
    load_ensembl_paralogs,
    load_paralogs,
    merge_all_features,
    parse_deletion_library_paralogs,
    read_coding_genes,
    read_deletion_library_paralogs,
)
from workflow.scripts.features.collect_dna_features import DnaConfig
from workflow.scripts.features.collect_rna_features import RnaConfig
from workflow.scripts.features.merge_features import MergeConfig


def test_get_ortholog_counts_counts_pipe_separated_entries(tmp_path):
    """A pipe-separated ortholog list counts entries; NONE maps to 0 via na_values."""
    f = tmp_path / "orthologs.txt"
    f.write_text("SPAC1002.01(name)\tOrthA|OrthB|OrthC\nSPAC1002.02(name)\tNONE\n")
    counts = get_ortholog_counts(f)
    assert counts.loc["SPAC1002.01"] == 3
    assert counts.loc["SPAC1002.02"] == 0


def _write_gene_meta(tmp_path: Path, ids: list[str]) -> Path:
    """Minimal PomBase gene_IDs_names_products.tsv listing `ids` as coding genes.

    `synonyms` must hold at least one real value: an all-empty column reads back as
    float64 and update_sysIDs()'s `.str.split(",")` then raises.
    """
    f = tmp_path / "gene_IDs_names_products.tsv"
    f.write_text(
        "gene_systematic_id\tgene_name\tsynonyms\tgene_type\n"
        + "".join(f"{g}\t{g}\t{g}-old\tprotein coding gene\n" for g in ids)
    )
    return f


def test_count_paralogs_distinguishes_none_from_one(tmp_path):
    """A blank-paralogue row means "none", not 1; unnamed genes are still counted."""
    meta = _write_gene_meta(tmp_path, ["g1", "g2", "g3", "g4"])
    export = pd.DataFrame({
        "Gene stable ID": ["g1", "g1", "g1", "g2", "g3", "g4"],
        "Gene name": ["a", "a", "a", "b", None, None],
        # g2: no paralogue at all (one row, blank paralogue column)
        # g3: 1 paralogue, but Ensembl ships no gene name for it
        "Schizosaccharomyces pombe paralogue gene stable ID": ["g9", "g10", None, None, "g8", "g7"],
    })
    tsv = tmp_path / "ensembl_paralogs.tsv"
    export.to_csv(tsv, sep="\t", index=False)

    pairs = load_ensembl_paralogs(tsv, meta)
    counts = count_paralogs(pairs, coding_genes=["g1", "g2", "g3"])

    assert counts.loc["g1", "paralog_count"] == 2
    assert counts.loc["g3", "paralog_count"] == 1, "unnamed genes must not be dropped"
    assert "g2" not in counts.index, "no paralogue -> absent, filled to 0 by the caller"
    assert "g4" not in counts.index, "genes outside coding_genes must be filtered out"


def test_parse_deletion_library_paralogs_explodes_none_and_normalises_case(tmp_path):
    """`|`-joined lists explode to long, NONE drops out, `.NNNC` normalises to `.NNNc`."""
    meta = _write_gene_meta(
        tmp_path,
        ["SPAC1002.13c", "SPAC1002.16c", "SPAC1002.12c", "SPBC2G2.17c", "SPAC1399.04c"],
    )
    xlsx = tmp_path / "deletion_library_categories.xlsx"
    pd.DataFrame({
        "Systematic ID": ["SPAC1002.13c", "SPAC1002.12c", "SPAC1002.16c"],
        "Paralogues": ["SPBC2G2.17C", "NONE", "SPAC1399.04C|SPAC1399.04c"],
    }).to_excel(xlsx, index=False)

    pairs = parse_deletion_library_paralogs(xlsx, meta)

    assert list(pairs.columns) == [
        "gene_systematic_id", "gene_name", "paralog_systematic_id", "paralog_name"
    ]
    assert set(pairs["gene_systematic_id"]) == {"SPAC1002.13c", "SPAC1002.16c"}
    assert set(pairs["paralog_systematic_id"]) == {"SPBC2G2.17c", "SPAC1399.04c"}
    assert len(pairs) == 2, "the two case-variants of one paralogue must collapse to one pair"


def test_parse_deletion_library_paralogs_falls_back_to_the_systematic_id(tmp_path):
    """An id PomBase has no name for (or does not know at all) keeps its own id."""
    # SPAC1002.13c is in the metadata with an EMPTY name -> falls back to its id.
    # SPBC2G2.17c is absent from the metadata entirely -> also falls back.
    meta = tmp_path / "gene_IDs_names_products.tsv"
    meta.write_text(
        "gene_systematic_id\tgene_name\tsynonyms\tgene_type\n"
        "SPAC1002.13c\t\tSPAC1002.13c-old\tprotein coding gene\n"
    )
    xlsx = tmp_path / "dl.xlsx"
    pd.DataFrame({
        "Systematic ID": ["SPAC1002.13c"],
        "Paralogues": ["SPBC2G2.17c"],
    }).to_excel(xlsx, index=False)

    pairs = parse_deletion_library_paralogs(xlsx, meta)

    assert pairs.loc[0, "gene_name"] == "SPAC1002.13c", "unnamed gene falls back to its id"
    assert pairs.loc[0, "paralog_name"] == "SPBC2G2.17c", "unknown id falls back to itself"


def test_deletion_library_paralogs_round_trip_through_parquet(tmp_path):
    """What build_deletion_library_paralogs writes is what the feature reads back."""
    meta = _write_gene_meta(tmp_path, ["g1", "g2"])
    xlsx = tmp_path / "dl.xlsx"
    pd.DataFrame({
        "Systematic ID": ["g1", "g2"],
        "Paralogues": ["g2", "NONE"],
    }).to_excel(xlsx, index=False)

    pairs = parse_deletion_library_paralogs(xlsx, meta)
    out = tmp_path / "deletion_library_paralogs.parquet"
    write_parquet(pairs, out)

    read_back = read_deletion_library_paralogs(out)
    assert read_back.columns.tolist() == pairs.columns.tolist()
    assert count_paralogs(read_back, coding_genes=["g1", "g2"]).loc["g1", "paralog_count"] == 1


def test_load_paralogs_rejects_an_unknown_source(tmp_path):
    """The config contract: only PARALOG_SOURCES are accepted."""
    meta = _write_gene_meta(tmp_path, ["g1"])
    with pytest.raises(ValueError, match="paralog source must be one of"):
        load_paralogs("pombase", tmp_path / "a", tmp_path / "b", meta)


def test_get_ortholog_counts_strips_parenthetical_gene_name(tmp_path):
    """The index's trailing (name) suffix is stripped before returning counts."""
    f = tmp_path / "orthologs.txt"
    f.write_text("SPAC1002.01(mrx11)\tOrthA\n")
    counts = get_ortholog_counts(f)
    assert "SPAC1002.01" in counts.index
    assert "SPAC1002.01(mrx11)" not in counts.index


def test_read_coding_genes_recovers_unique_gene_ids(tmp_path):
    """read_coding_genes returns the unique Gene_id set from a DNA-level parquet."""
    parquet_file = tmp_path / "dna_features.parquet"
    write_parquet(pd.DataFrame({"Gene_id": ["g1", "g1", "g2"], "Primary_candidate": [True, False, True]}), parquet_file)
    assert read_coding_genes(parquet_file) == ["g1", "g2"]


def test_dna_config_rejects_missing_pombase_dir(tmp_path):
    """DnaConfig.validate raises ValueError naming the missing PomBase dir."""
    cfg = DnaConfig(
        pombase_dir=tmp_path / "missing_pombase",
        genome_landmarks=tmp_path / "genome_landmarks.yaml",
        output_dna=tmp_path / "out" / "dna.parquet",
        output_codon_usage=tmp_path / "out" / "codon.tsv",
    )
    with pytest.raises(ValueError, match="does not exist"):
        cfg.validate()


def test_rna_config_gene_meta_file_property(tmp_path):
    """RnaConfig.gene_meta_file resolves under the PomBase Gene_metadata dir."""
    cfg = RnaConfig(
        pombase_dir=tmp_path / "pombase" / "2025-10-01",
        literature_dir=tmp_path / "lit",
        dna_features=tmp_path / "dna.parquet",
        output_rna=tmp_path / "out" / "rna.parquet",
    )
    assert cfg.gene_meta_file == tmp_path / "pombase" / "2025-10-01" / "Gene_metadata" / "gene_IDs_names_products.tsv"


def test_merge_all_features_fills_new_deletion_library_columns(tmp_path):
    """Sub_category/Growth_tier ride along with the other DeletionLibrary_* columns and get NA-filled."""
    dna_df = pd.DataFrame({"Gene_id": ["g1", "g2"], "Primary_candidate": [True, True]})
    empty_df = pd.DataFrame(index=pd.Index(["g1", "g2"], name="Gene_id"))
    network_df = pd.DataFrame(
        {"GO_term_richness": [1, 2], "PPI_degree": [0, 0], "GI_degree": [0, 0]},
        index=pd.Index(["g1", "g2"], name="Gene_id"),
    )
    phenotype_df = pd.DataFrame(
        {
            "DeletionLibrary_essentiality": ["V", None],
            "DeletionLibrary_category": ["WT-like", None],
            "RevisedDeletionLibrary_essentiality": ["V", None],
            "Sub_category": ["WT-like", None],
            "Growth_tier": [11, None],
        },
        index=pd.Index(["g1", "g2"], name="Gene_id"),
    )
    gene_meta = pd.DataFrame({"gene_systematic_id": ["g1", "g2"], "gene_name": ["a", "b"]})

    merged = merge_all_features(dna_df, empty_df, empty_df, empty_df, network_df, phenotype_df, gene_meta)

    assert merged.set_index("gene_systematic_id").loc["g2", "Sub_category"] == "Not_determined"
    assert merged.set_index("gene_systematic_id").loc["g2", "Growth_tier"] == 0
    assert merged.set_index("gene_systematic_id").loc["g1", "Growth_tier"] == 11
    assert merged["Growth_tier"].dtype == int


def test_merge_config_rejects_missing_level_parquet(tmp_path):
    """MergeConfig.validate raises when a per-level parquet input is absent."""
    real = tmp_path / "real"
    real.mkdir()
    (real / "dna.parquet").write_bytes(b"")
    cfg = MergeConfig(
        pombase_dir=real,
        dna_features=real / "dna.parquet",
        rna_features=real / "missing_rna.parquet",
        protein_features=real / "dna.parquet",
        evolutionary_features=real / "dna.parquet",
        network_features=real / "dna.parquet",
        phenotype_features=real / "dna.parquet",
        output_features=tmp_path / "out" / "features.tsv",
    )
    with pytest.raises(ValueError, match="does not exist"):
        cfg.validate()
