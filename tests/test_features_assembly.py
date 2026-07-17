"""Tests for the split feature-collection layout: assembly helpers + per-level driver configs."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.features.assembly import get_ortholog_counts, read_coding_genes
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


def test_get_ortholog_counts_strips_parenthetical_gene_name(tmp_path):
    """The index's trailing (name) suffix is stripped before returning counts."""
    f = tmp_path / "orthologs.txt"
    f.write_text("SPAC1002.01(mrx11)\tOrthA\n")
    counts = get_ortholog_counts(f)
    assert "SPAC1002.01" in counts.index
    assert "SPAC1002.01(mrx11)" not in counts.index


def test_read_coding_genes_recovers_unique_gene_ids(tmp_path):
    """read_coding_genes returns the unique Gene_id set from a DNA-level pickle."""
    pkl = tmp_path / "dna_features.pkl"
    pd.DataFrame({"Gene_id": ["g1", "g1", "g2"], "Primary_candidate": [True, False, True]}).to_pickle(pkl)
    assert read_coding_genes(pkl) == ["g1", "g2"]


def test_dna_config_rejects_missing_pombase_dir(tmp_path):
    """DnaConfig.validate raises ValueError naming the missing PomBase dir."""
    cfg = DnaConfig(
        pombase_dir=tmp_path / "missing_pombase",
        genome_landmarks=tmp_path / "genome_landmarks.yaml",
        output_dna=tmp_path / "out" / "dna.pkl",
        output_codon_usage=tmp_path / "out" / "codon.tsv",
    )
    with pytest.raises(ValueError, match="does not exist"):
        cfg.validate()


def test_rna_config_gene_meta_file_property(tmp_path):
    """RnaConfig.gene_meta_file resolves under the PomBase Gene_metadata dir."""
    cfg = RnaConfig(
        pombase_dir=tmp_path / "pombase" / "2025-10-01",
        literature_dir=tmp_path / "lit",
        dna_features=tmp_path / "dna.pkl",
        output_rna=tmp_path / "out" / "rna.pkl",
    )
    assert cfg.gene_meta_file == tmp_path / "pombase" / "2025-10-01" / "Gene_metadata" / "gene_IDs_names_products.tsv"


def test_merge_config_rejects_missing_level_pickle(tmp_path):
    """MergeConfig.validate raises when a per-level pickle input is absent."""
    real = tmp_path / "real"
    real.mkdir()
    (real / "dna.pkl").write_bytes(b"")
    cfg = MergeConfig(
        pombase_dir=real,
        dna_features=real / "dna.pkl",
        rna_features=real / "missing_rna.pkl",
        protein_features=real / "dna.pkl",
        evolutionary_features=real / "dna.pkl",
        network_features=real / "dna.pkl",
        phenotype_features=real / "dna.pkl",
        output_features=tmp_path / "out" / "features.tsv",
    )
    with pytest.raises(ValueError, match="does not exist"):
        cfg.validate()
