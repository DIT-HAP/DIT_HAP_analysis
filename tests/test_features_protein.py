"""Tests for workflow/src/features/protein.py — peptide + pLDDT feature extraction."""

import gzip
import pytest
from pathlib import Path
from workflow.src.features.protein import (
    calculate_aliphatic_index_biopython,
    extract_protein_features_from_peptide_sequence,
    extract_pLDDT_pdb_gz,
    pLDDT_statistics_report,
)

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")

# A real, minimal-but-valid pdb.gz record (single CA atom) for pLDDT extraction tests.
_MINIMAL_PDB = (
    "HEADER    TEST\n"
    "ATOM      1  CA  ALA A   1      11.104  13.207   2.000  1.00 87.50           C\n"
    "TER\n"
    "END\n"
)


def test_calculate_aliphatic_index_known_sequence():
    """Aliphatic index of a single Ala residue (100% A) matches the Ikai formula by hand."""
    assert calculate_aliphatic_index_biopython("A") == 10000.0


def test_aliphatic_index_all_glycine_is_zero():
    """A sequence with none of A/V/L/I contributes zero to the aliphatic index."""
    assert calculate_aliphatic_index_biopython("GGGG") == 0.0


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_extract_protein_features_from_real_peptide_fasta():
    """Feature extraction from the real peptide.fa produces the known SPAC1002.01 row."""
    peptide_fasta = POMBASE_DIR / "genome_sequence_and_features" / "peptide.fa"
    df = extract_protein_features_from_peptide_sequence(peptide_fasta)
    row = df.set_index("Gene_id").loc["SPAC1002.01"]
    assert "aromaticity" in df.columns
    assert "aa_percent_Ala" in df.columns
    assert 0 <= row["aromaticity"] <= 1


def test_extract_pLDDT_pdb_gz_reads_bfactor_as_plddt(tmp_path):
    """extract_pLDDT_pdb_gz reads the bfactor column as the pLDDT score."""
    pdb_gz = tmp_path / "AF-P00000-F1-model_v6.pdb.gz"
    with gzip.open(pdb_gz, "wt") as f:
        f.write(_MINIMAL_PDB)
    pLDDT = extract_pLDDT_pdb_gz(pdb_gz)
    assert pLDDT == [87.50]


def test_pLDDT_statistics_report_computes_disorder_fraction(tmp_path):
    """A structure with one high-confidence residue (>50) has disorder_fraction 0.0."""
    pdb_gz = tmp_path / "AF-P00000-F1-model_v6.pdb.gz"
    with gzip.open(pdb_gz, "wt") as f:
        f.write(_MINIMAL_PDB)
    report = pLDDT_statistics_report(tmp_path, structure_format="pdb.gz")
    row = report.set_index("uniprot_id").loc["P00000"]
    assert row["disorder_fraction"] == 0.0
    assert row["mean_pLDDT"] == 87.5
