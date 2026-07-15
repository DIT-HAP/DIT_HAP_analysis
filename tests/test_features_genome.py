"""Tests for workflow/src/features/genome.py — DNA-level feature extraction."""

import pytest
from pathlib import Path
from workflow.src.features.genome import (
    determine_primary_candidate,
    CHROMOSOME_END,
    CENTROMERE_POSITIONS,
)

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")


def test_chromosome_end_has_three_chromosomes():
    """CHROMOSOME_END covers chromosomes I, II, III on both arms."""
    assert set(CHROMOSOME_END["left"]) == {"I", "II", "III"}
    assert set(CHROMOSOME_END["right"]) == {"I", "II", "III"}


def test_centromere_positions_are_start_less_than_end():
    """Each centromere interval has start < end."""
    for chrom, (start, end) in CENTROMERE_POSITIONS.items():
        assert start < end, chrom


def test_determine_primary_candidate_hardcoded_exception_sPBC119_04():
    """SPBC119.04's .1 transcript is always primary regardless of peptide length match."""
    assert determine_primary_candidate("SPBC119.04", "SPBC119.04.1", 100, 999) is True
    assert determine_primary_candidate("SPBC119.04", "SPBC119.04.2", 100, 999) is False


def test_determine_primary_candidate_falls_back_to_length_match():
    """For genes without a hardcoded exception, primary candidacy is peptide-length equality."""
    assert determine_primary_candidate("SPAC1002.01", "SPAC1002.01.1", 162, 162) is True
    assert determine_primary_candidate("SPAC1002.01", "SPAC1002.01.2", 100, 162) is False


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_pombase_genome_config_loads_real_data():
    """PombaseGenomeConfig.from_pombase_dir loads the real genome FASTA and peptide_stats.tsv."""
    from workflow.src.features.genome import PombaseGenomeConfig

    cfg = PombaseGenomeConfig.from_pombase_dir(POMBASE_DIR)
    assert "I" in cfg.genome_length_dict
    assert cfg.primary_peptide_length["SPAC1002.01"] == 162


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_dna_level_features_for_one_known_gene(tmp_path):
    """DNA_level_features.from_gffutils_feature reproduces the known SPAC1002.01 record."""
    import gffutils
    from workflow.src.features.genome import PombaseGenomeConfig, DNA_level_features

    cfg = PombaseGenomeConfig.from_pombase_dir(POMBASE_DIR)
    db_path = tmp_path / "test.db"
    db = gffutils.create_db(cfg.gff3_file, str(db_path), force=True, merge_strategy="create_unique")

    mRNA = db["SPAC1002.01.1"]
    feat = DNA_level_features.from_gffutils_feature(mRNA, db, cfg)

    assert feat.Gene_id == "SPAC1002.01"
    assert feat.Chromosome == "I"
    assert feat.Strand == "+"
    assert feat.Peptide_length == 162
    assert feat.Primary_candidate is True
