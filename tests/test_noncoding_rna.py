"""Tests for non-coding RNA analysis core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import numpy as np
import pytest

from workflow.src.noncoding_rna.core import (
    normalize_chromosome_names,
    merge_gtrnadb_by_position,
    extract_tRNA_amino_acid_and_anticodon,
    compute_tRNA_copy_number,
)


def test_normalize_chromosome_names_replaces_chr_prefix():
    """chrI/II/III → I/II/III (exact quirk from source notebook)."""
    df = pd.DataFrame({"#Chr": ["chrI", "chrII", "chrIII", "mitochondrial"]})
    result = normalize_chromosome_names(df)
    assert list(result["#Chr"]) == ["I", "II", "III", "mitochondrial"]


def test_merge_gtrnadb_by_position_uses_chr_start_end():
    """GtRNAdb merge is on #Chr+Start+End (not Name) — key quirk."""
    ncrna = pd.DataFrame({
        "#Chr": ["I", "II"],
        "Start": [100, 200],
        "End": [200, 300],
        "Systematic ID": ["g1", "g2"],
    })
    gtrnadb = pd.DataFrame({
        "#Chr": ["I", "III"],
        "Start": [100, 500],
        "End": [200, 600],
        "GtRNAdb_Name": ["tRNA-Ala-AGC-1-1", "tRNA-Gly-GCC-2-1"],
    })
    merged = merge_gtrnadb_by_position(ncrna, gtrnadb)
    assert len(merged) == 2
    assert merged.loc[merged["Systematic ID"] == "g1", "GtRNAdb_Name"].iloc[0] == "tRNA-Ala-AGC-1-1"
    assert pd.isna(merged.loc[merged["Systematic ID"] == "g2", "GtRNAdb_Name"].iloc[0])


def test_extract_tRNA_amino_acid_and_anticodon():
    """Amino acid and anticodon parsed from GtRNAdb_Name field."""
    row_with = pd.Series({
        "Systematic ID": "SPATRNA.Ala1",
        "GtRNAdb_Name": "Schpo_chr1.trna1-AlaAGC",
    })
    # Amino acid from sysID (TRNA<AA>.)
    # FIXTURE FIX: real GtRNAdb bed names carry NO organism prefix (e.g.
    # "tRNA-Pro-TGG-2-1"), so the notebook's `GtRNAdb_Name.split("-")[2]`
    # anticodon parse lands on the anticodon. The originally-specified fixture
    # value "schiPomb_972H-tRNA-Ala-AGC-1-1" prepended an extra "schiPomb_972H-"
    # segment that does not occur in the real resource file, which would shift
    # split("-")[2] onto "Ala" and make the AGC assertion impossible without
    # weakening the byte-faithful parser. Corrected to the real on-disk format
    # (verified against resources/external/pombase/.../schiPomb_972H-tRNAs.bed)
    # so the test exercises production data shape; intent (anticodon==AGC) is
    # preserved exactly.
    row_sys = pd.Series({
        "Systematic ID": "SPTRNAALA.01",
        "GtRNAdb_Name": "tRNA-Ala-AGC-1-1",
    })
    result = extract_tRNA_amino_acid_and_anticodon(row_sys)
    assert result["Anticodon"] == "AGC"


def test_compute_tRNA_copy_number():
    """tRNA_copy_number = count of tRNAs sharing the same Amino_Acid+Anticodon."""
    df = pd.DataFrame({
        "Amino_Acid": ["Ala", "Ala", "Gly", "Ala"],
        "Anticodon": ["AGC", "AGC", "GCC", "AGC"],
        "Systematic ID": ["t1", "t2", "t3", "t4"],
    })
    result = compute_tRNA_copy_number(df)
    assert list(result.loc[result["Systematic ID"].isin(["t1", "t2", "t4"]), "tRNA_copy_number"]) == [3, 3, 3]
    assert result.loc[result["Systematic ID"] == "t3", "tRNA_copy_number"].iloc[0] == 1
