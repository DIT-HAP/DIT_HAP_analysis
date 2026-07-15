"""Tests for workflow/scripts/features/collect_pombe_features.py helper functions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.scripts.features.collect_pombe_features import get_ortholog_counts, InputOutputConfig
import pandas as pd
import pytest


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


def test_input_output_config_rejects_missing_input(tmp_path):
    """InputOutputConfig.__post_init__ raises ValueError naming the missing path."""
    real_dir = tmp_path / "real"
    real_dir.mkdir()
    with pytest.raises(ValueError, match="does not exist"):
        InputOutputConfig(
            pombase_dir=tmp_path / "missing_pombase",
            alphafold_dir=real_dir,
            literature_dir=real_dir,
            deletion_library_xlsx=real_dir / "x.xlsx",
            essentiality_verification_csv=real_dir / "x.csv",
            biogrid_tsv=real_dir / "x.tsv",
            ensembl_paralogs_tsv=real_dir / "x.tsv",
            output_features=tmp_path / "out" / "features.tsv",
            output_codon_usage=tmp_path / "out" / "codon.tsv",
        )
