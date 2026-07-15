"""Tests for workflow/src/io.py file-extension dispatch."""

import pytest
import pandas as pd
from pathlib import Path
from workflow.src.io import read_file


def test_read_tsv(tmp_path):
    """A .tsv file is read with tab separator."""
    f = tmp_path / "data.tsv"
    f.write_text("a\tb\n1\t2\n")
    df = read_file(f)
    assert list(df.columns) == ["a", "b"]
    assert df.iloc[0]["a"] == 1


def test_read_csv(tmp_path):
    """A .csv file is read with comma separator."""
    f = tmp_path / "data.csv"
    f.write_text("a,b\n1,2\n")
    df = read_file(f)
    assert list(df.columns) == ["a", "b"]


def test_read_bed_uses_tab_separator(tmp_path):
    """A .bed file is read with tab separator like tsv."""
    f = tmp_path / "regions.bed"
    f.write_text("chr1\t0\t100\n")
    df = read_file(f, header=None)
    assert df.shape == (1, 3)


def test_unsupported_extension_raises(tmp_path):
    """An unrecognized extension raises ValueError naming the file."""
    f = tmp_path / "data.json"
    f.write_text("{}")
    with pytest.raises(ValueError, match="Unsupported file type: data.json"):
        read_file(f)
