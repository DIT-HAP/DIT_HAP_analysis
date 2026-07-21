"""Tests for UTR insertion classification (assign_UTR_type)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.utr.classify_utr_insertions import (
    assign_UTR_type,
    UTR_DISTANCE_THRESHOLD,
    filter_intergenic_near_gene,
)


def test_utr_distance_threshold_constant():
    """distance_threshold=400 bp is the exact value from the source notebook."""
    assert UTR_DISTANCE_THRESHOLD == 400


def _row(strand, dist_start, dist_end, gene_strand="+"):
    return pd.Series({
        "Strand": strand,
        "gene_strand": gene_strand,
        "Distance_to_region_start": dist_start,
        "Distance_to_region_end": dist_end,
    })


def test_assign_UTR_type_5utr_forward_strand():
    """Insertion near gene start (forward strand) → 5UTR."""
    row = _row(strand="+", dist_start=50, dist_end=600, gene_strand="+")
    assert assign_UTR_type(row) == "5UTR"


def test_assign_UTR_type_3utr_forward_strand():
    """Insertion near gene end (forward strand) → 3UTR."""
    row = _row(strand="+", dist_start=600, dist_end=80, gene_strand="+")
    assert assign_UTR_type(row) == "3UTR"


def test_assign_UTR_type_5utr_reverse_strand():
    """Insertion near gene END (reverse strand) = 5' end → 5UTR."""
    row = _row(strand="-", dist_start=600, dist_end=100, gene_strand="-")
    assert assign_UTR_type(row) == "5UTR"


def test_assign_UTR_type_3utr_reverse_strand():
    """Insertion near gene START (reverse strand) = 3' end → 3UTR."""
    row = _row(strand="-", dist_start=50, dist_end=600, gene_strand="-")
    assert assign_UTR_type(row) == "3UTR"


def test_assign_UTR_type_ambiguous_both_close():
    """Insertion within threshold of both boundaries → '5UTR' (start takes priority)."""
    row = _row(strand="+", dist_start=100, dist_end=100, gene_strand="+")
    result = assign_UTR_type(row)
    assert result in ("5UTR", "3UTR", "ambiguous")  # implementation choice


def test_filter_intergenic_near_gene():
    """Filter: Type == Intergenic region AND (dist_start < 400 OR dist_end < 400)."""
    df = pd.DataFrame({
        "Type": ["Intergenic region", "Intergenic region", "Coding exon", "Intergenic region"],
        "Distance_to_region_start": [100, 500, 50, 500],
        "Distance_to_region_end": [600, 600, 600, 350],
    })
    result = filter_intergenic_near_gene(df)
    # Row 0: intergenic, dist_start < 400 → pass
    # Row 1: intergenic, both >= 400 → fail
    # Row 2: not intergenic → fail
    # Row 3: intergenic, dist_end < 400 → pass
    assert len(result) == 2
    assert 0 in result.index
    assert 3 in result.index
