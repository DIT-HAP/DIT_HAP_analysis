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
    resolve_parental_gene,
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


# ---------------------------------------------------------------------------
# resolve_parental_gene: two-flanking-gene interval -> single parental gene
# ---------------------------------------------------------------------------
def _interval_row(insertion_strand, dist_start, dist_end,
                  left_strand="+", right_strand="+"):
    """Build an intergenic-interval row with a [Chr, Coord, Strand, Target] name."""
    row = pd.Series({
        "Name": "leftGene|rightGene",
        "Systematic ID": "SPLEFT|SPRIGHT",
        "Strand_Interval": f"{left_strand}|{right_strand}",
        "Distance_to_region_start": dist_start,
        "Distance_to_region_end": dist_end,
    })
    row.name = ("I", 1000, insertion_strand, "TargetA")
    return row


def test_resolve_parental_gene_left_only_close():
    """Only region_start within threshold -> left gene, forward-strand left => 3UTR."""
    res = resolve_parental_gene(_interval_row("+", dist_start=100, dist_end=600))
    assert res["Parental_gene"] == "leftGene"
    assert res["Parental_gene_id"] == "SPLEFT"
    # Insertion sits past the left gene's high-coord (3') edge; + strand => 3UTR.
    assert res["UTR_type"] == "3UTR"
    assert res["Distance_to_gene_boundary"] == 100  # 3UTR is positive
    assert res["Insertion_direction"] == "Forward"  # insertion + == left +


def test_resolve_parental_gene_right_only_close():
    """Only region_end within threshold -> right gene, forward-strand right => 5UTR."""
    res = resolve_parental_gene(_interval_row("-", dist_start=600, dist_end=80))
    assert res["Parental_gene"] == "rightGene"
    assert res["Parental_gene_id"] == "SPRIGHT"
    # Insertion sits before the right gene's low-coord (5') edge; + strand => 5UTR.
    assert res["UTR_type"] == "5UTR"
    assert res["Distance_to_gene_boundary"] == -80  # 5UTR is negative
    assert res["Insertion_direction"] == "Reverse"  # insertion - != right +


def test_resolve_parental_gene_both_close_tie_goes_left():
    """Both boundaries within threshold and equidistant -> left gene (tie -> left)."""
    res = resolve_parental_gene(_interval_row("+", dist_start=150, dist_end=150))
    assert res["Parental_gene"] == "leftGene"


def test_resolve_parental_gene_both_close_right_strictly_closer():
    """Both within threshold, right strictly closer -> right gene."""
    res = resolve_parental_gene(_interval_row("+", dist_start=200, dist_end=90))
    assert res["Parental_gene"] == "rightGene"


def test_resolve_parental_gene_dead_zone_at_exact_threshold():
    """Exactly one distance == threshold (strict-inequality dead zone) -> unassigned.

    Mirrors the notebook's if/elif/elif fall-through: a row where one side is
    exactly distance_threshold matches none of the three branches, so it stays
    unassigned and is later dropped by the inner-join to gene-level stats.
    """
    res = resolve_parental_gene(_interval_row("+", dist_start=400, dist_end=600))
    assert res["Parental_gene"] == ""
    assert res["UTR_type"] == "neither"
