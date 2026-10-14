"""Tests for the intra-gene DR segmentation stage (4a_domain_differences).

The stage turns the release tables into one positioned DR profile per gene
(workflow/src/domain_differences/core.py) and cuts each profile with a tuned TGUH
(workflow/src/domain_differences/segmentation.py). The primitives below are what
the three compute scripts are built from; the algorithm's own calibration — MAD
consistency, the degenerate cases, the permutation p's — lives in
segmentation._check() and is re-pinned here only where a script's contract
depends on it.

The two failures worth remembering are both in core.py: a count that sliced by
label instead of position, and a region string that packs two intervals.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.domain_differences.core import (
    IN_GENE_Q,
    collapse_sites,
    hand_region_map,
    hand_terms,
    read_gene_table,
    read_insertions,
    segment_counts,
)
from workflow.src.domain_differences.segmentation import (
    merge_threshold,
    noise_sigma,
    region_contrasts,
    tguh_tuned,
    universe_floor,
)


# ---------------------------------------------------------------------------
# fixture builders — the release tables' key columns, in their real order
# ---------------------------------------------------------------------------
_ANN_COLUMNS = [
    "Chr", "Coordinate", "Strand", "Target", "Chr_Interval", "Start_Interval",
    "End_Interval", "Transcript", "Length", "Strand_Interval", "Feature",
    "Systematic ID", "Type", "Accumulated_CDS_bases", "Name", "FYPOviability",
    "DeletionLibrary_essentiality", "ParentalRegion_start", "ParentalRegion_end",
    "ParentalRegion_length", "Distance_to_region_start", "Distance_to_region_end",
    "Fraction_to_region_start", "Fraction_to_region_end", "Distance_to_start_codon",
    "Distance_to_stop_codon", "Fraction_to_start_codon", "Fraction_to_stop_codon",
    "Residue_affected", "Residue_frame", "Insertion_direction",
]
_FIT_COLUMNS = ["Chr", "Coordinate", "Strand", "Target", "A", "DR", "DL"]

_KEY = ("Chr", "Coordinate", "Strand", "Target")


def _write_release(tmp_path: Path, annotations: list[dict], fitting: list[dict]) -> tuple:
    """Write the two release tables the stage reads and return their paths."""
    ann_path = tmp_path / "annotations.tsv.gz"
    pd.DataFrame(annotations, columns=_ANN_COLUMNS).to_csv(
        ann_path, sep="\t", index=False, compression="gzip")
    fit_path = tmp_path / "fitting_results.tsv"
    pd.DataFrame(fitting, columns=_FIT_COLUMNS).to_csv(fit_path, sep="\t", index=False)
    genes_path = tmp_path / "gene_level.tsv"
    sids = sorted({a["Systematic ID"] for a in annotations if isinstance(a["Systematic ID"], str)})
    pd.DataFrame({"Systematic ID": sids, "Name": ["named_gene"] * len(sids),
                  "DR": [-0.75] * len(sids)}).to_csv(genes_path, sep="\t", index=False)
    return fit_path, ann_path, genes_path


def _ann(chr_, coord, strand, target, type_, stop, start, sid, residue, direction,
         feature="CDS"):
    """One annotation row with every column the stage reads; the rest are blanks."""
    row = dict.fromkeys(_ANN_COLUMNS, np.nan)
    row.update({
        "Chr": chr_, "Coordinate": coord, "Strand": strand, "Target": target,
        "Feature": feature, "Systematic ID": sid, "Type": type_,
        "Distance_to_start_codon": start, "Distance_to_stop_codon": stop,
        "Residue_affected": residue, "Residue_frame": residue,
        "Insertion_direction": direction,
    })
    return row


def _fit(chr_, coord, strand, target, dr, dl=3.0):
    """One fitting row; A/DL/t's are not read by the stage beyond their presence."""
    return {"Chr": chr_, "Coordinate": coord, "Strand": strand, "Target": target,
            "A": 1.0, "DR": dr, "DL": dl}


# ---------------------------------------------------------------------------
# IN_GENE_Q — the universe definition
# ---------------------------------------------------------------------------
def test_in_gene_filter_is_the_coding_gene_one():
    """The scan is defined on coding-gene insertions, not on every non-intergenic one.

    Byte-faithful to compute_coverage_stats.py. release_schema.IN_GENE_FILTER keeps
    any non-intergenic Feature (CDS + intron + UTR), which is a different and wider
    universe — swapping them silently changes every number downstream.
    """
    assert IN_GENE_Q == 'Type == "Coding gene" and Distance_to_stop_codon > 4'


# ---------------------------------------------------------------------------
# read_insertions — the universe and the CDS positions
# ---------------------------------------------------------------------------
def test_read_insertions_keeps_only_in_gene_rows(tmp_path):
    """Intergenic insertions and those within 4 bp of the stop codon are dropped."""
    annotations = [
        _ann("I", 100, "+", "T1", "Coding gene", 80, 20, "g1", 7, "Forward"),
        # stop codon distance 4 -> NOT > 4, excluded
        _ann("I", 200, "+", "T2", "Coding gene", 4, 60, "g1", 21, "Forward"),
        _ann("I", 300, "+", "T3", "Intergenic region", 0, np.nan, np.nan, np.nan, "Forward"),
    ]
    fitting = [_fit("I", 100, "+", "T1", -0.5), _fit("I", 200, "+", "T2", -0.2),
               _fit("I", 300, "+", "T3", 0.1)]
    paths = _write_release(tmp_path, annotations, fitting)

    out = read_insertions(*paths)
    assert out["Coordinate"].tolist() == [100]


def test_read_insertions_positions_on_the_cds(tmp_path):
    """frac = Distance_to_start / (start + stop): 0 at the start codon, 1 at the stop."""
    annotations = [
        _ann("I", 100, "+", "T1", "Coding gene", 75, 25, "g1", 9, "Forward"),
        _ann("I", 200, "+", "T2", "Coding gene", 10, 90, "g1", 30, "Forward"),
    ]
    fitting = [_fit("I", 100, "+", "T1", -0.5), _fit("I", 200, "+", "T2", -0.4)]
    out = read_insertions(*_write_release(tmp_path, annotations, fitting)).set_index("Coordinate")
    assert out.loc[100, "frac"] == pytest.approx(0.25)
    assert out.loc[200, "frac"] == pytest.approx(0.9)


def test_read_insertions_duplicate_annotation_keeps_the_in_gene_row(tmp_path):
    """A coordinate overlapping two Features must yield one insertion, not a cross product.

    Annotations may carry duplicate index entries; the collapse keeps whichever row
    passes the in-gene filter, so an insertion whose intergenic and CDS annotations
    share a key is still counted once (as the CDS one).
    """
    annotations = [
        _ann("I", 100, "+", "T1", "Intergenic region", 0, np.nan, np.nan, np.nan, "Forward"),
        _ann("I", 100, "+", "T1", "Coding gene", 50, 30, "g1", 11, "Forward"),
    ]
    fitting = [_fit("I", 100, "+", "T1", -0.5)]
    paths = _write_release(tmp_path, annotations, fitting)

    out = read_insertions(*paths)
    assert len(out) == 1
    assert out["Systematic ID"].iloc[0] == "g1"


def test_read_insertions_rejects_a_duplicated_fitting_index(tmp_path):
    """A duplicated insertion key is an error, not a silent fan-out into a cross product."""
    annotations = [_ann("I", 100, "+", "T1", "Coding gene", 50, 30, "g1", 11, "Forward")]
    fitting = [_fit("I", 100, "+", "T1", -0.5), _fit("I", 100, "+", "T1", -0.5)]
    paths = _write_release(tmp_path, annotations, fitting)

    with pytest.raises(ValueError, match="duplicate index"):
        read_insertions(*paths)


# ---------------------------------------------------------------------------
# collapse_sites — the by-site unit
# ---------------------------------------------------------------------------
def _insertions() -> pd.DataFrame:
    """Two sites: one with both directions, one Forward-only, at ordered coordinates."""
    return pd.DataFrame({
        "Systematic ID": ["g1"] * 3,
        "name": ["named"] * 3,
        "gene_DR": [-0.75] * 3,
        "Chr": ["I"] * 3,
        "Coordinate": [200, 100, 100],
        "Strand": ["+"] * 3,
        "Insertion_direction": ["Forward", "Forward", "Reverse"],
        "Residue_affected": [30, 10, 10],
        "Residue_frame": [30, 10, 10],
        "frac": [0.6, 0.2, 0.2],
        "DR": [-0.40, -1.00, -0.60],
        "DL": [3.0, 3.0, 3.0],
    })


def test_collapse_sites_averages_the_two_directions():
    """One row per Coordinate, DR = mean of the pair; the pair is a within-site replicate."""
    out = collapse_sites(_insertions()).set_index("Coordinate")
    assert len(out) == 2
    assert out.loc[100, "DR"] == pytest.approx(-0.8)   # mean(-1.0, -0.6)
    assert out.loc[200, "DR"] == pytest.approx(-0.4)   # Forward only -> itself


def test_collapse_sites_keeps_both_channels_and_the_evidence_count():
    """DR_fwd / DR_rev carry the pair for the plot, n_insertions the readings behind the site."""
    out = collapse_sites(_insertions()).set_index("Coordinate")
    assert out.loc[100, "DR_fwd"] == pytest.approx(-1.0)
    assert out.loc[100, "DR_rev"] == pytest.approx(-0.6)
    assert out.loc[200, "n_dir"] == 1 and out.loc[100, "n_dir"] == 2
    assert out.loc[100, "n_insertions"] == 2 and out.loc[200, "n_insertions"] == 1


def test_collapse_sites_orders_along_the_cds():
    """Rows come back sorted by frac — the order the segmentation reads as 'along the gene'."""
    out = collapse_sites(_insertions())
    assert out["Coordinate"].tolist() == [100, 200]


# ---------------------------------------------------------------------------
# segment_counts — n_insertions vs n_sites, and why the slice must be positional
# ---------------------------------------------------------------------------
def test_segment_counts_counts_sites_not_insertions():
    """n_sites counts runs of equal Coordinate; n_insertions sums the insertions behind them."""
    frame = pd.DataFrame({
        "Coordinate": [10, 10, 20, 30, 30],
        "n_insertions": [2, 2, 1, 2, 2],
        "DR": [-1.0] * 5,
    })
    # positions 0:3 cover coordinates 10,10,20 -> 2 sites / 5 insertions
    assert segment_counts(frame, [(0, 3), (3, 5)]) == [(5, 2), (4, 1)]


def test_segment_counts_is_positional_not_label_based():
    """Segments are (lo, hi) index RANGES into the profile, so the frame's own labels are irrelevant.

    A profile arrives with a reset index, but the count used to be taken off a
    Series indexed like the frame, which turned ``sites[lo:hi]`` into a LABEL
    slice: the same call on a label-shifted frame counted the wrong points (and
    silently, since the counts stayed plausible).
    """
    frame = pd.DataFrame({
        "Coordinate": [10, 20, 20, 30],
        "n_insertions": [1, 2, 2, 1],
        "DR": [-1.0] * 4,
    }, index=[10, 11, 12, 13])
    assert segment_counts(frame, [(0, 2), (2, 4)]) == [(3, 2), (3, 1)]


def test_segment_counts_without_evidence_column_counts_rows():
    """In insertion mode each row is one insertion, so n_insertions is the segment size."""
    frame = pd.DataFrame({"Coordinate": [10, 10, 20], "DR": [-1.0] * 3})
    assert segment_counts(frame, [(0, 2), (2, 3)]) == [(2, 1), (1, 1)]


# ---------------------------------------------------------------------------
# hand_region_map — the curated regions the gallery shades
# ---------------------------------------------------------------------------
def _hand(rows: list[dict]) -> pd.DataFrame:
    """The columns of domain_hand_regions.tsv that hand_region_map reads."""
    return pd.DataFrame(rows, columns=["Gene", "sid", "term", "regions", "n_regions"])


def test_hand_region_map_parses_a_single_interval():
    """``a-b`` becomes (a, b)."""
    hand = _hand([{"Gene": "mcm4", "sid": "SP1", "term": "C", "regions": "740-870", "n_regions": 1}])
    assert hand_region_map(hand) == {"SP1": [(740, 870)]}


def test_hand_region_map_parses_both_termini():
    """A gene annotated at both ends packs two intervals into ``a-b|c-d``.

    Splitting the whole string on '-' alone reads the junction as a number and
    raises; the '|' is the interval separator and has to come off first.
    """
    hand = _hand([{"Gene": "x", "sid": "SP2", "term": "N+C",
                   "regions": "1-24|740-870", "n_regions": 2}])
    assert hand_region_map(hand) == {"SP2": [(1, 24), (740, 870)]}


def test_hand_region_map_drops_rows_without_a_usable_region():
    """Unmapped rows and rows whose terminus could not be resolved are not shadable."""
    hand = _hand([
        {"Gene": "unknown", "sid": np.nan, "term": "C", "regions": "", "n_regions": 0},
        {"Gene": "y", "sid": "SP3", "term": "C", "regions": "", "n_regions": 0},
        {"Gene": "z", "sid": "SP4", "term": "C", "regions": "10-20", "n_regions": 1},
    ])
    assert hand_region_map(hand) == {"SP4": [(10, 20)]}


# ---------------------------------------------------------------------------
# the threshold and the universe floor
# ---------------------------------------------------------------------------
def test_noise_sigma_is_the_adjacent_difference_scale():
    """sigma = 1.4826 * median|diff| / sqrt(2); a flat profile has no noise to scale by."""
    assert noise_sigma(np.array([-0.5, -0.5, -0.5])) == 0.0
    assert noise_sigma(np.array([-0.5])) == 0.0
    # |diff| constant at d -> median = d, so sigma = 1.0486 * d
    assert noise_sigma(np.array([0.0, 0.1, 0.2])) == pytest.approx(1.4826 * 0.1 / np.sqrt(2))


def test_merge_threshold_is_affine_in_the_genes_own_noise():
    """lam = z*sigma + sigma_add, floored at 0 so a negative z cannot invert the rule."""
    x = np.array([0.0, 0.1, 0.2])
    assert merge_threshold(x, z=1.5, sigma_add=0.1) == pytest.approx(1.5 * noise_sigma(x) + 0.1)
    assert merge_threshold(x, z=-5.0, sigma_add=0.1) == 0.0


def test_universe_floor_is_two_segments_plus_a_margin():
    """2*min_seg is the structural minimum; the extra 2 keeps 2- and 3-point genes out."""
    assert universe_floor(1) == 4
    assert universe_floor(2) == 6


def test_threshold_without_a_floor_merges_everything():
    """z=0, sigma_add=0 means 'merge everything', so a flat series is one segment, not n."""
    assert tguh_tuned(np.full(12, -0.5), 1, z=0.0, sigma_add=0.0) == [(0, 12)]


def test_a_clean_step_is_found_and_a_two_point_gene_never_is():
    """Two flat levels at constant noise split at the change point; a two-point gene cannot.

    The two-point case is the reason the universe floor exists: with one
    difference there is nothing to average, so sigma is that difference and
    lam = 1.5 * 1.0486 * |d| + 0.1 always exceeds |d|.
    """
    rng = np.random.default_rng(0)
    x = np.concatenate([np.full(12, -0.9), np.full(12, -0.1)]) + rng.normal(scale=0.05, size=24)
    assert tguh_tuned(x, 1, z=1.5, sigma_add=0.1) == [(0, 12), (12, 24)]

    two_point = np.array([-0.9, 0.3])
    assert len(tguh_tuned(two_point, 1, z=1.5, sigma_add=0.1)) == 1


# ---------------------------------------------------------------------------
# read_gene_table — the release's gene-level row (name + DR)
# ---------------------------------------------------------------------------
def test_read_gene_table_normalizes_legacy_metric_headers(tmp_path):
    """Older releases ship the metric as `um`; the summary's gene_DR must not come out empty."""
    path = tmp_path / "gene_level.tsv"
    pd.DataFrame({"Systematic ID": ["SP1", "SP2"], "Name": ["mcm4", np.nan],
                  "um": [-0.9, -0.1]}).to_csv(path, sep="\t", index=False)

    table = read_gene_table(path)
    assert table.loc["SP1", "gene_DR"] == pytest.approx(-0.9)
    assert table.loc["SP1", "name"] == "mcm4"
    assert table.loc["SP2", "name"] == "SP2"   # unnamed genes fall back to the ID


def test_read_insertions_carries_gene_dr_onto_every_insertion(tmp_path):
    """gene_DR reaches the summary by riding through the parquet, so it is the release's value."""
    annotations = [_ann("I", 100, "+", "T1", "Coding gene", 80, 20, "g1", 7, "Forward")]
    out = read_insertions(*_write_release(tmp_path, annotations, [_fit("I", 100, "+", "T1", -0.5)]))
    assert out["gene_DR"].tolist() == [-0.75]


# ---------------------------------------------------------------------------
# hand_terms — the curated call the summary scores itself against
# ---------------------------------------------------------------------------
def test_hand_terms_reads_only_usable_rows():
    """A row with no mapped gene or no resolvable interval carries no call."""
    hand = pd.concat([
        _hand([{"Gene": "mcm4", "sid": "SP1", "term": "C", "regions": "740-870", "n_regions": 1}]),
        _hand([{"Gene": "x", "sid": "SP2", "term": "N+C", "regions": "1-24|740-870", "n_regions": 2}]),
        _hand([{"Gene": "y", "sid": np.nan, "term": "C", "regions": "", "n_regions": 0}]),
    ], ignore_index=True)
    assert hand_terms(hand) == {"SP1": "C", "SP2": "N+C"}


# ---------------------------------------------------------------------------
# region_contrasts — the statistic the classification is built on
# ---------------------------------------------------------------------------
def test_region_contrasts_are_immune_to_a_gene_level_offset():
    """A constant added to a whole profile cancels: DR's gene-level offset cannot bias the call."""
    means = np.array([-0.9, -0.3, -1.1])
    weights = np.array([5.0, 4.0, 6.0])
    assert region_contrasts(means + 7.0, weights) == pytest.approx(region_contrasts(means, weights))


def test_region_contrasts_compare_each_segment_with_the_weighted_rest():
    """The rest mean is weighted, so a one-site segment cannot outvote a large one."""
    means = np.array([-0.2, -1.0, -1.0])
    weights = np.array([1.0, 10.0, 10.0])
    # segment 0: -0.2 vs (10*-1 + 10*-1)/20 = -1.0 -> +0.8
    assert region_contrasts(means, weights)[0] == pytest.approx(0.8)
    # segment 1: -1.0 vs (1*-0.2 + 10*-1)/11 = -0.927 -> -0.073
    assert region_contrasts(means, weights)[1] == pytest.approx(-0.0727, abs=1e-3)


def test_a_single_segment_has_no_contrast():
    """Nothing to compare against: a one-segment profile gets 0, not a divide-by-zero."""
    assert region_contrasts(np.array([-0.5]), np.array([3.0])).tolist() == [0.0]
