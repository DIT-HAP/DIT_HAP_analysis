#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UTR Insertion Classification
============================

Per-dataset: classifies intergenic insertions that fall within
UTR_DISTANCE_THRESHOLD (400 bp) of a flanking gene boundary as 5'UTR or 3'UTR
insertions (strand-aware), merges them with insertion- and gene-level
curve-fitting statistics, and reports per-insertion um_ratio (insertion DR /
gene DR) and A_ratio (insertion A / gene A). Ported from the deterministic
section of DIT_HAP_pipeline/workflow/notebooks/upstream_and_downstream_analysis.ipynb.

The downstream human-review / plotting notebook lives at
notebooks/domain_analysis/review_utr_insertions.ipynb (Task 9) — this script
only produces the deterministic per-insertion table.

Input
-----
- Insertion-level fitting_results.tsv (MultiIndex [Chr, Coordinate, Strand,
  Target]) — carries per-insertion A / DR. Legacy releases ship the pre-rename
  um/lam headers; normalized to DR/DL on load.
- Insertion-level annotations.tsv(.gz) (same MultiIndex). For intergenic rows,
  Name / Systematic ID / Strand_Interval are pipe-separated "left|right" pairs
  describing the two genes flanking the intergenic interval, and
  Distance_to_region_start / _end are the distances to the LEFT gene's 3' end
  and the RIGHT gene's 5' end respectively.
- Gene-level fitting_results.tsv (Systematic ID, Name, DeletionLibrary_essentiality,
  A, DR, ...). Same legacy um/lam -> DR/DL normalization. Joined to insertions
  on the short gene Name (the annotations' Name field carries short names).

Output
------
- utr_insertion_stats.tsv: one row per UTR insertion (an intergenic insertion
  assigned to a flanking gene), with Parental_gene, UTR_type (5UTR/3UTR),
  Insertion_direction, Distance_to_gene_boundary, insertion/gene A + DR, and the
  derived um_ratio + A_ratio.

Usage
-----
    python classify_utr_insertions.py \\
        --fitting-results .../release/insertion_level/fitting_results.tsv \\
        --annotations .../release/insertion_level/annotations.tsv.gz \\
        --gene-level .../release/gene_level/fitting_results.tsv \\
        --distance-threshold 400 \\
        --output-stats results/utr/{dataset}/utr_insertion_stats.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-20
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Distance (bp) from a gene boundary within which an intergenic insertion is
# treated as a UTR insertion. Byte-faithful to the source notebook's
# Config.distance_threshold.
UTR_DISTANCE_THRESHOLD = 400

# Legacy -> current metric column names (same quirk as
# workflow/scripts/coverage/compute_coverage_stats.py and
# workflow/src/clustering/candidates.py): some datasets' fitting_results.tsv
# still ship the pre-rename um/lam headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# The intergenic-region label used by the annotations table's Type column.
_INTERGENIC = "Intergenic region"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class UTRConfig:
    """Inputs, outputs, and threshold for the UTR insertion classification."""
    fitting_results: Path
    annotations: Path
    gene_level: Path
    output_stats: Path
    distance_threshold: int = UTR_DISTANCE_THRESHOLD

    def validate(self) -> None:
        """Raise ValueError if any input is missing or the threshold is invalid, then ensure the output dir exists."""
        for path in [self.fitting_results, self.annotations, self.gene_level]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if self.distance_threshold <= 0:
            raise ValueError(f"distance_threshold must be positive, got {self.distance_threshold}")
        self.output_stats.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def _normalize_legacy_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Rename legacy um/lam metric columns to DR/DL when the new names are absent."""
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in df.columns and new not in df.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        df = df.rename(columns=rename)
    return df


def load_gene_level(gene_level_path: Path) -> pd.DataFrame:
    """Load gene-level fitting statistics, normalizing legacy um/lam -> DR/DL columns."""
    gene_result = pd.read_csv(gene_level_path, sep="\t")
    return _normalize_legacy_metrics(gene_result)


def load_insertion_level(fitting_results_path: Path, annotations_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load insertion-level fitting results + annotations, both indexed by [Chr, Coordinate, Strand, Target].

    Intergenic annotation rows are uniquely indexed in the release/ files (a
    coordinate falls in at most one intergenic interval), so — unlike the
    coding-annotation dedup that compute_coverage_stats.py needs — no
    duplicate collapsing is required here. fitting_results carries the
    per-insertion A/DR joined in downstream.
    """
    fitting_results = pd.read_csv(fitting_results_path, sep="\t", index_col=[0, 1, 2, 3])
    fitting_results = _normalize_legacy_metrics(fitting_results)
    annotations = pd.read_csv(annotations_path, sep="\t", index_col=[0, 1, 2, 3])
    return fitting_results, annotations


# =============================================================================
# CORE LOGIC — UTR classification primitives (unit-tested)
# =============================================================================
def filter_intergenic_near_gene(df: pd.DataFrame, distance_threshold: int = UTR_DISTANCE_THRESHOLD) -> pd.DataFrame:
    """Keep intergenic rows within ``distance_threshold`` of either flanking gene boundary.

    Byte-faithful to the source notebook's UTR_insertions query:
    ``Type == 'Intergenic region' and ((Distance_to_region_start < thr) or
    (Distance_to_region_end < thr))``. The original index is preserved so
    callers can join back onto insertion-level tables.
    """
    mask = (df["Type"] == _INTERGENIC) & (
        (df["Distance_to_region_start"] < distance_threshold)
        | (df["Distance_to_region_end"] < distance_threshold)
    )
    return df[mask]


def assign_UTR_type(row: pd.Series, distance_threshold: int = UTR_DISTANCE_THRESHOLD) -> str:
    """Classify a single insertion as ``"5UTR"``, ``"3UTR"``, or ``"neither"`` for one gene.

    This is the strand-aware primitive extracted from the notebook's
    assign_UTR_type. It reasons about ONE gene at a time using the gene's
    coding strand (``gene_strand``) and the genomic distances to the
    intergenic region's start/end:

    - ``Distance_to_region_start`` is the distance to the LOW-coordinate
      boundary; ``Distance_to_region_end`` the HIGH-coordinate boundary.
    - For a ``+`` strand gene, its start codon sits at the low-coordinate end,
      so being close to region_start (< threshold) means 5'UTR and close to
      region_end means 3'UTR.
    - For a ``-`` strand gene the orientation flips: close to region_end is the
      5'UTR, close to region_start is the 3'UTR.

    Ambiguity resolution (both boundaries within threshold): the start-boundary
    check is evaluated FIRST, so a ``+`` gene equidistant from both boundaries
    is called ``5UTR`` and a ``-`` gene ``5UTR`` as well (each toward its own 5'
    end). This matches the test's documented "start takes priority" hint.

    When NEITHER boundary is within threshold (which filter_intergenic_near_gene
    normally excludes, but is possible if called directly) the function returns
    ``"neither"`` rather than raising.
    """
    gene_strand = row["gene_strand"]
    dist_start = row["Distance_to_region_start"]
    dist_end = row["Distance_to_region_end"]

    if gene_strand == "+":
        if dist_start < distance_threshold:
            return "5UTR"
        if dist_end < distance_threshold:
            return "3UTR"
    else:  # "-" strand gene: 5' end is at the high-coordinate boundary
        if dist_end < distance_threshold:
            return "5UTR"
        if dist_start < distance_threshold:
            return "3UTR"
    return "neither"


# =============================================================================
# CORE LOGIC — intergenic interval -> parental gene resolution
# =============================================================================
def resolve_parental_gene(row: pd.Series, distance_threshold: int = UTR_DISTANCE_THRESHOLD) -> pd.Series:
    """Map a two-gene intergenic interval row to a single parental gene, then classify it.

    The release/ annotations describe each intergenic interval by the pair of
    genes flanking it: ``Name``/``Systematic ID``/``Strand_Interval`` are
    pipe-separated ``left|right`` values, and the distances are relative to the
    INTERGENIC INTERVAL, not to a single gene:
    ``Distance_to_region_start`` is the distance to the LOW-coordinate end of
    the interval (i.e. the LEFT gene's near boundary) and
    ``Distance_to_region_end`` the distance to the HIGH-coordinate end (the
    RIGHT gene's near boundary).

    Byte-faithful to the notebook's assign_UTR_type, this picks the closer
    flanking gene (left when only region_start is within threshold, right when
    only region_end is, and the strictly closer one — ties to the left — when
    both are). It then delegates the strand-aware 5'/3' call to the tested
    ``assign_UTR_type`` primitive, which expects SINGLE-GENE coordinates. The
    coordinate systems differ, so the interval distances are remapped:

    - LEFT gene: the insertion sits just past that gene's HIGH-coordinate (3')
      edge, so its single-gene ``Distance_to_region_end`` = interval
      ``Distance_to_region_start`` (and its ``Distance_to_region_start`` = inf).
    - RIGHT gene: the insertion sits just before that gene's LOW-coordinate (5')
      edge, so its single-gene ``Distance_to_region_start`` = interval
      ``Distance_to_region_end`` (and its ``Distance_to_region_end`` = inf).

    This remap reproduces the notebook's per-branch 5'/3' assignment exactly,
    including its three-way ``elif`` chain: ``start < thr and end > thr`` (left),
    ``start < thr and end < thr`` (closer of the two, tie -> left), ``start > thr
    and end < thr`` (right). Because each comparison is strict, a row where
    EXACTLY ONE distance equals ``distance_threshold`` falls through all three
    branches in the notebook (it is not excluded by the pre-filter's `<thr`
    disjunction, but no `elif` matches either), leaving the notebook's
    Parental_gene empty — that row is then silently dropped by the inner-join
    to gene-level stats. This function reproduces that dead zone explicitly
    below rather than assigning the row to a gene, so that the total output
    row count matches the notebook byte-for-byte (verified: 30,537 rows on the
    HD_DIT_HAP release set, identical to the notebook).

    Returns a Series with Parental_gene, Parental_gene_id, UTR_type,
    Insertion_direction, and Distance_to_gene_boundary. The signed
    Distance_to_gene_boundary places upstream (5'UTR) insertions at negative
    coordinates and downstream (3'UTR) at positive, matching the notebook.
    """
    left_name, right_name = str(row["Name"]).split("|")
    left_id, right_id = str(row["Systematic ID"]).split("|")
    left_strand, right_strand = str(row["Strand_Interval"]).split("|")
    insertion_strand = row.name[2]  # Strand is the 3rd MultiIndex level
    dist_start = row["Distance_to_region_start"]
    dist_end = row["Distance_to_region_end"]

    start_close = dist_start < distance_threshold
    end_close = dist_end < distance_threshold
    start_far = dist_start > distance_threshold
    end_far = dist_end > distance_threshold

    # Mirror the notebook's exact elif chain (see docstring for the dead zone
    # this leaves at dist == distance_threshold exactly).
    if start_close and end_far:
        gene_name, gene_id, gene_strand, near = left_name, left_id, left_strand, dist_start
        is_left = True
    elif start_close and end_close:
        if dist_start <= dist_end:
            gene_name, gene_id, gene_strand, near = left_name, left_id, left_strand, dist_start
            is_left = True
        else:
            gene_name, gene_id, gene_strand, near = right_name, right_id, right_strand, dist_end
            is_left = False
    elif start_far and end_close:
        gene_name, gene_id, gene_strand, near = right_name, right_id, right_strand, dist_end
        is_left = False
    else:
        # Dead zone (exactly one distance == distance_threshold) or neither
        # boundary within threshold (only reachable if called without the
        # pre-filter). The notebook assigns nothing in either case.
        return pd.Series({
            "Parental_gene": "",
            "Parental_gene_id": "",
            "UTR_type": "neither",
            "Insertion_direction": "",
            "Distance_to_gene_boundary": 0,
        })

    # Remap interval distances to the single-gene convention the primitive uses
    # (see docstring). float("inf") disables the far boundary's branch.
    if is_left:
        primitive_row = pd.Series({
            "gene_strand": gene_strand,
            "Distance_to_region_start": float("inf"),
            "Distance_to_region_end": near,
        })
    else:
        primitive_row = pd.Series({
            "gene_strand": gene_strand,
            "Distance_to_region_start": near,
            "Distance_to_region_end": float("inf"),
        })
    utr_type = assign_UTR_type(primitive_row, distance_threshold)

    # Signed distance to the gene boundary: negative upstream (5'UTR), positive
    # downstream (3'UTR), matching the notebook's Distance_to_gene_boundary.
    signed_distance = near if utr_type == "3UTR" else -near

    insertion_direction = "Forward" if insertion_strand == gene_strand else "Reverse"

    return pd.Series({
        "Parental_gene": gene_name,
        "Parental_gene_id": gene_id,
        "UTR_type": utr_type,
        "Insertion_direction": insertion_direction,
        "Distance_to_gene_boundary": signed_distance,
    })


def classify_utr_insertions(
    fitting_results: pd.DataFrame,
    annotations: pd.DataFrame,
    gene_result: pd.DataFrame,
    distance_threshold: int = UTR_DISTANCE_THRESHOLD,
) -> pd.DataFrame:
    """Assemble the per-insertion UTR table with um_ratio / A_ratio.

    Mirrors the notebook's UTR_insertion_LFCs construction:
      1. filter_intergenic_near_gene on the annotations,
      2. resolve each interval to a parental gene + UTR type,
      3. inner-join to the insertion-level fitting stats (per-insertion A/DR),
      4. inner-join to gene-level stats on the parental gene's short Name,
      5. compute A_ratio = A_insertion / A_gene and um_ratio = DR_insertion / DR_gene.
    """
    utr = filter_intergenic_near_gene(annotations, distance_threshold).copy()
    logger.info(f"Intergenic insertions within {distance_threshold}bp of a gene: {len(utr):,}")

    assigned = utr.apply(lambda row: resolve_parental_gene(row, distance_threshold), axis=1)
    # annotations.tsv.gz already carries an (all-NaN, for intergenic rows)
    # Insertion_direction column; drop it so the resolved Forward/Reverse
    # value below is the only one that survives (a naive concat would create
    # a duplicate `Insertion_direction`/`Insertion_direction.1` pair, and any
    # downstream `df["Insertion_direction"]` lookup would silently see NaN).
    utr = utr.drop(columns=[c for c in assigned.columns if c in utr.columns])
    utr = pd.concat([utr, assigned], axis=1)

    # Attach per-insertion fitting stats (A/DR) via the shared MultiIndex.
    utr = utr.join(fitting_results[["A", "DR"]], how="inner", rsuffix="_fit")
    logger.info(f"UTR insertions with insertion-level fitting stats: {len(utr):,}")

    utr = utr.reset_index()

    # Join gene-level stats on the parental gene's short Name.
    gene_cols = ["Name", "Systematic ID", "DeletionLibrary_essentiality", "A", "DR"]
    gene_cols = [c for c in gene_cols if c in gene_result.columns]
    merged = pd.merge(
        utr,
        gene_result[gene_cols],
        left_on="Parental_gene",
        right_on="Name",
        how="inner",
        suffixes=("_insertion", "_gene"),
    )
    logger.info(f"UTR insertions with matched gene-level stats: {len(merged):,}")

    merged["A_ratio"] = merged["A_insertion"] / merged["A_gene"]
    merged["um_ratio"] = merged["DR_insertion"] / merged["DR_gene"]

    return merged


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: UTRConfig) -> None:
    """Load -> classify UTR insertions -> save per-insertion TSV."""
    config.validate()

    gene_result = load_gene_level(config.gene_level)
    fitting_results, annotations = load_insertion_level(config.fitting_results, config.annotations)

    stats = classify_utr_insertions(fitting_results, annotations, gene_result, config.distance_threshold)
    stats.to_csv(config.output_stats, sep="\t", index=False)

    n5 = int((stats["UTR_type"] == "5UTR").sum())
    n3 = int((stats["UTR_type"] == "3UTR").sum())
    logger.success(
        f"UTR insertions classified: {len(stats):,} total "
        f"({n5:,} 5'UTR, {n3:,} 3'UTR) -> {config.output_stats}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Classify intergenic insertions near gene boundaries as 5'/3' UTR insertions")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Insertion-level fitting_results.tsv")
    parser.add_argument("--annotations", type=Path, required=True, help="Insertion-level annotations.tsv(.gz)")
    parser.add_argument("--gene-level", type=Path, required=True, help="Gene-level fitting_results.tsv")
    parser.add_argument("--distance-threshold", type=int, default=UTR_DISTANCE_THRESHOLD, help="UTR distance threshold in bp (default: 400)")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output per-insertion UTR stats TSV")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = UTRConfig(
            fitting_results=args.fitting_results,
            annotations=args.annotations,
            gene_level=args.gene_level,
            output_stats=args.output_stats,
            distance_threshold=args.distance_threshold,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
