"""
Intra-gene DR Segmentation — Core Logic
=======================================

Library module (no CLI, no file output) behind the ``4a_domain_differences``
stage. Turns the upstream release tables into one positioned DR profile per gene,
and turns the hand-curated terminus spreadsheet into explicit residue intervals.
The segmentation itself lives in :mod:`domain_differences.segmentation`.

The analysis is an *unbiased* scan: it uses no domain annotation at all, only the
position of each in-gene insertion along the spliced CDS, and calls a boundary
wherever two adjacent segments' mean DR differ by more than the merge threshold.
See ``tmp/domain_scan/README.md`` for the measurements behind every choice here,
and ``tmp/domain_scan/explainer/`` for the algorithm walkthrough.

Input
-----
- ``insertion_level/fitting_results.tsv`` — defines the insertion universe
  (unique per Chr/Coordinate/Strand/Target).
- ``insertion_level/annotations.tsv.gz`` — the in-gene call and the CDS positions.
- ``gene_level/fitting_results.tsv`` — gene names and the gene-level DR.
- ``results/1c_annotation/.../gene_annotation_reference.protein.parquet`` — the
  per-gene context columns (FYPOviability, gRNA DR).
- ``resources/curated/non_esential_domain_candidates.xlsx`` — the hand annotation.
- ``resources/external/pombase/{pombase_version}/Protein_features/peptide_stats.tsv``
  — protein lengths, for turning an "affected C-terminal residue" into an interval.

Output
------
Pure functions returning DataFrames; callers own the disk. See
``workflow/scripts/domain_differences/``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Local Imports
from io_table import read_parquet
from release_schema import read_gene_level

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# In-gene insertion criterion, byte-faithful to the source notebook and to
# workflow/scripts/coverage/compute_coverage_stats.py. Deliberately NOT
# release_schema.IN_GENE_FILTER: that one keeps any non-intergenic Feature
# (CDS + intron + UTR), while the scan is defined on coding-gene insertions
# only and at least 5 bp upstream of the stop codon. Switching to the other
# filter changes the universe and every number downstream.
IN_GENE_Q = 'Type == "Coding gene" and Distance_to_stop_codon > 4'

# Release columns the stage carries through the pipeline. ``Coordinate`` is the
# insertion site key (one site = one TTAA, reported once per direction).
INSERTION_COLUMNS = [
    "Systematic ID", "name", "gene_DR", "Chr", "Coordinate", "Strand",
    "Insertion_direction", "Residue_affected", "Residue_frame", "frac", "DR", "DL",
]

# The three sheets of the hand spreadsheet, in the order the plot groups them.
XLSX_SHEETS = [
    ("Non-essential C terminus", "C"),
    ("Non-essential N terminus", "N"),
    ("Non-essential N and C termini", "N+C"),
]

# Terminus groups in gallery order (C first, then N, then both, then unannotated).
TERM_ORDER = ("C", "N", "N+C")
UNANNOTATED = "unannotated"

# Two spreadsheet rows need repairing before the gene-name lookup:
#   "nes2" is not a PomBase gene name anywhere — it is nse2 with the letters
#   transposed, and nse2 (SPAC16A10.06c, 250 aa) matches the recorded 247/15bp
#   exactly (its last two insertions sit at that residue, DR ~ 0).
#   One Gene cell was auto-converted to a date by Excel: "oct1" -> 2025-10-01.
#   oct1 (SPAC1F3.10c, 762 aa) matches the recorded 756/23bp the same way.
XLSX_ALIASES = {"nes2": "nse2", "2025-10-01 00:00:00": "oct1"}


# =============================================================================
# CORE LOGIC — release tables
# =============================================================================
def read_insertions(fitting_results: Path, annotations: Path, gene_level: Path) -> pd.DataFrame:
    """Every in-gene insertion of every gene, positioned along the CDS.

    The insertion universe is ``fitting_results`` (one row per real insertion);
    ``annotations`` supplies the in-gene call and the CDS distances and may carry
    duplicate index entries (a coordinate overlapping two Features). Duplicates
    are collapsed keeping the row that passes the in-gene filter, so one insertion
    is always one point — the same universe the source notebook's merge produced.
    """
    fitting = pd.read_csv(fitting_results, sep="\t", index_col=[0, 1, 2, 3])
    # reindex() and the join below both require a unique key; a duplicated
    # Coordinate/Strand/Target would silently fan out into a cross product.
    if not fitting.index.is_unique:
        raise ValueError(f"insertion fitting_results has {int(fitting.index.duplicated().sum())} duplicate index rows")

    ann = pd.read_csv(annotations, sep="\t", index_col=[0, 1, 2, 3])
    ann = (
        ann.assign(_in_gene=ann.eval(IN_GENE_Q))
        .sort_values("_in_gene", ascending=False, kind="stable")
        .loc[lambda d: ~d.index.duplicated(keep="first")]
        .drop(columns="_in_gene")
    )
    merged = ann.reindex(fitting.index).join(fitting[["DR", "DL"]], how="inner")
    in_gene = merged.query(IN_GENE_Q).dropna(subset=["DR"]).copy()

    span = in_gene["Distance_to_start_codon"] + in_gene["Distance_to_stop_codon"]
    # nucleotide-linear fraction along the CDS, 0 at the start codon and 1 at the
    # stop codon; a zero-length span (both distances 0) stays NaN rather than
    # dividing by zero. Kept alongside Residue_affected because it is the column
    # the prototype's tables were built on.
    in_gene["frac"] = (in_gene["Distance_to_start_codon"] / span.where(span != 0)).clip(lower=0, upper=1)

    genes = read_gene_table(gene_level)
    out = in_gene.reset_index()
    out["name"] = out["Systematic ID"].map(genes["name"]).fillna(out["Systematic ID"])
    out["gene_DR"] = out["Systematic ID"].map(genes["gene_DR"])
    return (
        out[INSERTION_COLUMNS]
        .sort_values(["Systematic ID", "frac"], kind="stable")
        .reset_index(drop=True)
    )


def read_gene_table(gene_level: Path) -> pd.DataFrame:
    """Release gene-level table indexed by Systematic ID: ``name`` and ``gene_DR``.

    A duplicated Systematic ID would turn a caller's ``Series.map`` into a cross
    product, so the last row wins — the rule the previous dict comprehension
    applied silently.
    """
    g = read_gene_level(gene_level)
    table = (g.assign(name=g["Name"].fillna(g["Systematic ID"]))
             .set_index("Systematic ID")[["name", "DR"]]
             .rename(columns={"DR": "gene_DR"}))
    return table[~table.index.duplicated(keep="last")]


def read_protein_lengths(peptide_stats: Path) -> dict[str, int]:
    """{Systematic ID: protein length in residues} from PomBase's peptide_stats.tsv."""
    p = pd.read_csv(peptide_stats, sep="\t")
    return {
        sid: int(n) for sid, n in zip(p["Systematic_ID"], p["Residues"]) if pd.notna(n)
    }


def read_gene_context(annotation: Path) -> tuple[dict[str, str], dict[str, float]]:
    """(FYPOviability, gRNA DR) per gene, from 1c's annotation reference.

    The reference is the project's annotation layer, keyed by systematic id, and
    it is where ``FYPOviability`` and the curated ``gRNA_DR`` already meet — the
    viability call lives in the release table and the gRNA fit in a curated
    spreadsheet, and joining them here would mean this stage re-deriving a join
    that already exists two chapters upstream.
    """
    reference = read_parquet(annotation)
    return (
        {sid: str(v) for sid, v in reference["FYPOviability"].items()},
        {sid: float(v) for sid, v in reference["gRNA_DR"].dropna().items()},
    )


def collapse_sites(insertions: pd.DataFrame) -> pd.DataFrame:
    """One row per insertion SITE, the Forward/Reverse pair averaged into one DR.

    Every site carries the same TTAA insertion in both orientations, so the pair
    is a within-site replicate: over the n>=6 universe 15,436 of 15,844 sites hold
    exactly two rows, and the Forward-minus-Reverse difference has median 0.000 —
    averaging is unbiased and removes that replicate's noise. It also stops the
    segmentation from cutting *between* the two copies of one site.

    Two things that bite later:

    - The grouping key is ``Coordinate``, not ``Residue_affected``: two insertions
      in the same codon at different nucleotides are two genuinely different
      sites that still share a displayed residue, so adjacent segments sharing an
      endpoint residue is not a bug.
    - Every count downstream (``n_sites``, ``min_seg``, the universe floor) is
      counted in SITES when this is used, not insertions.
    """
    d = insertions.assign(
        _fwd=insertions["DR"].where(insertions["Insertion_direction"].eq("Forward")),
        _rev=insertions["DR"].where(insertions["Insertion_direction"].eq("Reverse")),
    )
    out = d.groupby("Coordinate", sort=False).agg(
        **{
            "Systematic ID": ("Systematic ID", "first"),
            "name": ("name", "first"),
            "gene_DR": ("gene_DR", "first"),
            "Chr": ("Chr", "first"),
            "Strand": ("Strand", "first"),
            "Residue_affected": ("Residue_affected", "first"),
            "Residue_frame": ("Residue_frame", "first"),
            "frac": ("frac", "first"),
            "DR": ("DR", "mean"),
            "DR_fwd": ("_fwd", "mean"),
            "DR_rev": ("_rev", "mean"),
            "n_insertions": ("DR", "size"),
        }
    )
    out["n_dir"] = out[["DR_fwd", "DR_rev"]].notna().sum(axis=1)
    # stable sort: tied fractions (the two directions of one site already
    # collapsed, so only same-codon sites) must not reorder run to run
    return out.reset_index().sort_values("frac", kind="stable").reset_index(drop=True)


def ordered_genes(insertions: pd.DataFrame, *, by_site: bool) -> dict[str, pd.DataFrame]:
    """{Systematic ID: that gene's profile ordered along the CDS}, optionally collapsed to sites."""
    return {
        sid: (collapse_sites(d) if by_site else d.reset_index(drop=True))
        for sid, d in insertions.groupby("Systematic ID", sort=True)
    }


def segment_counts(frame: pd.DataFrame, segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """(insertions, sites) behind each segment of an already-ordered profile.

    The two differ by construction in both modes, and both are wanted: ``n_insertions`` is
    the evidence — how many experimental readings a segment rests on — while
    ``n_sites`` is the statistic — how many independent insertion sites the
    segmentation actually used. Counting by segment membership (not by residue
    interval) matters because adjacent segments can share an endpoint residue.

    In a collapsed profile each row IS a site, so the raw insertions behind it are
    the ``n_insertions`` column ``collapse_sites()`` kept; in an insertion profile it is
    the other way round, and the sites are the runs of equal Coordinate.
    """
    coord = frame["Coordinate"].to_numpy()
    runs = np.ones(len(coord), dtype=bool)
    if len(coord) > 1:
        runs[1:] = coord[1:] != coord[:-1]
    per_row_insertions = frame["n_insertions"].to_numpy() if "n_insertions" in frame.columns else None
    counts = []
    for lo, hi in segments:
        n_sites = int(runs[lo:hi].sum())
        n_insertions = int(per_row_insertions[lo:hi].sum()) if per_row_insertions is not None else hi - lo
        counts.append((n_insertions, n_sites))
    return counts


# =============================================================================
# CORE LOGIC — hand-curated non-essential termini
# =============================================================================
def gene_name_to_sid(gene_level: Path) -> dict[str, str]:
    """{lowercased gene name: Systematic ID}, from the same release the scan reads."""
    genes = read_gene_table(gene_level)
    return {str(n).lower(): s for s, n in genes["name"].items()}


def read_hand_regions(
    xlsx: Path, peptide_stats: Path, name_to_sid: dict[str, str]
) -> pd.DataFrame:
    """The hand annotation as explicit per-gene residue intervals.

    A row of the spreadsheet says "the last X residues of this protein are not
    essential" (C terminus sheet) or "the first X are" (N terminus sheet), by
    giving the first/last affected residue. Protein length comes from PomBase
    ``peptide_stats.tsv``; the spreadsheet also records the bp distance from that
    residue to the stop codon, so the reading is checkable:
    ``(protein_len - affected_C + 1) * 3`` must reproduce it to within 3 bp.
    """
    lengths = read_protein_lengths(peptide_stats)
    rows = []
    for sheet, term in XLSX_SHEETS:
        for _, r in pd.read_excel(xlsx, sheet_name=sheet).iterrows():
            gene = str(r["Gene"])
            sid = name_to_sid.get(XLSX_ALIASES.get(gene, gene).lower())
            length = lengths.get(sid)
            n_res = r.get("Affected N terminal residue (aa)")
            c_res = r.get("Affected C terminal residue (aa)")
            regions = []
            if term.startswith("N") and pd.notna(n_res):
                regions.append((1, int(n_res)))
            if term.endswith("C") and pd.notna(c_res) and pd.notna(length):
                regions.append((int(c_res), int(length)))
            dist_stop = r.get("Distance to stop codon (bp)")
            comment = r.get("Comment")
            comment = "" if pd.isna(comment) else str(comment)
            bp_check = ""
            if pd.notna(c_res) and pd.notna(length) and pd.notna(dist_stop):
                implied = (length - int(c_res) + 1) * 3
                bp_check = ("ok" if abs(implied - dist_stop) <= 3
                            else f"off({implied:g} vs {dist_stop:g})")
            rows.append({
                "Gene": gene, "sid": sid, "term": term, "protein_len": length,
                "regions": "|".join(f"{a}-{b}" for a, b in regions),
                "n_regions": len(regions),
                "affected_N": n_res, "affected_C": c_res,
                "dist_start_bp": r.get("Distance to start codon (bp)"),
                "dist_stop_bp": dist_stop, "comment": comment,
                "bp_check": bp_check,
            })
    return pd.DataFrame(rows)


def hand_terms(hand: pd.DataFrame) -> dict[str, str]:
    """{Systematic ID: curated terminus call} — ``C``, ``N`` or ``N+C``; unusable rows dropped.

    No gene appears in two sheets of the spreadsheet (checked: 0 of the 235 usable
    rows share a Systematic ID), so one lookup is enough; a future double entry
    would keep the last.
    """
    usable = hand[hand["sid"].notna() & hand["n_regions"].gt(0)]
    return dict(zip(usable["sid"], usable["term"]))


def hand_region_map(hand: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    """{Systematic ID: [(start, end), ...]} over the genes that have a usable interval.

    ``regions`` packs intervals as ``a-b`` joined by ``|``, so a gene annotated at
    both termini carries two of them.
    """
    usable = hand[hand["sid"].notna() & hand["n_regions"].gt(0)]
    out: dict[str, list[tuple[int, int]]] = {}
    for sid, group in usable.groupby("sid", sort=False):
        out[sid] = [
            tuple(int(v) for v in interval.split("-"))
            for regions in group["regions"]
            for interval in str(regions).split("|")
        ]
    return out
