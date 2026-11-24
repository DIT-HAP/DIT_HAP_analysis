#!/usr/bin/env python3

"""
Figure 3 — Coherence of annotated gene groups in the DR-DL fitness space
========================================================================

Panel A is the coherence method schematic, drawn outside the pipeline and
supplied as a one-page vector PDF: the observed median pairwise distance in
(DR, DL) space, the permutation null, and the resulting coherence z-score.

Panel B sits to the right of A and asks how significance relates to that
z-score across every annotated term: one point per term in the plane of
coherence z (more negative = tighter than random) against the BH-adjusted
p (FDR) on a linear scale. The same two axes carry marginal histograms along
the top and the right, split into all terms (grey) and the coherent cohort
(blue), so the shape of each marginal can be read without leaving the panel.
Dashed guides mark the cohort cut (z < -2, FDR <= 0.05); the coherent cohort
hugs the bottom left.

Panel C drops the coherent terms onto the fitness axes themselves: one point per
term at the mean DR and mean DL of the genes that were scored for its coherence
test. Marker area is the term's scored size, so a small term and a 200-gene
complex can be told apart. The means are computed here from the release's
per-gene DR/DL table over exactly the members listed in the term's
``scored_member_names`` cell; this is the term-level centre of the same gene set
whose median pairwise distance panel B scores, not the geometric median the
coherence stage publishes (that statistic weights genes differently, and the
panel asks for a plain mean).

Panel D returns to the gene level and makes the grid partition visible: every
fitted gene is scattered in raw DR-DL space, two DR cuts and two DL cuts divide
the plane into nine cells, and each cell is labelled with its cluster id and
gene count. The cuts are the grid9 variant's (config/analysis.yaml) translated
back to raw coordinates: DR is not rescaled before clustering (it is only
clamped below -1.3, which cannot move a gene across a cut), while DL is divided
by 10, so the scaled 0.05 / 0.3 cuts are drawn at DL 0.5 / 3.0. Cell membership
is checked against the table's ``cluster`` column, so a table from any other
partition fails instead of drawing a grid that does not describe it.

Panel E shows the same nine clusters' GO biological-process slim enrichment as
a term-by-cluster dot matrix. Rows are the union of each cluster's three most
significant BP slim terms, ordered by the cluster where the term peaks; one dot
marks each significant (term, cluster) cell, its colour is -log10(BH-adjusted
p) and its area is the gene ratio. Cluster columns are the panel D clusters, so
the two panels share one cluster order. The panel is a narrow column: term names
are word-wrapped to a 60 pt line budget, the left and bottom axes/ticks stay
visible, and the gene-ratio size key plus a shortened colour bar sit in a
compact strip to the right of the matrix.

The figure is a two-row block. A (135 x 135 pt) sits above C (150 x 135 pt) in
the first column; B (150 x 135 pt) sits above D (135 x 135 pt) in the second;
E (150 x 277.2 pt) spans both rows to their right. C and D are multipanel
`below=` subtrees of A and B, so they stay under their parents while E is one
tall root panel. E's 277.2 pt axes = 2 x 135 (the two matrix rows) + 7.2 (the
second row's panel-letter band), which makes E's top and bottom align with the
A/C axes. A is embedded as vector content, scaled from its 184.56 pt page box
to its 135 pt square, the same pypdf overlay Figure 1 uses for its schematic:
matplotlib cannot embed a PDF without rasterizing it. The review PNG is
rasterized from the merged PDF rather than from the canvas.

Panels B and C come from the de-duplicated representative set (the stage's
paper-facing display layer), so one term contributes one point and the coherent
highlight in B is exactly the set drawn in C. Panels D and E instead read the
grid9 variant's final clusters and that variant's slim enrichment, because
their subject is the grid partition, not the coherence cohorts.

Every drawing decision lives in this file. The script imports no stage module
from ``workflow/src`` and calls no upstream renderer, so changing Figure 3
cannot change another stage's output.

Input
-----
- results/3a_coherence/<dataset>/dedup/coherence_terms_representatives.tsv
  Columns: source, group_id, group_name, n_scored_members,
  scored_member_names, median_pairwise_distance_z, q_value_coherent.
- <HD_DIT_HAP release>/gene_level/fitting_results.tsv
  Columns: Name, DR, DL (one row per fitted gene).
- results/3b_clustering/<dataset>/grid9/final_clusters.tsv
  Columns: DR, DL, cluster (one row per fitted gene; raw DR/DL plus the
  grid-applied, DR-ascending 1..9 cluster id).
- results/3c_enrichment/<dataset>/grid9/<pombase_version>/go_enrichment_slim.tsv
  Columns: Cluster, namespace, term_id, term, p_fdr, gene_ratio.
- resources/figures/figure3a.pdf
  One page whose page box is its content box (184.56 x 184.56 pt); it is
  scaled to the 135 x 135 pt panel square at assembly time.

Output
------
- results/7a_figures/figure3.pdf                content mounted on an A4 page
- results/7a_figures/figure3.review.png         360 dpi review raster
- results/7a_figures/_work/figure3_canvas.pdf   canvas before the A4 mount (debug)

Usage
-----
    python workflow/scripts/figures/figure3.py \\
        --representatives \\
            results/3a_coherence/HD_DIT_HAP/dedup/coherence_terms_representatives.tsv \\
        --fitting-results <HD_DIT_HAP release>/gene_level/fitting_results.tsv \\
        --grid-clusters results/3b_clustering/HD_DIT_HAP/grid9/final_clusters.tsv \\
        --go-slim \\
            results/3c_enrichment/HD_DIT_HAP/grid9/2026-06-01/go_enrichment_slim.tsv \\
        --schematic resources/figures/figure3a.pdf \\
        --output-pdf results/7a_figures/figure3.pdf \\
        --output-review results/7a_figures/figure3.review.png \\
        --canvas results/7a_figures/_work/figure3_canvas.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-16
Version:  1.1.0
"""

# =============================================================================
# 1. IMPORTS (Standard Library -> Third-Party)
# =============================================================================
import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys

import cnsplots as cns
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib.textpath import TextToPath
import numpy as np
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

# =============================================================================
# 2. CONSTANTS (house style, layout, data vocabulary)
# =============================================================================
# --- House style (self-contained copy; see docs/figures.md) -----------------
HOUSE_PALETTE = "Cell"
FONT_SANS_SERIF = (
    "Arial",
    "Helvetica",
    "Nimbus Sans",
    "Liberation Sans",
    "DejaVu Sans",
)
HOUSE_FIGURE_DPI = 360
FURNITURE_COLOR = cns.GRAY

# --- Canvas and panel geometry, in points -----------------------------------
# The manuscript page is A4 (210 x 297 mm). multipanel sizes the figure by
# width only, so the content canvas is A4 wide minus the page margins and the
# assembly mounts it PAGE_MARGIN_PT in from the page's top-left corner; the
# rest of the page stays blank.
MM_PER_INCH = 25.4
A4_WIDTH_PT = 210.0 / MM_PER_INCH * 72.0  # 595.28
A4_HEIGHT_PT = 297.0 / MM_PER_INCH * 72.0  # 841.89
PAGE_MARGIN_PT = 40.0  # white border at the left, top and right
CANVAS_WIDTH_PT = A4_WIDTH_PT - 2 * PAGE_MARGIN_PT  # 515.28

# All three panels share one row and one 135 pt height: A is the square
# schematic box, B and C are 150 pt wide. The wide-but-short B/C keep their
# plotting areas closer to the data's own proportions (a wide volcano, a
# DR-DL cloud with room for DL but no dead vertical space). A's page box is
# scaled to fit its square (vector, via pypdf); B and C lay their insets out in
# their own pt coordinates.
PANEL_A_SIZE_PT = 135.0
PANEL_B_WIDTH_PT = 150.0
PANEL_B_HEIGHT_PT = 135.0
PANEL_B_LEFT_DECORATION_PT = 28.0  # y tick labels + y label
PANEL_B_BOTTOM_DECORATION_PT = 22.0  # x tick labels + x label
PANEL_B_TOP_HISTOGRAM_HEIGHT_PT = 15.0
PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT = 16.0
PANEL_B_HISTOGRAM_GAP_PT = 3.0
PANEL_B_RIGHT_MARGIN_PT = 4.0
PANEL_B_PLOT_WIDTH_PT = (
    PANEL_B_WIDTH_PT
    - PANEL_B_LEFT_DECORATION_PT
    - PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT
    - PANEL_B_HISTOGRAM_GAP_PT
    - PANEL_B_RIGHT_MARGIN_PT
)  # 99

# Panel C is 150 x 135 pt: the plot plus a narrow size-legend band on its right.
PANEL_C_WIDTH_PT = 150.0
PANEL_C_HEIGHT_PT = 135.0
PANEL_C_LEFT_DECORATION_PT = 28.0
PANEL_C_BOTTOM_DECORATION_PT = 22.0
PANEL_C_TOP_MARGIN_PT = 4.0
PANEL_C_PLOT_WIDTH_PT = 92.0
PANEL_C_LEGEND_GAP_PT = 4.0

# Layout: A sits above C (first column), B above D (second column), and E
# spans both rows to their right. C/D are multipanel subtrees (`below=`), not
# second-row flow panels, so they stay under their parents while E remains one
# tall root panel. E's axes are 277.2 pt = 2 x 135 + 7.2 (the C/D letter band),
# so its top and bottom align with A/C's axes; the nominal "two 135 pt rows"
# total 270 pt is the matrix content, the extra 7.2 pt is the letter band.
PANEL_D_WIDTH_PT = 135.0
PANEL_D_HEIGHT_PT = 135.0
PANEL_D_LEFT_DECORATION_PT = 28.0
PANEL_D_BOTTOM_DECORATION_PT = 22.0
PANEL_D_TOP_MARGIN_PT = 8.0
PANEL_D_PLOT_WIDTH_PT = 107.0
PANEL_D_PLOT_HEIGHT_PT = PANEL_D_HEIGHT_PT - PANEL_D_BOTTOM_DECORATION_PT - PANEL_D_TOP_MARGIN_PT

PANEL_E_WIDTH_PT = 150.0
PANEL_E_HEIGHT_PT = 277.2
PANEL_E_LEFT_LABEL_PT = 62.0
PANEL_E_BOTTOM_DECORATION_PT = 22.0
PANEL_E_TOP_TITLE_PT = 9.0
PANEL_E_PLOT_WIDTH_PT = 66.0
PANEL_E_LEGEND_GAP_PT = 3.0
PANEL_E_LEGEND_WIDTH_PT = 18.0
PANEL_E_LEGEND_LEFT_PT = PANEL_E_LEFT_LABEL_PT + PANEL_E_PLOT_WIDTH_PT + PANEL_E_LEGEND_GAP_PT
PANEL_E_PLOT_HEIGHT_PT = PANEL_E_HEIGHT_PT - PANEL_E_BOTTOM_DECORATION_PT - PANEL_E_TOP_TITLE_PT
PANEL_E_COLORBAR_WIDTH_PT = 4.0
PANEL_E_COLORBAR_HEIGHT_PT = 90.0
PANEL_E_COLORBAR_TOP_PT = 234.0
PANEL_E_SIZE_LEGEND_TOP_PT = 268.0
PANEL_E_SIZE_LEGEND_HEIGHT_PT = 30.0

# Panel B draws the marginal histograms on the same bin grid as the main plot;
# 40 bins keep the count profile smooth at this size without inventing detail.
PANEL_B_BINS = 40

# --- Data vocabulary --------------------------------------------------------
# The de-duplicated representative table is the stage's display layer: one row
# per non-redundant term, with the pooled BH q carried through from combined/.
REPRESENTATIVE_COLUMNS = (
    "source",
    "group_id",
    "group_name",
    "n_scored_members",
    "scored_member_names",
    "median_pairwise_distance_z",
    "q_value_coherent",
)
FITTING_RESULT_COLUMNS = ("Name", "DR", "DL")
SOURCE_COLUMN = "source"
GROUP_ID_COLUMN = "group_id"
GROUP_NAME_COLUMN = "group_name"
SCORED_SIZE_COLUMN = "n_scored_members"
MEMBER_COLUMN = "scored_member_names"
Z_COLUMN = "median_pairwise_distance_z"
Q_COLUMN = "q_value_coherent"
GENE_COLUMN = "Name"
DR_COLUMN = "DR"
DL_COLUMN = "DL"
MEMBER_SEPARATOR = ","

# Panel D reads the grid variant's final clusters; panel E the matching slim table.
CLUSTER_COLUMN = "cluster"
SLIM_CLUSTER_COLUMN = "Cluster"
ONTOLOGY_COLUMN = "namespace"
TERM_ID_COLUMN = "term_id"
TERM_COLUMN = "term"
FDR_COLUMN = "p_fdr"
GENE_RATIO_COLUMN = "gene_ratio"
GRID_CLUSTER_COLUMNS = (DR_COLUMN, DL_COLUMN, CLUSTER_COLUMN)
SLIM_COLUMNS = (
    SLIM_CLUSTER_COLUMN,
    ONTOLOGY_COLUMN,
    TERM_ID_COLUMN,
    TERM_COLUMN,
    FDR_COLUMN,
    GENE_RATIO_COLUMN,
)

# The grid9 cuts in RAW DR/DL space. config/analysis.yaml holds them in scaled
# space: DR is not rescaled (only clamped below -1.3, which cannot move a gene
# across a cut), while DL is divided by 10 before clustering, so the scaled
# 0.05 / 0.3 cuts become DL 0.5 / 3.0 here.
GRID_DR_CUTS = (-0.8, -0.4)
GRID_DL_CUTS = (0.5, 3.0)
SLIM_NAMESPACE = "BP"
SLIM_TERMS_PER_CLUSTER = 3

# The cohort definition is config.coherence's; the values are fixed here so the
# figure draws the cut the table was thresholded on.
COHERENT_Q_MAX = 0.05
COHERENT_Z_CUT = -2.0

# --- Panel vocabulary -------------------------------------------------------
COHERENT_COLOR = cns.get_hexcolors_from_apalette([3], HOUSE_PALETTE)[0]  # blue-grey
BACKGROUND_COLOR = "#A3A3A3"  # house grey, the non-coherent cloud
PLOT_TICK_FONT_PT = 6.0
PLOT_LABEL_FONT_PT = 6.5
LEGEND_FONT_PT = 5.0
LEGEND_TITLE_FONT_PT = 5.5
PANEL_B_POINT_SIZE = 4.0
PANEL_B_COHERENT_SIZE = 6.0
PANEL_B_POINT_ALPHA = 0.50
PANEL_B_COHERENT_ALPHA = 0.80
CUTOFF_LINE_WIDTH_PT = 0.5
CUTOFF_LABEL_FONT_PT = 4.5
X_PADDING = 0.03
FDR_HEADROOM = 1.03
COHERENCE_TICK_STEP = 2.0
FDR_TICK_STEP = 0.2
MEAN_DR_TICK_STEP = 0.2
MEAN_DL_TICK_STEP = 1.0
MEAN_DR_PADDING = 0.04
MEAN_DL_HEADROOM = 1.06
PANEL_B_X_LABEL = "Coherence z-score"
PANEL_B_Y_LABEL = r"FDR (BH-adjusted $p$)"
PANEL_C_X_LABEL = "Mean DR"
PANEL_C_Y_LABEL = "Mean DL"
TERM_SIZE_LEGEND_TITLE = "Genes"
TERM_SIZE_LEGEND_VALUES = (5, 20, 50, 100, 200)
TERM_SIZE_SCALE = 0.35  # marker area in pt^2 per scored gene
TERM_POINT_ALPHA = 0.80
TERM_POINT_EDGE_PT = 0.4
TERM_POINT_EDGE_COLOR = "white"

# --- Panel D: every gene in the DR-DL plane, cut by the grid9 partition ------
GENE_POINT_SIZE = 2.0
GENE_POINT_ALPHA = 0.45
GRID_LINE_WIDTH_PT = 0.6
GRID_CELL_FONT_PT = 4.2
GRID_CUT_FONT_PT = 4.0
GRID_LABEL_COLOR = "#4D4D4D"
GRID_LABEL_X_INSET_PT = 1.5
GRID_LABEL_Y_INSET_PT = 0.6
GRID_LABEL_STROKE_PT = 1.2
GENE_DR_TICK_STEP = 0.5
GENE_DL_TICK_STEP = 2.0
PANEL_D_X_LABEL = "DR"
PANEL_D_Y_LABEL = "DL"

# --- Panel E: GO BP slim term x cluster dot matrix --------------------------
SLIM_COLOR_MAP = "viridis"
SLIM_DOT_EDGE_COLOR = "white"
SLIM_DOT_EDGE_PT = 0.3
SLIM_SIZE_BASE = 2.0
SLIM_SIZE_SCALE = 18.0
SLIM_SIZE_LEGEND_VALUES = (0.02, 0.1, 0.2, 0.3)
SLIM_SIZE_LEGEND_TITLE = "Gene ratio"
SLIM_TITLE = "GO BP slim enrichment"
SLIM_X_LABEL = "Cluster  (dot size = gene ratio)"
SLIM_CBAR_LABEL = r"$-\log_{10}$(FDR)"
SLIM_LABEL_LINE_PT = 60.0
SLIM_TITLE_FONT_PT = 5.5
SLIM_ROW_FONT_PT = 4.5
SLIM_TICK_FONT_PT = 5.0
SLIM_X_TICK_PAD_PT = 2.0
SLIM_X_LABEL_FONT_PT = 5.5
SLIM_CBAR_FONT_PT = 3.4
SLIM_SIZE_LEGEND_FONT_PT = 3.4
SLIM_SIZE_LEGEND_TITLE_FONT_PT = 3.6
SLIM_CBAR_STEP_LARGE = 10.0
SLIM_CBAR_STEP_SMALL = 5.0
SLIM_CBAR_LARGE_CUT = 20.0
SLIM_CBAR_SMALL_CUT = 5.0


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class Figure3Config:
    """Resolved input/output paths for Figure 3."""

    representatives: Path
    fitting_results: Path
    grid_clusters: Path
    go_slim: Path
    schematic: Path
    output_pdf: Path
    output_review: Path
    canvas_pdf: Path

    def validate(self) -> None:
        """Raise if a required input is missing, then create the output dirs."""
        for required in [
            self.representatives,
            self.fitting_results,
            self.grid_clusters,
            self.go_slim,
            self.schematic,
        ]:
            if not required.exists():
                raise ValueError(f"Required input not found: {required}")
        for output in [self.output_pdf, self.output_review, self.canvas_pdf]:
            output.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class TermCloud:
    """Every representative term's z-score and FDR (q), plus the cohort mask."""

    z: np.ndarray
    fdr: np.ndarray
    coherent: np.ndarray


@dataclass(kw_only=True, slots=True, frozen=True)
class CoherentTerm:
    """One coherent term's mean fitness coordinates and scored size."""

    source: str
    group_name: str
    n_scored: int
    mean_dr: float
    mean_dl: float


@dataclass(kw_only=True, slots=True, frozen=True)
class GeneCloud:
    """Every fitted gene's raw DR/DL and its grid-applied cluster id."""

    dr: np.ndarray
    dl: np.ndarray
    cluster: np.ndarray


@dataclass(kw_only=True, slots=True, frozen=True)
class GridCell:
    """One cluster's grid cell: bin indices on each axis and its gene count."""

    cluster: int
    dr_bin: int
    dl_bin: int
    n_genes: int


@dataclass(kw_only=True, slots=True, frozen=True)
class SlimDot:
    """One significant (term, cluster) cell of the BP slim matrix."""

    term_id: str
    cluster: int
    significance: float
    gene_ratio: float


@dataclass(kw_only=True, slots=True, frozen=True)
class SlimPanel:
    """The BP slim matrix: row order, cluster order, dots and the colour ceiling."""

    term_ids: tuple[str, ...]
    terms: tuple[str, ...]
    clusters: tuple[int, ...]
    dots: tuple[SlimDot, ...]
    color_max: float
    terms_per_cluster: int


# =============================================================================
# 4. DATA PREPARATION (pure; schemas are validated at the I/O boundary)
# =============================================================================
def coherent_mask(table: pd.DataFrame) -> np.ndarray:
    """Return the cohort mask defined on the table's q and z columns."""
    is_significant = table[Q_COLUMN].to_numpy(dtype=float) <= COHERENT_Q_MAX
    is_tight = table[Z_COLUMN].to_numpy(dtype=float) < COHERENT_Z_CUT
    return is_significant & is_tight


def build_term_cloud(table: pd.DataFrame) -> TermCloud:
    """Return the panel B arrays for every representative term."""
    z = table[Z_COLUMN].to_numpy(dtype=float)
    q = table[Q_COLUMN].to_numpy(dtype=float)
    if not np.isfinite(z).all() or not np.isfinite(q).all():
        raise ValueError("coherence metrics carry non-finite z or q values")
    if not ((q >= 0) & (q <= 1)).all():
        raise ValueError("q_value_coherent must lie in [0, 1]")
    return TermCloud(
        z=z,
        fdr=q,
        coherent=coherent_mask(table),
    )


def parse_member_names(cell: str) -> tuple[str, ...]:
    """Split a flattened scored-member cell into gene names."""
    return tuple(name.strip() for name in cell.split(MEMBER_SEPARATOR) if name.strip())


def member_frame(table: pd.DataFrame) -> pd.DataFrame:
    """Explode the coherent rows into one (term, member) row per scored gene."""
    coherent = table.loc[coherent_mask(table)].copy()
    if coherent.empty:
        raise ValueError("the representative table has no coherent terms")
    exploded = coherent[
        [SOURCE_COLUMN, GROUP_ID_COLUMN, GROUP_NAME_COLUMN, SCORED_SIZE_COLUMN]
    ].copy()
    exploded[MEMBER_COLUMN] = coherent[MEMBER_COLUMN].map(parse_member_names)
    return exploded.explode(MEMBER_COLUMN).reset_index(drop=True)


def build_coherent_terms(table: pd.DataFrame, fitting: pd.DataFrame) -> tuple[CoherentTerm, ...]:
    """Aggregate each coherent term's scored members to mean DR and mean DL."""
    coherent = table.loc[coherent_mask(table)].copy()
    members = member_frame(table)
    term_keys = [SOURCE_COLUMN, GROUP_ID_COLUMN]

    # The parsed member list must reproduce the size the coherence test used;
    # a mismatch means the flattened cell was truncated or re-formatted.
    expected = coherent.set_index(term_keys, verify_integrity=True)[SCORED_SIZE_COLUMN]
    parsed_counts = members.groupby(term_keys, sort=False).size()
    mismatch = parsed_counts.reindex(expected.index) != expected
    if mismatch.any():
        raise ValueError(
            f"{int(mismatch.sum())} coherent terms list a member count that differs from "
            f"n_scored_members"
        )

    scores = fitting.set_index(GENE_COLUMN)[[DR_COLUMN, DL_COLUMN]]
    if scores.index.duplicated().any():
        raise ValueError("fitting_results.tsv has duplicate gene names")

    joined = members.join(scores, on=MEMBER_COLUMN, how="left")
    missing = joined[DR_COLUMN].isna()
    if missing.any():
        examples = sorted(joined.loc[missing, MEMBER_COLUMN].unique())[:5]
        raise ValueError(f"scored members without a fitted DR/DL: {examples}")

    means = joined.groupby(term_keys, sort=False)[[DR_COLUMN, DL_COLUMN]].mean()
    term_rows = means.reset_index().rename(columns={DR_COLUMN: "mean_dr", DL_COLUMN: "mean_dl"})
    term_rows = term_rows.merge(
        coherent[[SOURCE_COLUMN, GROUP_ID_COLUMN, GROUP_NAME_COLUMN, SCORED_SIZE_COLUMN]],
        on=term_keys,
        how="inner",
        validate="one_to_one",
    )
    return tuple(
        CoherentTerm(
            source=str(row[SOURCE_COLUMN]),
            group_name=str(row[GROUP_NAME_COLUMN]),
            n_scored=int(row[SCORED_SIZE_COLUMN]),
            mean_dr=float(row["mean_dr"]),
            mean_dl=float(row["mean_dl"]),
        )
        for _, row in term_rows.iterrows()
    )


def build_gene_cloud(table: pd.DataFrame) -> GeneCloud:
    """Return every fitted gene's raw DR/DL and validated 1-based cluster id."""
    dr = table[DR_COLUMN].to_numpy(dtype=float)
    dl = table[DL_COLUMN].to_numpy(dtype=float)
    raw_cluster = table[CLUSTER_COLUMN].to_numpy(dtype=float)
    if not np.isfinite(dr).all() or not np.isfinite(dl).all():
        raise ValueError("the grid cluster table carries non-finite DR or DL")
    if not np.equal(raw_cluster, np.floor(raw_cluster)).all():
        raise ValueError("grid cluster ids must be whole numbers")
    cluster = raw_cluster.astype(int)
    if cluster.min() < 1:
        raise ValueError("grid cluster ids must be 1-based")
    return GeneCloud(dr=dr, dl=dl, cluster=cluster)


def build_grid_cells(cloud: GeneCloud) -> tuple[GridCell, ...]:
    """Return each cluster's single grid cell and gene count, failing on a straddle."""
    frame = pd.DataFrame(
        {
            CLUSTER_COLUMN: cloud.cluster,
            "dr_bin": np.digitize(cloud.dr, GRID_DR_CUTS),
            "dl_bin": np.digitize(cloud.dl, GRID_DL_CUTS),
        }
    )
    grouped = frame.groupby(CLUSTER_COLUMN, sort=True)
    unique_bins = grouped[["dr_bin", "dl_bin"]].nunique()
    if (unique_bins > 1).any(axis=None):
        raise ValueError(
            "a cluster's genes straddle a grid cut; the table is not the grid9 partition"
        )
    cell_bins = grouped[["dr_bin", "dl_bin"]].first()
    counts = grouped.size()
    return tuple(
        GridCell(
            cluster=int(cluster),
            dr_bin=int(cell_bins.loc[cluster, "dr_bin"]),
            dl_bin=int(cell_bins.loc[cluster, "dl_bin"]),
            n_genes=int(counts.loc[cluster]),
        )
        for cluster in cell_bins.index.to_numpy()
    )


def slim_candidate_terms(table: pd.DataFrame, top_per_cluster: int) -> pd.DataFrame:
    """Return the union of each cluster's top BP slim terms, ordered by home cluster."""
    bp = table.loc[table[ONTOLOGY_COLUMN] == SLIM_NAMESPACE]
    if bp.empty:
        raise ValueError(f"the slim table has no {SLIM_NAMESPACE} rows")
    ranked = bp.sort_values([SLIM_CLUSTER_COLUMN, FDR_COLUMN, TERM_COLUMN])
    top = ranked.groupby(SLIM_CLUSTER_COLUMN, sort=True).head(top_per_cluster)
    strongest = top.loc[top.groupby(TERM_ID_COLUMN)[FDR_COLUMN].idxmin()]
    strongest = strongest.set_index(TERM_ID_COLUMN)

    terms = top[[TERM_ID_COLUMN, TERM_COLUMN]].drop_duplicates(TERM_ID_COLUMN).copy()
    terms["home_cluster"] = terms[TERM_ID_COLUMN].map(strongest[SLIM_CLUSTER_COLUMN])
    terms["best_fdr"] = terms[TERM_ID_COLUMN].map(strongest[FDR_COLUMN])
    return terms.sort_values(["home_cluster", "best_fdr", TERM_ID_COLUMN])


def build_slim_panel(
    table: pd.DataFrame, clusters: tuple[int, ...], top_per_cluster: int
) -> SlimPanel:
    """Return the BP slim dot matrix over the same cluster order as panel D."""
    terms = slim_candidate_terms(table, top_per_cluster)
    is_selected = (table[ONTOLOGY_COLUMN] == SLIM_NAMESPACE) & table[TERM_ID_COLUMN].isin(
        terms[TERM_ID_COLUMN]
    )
    selected = table.loc[is_selected]
    if not ((selected[FDR_COLUMN] > 0) & (selected[FDR_COLUMN] <= 1)).all():
        raise ValueError("p_fdr must lie in (0, 1] for the -log10 colour scale")

    unknown = set(selected[SLIM_CLUSTER_COLUMN].astype(int)) - set(clusters)
    if unknown:
        raise ValueError(f"slim terms carry clusters absent from the grid table: {sorted(unknown)}")

    dots = tuple(
        SlimDot(
            term_id=str(term_id),
            cluster=int(cluster),
            significance=float(-np.log10(fdr)),
            gene_ratio=float(gene_ratio),
        )
        for term_id, cluster, fdr, gene_ratio in selected[
            [TERM_ID_COLUMN, SLIM_CLUSTER_COLUMN, FDR_COLUMN, GENE_RATIO_COLUMN]
        ].itertuples(index=False, name=None)
    )
    return SlimPanel(
        term_ids=tuple(str(term_id) for term_id in terms[TERM_ID_COLUMN]),
        terms=tuple(str(term) for term in terms[TERM_COLUMN]),
        clusters=clusters,
        dots=dots,
        color_max=float(np.ceil(max(dot.significance for dot in dots))),
        terms_per_cluster=top_per_cluster,
    )


def padded_limits(low: float, high: float, padding: float) -> tuple[float, float]:
    """Pad a data range by a fraction of its span."""
    span = high - low
    return low - padding * span, high + padding * span


def tick_values(low: float, high: float, step: float) -> np.ndarray:
    """Tick the fixed grid that falls inside the plotted range."""
    first = np.ceil(low / step) * step
    last = np.floor(high / step) * step
    return np.arange(first, last + step / 2, step)


def histogram_edges(low: float, high: float, bins: int) -> np.ndarray:
    """Return bin edges over [low, high]."""
    return np.linspace(low, high, bins + 1)


def stacked_counts(
    values: np.ndarray, mask: np.ndarray, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the background and highlighted counts per bin."""
    all_counts, _ = np.histogram(values, bins=edges)
    coherent_counts, _ = np.histogram(values[mask], bins=edges)
    return all_counts - coherent_counts, coherent_counts


# =============================================================================
# 5. RENDERING
# =============================================================================
def apply_house_style() -> None:
    """Configure cnsplots with the project house style (Cell palette, Arial first)."""
    cns.settings.font_sans_serif = FONT_SANS_SERIF
    cns.settings.panel_label_fontname = None
    cns.settings.pdf_fonttype = 42  # TrueType embedding, editable in Illustrator
    cns.settings.svg_fonttype = "none"
    cns.settings.figure_dpi = HOUSE_FIGURE_DPI
    cns.settings.savefig_dpi = HOUSE_FIGURE_DPI
    cns.settings.savefig_transparent = False
    cns.settings.palette_qual = HOUSE_PALETTE
    cns.settings.palette_seq = "viridis"
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE, color_map="viridis")


def draw_figure(
    cloud: TermCloud,
    terms: tuple[CoherentTerm, ...],
    genes: GeneCloud,
    cells: tuple[GridCell, ...],
    slim: SlimPanel,
) -> tuple[cns.multipanel, Axes]:
    """Draw A above C, B above D, and E spanning both rows to their right.

    A, B and E are the row's root panels; C and D are `below=` subtrees of A
    and B. Returns the multipanel itself, not just its figure: cnsplots keeps
    the measured-layout draw callback as a weak reference, so the caller must
    hold the object until the file is written or the panels stay unmeasured.
    Panel A's axes is returned so the caller can read its settled box for the
    pypdf overlay.
    """
    mp = cns.multipanel(max_width=CANVAS_WIDTH_PT, unit="pt")

    schematic_ax = mp.panel(
        "A",
        width=PANEL_A_SIZE_PT,
        height=PANEL_A_SIZE_PT,
        margin_bottom=0,
    )
    blank_container(schematic_ax)

    container_b = mp.panel(
        "B",
        width=PANEL_B_WIDTH_PT,
        height=PANEL_B_HEIGHT_PT,
        margin_bottom=0,
    )
    blank_container(container_b)
    render_panel_b(container_b, cloud, PANEL_B_HEIGHT_PT)

    container_e = mp.panel(
        "E",
        width=PANEL_E_WIDTH_PT,
        height=PANEL_E_HEIGHT_PT,
        margin_right=0,
        margin_bottom=0,
    )
    blank_container(container_e)
    render_panel_e(container_e, slim)

    container_c = mp.panel(
        "C",
        width=PANEL_C_WIDTH_PT,
        height=PANEL_C_HEIGHT_PT,
        below="A",
        margin_bottom=0,
    )
    blank_container(container_c)
    render_panel_c(container_c, terms, PANEL_C_HEIGHT_PT)

    container_d = mp.panel(
        "D",
        width=PANEL_D_WIDTH_PT,
        height=PANEL_D_HEIGHT_PT,
        below="B",
        margin_bottom=0,
    )
    blank_container(container_d)
    render_panel_d(container_d, genes, cells)
    return mp, schematic_ax


def blank_container(container: Axes) -> None:
    """Clear the ticks, turn the axis off and give the container a pt coordinate system."""
    # Clear the ticks before turning the axis off: multipanel measures y-axis
    # decoration even on an axis-off axes, and placeholder ticks would widen
    # the panel's left reserve.
    container.set_xticks([])
    container.set_yticks([])
    container.set_axis_off()


def render_panel_b(container: Axes, cloud: TermCloud, height_pt: float) -> None:
    """Draw the z-vs-FDR cloud and its two marginal histograms."""
    container.set_xlim(0, PANEL_B_WIDTH_PT)
    container.set_ylim(0, height_pt)
    plot_height_pt = (
        height_pt
        - PANEL_B_BOTTOM_DECORATION_PT
        - PANEL_B_TOP_HISTOGRAM_HEIGHT_PT
        - PANEL_B_HISTOGRAM_GAP_PT
    )
    x_limits = padded_limits(float(cloud.z.min()), float(cloud.z.max()), X_PADDING)
    y_limits = (0.0, float(cloud.fdr.max()) * FDR_HEADROOM)

    plot_ax = container.inset_axes(
        (
            PANEL_B_LEFT_DECORATION_PT / PANEL_B_WIDTH_PT,
            PANEL_B_BOTTOM_DECORATION_PT / height_pt,
            PANEL_B_PLOT_WIDTH_PT / PANEL_B_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    top_ax = container.inset_axes(
        (
            PANEL_B_LEFT_DECORATION_PT / PANEL_B_WIDTH_PT,
            (PANEL_B_BOTTOM_DECORATION_PT + plot_height_pt + PANEL_B_HISTOGRAM_GAP_PT) / height_pt,
            PANEL_B_PLOT_WIDTH_PT / PANEL_B_WIDTH_PT,
            PANEL_B_TOP_HISTOGRAM_HEIGHT_PT / height_pt,
        )
    )
    right_ax = container.inset_axes(
        (
            (PANEL_B_LEFT_DECORATION_PT + PANEL_B_PLOT_WIDTH_PT + PANEL_B_HISTOGRAM_GAP_PT)
            / PANEL_B_WIDTH_PT,
            PANEL_B_BOTTOM_DECORATION_PT / height_pt,
            PANEL_B_RIGHT_HISTOGRAM_WIDTH_PT / PANEL_B_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    render_panel_b_plot(plot_ax, cloud, x_limits, y_limits)
    render_panel_b_marginals(top_ax, right_ax, plot_ax, cloud, x_limits, y_limits)


def render_panel_b_plot(
    ax: Axes,
    cloud: TermCloud,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    """Scatter every term and mark the cohort cut on both axes."""
    is_background = ~cloud.coherent
    ax.scatter(
        cloud.z[is_background],
        cloud.fdr[is_background],
        s=PANEL_B_POINT_SIZE,
        color=BACKGROUND_COLOR,
        alpha=PANEL_B_POINT_ALPHA,
        linewidths=0,
    )
    ax.scatter(
        cloud.z[cloud.coherent],
        cloud.fdr[cloud.coherent],
        s=PANEL_B_COHERENT_SIZE,
        color=COHERENT_COLOR,
        alpha=PANEL_B_COHERENT_ALPHA,
        linewidths=0,
    )
    ax.axvline(
        COHERENT_Z_CUT,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
        zorder=0,
    )
    ax.axhline(
        COHERENT_Q_MAX,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
        zorder=0,
    )
    ax.text(
        x_limits[1],
        COHERENT_Q_MAX,
        f"FDR = {COHERENT_Q_MAX:g}",
        ha="right",
        va="bottom",
        fontsize=CUTOFF_LABEL_FONT_PT,
        color=FURNITURE_COLOR,
    )
    ax.text(
        COHERENT_Z_CUT,
        y_limits[1],
        f"$z = {COHERENT_Z_CUT:g}$",
        ha="left",
        va="top",
        rotation=90,
        fontsize=CUTOFF_LABEL_FONT_PT,
        color=FURNITURE_COLOR,
    )
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xticks(tick_values(*x_limits, COHERENCE_TICK_STEP))
    ax.set_yticks(tick_values(y_limits[0], y_limits[1], FDR_TICK_STEP))
    ax.set_xlabel(PANEL_B_X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(PANEL_B_Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    add_panel_b_legend(ax, cloud)


def add_panel_b_legend(ax: Axes, cloud: TermCloud) -> None:
    """Name the two point styles in the empty upper-left corner."""
    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=3.0, color=BACKGROUND_COLOR),
        Line2D([], [], marker="o", linestyle="", markersize=3.5, color=COHERENT_COLOR),
    ]
    labels = [
        f"All terms (n = {len(cloud.z):,})",
        f"Coherent (n = {int(cloud.coherent.sum()):,})",
    ]
    ax.legend(
        handles,
        labels,
        loc="upper left",
        frameon=False,
        fontsize=LEGEND_FONT_PT,
        handlelength=0.8,
        handletextpad=0.4,
        labelspacing=0.3,
        borderaxespad=0.2,
    )


def render_panel_b_marginals(
    top_ax: Axes,
    right_ax: Axes,
    plot_ax: Axes,
    cloud: TermCloud,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    """Draw the z and FDR distributions along their own axes."""
    z_edges = histogram_edges(x_limits[0], x_limits[1], PANEL_B_BINS)
    fdr_edges = histogram_edges(y_limits[0], y_limits[1], PANEL_B_BINS)
    z_background, z_coherent = stacked_counts(cloud.z, cloud.coherent, z_edges)
    q_background, q_coherent = stacked_counts(cloud.fdr, cloud.coherent, fdr_edges)

    # Empty bins are skipped: a zero-height bar still renders as a hairline at
    # the baseline, which would draw a false rule across the histogram.
    z_occupied = (z_background + z_coherent) > 0
    z_highlighted = z_coherent > 0
    top_ax.bar(
        z_edges[:-1][z_occupied],
        z_background[z_occupied],
        width=np.diff(z_edges)[z_occupied],
        align="edge",
        color=BACKGROUND_COLOR,
        linewidth=0,
    )
    top_ax.bar(
        z_edges[:-1][z_highlighted],
        z_coherent[z_highlighted],
        bottom=z_background[z_highlighted],
        width=np.diff(z_edges)[z_highlighted],
        align="edge",
        color=COHERENT_COLOR,
        linewidth=0,
    )
    q_occupied = (q_background + q_coherent) > 0
    q_highlighted = q_coherent > 0
    right_ax.barh(
        fdr_edges[:-1][q_occupied],
        q_background[q_occupied],
        height=np.diff(fdr_edges)[q_occupied],
        align="edge",
        color=BACKGROUND_COLOR,
        linewidth=0,
    )
    right_ax.barh(
        fdr_edges[:-1][q_highlighted],
        q_coherent[q_highlighted],
        left=q_background[q_highlighted],
        height=np.diff(fdr_edges)[q_highlighted],
        align="edge",
        color=COHERENT_COLOR,
        linewidth=0,
    )

    top_ax.axvline(
        COHERENT_Z_CUT,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
    )
    right_ax.axhline(
        COHERENT_Q_MAX,
        color=FURNITURE_COLOR,
        linestyle="--",
        linewidth=CUTOFF_LINE_WIDTH_PT,
    )
    top_ax.sharex(plot_ax)
    right_ax.sharey(plot_ax)
    top_ax.set_ylim(0, float((z_background + z_coherent).max() * 1.08))
    right_ax.set_xlim(0, float((q_background + q_coherent).max() * 1.08))
    hide_frame(top_ax)
    hide_frame(right_ax)


def hide_frame(ax: Axes) -> None:
    """Remove every frame element from a marginal axes, sharing axes included."""
    ax.set_axis_off()
    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_panel_c(container: Axes, terms: tuple[CoherentTerm, ...], height_pt: float) -> None:
    """Draw one mean-DR/mean-DL bubble per coherent term."""
    container.set_xlim(0, PANEL_C_WIDTH_PT)
    container.set_ylim(0, height_pt)
    plot_height_pt = height_pt - PANEL_C_BOTTOM_DECORATION_PT - PANEL_C_TOP_MARGIN_PT
    ax = container.inset_axes(
        (
            PANEL_C_LEFT_DECORATION_PT / PANEL_C_WIDTH_PT,
            PANEL_C_BOTTOM_DECORATION_PT / height_pt,
            PANEL_C_PLOT_WIDTH_PT / PANEL_C_WIDTH_PT,
            plot_height_pt / height_pt,
        )
    )
    mean_dr = np.array([term.mean_dr for term in terms])
    mean_dl = np.array([term.mean_dl for term in terms])
    sizes = np.array([term.n_scored for term in terms], dtype=float)
    x_limits = padded_limits(float(mean_dr.min()), float(mean_dr.max()), MEAN_DR_PADDING)
    y_limits = (0.0, float(mean_dl.max()) * MEAN_DL_HEADROOM)

    ax.scatter(
        mean_dr,
        mean_dl,
        s=sizes * TERM_SIZE_SCALE,
        color=COHERENT_COLOR,
        alpha=TERM_POINT_ALPHA,
        edgecolors=TERM_POINT_EDGE_COLOR,
        linewidths=TERM_POINT_EDGE_PT,
    )
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xticks(tick_values(*x_limits, MEAN_DR_TICK_STEP))
    ax.set_yticks(tick_values(y_limits[0], y_limits[1], MEAN_DL_TICK_STEP))
    ax.set_xlabel(PANEL_C_X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(PANEL_C_Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    add_panel_c_size_legend(ax)


def add_panel_c_size_legend(ax: Axes) -> None:
    """Draw the marker-area key to the right of the panel C plot."""
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor=COHERENT_COLOR,
            markeredgecolor=TERM_POINT_EDGE_COLOR,
            markeredgewidth=TERM_POINT_EDGE_PT,
            markersize=np.sqrt(size * TERM_SIZE_SCALE),
        )
        for size in TERM_SIZE_LEGEND_VALUES
    ]
    labels = [f"{size:,}" for size in TERM_SIZE_LEGEND_VALUES]
    legend_x = (PANEL_C_PLOT_WIDTH_PT + PANEL_C_LEGEND_GAP_PT) / PANEL_C_PLOT_WIDTH_PT
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(legend_x, 1.0),
        frameon=False,
        fontsize=LEGEND_FONT_PT,
        title=TERM_SIZE_LEGEND_TITLE,
        title_fontsize=LEGEND_TITLE_FONT_PT,
        handletextpad=0.5,
        labelspacing=0.6,
        borderaxespad=0.0,
    )


def render_panel_d(container: Axes, cloud: GeneCloud, cells: tuple[GridCell, ...]) -> None:
    """Scatter every gene in DR-DL space and draw the nine-cell grid with its counts."""
    container.set_xlim(0, PANEL_D_WIDTH_PT)
    container.set_ylim(0, PANEL_D_HEIGHT_PT)
    ax = container.inset_axes(
        (
            PANEL_D_LEFT_DECORATION_PT / PANEL_D_WIDTH_PT,
            PANEL_D_BOTTOM_DECORATION_PT / PANEL_D_HEIGHT_PT,
            PANEL_D_PLOT_WIDTH_PT / PANEL_D_WIDTH_PT,
            PANEL_D_PLOT_HEIGHT_PT / PANEL_D_HEIGHT_PT,
        )
    )
    x_limits = padded_limits(float(cloud.dr.min()), float(cloud.dr.max()), X_PADDING)
    y_limits = (0.0, float(cloud.dl.max()) * MEAN_DL_HEADROOM)

    ax.scatter(
        cloud.dr,
        cloud.dl,
        s=GENE_POINT_SIZE,
        color=BACKGROUND_COLOR,
        alpha=GENE_POINT_ALPHA,
        linewidths=0,
    )
    draw_grid_cuts(ax, x_limits, y_limits)
    add_grid_cell_labels(ax, cells, x_limits, y_limits)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xticks(tick_values(*x_limits, GENE_DR_TICK_STEP))
    ax.set_yticks(tick_values(*y_limits, GENE_DL_TICK_STEP))
    ax.set_xlabel(PANEL_D_X_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.set_ylabel(PANEL_D_Y_LABEL, fontsize=PLOT_LABEL_FONT_PT)
    ax.tick_params(axis="both", labelsize=PLOT_TICK_FONT_PT)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_grid_cuts(ax: Axes, x_limits: tuple[float, float], y_limits: tuple[float, float]) -> None:
    """Draw the two DR and two DL cuts, with their values at the top and right edges."""
    for cut in GRID_DR_CUTS:
        ax.axvline(cut, color=FURNITURE_COLOR, linewidth=GRID_LINE_WIDTH_PT, zorder=2)
        ax.text(
            cut,
            y_limits[1],
            f"{cut:g}",
            ha="right",
            va="top",
            fontsize=GRID_CUT_FONT_PT,
            color=FURNITURE_COLOR,
            zorder=4,
            path_effects=label_path_effects(),
        )
    for cut in GRID_DL_CUTS:
        ax.axhline(cut, color=FURNITURE_COLOR, linewidth=GRID_LINE_WIDTH_PT, zorder=2)
        ax.text(
            x_limits[1],
            cut,
            f"{cut:g}",
            ha="right",
            va="bottom",
            fontsize=GRID_CUT_FONT_PT,
            color=FURNITURE_COLOR,
            zorder=4,
            path_effects=label_path_effects(),
        )


def add_grid_cell_labels(
    ax: Axes,
    cells: tuple[GridCell, ...],
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    """Write each cell's cluster id and gene count on the cell's left edge."""
    dr_edges = (x_limits[0], *GRID_DR_CUTS, x_limits[1])
    dl_edges = (y_limits[0], *GRID_DL_CUTS, y_limits[1])
    x_offset = (x_limits[1] - x_limits[0]) * GRID_LABEL_X_INSET_PT / PANEL_D_PLOT_WIDTH_PT
    y_offset = (y_limits[1] - y_limits[0]) * GRID_LABEL_Y_INSET_PT / PANEL_D_PLOT_HEIGHT_PT

    for cell in cells:
        x = dr_edges[cell.dr_bin] + x_offset
        # The bottom strip holds the DL = 0 pile-up: anchor its label at the
        # strip's top so the text stays above the piled-up genes; the taller
        # cells are labelled at mid-height where their cloud is thin.
        y = (
            dl_edges[1] - y_offset
            if cell.dl_bin == 0
            else (dl_edges[cell.dl_bin] + dl_edges[cell.dl_bin + 1]) / 2
        )
        vertical = "top" if cell.dl_bin == 0 else "center"
        ax.text(
            x,
            y,
            f"{cell.cluster}: {cell.n_genes:,}",
            ha="left",
            va=vertical,
            fontsize=GRID_CELL_FONT_PT,
            color=GRID_LABEL_COLOR,
            zorder=4,
            path_effects=label_path_effects(),
        )


def label_path_effects() -> list[pe.AbstractPathEffect]:
    """Return the white halo that keeps panel D's grid labels legible over the cloud."""
    return [pe.withStroke(linewidth=GRID_LABEL_STROKE_PT, foreground="white")]


def wrap_slim_terms(terms: tuple[str, ...]) -> list[str]:
    """Wrap every term label to the label column at the row font size."""
    font_path = findfont(FontProperties(family=list(FONT_SANS_SERIF)))
    font = FontProperties(fname=font_path, size=SLIM_ROW_FONT_PT)
    engine = TextToPath()
    return [_wrap_words(term, font, engine) for term in terms]


def _wrap_words(text: str, font: FontProperties, engine: TextToPath) -> str:
    """Greedy word-wrap one label to SLIM_LABEL_LINE_PT."""
    lines: list[str] = []
    current = ""
    for word in text.split():
        candidate = f"{current} {word}".strip()
        width = engine.get_text_width_height_descent(candidate, font, False)[0]
        if current and width > SLIM_LABEL_LINE_PT:
            lines.append(current)
            current = word
        else:
            current = candidate
    lines.append(current)
    return "\n".join(lines)


def render_panel_e(container: Axes, panel: SlimPanel) -> None:
    """Draw the BP slim term x cluster matrix with a size key and a short colour bar."""
    container.set_xlim(0, PANEL_E_WIDTH_PT)
    container.set_ylim(0, PANEL_E_HEIGHT_PT)
    ax = container.inset_axes(
        (
            PANEL_E_LEFT_LABEL_PT / PANEL_E_WIDTH_PT,
            PANEL_E_BOTTOM_DECORATION_PT / PANEL_E_HEIGHT_PT,
            PANEL_E_PLOT_WIDTH_PT / PANEL_E_WIDTH_PT,
            PANEL_E_PLOT_HEIGHT_PT / PANEL_E_HEIGHT_PT,
        )
    )
    x_of = {cluster: index for index, cluster in enumerate(panel.clusters)}
    row_of = {term_id: index for index, term_id in enumerate(panel.term_ids)}
    x = np.array([x_of[dot.cluster] for dot in panel.dots])
    y = np.array([len(panel.term_ids) - 1 - row_of[dot.term_id] for dot in panel.dots])
    significance = np.array([dot.significance for dot in panel.dots])
    gene_ratio = np.array([dot.gene_ratio for dot in panel.dots])
    normalizer = Normalize(vmin=0.0, vmax=panel.color_max)

    ax.scatter(
        x,
        y,
        s=SLIM_SIZE_BASE + gene_ratio * SLIM_SIZE_SCALE,
        c=significance,
        cmap=SLIM_COLOR_MAP,
        norm=normalizer,
        edgecolors=SLIM_DOT_EDGE_COLOR,
        linewidths=SLIM_DOT_EDGE_PT,
        zorder=2,
    )
    ax.set_xlim(-0.5, len(panel.clusters) - 0.5)
    ax.set_ylim(-0.5, len(panel.term_ids) - 0.5)
    ax.set_xticks(range(len(panel.clusters)))
    ax.set_xticklabels([str(cluster) for cluster in panel.clusters])
    ax.set_yticks(range(len(panel.term_ids)))
    ax.set_yticklabels(list(reversed(wrap_slim_terms(panel.terms))))
    ax.set_xlabel(SLIM_X_LABEL, fontsize=SLIM_X_LABEL_FONT_PT, labelpad=2.0)
    ax.tick_params(
        axis="x", labelsize=SLIM_TICK_FONT_PT, length=2.0, width=0.4, pad=SLIM_X_TICK_PAD_PT
    )
    ax.tick_params(axis="y", labelsize=SLIM_ROW_FONT_PT, length=1.5, width=0.4, pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.5)

    panel_title = f"{SLIM_TITLE}  ·  top {panel.terms_per_cluster}/cluster"
    container.text(
        0,
        PANEL_E_HEIGHT_PT - 1,
        panel_title,
        ha="left",
        va="top",
        fontsize=SLIM_TITLE_FONT_PT,
    )
    add_slim_size_legend(container)
    add_slim_colorbar(container, normalizer, panel.color_max)


def add_slim_size_legend(container: Axes) -> None:
    """Draw the gene-ratio dot-size key above the colour bar."""
    height = PANEL_E_SIZE_LEGEND_HEIGHT_PT
    legend_ax = container.inset_axes(
        (
            PANEL_E_LEGEND_LEFT_PT / PANEL_E_WIDTH_PT,
            (PANEL_E_SIZE_LEGEND_TOP_PT - height) / PANEL_E_HEIGHT_PT,
            PANEL_E_LEGEND_WIDTH_PT / PANEL_E_WIDTH_PT,
            height / PANEL_E_HEIGHT_PT,
        )
    )
    legend_ax.set_xlim(0, PANEL_E_LEGEND_WIDTH_PT)
    legend_ax.set_ylim(0, height)
    legend_ax.set_axis_off()
    legend_ax.text(
        0,
        height - 1,
        SLIM_SIZE_LEGEND_TITLE,
        ha="left",
        va="top",
        fontsize=SLIM_SIZE_LEGEND_TITLE_FONT_PT,
    )
    for index, value in enumerate(SLIM_SIZE_LEGEND_VALUES):
        y = height - 9.0 - index * 6.5
        legend_ax.scatter(
            [3.0],
            [y],
            s=SLIM_SIZE_BASE + value * SLIM_SIZE_SCALE,
            color=GRID_LABEL_COLOR,
            linewidths=0,
        )
        legend_ax.text(
            6.0,
            y,
            f"{value:g}",
            ha="left",
            va="center",
            fontsize=SLIM_SIZE_LEGEND_FONT_PT,
        )


def add_slim_colorbar(container: Axes, normalizer: Normalize, color_max: float) -> None:
    """Draw a short -log10(FDR) colour bar below the size key."""
    cax = container.inset_axes(
        (
            PANEL_E_LEGEND_LEFT_PT / PANEL_E_WIDTH_PT,
            (PANEL_E_COLORBAR_TOP_PT - PANEL_E_COLORBAR_HEIGHT_PT) / PANEL_E_HEIGHT_PT,
            PANEL_E_COLORBAR_WIDTH_PT / PANEL_E_WIDTH_PT,
            PANEL_E_COLORBAR_HEIGHT_PT / PANEL_E_HEIGHT_PT,
        )
    )
    ColorbarBase(cax, cmap=SLIM_COLOR_MAP, norm=normalizer, orientation="vertical")
    ticks = significance_ticks(color_max)
    cax.set_yticks(ticks)
    cax.set_yticklabels([f"{tick:g}" for tick in ticks])
    cax.tick_params(labelsize=SLIM_CBAR_FONT_PT, length=1.0, width=0.3, pad=1.0)
    cax.yaxis.set_label_position("right")
    cax.set_ylabel(SLIM_CBAR_LABEL, fontsize=SLIM_CBAR_FONT_PT, rotation=270, labelpad=4)


def significance_ticks(color_max: float) -> np.ndarray:
    """Tick the colour scale at 10, 5 or 1 unit intervals, depending on its ceiling."""
    if color_max >= SLIM_CBAR_LARGE_CUT:
        step = SLIM_CBAR_STEP_LARGE
    elif color_max >= SLIM_CBAR_SMALL_CUT:
        step = SLIM_CBAR_STEP_SMALL
    else:
        step = 1.0
    return tick_values(0.0, color_max, step)


# =============================================================================
# 6. ASSEMBLY
# =============================================================================
def axes_box_pt(ax: Axes, figure: Figure) -> tuple[float, float, float, float]:
    """Return an axes box as (x, y, width, height) in points from the page's bottom-left."""
    width_in, height_in = figure.get_size_inches()
    width_pt, height_pt = width_in * 72, height_in * 72
    box = ax.get_position()
    return box.x0 * width_pt, box.y0 * height_pt, box.width * width_pt, box.height * height_pt


def assemble_page(
    canvas_pdf: Path,
    schematic_pdf: Path,
    output_pdf: Path,
    box_pt: tuple[float, float, float, float],
) -> None:
    """Overlay the schematic on the canvas, then mount the content atop an A4 page."""
    # The box is in PDF points already: bbox_inches=None keeps the page equal to
    # the figure, so the figure fractions flip straight into the page's origin.
    x_pt, y_pt, box_width_pt, box_height_pt = box_pt
    # The canvas page must belong to a writer before its content can grow: pypdf
    # 7 refuses content replacements on reader-owned pages.
    writer = PdfWriter(clone_from=str(canvas_pdf))
    canvas_page = writer.pages[0]
    schematic_page = PdfReader(str(schematic_pdf)).pages[0]

    source_width = float(schematic_page.mediabox.width)
    source_height = float(schematic_page.mediabox.height)
    scale = min(box_width_pt / source_width, box_height_pt / source_height)

    x_pt += (box_width_pt - source_width * scale) / 2
    y_pt += (box_height_pt - source_height * scale) / 2

    canvas_page.merge_transformed_page(
        schematic_page, Transformation().scale(scale).translate(x_pt, y_pt)
    )

    # A4 page: the canvas keeps its designed size, inset from the page's
    # top-left corner; the rest of the page stays blank.
    content_height_pt = float(canvas_page.mediabox.height)
    bottom_pt = A4_HEIGHT_PT - PAGE_MARGIN_PT - content_height_pt
    page = writer.add_blank_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
    page.merge_transformed_page(canvas_page, Transformation().translate(PAGE_MARGIN_PT, bottom_pt))
    writer.remove_page(0)  # the content now lives on the A4 page
    writer.write(str(output_pdf))


def rasterize_review(pdf: Path, review_png: Path) -> None:
    """Rasterize the assembled PDF at the house review resolution."""
    subprocess.run(
        [
            "pdftocairo",
            "-png",
            "-r",
            str(HOUSE_FIGURE_DPI),
            "-singlefile",
            str(pdf),
            str(review_png.with_suffix("")),
        ],
        check=True,
    )


# =============================================================================
# 7. IMPERATIVE SHELL (boundary validation, orchestration, CLI)
# =============================================================================
def read_tsv_columns(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    """Read only the required TSV columns and fail once if any is missing."""
    table = pd.read_csv(path, sep="\t", usecols=lambda column: column in required)
    missing_columns = set(required) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {sorted(missing_columns)}")
    return table


def run(config: Figure3Config) -> tuple[float, float]:
    """Build Figure 3 end to end and return the page's (width, height) in points."""
    config.validate()
    apply_house_style()

    representatives = read_tsv_columns(config.representatives, REPRESENTATIVE_COLUMNS)
    fitting = read_tsv_columns(config.fitting_results, FITTING_RESULT_COLUMNS)
    grid_clusters = read_tsv_columns(config.grid_clusters, GRID_CLUSTER_COLUMNS)
    go_slim = read_tsv_columns(config.go_slim, SLIM_COLUMNS)
    cloud = build_term_cloud(representatives)
    terms = build_coherent_terms(representatives, fitting)
    genes = build_gene_cloud(grid_clusters)
    cells = build_grid_cells(genes)
    gene_clusters = tuple(int(cluster) for cluster in np.unique(genes.cluster))
    slim = build_slim_panel(go_slim, gene_clusters, SLIM_TERMS_PER_CLUSTER)
    logger.info(
        f"Panel B: {len(cloud.z):,} terms ({int(cloud.coherent.sum()):,} coherent); "
        f"panel C: {len(terms):,} coherent terms"
    )
    logger.info(
        f"Panel D: {len(genes.dr):,} genes over {len(cells)} grid cells; "
        f"panel E: {len(slim.term_ids)} BP slim terms, {len(slim.dots)} dots"
    )

    mp, schematic_ax = draw_figure(cloud, terms, genes, cells, slim)
    assert mp.fig is not None
    # cnsplots sets savefig.bbox to "tight", which would crop the page to the
    # drawn artists and shift panel A's box away from the overlay below; keep
    # the designed canvas. The save draws the figure, which is also what
    # settles multipanel's measured layout before the box is read back.
    cns.savefig(config.canvas_pdf, fig=mp.fig, bbox_inches=None)
    assemble_page(
        config.canvas_pdf,
        config.schematic,
        config.output_pdf,
        axes_box_pt(schematic_ax, mp.fig),
    )
    rasterize_review(config.output_pdf, config.output_review)

    logger.success(
        f"Wrote {config.output_pdf} (A4 {A4_WIDTH_PT:.2f} x {A4_HEIGHT_PT:.2f} pt) "
        f"and {config.output_review}"
    )
    return A4_WIDTH_PT, A4_HEIGHT_PT


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build the composite Figure 3")
    parser.add_argument(
        "--representatives",
        type=Path,
        required=True,
        help="De-duplicated representative terms TSV",
    )
    parser.add_argument(
        "--fitting-results",
        type=Path,
        required=True,
        help="Release per-gene fitting results (Name + DR + DL)",
    )
    parser.add_argument(
        "--grid-clusters",
        type=Path,
        required=True,
        help="grid9 final_clusters.tsv (DR + DL + cluster)",
    )
    parser.add_argument(
        "--go-slim",
        type=Path,
        required=True,
        help="GO slim enrichment TSV (BP rows are selected)",
    )
    parser.add_argument("--schematic", type=Path, required=True, help="Panel A schematic PDF")
    parser.add_argument("--output-pdf", type=Path, required=True, help="Final figure PDF")
    parser.add_argument("--output-review", type=Path, required=True, help="Review PNG")
    parser.add_argument("--canvas", type=Path, required=True, help="Canvas PDF before the A4 mount")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging"
    )
    return parser.parse_args()


def setup_logger(verbose: bool) -> None:
    """Configure the loguru sink at INFO or DEBUG."""
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO")


def main() -> int:
    """Build the config, render the figure, and report the outcome."""
    args = parse_args()
    setup_logger(args.verbose)
    config = Figure3Config(
        representatives=args.representatives,
        fitting_results=args.fitting_results,
        grid_clusters=args.grid_clusters,
        go_slim=args.go_slim,
        schematic=args.schematic,
        output_pdf=args.output_pdf,
        output_review=args.output_review,
        canvas_pdf=args.canvas,
    )
    try:
        run(config)
    except Exception as e:
        logger.exception(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
