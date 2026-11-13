"""Tests for the final-figure stage: layout, frames, and PDF assembly per figure."""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cnsplots")
pytest.importorskip("pypdf")

import cnsplots as cns
from matplotlib import pyplot as plt
from matplotlib.transforms import Bbox
from pypdf import PdfReader
from workflow.scripts.figures import figure_s1
from workflow.scripts.figures.figure1 import (
    A4_HEIGHT_PT,
    A4_WIDTH_PT,
    AGG_PANEL_HEIGHT_PT,
    AGG_PANEL_WIDTH_PT,
    AGGREGATION_GENE_ID,
    CANVAS_WIDTH_PT,
    COVERAGE_CLASSES,
    COVERAGE_PANEL_WIDTH_PT,
    DIT_HAP_SERIES,
    FACS_CONTAINER_HEIGHT_PT,
    FACS_CONTAINER_WIDTH_PT,
    GENERATION_COLUMNS,
    GENERATIONS,
    ROW_PANEL_SIZE_PT,
    SPIKEIN_POINT_SIZE,
    Figure1Config,
    apply_house_style,
    axes_box_pt,
    build_aggregation_data,
    build_coverage_frame,
    build_facs_comparison_frame,
    build_spikein_frame,
    build_trajectory_frame,
    draw_figure,
    render_facs_panel,
    render_lfc_trajectory_panel,
    render_spikein_panel,
    run,
)


@pytest.fixture(autouse=True)
def _house_style() -> None:
    """Draw with the production style: layout constants are calibrated at 360 dpi."""
    apply_house_style()


# =============================================================================
# Fixtures
# =============================================================================
def _coverage_table() -> pd.DataFrame:
    """One row per (dataset, FYPO class), the shape compare_lfc_coverage writes."""
    rows = [("LD_haploid", cls, 0, 0, 100, 0.0, 0.50) for cls in COVERAGE_CLASSES] + [
        ("HD_DIT_HAP", cls, 0, 0, 200, 0.0, 0.80) for cls in COVERAGE_CLASSES
    ]
    return pd.DataFrame(
        rows,
        columns=["dataset", "fypo_class", "n_total", "n_dr", "n_lfc", "pct_dr", "pct_lfc"],
    )


def _spikein_table() -> pd.DataFrame:
    """Two insertion sites across three non-zero dilutions, plus the reference."""
    records = [("Spikein0", "Spike-in Insertion 1", -9.4, float("-inf"))]
    for site in (1, 2):
        for dilution, reading in ((1, -8.0), (2, -6.0), (3, -4.0)):
            records.append(
                (
                    f"Spikein{dilution}",
                    f"Spike-in Insertion {site}",
                    float(dilution),
                    reading - site,
                )
            )
    return pd.DataFrame(
        records,
        columns=["Sample", "Name", "Relative_Dilution_Ratio", "Relative_Read_Ratio"],
    )


def _write_schematic(path: Path) -> None:
    """Write a one-page stand-in for panel A whose page box is exactly 200 x 50 pt."""
    fig = cns.figure(width=200, height=50, unit="pt")
    fig.text(0.5, 0.5, "SCHEMATIC", ha="center", va="center")
    # cnsplots sets savefig.bbox="tight"; pin the page to the canvas so the
    # assembly sees the size this fixture advertises.
    fig.savefig(path, format="pdf", bbox_inches=Bbox.from_bounds(0.0, 0.0, *fig.get_size_inches()))


def _replicate_table(x: str, y: str) -> pd.DataFrame:
    """A strictly positive pair of read-count columns plus one unplottable row.

    Independent draws, because the density estimate needs a two-dimensional
    cloud: collinear points give gaussian_kde a singular covariance.
    """
    rng = np.random.default_rng(0)
    table = pd.DataFrame({x: rng.integers(1, 1000, size=40), y: rng.integers(1, 1000, size=40)})
    # Zero on either side cannot appear on a log axis; the panel must drop it.
    return pd.concat([table, pd.DataFrame({x: [0], y: [5]})], ignore_index=True)


def _gene_lfc_table() -> pd.DataFrame:
    """Two inviable and two viable genes over the release's five generations."""
    return pd.DataFrame(
        {
            "Systematic ID": ["G1", "G2", "G3", "G4", "G5"],
            "Name": ["cut8", "fal1", "leu1", "rpl3601", "g5"],
            "FYPOviability": ["inviable", "inviable", "viable", "viable", "unknown"],
            "YES0": [0.0, 0.0, 0.0, 0.0, 0.0],
            "YES1": [-1.0, -2.0, -0.1, 0.1, -9.0],
            "YES2": [-3.0, -4.0, -0.1, 0.1, -9.0],
            "YES3": [-5.0, -6.0, -0.2, 0.0, -9.0],
            "YES4": [-6.0, -7.0, -0.2, -0.1, -9.0],
        }
    )


def _gene_lfc_with_aggregation() -> pd.DataFrame:
    """The panel E gene table plus the worked gene panel D pins by systematic ID."""
    worked_gene = pd.DataFrame(
        {
            "Systematic ID": [AGGREGATION_GENE_ID],
            "Name": ["cdc7"],
            "FYPOviability": ["inviable"],
            "YES0": [0.0],
            "YES1": [-1.4],
            "YES2": [-2.32],
            "YES3": [-3.56],
            "YES4": [-4.2],
        }
    )
    return pd.concat([_gene_lfc_table(), worked_gene], ignore_index=True)


def _insertion_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Two TTAA sites on one plus-strand gene; site 1,200 carries both orientations."""
    keys = ["Chr", "Coordinate", "Strand", "Target"]
    lfc = pd.DataFrame(
        [
            ("I", 1200, "+", "TTAA", 0.0, -2.0, -4.0, -6.0, -7.0),
            ("I", 1200, "-", "TTAA", 0.0, -4.0, -4.0, -6.0, -7.0),
            ("I", 1800, "+", "TTAA", 0.0, 0.5, 0.2, 0.1, 0.0),
        ],
        columns=[*keys, *GENERATION_COLUMNS],
    )
    padj = pd.DataFrame(
        [
            ("I", 1200, "+", "TTAA", 1.0, 1e-4, 1e-5, 1e-6, 1e-6),
            ("I", 1200, "-", "TTAA", 1.0, 1e-5, 1e-5, 1e-6, 1e-6),
            ("I", 1800, "+", "TTAA", 1.0, 0.8, 0.5, 0.3, 0.9),
        ],
        columns=[*keys, *GENERATION_COLUMNS],
    )
    weights = pd.DataFrame(
        [
            ("I", coordinate, strand, "TTAA", timepoint, weight, AGGREGATION_GENE_ID)
            for coordinate, strand, weight in ((1200, "+", 0.4), (1200, "-", 0.2), (1800, "+", 0.4))
            for timepoint in GENERATION_COLUMNS
        ],
        columns=[*keys, "Timepoint", "Weight", "Systematic ID"],
    )
    return lfc, padj, weights


def _gene_intervals() -> pd.DataFrame:
    """A 1,000 bp ORF with one intron, in PomBase's primary-transcript BED shape."""
    segments = (("CDS", 1000, 1400), ("intron", 1400, 1500), ("CDS", 1500, 2000))
    return pd.DataFrame(
        [
            {
                "Chr": "I",
                "Systematic ID": AGGREGATION_GENE_ID,
                "Feature": feature,
                "Start": start,
                "End": end,
                "Strand": "+",
                "ParentalRegion_start": 1000,
                "ParentalRegion_end": 2000,
            }
            for feature, start, end in segments
        ]
    )


def _facs_table() -> pd.DataFrame:
    """The panel's four genes over three timepoints; cut8 carries both strains."""
    blocks = (
        ("cut8", "Strain#1", "LV1"),
        ("cut8", "Strain#2", "LV2"),
        ("fal1", "Strain#1", "LV3"),
        ("leu1", "Strain#1", "LV4"),
        ("rpl3601", "Strain#1", "LV5"),
    )
    records = [
        (f"SYS_{gene}", gene, strain_id, strain, f"YES{index}", generation, lfc)
        for gene, strain, strain_id in blocks
        for index, (generation, lfc) in enumerate(((0.0, 0.0), (2.0, 1.0), (4.0, 2.0)))
    ]
    return pd.DataFrame(
        records,
        columns=[
            "Systematic ID",
            "Name",
            "Strain",
            "Strain_alias",
            "Timepoints",
            "Generations",
            "LFC",
        ],
    )


def _annotation_table() -> pd.DataFrame:
    """Two genes share three in-gene insertions; intergenic rows must not count."""
    return pd.DataFrame(
        {
            "Type": ["Coding gene", "Coding gene", "Coding gene", "Intergenic region"],
            "Systematic ID": ["G1", "G1", "G2", "G1"],
        }
    )


@pytest.fixture
def figure_s1_inputs(tmp_path: Path) -> figure_s1.FigureS1Config:
    """Write fixture parquets and return a Figure S1 config pointing at them."""
    pbl_pbr = tmp_path / "pbl_pbr.parquet"
    tech = tmp_path / "tech.parquet"
    bio = tmp_path / "bio.parquet"
    annotations = tmp_path / "annotations.tsv"
    _replicate_table("PBL", "PBR").to_parquet(pbl_pbr)
    _replicate_table("Reads_1", "Reads_2").to_parquet(tech)
    _replicate_table("Reads_1", "Reads_2").to_parquet(bio)
    _annotation_table().to_csv(annotations, sep="\t", index=False)
    return figure_s1.FigureS1Config(
        pbl_pbr=pbl_pbr,
        tech=tech,
        bio=bio,
        annotations=annotations,
        output_pdf=tmp_path / "out" / "figure_s1.pdf",
        output_review=tmp_path / "out" / "figure_s1.review.png",
        canvas_pdf=tmp_path / "out" / "_work" / "figure_s1_canvas.pdf",
    )


@pytest.fixture
def figure1_inputs(tmp_path: Path) -> Figure1Config:
    """Write fixture inputs and return a config pointing at them."""
    coverage = tmp_path / "lfc_coverage_comparison.tsv"
    spikein = tmp_path / "spike_in_stats.tsv"
    gene_lfc = tmp_path / "gene_lfc.tsv"
    facs = tmp_path / "facs.xlsx"
    insertion_lfc = tmp_path / "insertion_lfc.tsv"
    insertion_padj = tmp_path / "insertion_padj.tsv"
    insertion_weights = tmp_path / "insertion_weights.tsv"
    gene_intervals = tmp_path / "coding_gene_primary_transcripts.bed"
    schematic = tmp_path / "figure1a.pdf"

    _coverage_table().to_csv(coverage, sep="\t", index=False)
    _spikein_table().to_csv(spikein, sep="\t", index=False)
    _gene_lfc_with_aggregation().to_csv(gene_lfc, sep="\t", index=False)
    _facs_table().to_excel(facs, index=False)
    site_lfc, site_padj, site_weights = _insertion_tables()
    site_lfc.to_csv(insertion_lfc, sep="\t", index=False)
    site_padj.to_csv(insertion_padj, sep="\t", index=False)
    site_weights.to_csv(insertion_weights, sep="\t", index=False)
    intervals = _gene_intervals().rename(columns={"Chr": "#Chr"})
    intervals.to_csv(gene_intervals, sep="\t", index=False)
    _write_schematic(schematic)

    return Figure1Config(
        lfc_comparison=coverage,
        spikein=spikein,
        gene_lfc=gene_lfc,
        facs=facs,
        insertion_lfc=insertion_lfc,
        insertion_padj=insertion_padj,
        insertion_weights=insertion_weights,
        gene_intervals=gene_intervals,
        schematic=schematic,
        output_pdf=tmp_path / "out" / "figure1.pdf",
        output_review=tmp_path / "out" / "figure1.review.png",
        canvas_pdf=tmp_path / "out" / "_work" / "figure1_canvas.pdf",
    )


# =============================================================================
# Layout
# =============================================================================
def _panel_boxes(mp: cns.multipanel, labels: tuple[str, ...]) -> dict[str, tuple[float, ...]]:
    """Settle the flow layout and return each panel's point box by label."""
    mp.fig.canvas.draw()
    boxes = {}
    for label in labels:
        ax = mp.get_axes(label)
        assert ax is not None
        boxes[label] = axes_box_pt(ax, mp.fig)
    return boxes


def _panel_label_left_pt(mp: cns.multipanel, label: str) -> float:
    """Left edge of the panel letter multipanel drew for a panel, in points."""
    ax = mp.get_axes(label)
    assert ax is not None
    text = next(text for text in ax.texts if text.get_text() == label)
    return text.get_window_extent(mp.fig.canvas.get_renderer()).x0 * 72 / mp.fig.dpi


def test_draw_figure_places_the_two_column_grid():
    """A above; B/C/E stack on the left, D aggregation and F FACS on the right."""
    coverage = build_coverage_frame(_coverage_table())
    spikein = build_spikein_frame(_spikein_table())
    trajectory = build_trajectory_frame(_gene_lfc_table())
    facs_comparison = build_facs_comparison_frame(_facs_table(), _gene_lfc_table())
    insertion_lfc, insertion_padj, insertion_weights = _insertion_tables()
    aggregation = build_aggregation_data(
        _gene_lfc_with_aggregation(),
        insertion_lfc,
        insertion_padj,
        insertion_weights,
        _gene_intervals(),
    )
    mp, schematic_ax = draw_figure(
        coverage,
        spikein,
        trajectory,
        facs_comparison,
        aggregation,
        schematic_size_pt=(400.0, 120.0),
    )

    boxes = _panel_boxes(mp, ("B", "C", "D", "E", "F"))  # also settles the flow layout
    a_x, a_y, a_width, a_height = axes_box_pt(schematic_ax, mp.fig)
    b_x, b_y, b_width, b_height = boxes["B"]
    c_x, c_y, c_width, c_height = boxes["C"]
    d_x, d_y, d_width, d_height = boxes["D"]
    e_x, e_y, e_width, e_height = boxes["E"]
    f_x, f_y, f_width, f_height = boxes["F"]

    assert mp.fig.get_size_inches()[0] * 72 == pytest.approx(CANVAS_WIDTH_PT)
    assert 0 < a_x <= 10  # only the panel letter's column sits left of A
    assert (a_width, a_height) == pytest.approx((400.0, 120.0))  # actual size, no stretch
    assert (b_width, b_height) == pytest.approx((COVERAGE_PANEL_WIDTH_PT, ROW_PANEL_SIZE_PT))
    assert (c_width, c_height) == pytest.approx((ROW_PANEL_SIZE_PT, ROW_PANEL_SIZE_PT))
    assert (e_width, e_height) == pytest.approx((ROW_PANEL_SIZE_PT, ROW_PANEL_SIZE_PT))
    assert (d_width, d_height) == pytest.approx((AGG_PANEL_WIDTH_PT, AGG_PANEL_HEIGHT_PT))
    assert (f_width, f_height) == pytest.approx((FACS_CONTAINER_WIDTH_PT, FACS_CONTAINER_HEIGHT_PT))
    assert b_y + b_height < a_y  # the grid opens below A
    assert d_x > b_x + b_width  # D sits to the right of B
    assert c_y + c_height < b_y  # C sits below B
    assert e_x > c_x + c_width  # E sits to the right of C
    assert f_x > e_x + e_width  # F sits to the right of E
    assert e_y == pytest.approx(c_y, abs=1.5)  # C and E share the row
    assert f_y + f_height < d_y  # row 3 opens below D's row
    assert f_y + f_height > c_y + c_height  # F is the taller panel in row 3
    assert _panel_label_left_pt(mp, "A") == pytest.approx(_panel_label_left_pt(mp, "B"), abs=0.5)
    assert _panel_label_left_pt(mp, "B") == pytest.approx(_panel_label_left_pt(mp, "C"), abs=0.5)
    assert _panel_label_left_pt(mp, "E") > _panel_label_left_pt(mp, "C")
    assert _panel_label_left_pt(mp, "F") > _panel_label_left_pt(mp, "E")


# =============================================================================
# Frame builders
# =============================================================================
def test_build_coverage_frame_maps_labels_and_percent():
    """Only the two compared libraries survive, renamed and rescaled."""
    frame = build_coverage_frame(_coverage_table())

    assert len(frame) == 2 * len(COVERAGE_CLASSES)
    assert set(frame["dataset_label"]) == {"Traditional Tn-Seq", "DIT-HAP"}
    assert set(frame["class_label"]) == {
        "Inviable",
        "Viable",
        "Depends on\nconditions",
        "Unknown",
        "Total",
    }
    assert frame["coverage_pct"].max() == pytest.approx(80.0)


def test_build_coverage_frame_rejects_missing_class():
    """A table without every FYPO class fails at the boundary, not in the figure."""
    table = _coverage_table()
    table = table[table["fypo_class"] != "unknown"]

    with pytest.raises(ValueError, match="missing FYPO classes"):
        build_coverage_frame(table)


def test_trajectory_summary_keeps_only_the_two_classes_and_maps_generations():
    """Unknown stays out; each surviving gene contributes one row per generation."""
    frame = build_trajectory_frame(_gene_lfc_table())

    assert set(frame["FYPOviability"]) == {"inviable", "viable"}
    assert sorted(frame["generation"].unique()) == list(GENERATIONS)
    assert len(frame) == 4 * len(GENERATIONS)


def test_lfc_panel_draws_one_line_per_gene():
    """Each gene is its own curve; each class is one legend entry, not one curve."""
    frame = build_trajectory_frame(_gene_lfc_table())
    fig, ax = plt.subplots()

    render_lfc_trajectory_panel(ax, frame)

    gene_lines = ax.lines[:-1]  # the last artist is the zero reference line
    assert len(gene_lines) == 4  # two inviable + two viable; unknown stays out
    assert all(np.asarray(line.get_ydata()).shape == (len(GENERATIONS),) for line in gene_lines)
    legend = ax.get_legend()
    assert legend is not None
    assert [text.get_text() for text in legend.get_texts()] == [
        "Inviable (n = 2)",
        "Viable (n = 2)",
    ]
    fig.canvas.draw()
    legend_box = legend.get_window_extent()
    axes_box = ax.get_window_extent()
    assert legend_box.y0 < axes_box.y1  # the lower row overlaps the plot's top
    assert (legend_box.y0 + legend_box.y1) / 2 == pytest.approx(axes_box.y1)
    assert (legend_box.x0 + legend_box.x1) / 2 == pytest.approx((axes_box.x0 + axes_box.x1) / 2)
    plt.close(fig)


def test_facs_comparison_aligns_signs_and_keeps_both_strains():
    """FACS values flip onto the release's depletion-negative convention."""
    frame = build_facs_comparison_frame(_facs_table(), _gene_lfc_table())

    assert set(frame["gene"]) == {"cut8", "fal1", "leu1", "rpl3601"}
    facs = frame[frame["series"] != DIT_HAP_SERIES]
    assert (facs["lfc"] <= 0).all()  # the fixture's FACS LFC is depletion-positive
    assert set(facs["series"]) == {"FACS Strain#1", "FACS Strain#2"}
    assert frame["generation"].max() == pytest.approx(GENERATIONS[-1])


def test_facs_comparison_rejects_gene_absent_from_dit_hap():
    """A shown FACS gene with no release row would silently lose its DIT-HAP curve."""
    gene_lfc = _gene_lfc_table()
    gene_lfc.loc[gene_lfc["Name"] == "fal1", "Name"] = "not_in_lfc"

    with pytest.raises(ValueError, match="missing FACS genes"):
        build_facs_comparison_frame(_facs_table(), gene_lfc)


def test_facs_comparison_keeps_only_the_shown_genes():
    """Genes outside FACS_SHOWN_GENES stay out of the panel even when supplied."""
    facs = _facs_table()
    facs.loc[facs["Name"] == "fal1", "Name"] = "orc6"
    gene_lfc = _gene_lfc_table()
    gene_lfc.loc[gene_lfc["Name"] == "fal1", "Name"] = "orc6"

    frame = build_facs_comparison_frame(facs, gene_lfc)

    assert "orc6" not in set(frame["gene"])
    assert set(frame["gene"]) == {"cut8", "leu1", "rpl3601"}


def test_aggregation_data_reproduces_the_gene_level_curve():
    """Two orientations collapse to one site and the weighted sites sum to the release row."""
    insertion_lfc, insertion_padj, insertion_weights = _insertion_tables()
    data = build_aggregation_data(
        _gene_lfc_with_aggregation(),
        insertion_lfc,
        insertion_padj,
        insertion_weights,
        _gene_intervals(),
    )

    assert list(data.site_lfc.index) == [200.0, 800.0]  # bp from the start codon
    assert data.span_bp == pytest.approx(1000.0)
    # Site 1,200's two orientations (-2 and -4 at YES1) average by their weights.
    assert data.site_lfc.loc[200.0, "YES1"] == pytest.approx(-8.0 / 3.0)
    assert data.gene_curve.to_list() == pytest.approx([0.0, -1.4, -2.32, -3.56, -4.2])
    assert data.high_weight.to_list() == [True, False]
    assert (data.site_lfc * data.site_weight).sum().to_list() == pytest.approx(
        data.gene_curve.to_list()
    )


def test_aggregation_data_rejects_a_missing_gene_row():
    """Panel D is pinned to one gene; a gene-level table without it cannot draw."""
    insertion_lfc, insertion_padj, insertion_weights = _insertion_tables()

    with pytest.raises(ValueError, match="exactly one"):
        build_aggregation_data(
            _gene_lfc_table(),
            insertion_lfc,
            insertion_padj,
            insertion_weights,
            _gene_intervals(),
        )


def test_facs_panel_draws_one_mini_axes_per_gene_with_shared_limits():
    """Every gene gets its own axes; all share one y scale and the zero line."""
    frame = build_facs_comparison_frame(_facs_table(), _gene_lfc_table())
    fig, container = plt.subplots()

    render_facs_panel(container, frame)

    mini_axes = container.child_axes
    assert [ax.title.get_text() for ax in mini_axes] == ["cut8", "fal1", "leu1", "rpl3601"]
    assert len(mini_axes[0].lines) == 4  # DIT-HAP, two FACS strains, zero line
    assert len(mini_axes[1].lines) == 3  # fal1 has no Strain#2 block
    assert mini_axes[0].get_ylim() == mini_axes[1].get_ylim()
    legend = container.get_legend()
    assert legend is not None
    assert [text.get_text() for text in legend.get_texts()] == [
        "DIT-HAP",
        "FACS Strain#1",
        "FACS Strain#2",
    ]
    plt.close(fig)


def test_spikein_panel_uses_the_configured_point_size_and_legend_row():
    """Points take panel C's size constant; its legend is one row above the axes."""
    frame = build_spikein_frame(_spikein_table())
    fig, ax = plt.subplots()

    render_spikein_panel(ax, frame)

    point_sizes = {
        tuple(collection.get_sizes())
        for collection in ax.collections
        if len(collection.get_sizes())  # the fit's CI band carries no marker size
    }
    assert point_sizes == {(SPIKEIN_POINT_SIZE,)}
    legend = ax.get_legend()
    assert legend is not None
    fig.canvas.draw()
    assert legend.get_window_extent().y0 >= ax.get_window_extent().y1
    plt.close(fig)


def test_build_spikein_frame_drops_reference_and_relabels_sites():
    """The zero-dilution sample never reaches the fit; site names become #N."""
    frame = build_spikein_frame(_spikein_table())

    assert "Spikein0" not in set(frame["Sample"])
    assert set(frame["Name"]) == {"#1", "#2"}
    assert frame["Name"].is_monotonic_increasing


# =============================================================================
# Figure S1
# =============================================================================
def _figure_s1_frames() -> tuple[pd.DataFrame, ...]:
    """One strictly positive pair per panel, in PANEL_SPECS order."""
    pairs = (("PBL", "PBR"), ("Reads_1", "Reads_2"), ("Reads_1", "Reads_2"))
    tables = [_replicate_table(x, y) for x, y in pairs]
    return tuple(
        figure_s1.build_scatter_frame(table, spec)
        for table, spec in zip(tables, figure_s1.PANEL_SPECS, strict=True)
    )


def test_draw_s1_figure_left_aligns_panel_labels_across_rows():
    """A-C share one row; D opens the second row with its letter under A's."""
    histogram = figure_s1.build_count_histogram(
        figure_s1.build_insertion_counts(_annotation_table())
    )
    mp = figure_s1.draw_figure(_figure_s1_frames(), histogram)
    boxes = _panel_boxes(mp, ("A", "B", "C", "D"))

    a_y = boxes["A"][1]
    assert boxes["B"][1] == pytest.approx(a_y)  # one common row
    assert boxes["C"][1] == pytest.approx(a_y)
    assert boxes["A"][0] < boxes["B"][0] < boxes["C"][0]
    assert boxes["D"][1] + boxes["D"][3] < a_y  # D starts the second row
    assert mp.fig.get_size_inches()[0] * 72 == pytest.approx(figure_s1.CANVAS_WIDTH_PT)
    assert boxes["D"][0] + boxes["D"][2] < figure_s1.CANVAS_WIDTH_PT  # right side empty
    assert _panel_label_left_pt(mp, "A") == pytest.approx(_panel_label_left_pt(mp, "D"), abs=0.5)


def test_build_insertion_counts_s1_counts_only_in_gene_rows():
    """Intergenic rows never reach the count; a gene absent from the table has no row."""
    counts = figure_s1.build_insertion_counts(_annotation_table())

    assert counts.to_dict() == {"G1": 2, "G2": 1}


def test_build_count_histogram_s1_orders_bins_and_lumps_the_tail():
    """Bins run left to right in count order; the last bar holds every count >= 300."""
    counts = pd.Series([0, 10, 25, 299, 300, 1254])
    frame = figure_s1.build_count_histogram(counts)

    expected_labels = [str(edge) for edge in range(0, 300, 25)] + ["≥300"]
    assert frame["bin_label"].tolist() == expected_labels
    assert frame["n_genes"].tolist()[0] == 2
    assert frame["n_genes"].tolist()[1] == 1
    assert frame["n_genes"].tolist()[-1] == 2
    assert frame["n_genes"].sum() == len(counts)


def test_build_scatter_frame_s1_drops_non_positive_pairs():
    """A zero read count cannot be drawn on a log axis and is removed."""
    spec = figure_s1.PANEL_SPECS[0]
    frame = figure_s1.build_scatter_frame(_replicate_table("PBL", "PBR"), spec)

    assert (frame[["PBL", "PBR"]] > 0).all().all()
    assert len(frame) == 40


def test_read_panel_tables_s1_rejects_missing_column(figure_s1_inputs):
    """A tech table without Reads_2 fails at the boundary, not inside the KDE."""
    table = pd.read_parquet(figure_s1_inputs.tech)[["Reads_1"]]
    table.to_parquet(figure_s1_inputs.tech)

    with pytest.raises(ValueError, match="missing required columns"):
        figure_s1.read_panel_tables(figure_s1_inputs)


# =============================================================================
# End-to-end assembly
# =============================================================================
def test_run_merges_schematic_into_a_single_vector_page(figure1_inputs: Figure1Config):
    """The output is one page of the designed size carrying both layers' text."""
    if shutil.which("pdftocairo") is None:
        pytest.skip("pdftocairo (poppler) is not installed")

    width_pt, height_pt = run(figure1_inputs)

    reader = PdfReader(str(figure1_inputs.output_pdf))
    assert len(reader.pages) == 1

    page = reader.pages[0]
    assert width_pt == pytest.approx(A4_WIDTH_PT, abs=0.01)
    assert height_pt == pytest.approx(A4_HEIGHT_PT, abs=0.01)  # the manuscript page
    assert float(page.mediabox.width) == pytest.approx(width_pt, abs=0.01)
    assert float(page.mediabox.height) == pytest.approx(height_pt, abs=0.01)

    text = page.extract_text() or ""
    assert "SCHEMATIC" in text  # the overlaid external page survived
    assert "Covered genes (%)" in text  # the matplotlib canvas survived
    assert "FACS Strain#1" in text  # panel F's mini axes drew inside the canvas
    assert figure1_inputs.canvas_pdf.exists()
    assert figure1_inputs.output_review.exists()


def test_figure_s1_run_assembles_one_page(figure_s1_inputs: figure_s1.FigureS1Config):
    """Figure S1 is one page of the designed size carrying all three panels' labels."""
    if shutil.which("pdftocairo") is None:
        pytest.skip("pdftocairo (poppler) is not installed")

    width_pt, height_pt = figure_s1.run(figure_s1_inputs)

    reader = PdfReader(str(figure_s1_inputs.output_pdf))
    assert len(reader.pages) == 1

    page = reader.pages[0]
    assert width_pt == pytest.approx(figure_s1.A4_WIDTH_PT, abs=0.01)
    assert height_pt == pytest.approx(figure_s1.A4_HEIGHT_PT, abs=0.01)  # the manuscript page
    assert float(page.mediabox.width) == pytest.approx(width_pt, abs=0.01)
    assert float(page.mediabox.height) == pytest.approx(height_pt, abs=0.01)

    text = page.extract_text() or ""
    assert "PBL Reads" in text
    assert "Reads of Biological Replicate 1" in text
    assert "Insertions per gene" in text
    assert figure_s1_inputs.canvas_pdf.exists()
    assert figure_s1_inputs.output_review.exists()
