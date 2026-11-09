"""Tests for the final-figure stage: layout, frames, and PDF assembly per figure."""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("cnsplots")
pytest.importorskip("pypdf")

import cnsplots as cns
from matplotlib.transforms import Bbox
from pypdf import PdfReader
from workflow.scripts.figures import figure_s1
from workflow.scripts.figures.figure1 import (
    CANVAS_WIDTH_PT,
    COVERAGE_CLASSES,
    Figure1Config,
    build_coverage_frame,
    build_layout,
    build_spikein_frame,
    run,
)


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


@pytest.fixture
def figure_s1_inputs(tmp_path: Path) -> figure_s1.FigureS1Config:
    """Write fixture parquets and return a Figure S1 config pointing at them."""
    pbl_pbr = tmp_path / "pbl_pbr.parquet"
    tech = tmp_path / "tech.parquet"
    bio = tmp_path / "bio.parquet"
    _replicate_table("PBL", "PBR").to_parquet(pbl_pbr)
    _replicate_table("Reads_1", "Reads_2").to_parquet(tech)
    _replicate_table("Reads_1", "Reads_2").to_parquet(bio)
    return figure_s1.FigureS1Config(
        pbl_pbr=pbl_pbr,
        tech=tech,
        bio=bio,
        output_pdf=tmp_path / "out" / "figure_s1.pdf",
        output_review=tmp_path / "out" / "figure_s1.review.png",
    )


@pytest.fixture
def figure1_inputs(tmp_path: Path) -> Figure1Config:
    """Write fixture inputs and return a config pointing at them."""
    coverage = tmp_path / "lfc_coverage_comparison.tsv"
    spikein = tmp_path / "spike_in_stats.tsv"
    schematic = tmp_path / "figure1a.pdf"
    _coverage_table().to_csv(coverage, sep="\t", index=False)
    _spikein_table().to_csv(spikein, sep="\t", index=False)
    _write_schematic(schematic)
    return Figure1Config(
        lfc_comparison=coverage,
        spikein=spikein,
        schematic=schematic,
        output_pdf=tmp_path / "out" / "figure1.pdf",
        output_review=tmp_path / "out" / "figure1.review.png",
        canvas_pdf=tmp_path / "out" / "_work" / "figure1_canvas.pdf",
    )


# =============================================================================
# Layout
# =============================================================================
def test_build_layout_stacks_schematic_above_the_panel_row():
    """A spans the canvas width above a B/C row whose axes share one top edge."""
    layout = build_layout(453.72, 113.4)

    assert layout.schematic_height_pt == pytest.approx(493.0 * 113.4 / 453.72)
    assert layout.schematic_rect[0] == 0.0
    assert layout.schematic_rect[2] == pytest.approx(1.0)  # full canvas width
    assert layout.schematic_rect[1] > layout.panel_b_rect[1]  # A sits above B
    assert layout.panel_b_rect[0] < layout.panel_c_rect[0]
    assert layout.panel_b_rect[1] == layout.panel_c_rect[1]
    assert layout.panel_b_rect[3] == layout.panel_c_rect[3]
    assert layout.canvas_height_pt > layout.schematic_height_pt + 80.0


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


def test_build_spikein_frame_drops_reference_and_relabels_sites():
    """The zero-dilution sample never reaches the fit; site names become #N."""
    frame = build_spikein_frame(_spikein_table())

    assert "Spikein0" not in set(frame["Sample"])
    assert set(frame["Name"]) == {"#1", "#2"}
    assert frame["Name"].is_monotonic_increasing


# =============================================================================
# Figure S1
# =============================================================================
def test_build_layout_s1_places_three_panels_across_the_page():
    """Three equally sized panels share one top edge and stay inside the canvas."""
    layout = figure_s1.build_layout()

    assert len(layout.panel_rects) == len(figure_s1.PANEL_SPECS)
    assert layout.panel_xs_pt == tuple(sorted(layout.panel_xs_pt))
    assert layout.panel_rects[0][1] == layout.panel_rects[2][1]  # one common top edge
    assert layout.panel_rects[0][2] == layout.panel_rects[2][2]
    assert layout.canvas_height_pt == pytest.approx(
        figure_s1.PANEL_LETTER_ROOM_PT
        + figure_s1.PANEL_SIZE_PT
        + figure_s1.BELOW_DECORATION_PT
        + figure_s1.BOTTOM_MARGIN_PT
    )
    right_edge_pt = max(layout.panel_xs_pt) + figure_s1.PANEL_SIZE_PT
    assert right_edge_pt < figure_s1.CANVAS_WIDTH_PT


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

    run(figure1_inputs)

    reader = PdfReader(str(figure1_inputs.output_pdf))
    assert len(reader.pages) == 1

    layout = build_layout(200.0, 50.0)
    page = reader.pages[0]
    assert float(page.mediabox.width) == pytest.approx(CANVAS_WIDTH_PT, abs=0.01)
    assert float(page.mediabox.height) == pytest.approx(layout.canvas_height_pt, abs=0.01)

    text = page.extract_text() or ""
    assert "SCHEMATIC" in text  # the overlaid external page survived
    assert "Covered genes (%)" in text  # the matplotlib canvas survived
    assert figure1_inputs.canvas_pdf.exists()
    assert figure1_inputs.output_review.exists()


def test_figure_s1_run_assembles_one_page(figure_s1_inputs: figure_s1.FigureS1Config):
    """Figure S1 is one page of the designed size carrying all three panels' labels."""
    if shutil.which("pdftocairo") is None:
        pytest.skip("pdftocairo (poppler) is not installed")

    figure_s1.run(figure_s1_inputs)

    reader = PdfReader(str(figure_s1_inputs.output_pdf))
    assert len(reader.pages) == 1

    layout = figure_s1.build_layout()
    page = reader.pages[0]
    assert float(page.mediabox.width) == pytest.approx(figure_s1.CANVAS_WIDTH_PT, abs=0.01)
    assert float(page.mediabox.height) == pytest.approx(layout.canvas_height_pt, abs=0.01)

    text = page.extract_text() or ""
    assert "PBL Reads" in text
    assert "Reads of Biological Replicate 1" in text
    assert figure_s1_inputs.output_review.exists()
