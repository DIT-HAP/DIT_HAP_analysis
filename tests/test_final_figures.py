"""Tests for the final-figure stage: Figure 1's frames, layout and PDF assembly."""

from pathlib import Path
import shutil

import pandas as pd
import pytest

pytest.importorskip("cnsplots")
pytest.importorskip("pypdf")

import cnsplots as cns
from matplotlib.transforms import Bbox
from pypdf import PdfReader
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
