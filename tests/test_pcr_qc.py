"""Tests for the PCR / library-prep quality control stage (loader + generic plot + script config)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from workflow.src.data_config import MERGED_READS_SUBDIR, merged_reads_path
from workflow.src.plotting.generic import create_scatter_correlation_plot
from workflow.scripts.pcr_qc.prepare_pcr_qc_data import PCRQCConfig


# =============================================================================
# merged_reads_path (pre-release intermediate loader)
# =============================================================================
def test_merged_reads_path_resolves_for_registered_results_dir():
    """A dataset with results_dir yields results/8_merged/{sample}_{tp}_{cond}.tsv."""
    p = merged_reads_path("LD_DIT_HAP", "LD1328-7", "0h", "YES")
    assert p.name == "LD1328-7_0h_YES.tsv"
    assert p.parent.name == MERGED_READS_SUBDIR
    # ...and the subdir sits under the dataset's results_dir, not release/.
    assert "release" not in str(p)
    assert str(p).endswith("LD_DIT_HAP/results/8_merged/LD1328-7_0h_YES.tsv")


def test_merged_reads_path_raises_without_results_dir():
    """Release-only datasets (no results_dir) are deliberately unreachable here."""
    with pytest.raises(KeyError):
        merged_reads_path("HD_diploid", "x", "0h", "YES")


def test_merged_reads_path_raises_for_unknown_dataset():
    """An unregistered dataset name raises KeyError, not a silently-wrong path."""
    with pytest.raises(KeyError):
        merged_reads_path("NoSuchDataset", "x", "0h", "YES")


# =============================================================================
# create_scatter_correlation_plot (generic, domain-agnostic)
# =============================================================================
def test_scatter_correlation_returns_axes_with_stats_box():
    """The plot renders points and annotates a PCC/R² stats text box."""
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(0)
    x = rng.random(50) + 1
    y = x * 2 + rng.random(50) * 0.1
    fig, ax = plt.subplots()
    returned = create_scatter_correlation_plot(x, y, ax=ax)
    assert returned is ax
    texts = [t.get_text() for t in ax.texts]
    assert any("PCC" in t and "R²" in t for t in texts)
    plt.close(fig)


def test_scatter_correlation_log_filters_nonpositive():
    """xscale='log' drops non-positive x so log10 never sees a zero/negative."""
    import matplotlib.pyplot as plt

    # Two of the four points are non-positive in x; they must be excluded.
    x = np.array([0.0, -1.0, 10.0, 100.0])
    y = np.array([1.0, 2.0, 10.0, 100.0])
    fig, ax = plt.subplots()
    create_scatter_correlation_plot(x, y, ax=ax, xscale="log", yscale="log")
    stats = " ".join(t.get_text() for t in ax.texts)
    assert "Data points: 2" in stats  # only the two positive-x points survive
    assert ax.get_xscale() == "log"
    plt.close(fig)


# =============================================================================
# PCRQCConfig.validate (script config)
# =============================================================================
def _write_merged(path: Path) -> None:
    """Write a minimal merged reads TSV indexed by (Chr, Coordinate, Strand)."""
    df = pd.DataFrame(
        {
            "Chr": ["I", "I"],
            "Coordinate": [100, 200],
            "Strand": ["+", "-"],
            "PBL": [10, 20],
            "PBR": [15, 25],
            "Reads": [25, 45],
        }
    ).set_index(["Chr", "Coordinate", "Strand"])
    df.to_csv(path, sep="\t")


def _make_config(tmp_path: Path, *, missing: bool = False) -> PCRQCConfig:
    inputs = {}
    for key in ["pbl_pbr", "tech_rep_1", "tech_rep_2", "bio_rep_1", "bio_rep_2", "spikein"]:
        p = tmp_path / f"{key}.tsv"
        if not missing:
            _write_merged(p)
        inputs[key] = p
    return PCRQCConfig(output=tmp_path / "out" / "PCR_quality_control.pdf", **inputs)


def test_config_validate_creates_output_dir(tmp_path):
    """validate() makes the output parent dir when all inputs exist."""
    config = _make_config(tmp_path)
    config.validate()
    assert config.output.parent.is_dir()


def test_config_validate_raises_on_missing_input(tmp_path):
    """A missing input file surfaces as ValueError, not a downstream read error."""
    config = _make_config(tmp_path, missing=True)
    with pytest.raises(ValueError):
        config.validate()
