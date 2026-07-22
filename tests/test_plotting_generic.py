"""Tests for workflow/src/plotting/generic.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from workflow.src.plotting.generic import boxplot_with_violinplot


def test_boxplot_with_violinplot_ytick_count_and_labels():
    """One y-tick per label, each annotated with its sample size n=."""
    labels = ["WT", "small colonies", "spores"]
    values = [[0.1, 0.2, 0.3, 0.4], [0.5, 0.6], [0.9]]
    colors = ["#111111", "#222222", "#333333"]

    fig, ax = plt.subplots()
    boxplot_with_violinplot(labels, values, ax, colors)

    assert len(ax.get_yticks()) == 3
    tick_texts = [t.get_text() for t in ax.get_yticklabels()]
    assert tick_texts == ["WT (n=4)", "small colonies (n=2)", "spores (n=1)"]
    plt.close(fig)


def test_boxplot_with_violinplot_returns_axes():
    """Returns the same Axes it drew on (for composability)."""
    fig, ax = plt.subplots()
    out = boxplot_with_violinplot(["a"], [[1.0, 2.0, 3.0]], ax, ["#444444"])
    assert out is ax
    plt.close(fig)
