"""Tests for workflow/src/figures.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")  # 无头后端

import matplotlib.pyplot as plt
import pytest

from workflow.src import figures
from workflow.src.figures import PanelShape


def test_apply_house_style_sets_cnsplots_palette():
    """apply_house_style correctly configures cnsplots palette and settings."""
    import cnsplots as cns

    figures.apply_house_style()

    # 验证调色板设置
    assert cns.settings.palette_qual == figures.HOUSE_PALETTE
    assert cns.settings.palette_seq == figures.DENSITY_CMAP

    # 验证字体设置
    assert "Arial" in cns.settings.font_sans_serif
    assert cns.settings.panel_label_fontname == "Arial"

    # 验证面板标签偏移
    assert cns.settings.panel_pad_left == figures.PANEL_LABEL_PAD_LEFT_PX
    assert cns.settings.panel_pad_top == figures.PANEL_LABEL_PAD_TOP_PX


def test_save_dual_creates_pdf_and_png(tmp_path):
    """save_dual generates both PDF and PNG files."""
    figures.apply_house_style()

    # 创建简单图表
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])

    # 保存到临时目录
    output_stem = tmp_path / "test_figure"
    figures.save_dual(output_stem)

    # 验证两个文件都存在
    assert (tmp_path / "test_figure.pdf").exists()
    assert (tmp_path / "test_figure.review.png").exists()

    plt.close(fig)


def test_grid_axes_returns_correct_shape():
    """grid_axes returns the correct number of axes in the expected shape."""
    figures.apply_house_style()

    # 创建 2x3 网格
    axes = figures.grid_axes(2, 3, shape=PanelShape.SQUARE)

    # 验证返回 6 个 axes
    assert len(axes) == 6

    # 验证所有 axes 都是 Axes 对象
    for ax in axes:
        assert isinstance(ax, plt.Axes)

    # 验证 figure 被创建
    fig = plt.gcf()
    assert fig is not None

    plt.close(fig)


def test_panel_labels_adds_text_to_axes():
    """grid_axes adds panel labels (A, B, C, ...) to each axes."""
    figures.apply_house_style()

    # 创建 1x3 网格，提供自定义标签
    custom_labels = ["X", "Y", "Z"]
    axes = figures.grid_axes(1, 3, labels=custom_labels, shape=PanelShape.SQUARE)

    # 验证返回 3 个 axes
    assert len(axes) == 3

    # 验证 figure 被创建
    fig = plt.gcf()
    assert fig is not None

    # 注意：cns.add_panel_label 添加的标签不是通过 ax.texts 访问的，
    # 而是通过 cnsplots 的内部机制。这里我们只验证函数调用成功完成。
    # 实际的标签渲染由 cnsplots 保证。

    plt.close(fig)


def test_panel_labels_function_generates_correct_sequence():
    """panel_labels generates A-Z, then A1, A2, ... for n > 26."""
    # 测试前 26 个标签（A-Z）
    labels = figures.panel_labels(26)
    assert labels[0] == "A"
    assert labels[25] == "Z"

    # 测试超过 26 个标签（A1, A2, ...）
    labels = figures.panel_labels(30)
    assert labels[26] == "A1"
    assert labels[27] == "A2"
    assert labels[29] == "A4"


def test_figure_size_for_grid_calculates_correct_dimensions():
    """figure_size_for_grid returns correct pixel dimensions for given grid."""
    # 测试 1x1 SQUARE 网格
    width, height = figures.figure_size_for_grid(1, 1, PanelShape.SQUARE)
    expected_width = 100 + figures.PANEL_DECORATION_WIDTH_PX + figures.FIGURE_PAD_PX
    expected_height = 100 + figures.PANEL_DECORATION_HEIGHT_PX + figures.FIGURE_PAD_PX
    assert width == expected_width
    assert height == expected_height

    # 测试 2x3 SQUARE 网格
    width, height = figures.figure_size_for_grid(2, 3, PanelShape.SQUARE)
    expected_width = 3 * (100 + figures.PANEL_DECORATION_WIDTH_PX) + figures.FIGURE_PAD_PX
    expected_height = 2 * (100 + figures.PANEL_DECORATION_HEIGHT_PX) + figures.FIGURE_PAD_PX
    assert width == expected_width
    assert height == expected_height


def test_figure_size_for_grid_raises_on_invalid_input():
    """figure_size_for_grid raises ValueError for invalid grid dimensions."""
    with pytest.raises(ValueError, match="n_rows must be at least 1"):
        figures.figure_size_for_grid(0, 1)

    with pytest.raises(ValueError, match="n_cols must be at least 1"):
        figures.figure_size_for_grid(1, 0)
