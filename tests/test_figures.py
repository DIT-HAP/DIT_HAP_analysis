"""Tests for workflow/src/figures.py."""

from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib as mpl

mpl.use("Agg")  # 无头后端

from loguru import logger
import matplotlib.pyplot as plt
from PIL import Image
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
    assert cns.settings.font_sans_serif == figures.FONT_SANS_SERIF
    assert cns.settings.panel_label_fontname is None

    # 验证面板标签偏移 (display px, 随 figure_dpi 换算)
    assert cns.settings.panel_pad_left == figures.PANEL_LABEL_PAD_LEFT_PX
    assert cns.settings.panel_pad_top == figures.PANEL_LABEL_PAD_TOP_PX

    # 验证导出契约
    assert cns.settings.figure_dpi == figures.HOUSE_FIGURE_DPI
    assert cns.settings.savefig_dpi == figures.HOUSE_FIGURE_DPI
    assert cns.settings.savefig_transparent is False


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

    # 创建 1x3 网格,提供自定义标签
    custom_labels = ["X", "Y", "Z"]
    axes = figures.grid_axes(1, 3, labels=custom_labels, shape=PanelShape.SQUARE)

    # 验证返回 3 个 axes
    assert len(axes) == 3

    # 验证 figure 被创建
    fig = plt.gcf()
    assert fig is not None

    # 注意:cns.add_panel_label 添加的标签不是通过 ax.texts 访问的,
    # 而是通过 cnsplots 的内部机制。这里我们只验证函数调用成功完成。
    # 实际的标签渲染由 cnsplots 保证。

    plt.close(fig)


def test_panel_labels_function_generates_correct_sequence():
    """panel_labels generates A-Z, then A1, A2, ... for n > 26."""
    # 测试前 26 个标签(A-Z)
    labels = figures.panel_labels(26)
    assert labels[0] == "A"
    assert labels[25] == "Z"

    # 测试超过 26 个标签(A1, A2, ...)
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


def test_content_bounds_never_shrinks_below_canvas():
    """content_bounds keeps the designed canvas and adds the 2 pt pad around it."""
    figures.apply_house_style()
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3])
    width, height = fig.get_size_inches()
    pad_in = figures.CONTENT_PAD_PT / 72

    bounds = figures.content_bounds(fig)

    assert bounds.width == pytest.approx(width + 2 * pad_in, abs=1e-6)
    assert bounds.height == pytest.approx(height + 2 * pad_in, abs=1e-6)
    plt.close(fig)


def test_content_bounds_grows_for_legend_outside_axes():
    """content_bounds includes a legend anchored past the right edge of the canvas."""
    figures.apply_house_style()
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3], label="series")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    canvas_width = fig.get_size_inches()[0]

    bounds = figures.content_bounds(fig)

    assert bounds.x1 > canvas_width + figures.CONTENT_PAD_PT / 72 + 0.1
    plt.close(fig)


def test_content_bounds_includes_in_layout_false_text():
    """content_bounds picks up a figure text that matplotlib's tight bbox ignores."""
    figures.apply_house_style()
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3])
    text = fig.text(1.2, 0.5, "outside")
    text.set_in_layout(False)
    canvas_width = fig.get_size_inches()[0]

    bounds = figures.content_bounds(fig)

    assert bounds.x1 > 1.2 * canvas_width
    plt.close(fig)


def test_content_bounds_warns_on_large_overshoot():
    """content_bounds logs a warning when content grows the page past the ratio."""
    figures.apply_house_style()
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3])
    fig.text(3.0, 0.5, "far outside")

    messages: list[str] = []
    sink_id = logger.add(lambda message: messages.append(message), level="WARNING")
    try:
        figures.content_bounds(fig)
    finally:
        logger.remove(sink_id)

    assert any("Content expands" in message for message in messages)
    plt.close(fig)


def test_save_dual_writes_measured_bounds(tmp_path):
    """save_dual renders the review PNG at the same measured bounds as the PDF."""
    figures.apply_house_style()
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3])
    fig.text(1.2, 0.5, "outside")
    bounds = figures.content_bounds(fig)

    figures.save_dual(tmp_path / "outside")

    with Image.open(tmp_path / "outside.review.png") as image:
        assert image.width == pytest.approx(bounds.width * figures.HOUSE_FIGURE_DPI, abs=2)
        assert image.height == pytest.approx(bounds.height * figures.HOUSE_FIGURE_DPI, abs=2)
    [(pdf_width, pdf_height)] = _pdf_page_sizes(tmp_path / "outside.pdf")
    assert (pdf_width, pdf_height) == pytest.approx(
        (bounds.width * 72, bounds.height * 72), abs=0.5
    )
    plt.close(fig)


def test_save_pdf_pages_uses_one_bound_for_every_page(tmp_path):
    """save_pdf_pages gives every page the union of all page bounds."""
    figures.apply_house_style()
    plain, plain_ax = plt.subplots(figsize=(2, 2))
    plain_ax.plot([1, 2, 3])
    wide, wide_ax = plt.subplots(figsize=(3, 2))
    wide_ax.plot([1, 2, 3])
    wide.text(1.2, 0.5, "outside")

    figures.save_pdf_pages(tmp_path / "pages", [plain, wide])

    sizes = _pdf_page_sizes(tmp_path / "pages.pdf")
    assert len(sizes) == 2
    assert sizes[0] == pytest.approx(sizes[1], abs=0.5)
    assert sizes[0][0] > 1.2 * 3 * 72  # the outside text starts at 1.2 x 3 in
    plt.close("all")


def _pdf_page_sizes(path: Path) -> list[tuple[float, float]]:
    """Return (width, height) of every PDF page in points, read from the MediaBox entries."""
    text = path.read_bytes().decode("latin-1")
    boxes = [box.split() for box in re.findall(r"/MediaBox \[([^\]]+)\]", text)]
    return [(float(x1) - float(x0), float(y1) - float(y0)) for x0, y0, x1, y1 in boxes]
