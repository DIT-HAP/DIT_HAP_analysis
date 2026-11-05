#!/usr/bin/env python3

"""
Cluster Enrichment Dotplot (interactive)
========================================

The Altair counterpart of ``go_enrichment_dotplot.pdf``: the same cluster x term
matrix on one self-contained page, with hover detail the static figure cannot
carry — the FDR, the gene ratio, the population size and the enriched genes
behind each dot.

Rows use the same home-cluster ordering as the static figure (see
``enrichment.plotting.dotplot_terms``), so a screenshot of either matches the
other. Clicking a dot dims every other term, which is how a reader follows one
process across clusters; double-click clears the selection.

Input
-----
- go_enrichment_full_filtered.tsv (pop_count-capped GO table from finalize_enrichment)

Output
------
- go_enrichment_dotplot.html: one page per variant. vega-embed loads from a CDN,
  so opening it needs a network connection the first time.

Usage
-----
    python plot_enrichment_interactive.py \\
        --enrichment .../2026-06-01/go_enrichment_full_filtered.tsv \\
        --output .../2026-06-01/go_enrichment_dotplot.html \\
        --variant-label kmeans_direct9 \\
        --top-terms-per-cluster 5

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-05
Version:  1.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
# Standard library
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any

# Third-party
import altair as alt
from loguru import logger
import pandas as pd

# Local project (resolved relative to this script, never the current directory)
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from enrichment.plotting import dotplot_matrix, dotplot_terms  # noqa: E402
from enrichment.schema import ENRICHMENT_COLUMNS, EnrichmentCol  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# 2. CONSTANTS
# =============================================================================
TOOLTIP_COLUMNS = [
    EnrichmentCol.CLUSTER,
    EnrichmentCol.TERM_ID,
    EnrichmentCol.STUDY_COUNT,
    EnrichmentCol.STUDY_ITEMS,
    EnrichmentCol.POP_COUNT,
]
GENE_TRUNCATE = 12
PAGE_CSS = """
<style>
  body { font-family: "Helvetica Neue", Arial, sans-serif; background: #ffffff; }
  #vis { margin: 0 auto; }
</style>
"""


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class InteractiveConfig:
    """Input table, output page, and the row-selection budget."""

    enrichment: Path
    output: Path
    variant_label: str
    top_terms_per_cluster: int

    def validate(self) -> None:
        """Raise if the input table is missing, then create the output directory."""
        if not self.enrichment.exists():
            raise ValueError(f"Required input not found: {self.enrichment}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 4. CHART
# =============================================================================
def prepare_rows(
    enrichment: pd.DataFrame, top_terms_per_cluster: int
) -> tuple[pd.DataFrame, list[str]]:
    """Build the plot rows with term labels, counts, and a truncated gene list."""
    terms = dotplot_terms(enrichment, top_terms_per_cluster)
    matrix = dotplot_matrix(enrichment, terms)

    labels = terms.set_index(EnrichmentCol.TERM_ID)[[EnrichmentCol.TERM, EnrichmentCol.NAMESPACE]]
    counts = enrichment[TOOLTIP_COLUMNS].drop_duplicates(
        [EnrichmentCol.CLUSTER, EnrichmentCol.TERM_ID]
    )

    rows = matrix.join(labels, on=EnrichmentCol.TERM_ID)
    rows = rows.merge(
        counts,
        on=[EnrichmentCol.CLUSTER, EnrichmentCol.TERM_ID],
        how="left",
        validate="many_to_one",
    )
    rows["genes"] = rows[EnrichmentCol.STUDY_ITEMS].map(_truncate_genes)
    return rows, list(terms[EnrichmentCol.TERM])


def _truncate_genes(study_items: str) -> str:
    """Shorten a comma-separated gene list to the first few names plus a count."""
    genes = [gene for gene in str(study_items).split(",") if gene]
    if len(genes) <= GENE_TRUNCATE:
        return ", ".join(genes)
    return ", ".join(genes[:GENE_TRUNCATE]) + f" … (+{len(genes) - GENE_TRUNCATE} more)"


def build_chart(rows: pd.DataFrame, term_order: list[str], title: str) -> Any:
    """Assemble the interactive dotplot: hover detail, click to isolate one term."""
    term_pick = alt.selection_point(fields=["term"], on="click", clear="dblclick", name="term")
    height = max(260, 15 * len(term_order))

    # Altair's mark methods return a generic self-type that ty cannot resolve,
    # so the composed chart is held as Any instead of suppressing per line.
    base: Any = alt.Chart(rows).mark_circle(stroke="#4d4d4d", strokeWidth=0.4)
    chart = base.encode(
        x=alt.X("Cluster:O", title="Cluster", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("term:N", sort=term_order, title="", axis=alt.Axis(labelLimit=420)),
        size=alt.Size(
            "gene_ratio:Q",
            title="Gene ratio",
            scale=alt.Scale(range=[20, 420]),
            legend=alt.Legend(format=".2f"),
        ),
        color=alt.Color(
            "neg_log10_fdr:Q",
            title="-log10(FDR)",
            scale=alt.Scale(scheme="viridis"),
        ),
        opacity=alt.condition(term_pick, alt.value(1.0), alt.value(0.35)),
        tooltip=[
            alt.Tooltip("term:N", title="Term"),
            alt.Tooltip("namespace:N", title="Namespace"),
            alt.Tooltip("Cluster:O", title="Cluster"),
            alt.Tooltip("p_fdr:Q", title="FDR", format=".3g"),
            alt.Tooltip("gene_ratio:Q", title="Gene ratio", format=".2f"),
            alt.Tooltip("study_count:Q", title="Genes", format="d"),
            alt.Tooltip("pop_count:Q", title="Population", format="d"),
            alt.Tooltip("genes:N", title="Enriched genes"),
        ],
    )
    chart = chart.add_params(term_pick)
    chart = chart.properties(
        width=520,
        height=height,
        title=alt.TitleParams(text=title, anchor="start", fontSize=16),
    )
    return chart.configure_view(stroke=None)


def save_page(chart: alt.Chart, output: Path) -> None:
    """Write the page and inject the small house CSS block."""
    chart.save(output)
    html = output.read_text(encoding="utf-8")
    output.write_text(html.replace("</head>", f"{PAGE_CSS}</head>", 1), encoding="utf-8")


# =============================================================================
# 5. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Write one variant's interactive enrichment page")
    parser.add_argument("--enrichment", type=Path, required=True, help="Filtered GO enrichment TSV")
    parser.add_argument("--output", type=Path, required=True, help="Output HTML path")
    parser.add_argument("--variant-label", required=True, help="Variant label shown in the title")
    parser.add_argument(
        "--top-terms-per-cluster",
        type=int,
        required=True,
        help="Terms kept per cluster before the union",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then write the page."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = InteractiveConfig(
        enrichment=args.enrichment,
        output=args.output,
        variant_label=args.variant_label,
        top_terms_per_cluster=args.top_terms_per_cluster,
    )

    try:
        config.validate()
        enrichment = pd.read_csv(config.enrichment, sep="\t")
        missing_columns = set(ENRICHMENT_COLUMNS) - set(enrichment.columns)
        if missing_columns:
            logger.error(
                f"Missing required columns in {config.enrichment}: {sorted(missing_columns)}"
            )
            return 1
        rows, term_order = prepare_rows(enrichment, config.top_terms_per_cluster)
        chart = build_chart(rows, term_order, config.variant_label)
        save_page(chart, config.output)
        logger.success(f"Wrote interactive page -> {config.output}")
    except Exception as error:
        logger.exception(f"Interactive page failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
