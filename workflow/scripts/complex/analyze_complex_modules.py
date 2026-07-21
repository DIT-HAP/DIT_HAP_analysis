#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Named Complex Module Visualization
==================================

Per-dataset: visualizes a handful of named functional modules (cytoplasmic
translation, kinetochore, mitochondria, vesicle trafficking, vacuolar ATPase)
in the 2D DIT-HAP fitness (DR, DL) space, highlighting each module's member
genes against the genome-wide gene cloud. Ported from
DIT_HAP_pipeline/workflow/notebooks/complex_analysis.ipynb (section 4).

Module membership
-----------------
The source notebook hardcoded curated gene-symbol dicts per module. Here each
module is driven by the config `complex.modules` map:
  * a NON-EMPTY list is taken as an explicit list of gene SYMBOLS (Name) for
    that module (the notebook's curated behaviour);
  * an EMPTY list means "auto-resolve from the PomBase complex annotation" by
    matching the module's keyword(s) (see _MODULE_KEYWORDS) against
    GO_term_name, then taking the union of member genes. This lets the empty
    config placeholders resolve at runtime rather than requiring hand-curation.

Input
-----
- Curated final_clusters.tsv (Systematic ID, A, DR, DL, revised_cluster); the
  Batch-B human-curated cluster table. Legacy um/lam headers -> DR/DL on load.
- PomBase macromolecular_complex_annotation.tsv (complex_term_id, GO_term_name,
  systematic_id, symbol, ...).

Output
------
- complex_module_visualization.pdf: one feature-space subplot per module.
- module_visualization_done.flag: sentinel written on success.

Usage
-----
    python analyze_complex_modules.py \\
        --final-clusters resources/curated/final_clusters.tsv \\
        --complex-annotation .../macromolecular_complex_annotation.tsv \\
        --modules "{'cytoplasmic_translation': [], ...}" \\
        --output-flag results/complex/{dataset}/_work/module_visualization_done.flag \\
        --output-figure results/complex/{dataset}/complex_module_visualization.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-20
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import ast
import sys
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")  # headless: this script only writes a PDF, never displays
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.plotting.gene_level import plot_given_genes_on_feature_space  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Legacy -> current metric column names (same quirk as
# workflow/src/clustering/candidates.py's _LEGACY_METRIC_RENAME).
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# PomBase annotation column -> canonical name.
_ANNOTATION_RENAME = {"systematic_id": "Systematic ID", "symbol": "Name"}

# Case-insensitive substrings used to auto-resolve a module from the complex
# annotation's GO_term_name when its config gene list is empty. A complex is
# assigned to the module if its GO_term_name contains ANY of the keywords AND
# none of the module's excludes (e.g. cytoplasmic-translation ribosomes must
# not be mitochondrial). Keys match the config `complex.modules` map.
_MODULE_KEYWORDS: dict[str, dict[str, list[str]]] = {
    "cytoplasmic_translation": {
        "include": ["cytosolic large ribosomal", "cytosolic small ribosomal",
                    "eukaryotic translation", "aminoacyl-tRNA synthetase",
                    "translation release factor"],
        "exclude": ["mitochondrial"],
    },
    "kinetochore": {"include": ["kinetochore"], "exclude": []},
    "mitochondria": {"include": ["mitochondrial"], "exclude": []},
    "vesicle": {
        "include": ["SNARE", "exocyst", "TRAPP", "retromer", "COPI", "COPII", "GARP"],
        "exclude": [],
    },
    "vacuolar_ATPase": {
        "include": ["vacuolar proton-transporting V-type ATPase"], "exclude": [],
    },
}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ModuleConfig:
    """Inputs, outputs, and the module -> gene-symbol map for module visualization."""
    final_clusters: Path
    complex_annotation: Path
    output_flag: Path
    output_figure: Path
    modules: dict = field(default_factory=dict)

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.final_clusters, self.complex_annotation]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_flag, self.output_figure]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def parse_modules_arg(raw: str) -> dict[str, list[str]]:
    """Parse the --modules argument (a Python/JSON-ish dict literal) into a dict.

    Snakemake renders a config dict via str(), so it arrives as a Python dict
    literal (single quotes) rather than strict JSON — ast.literal_eval handles
    both. An empty / missing value yields an empty dict (all modules
    auto-resolve).
    """
    raw = (raw or "").strip()
    if not raw:
        return {}
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError) as exc:
        raise ValueError(f"Could not parse --modules as a dict literal: {exc}") from exc
    if not isinstance(parsed, dict):
        raise ValueError(f"--modules must be a dict, got {type(parsed).__name__}")
    return {str(k): list(v) if v else [] for k, v in parsed.items()}


def load_final_clusters(final_clusters_path: Path) -> pd.DataFrame:
    """Load curated final_clusters.tsv, normalizing legacy um/lam -> DR/DL columns."""
    clusters = pd.read_csv(final_clusters_path, sep="\t")
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in clusters.columns and new not in clusters.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        clusters = clusters.rename(columns=rename)
    for required in ["Systematic ID", "DR", "DL"]:
        if required not in clusters.columns:
            raise ValueError(f"final_clusters.tsv missing required column '{required}' (have: {list(clusters.columns)})")
    return clusters.dropna(subset=["DR", "DL"]).copy()


def load_complex_annotation(annotation_path: Path) -> pd.DataFrame:
    """Load PomBase macromolecular_complex_annotation.tsv (complex -> member genes)."""
    annotation = pd.read_csv(annotation_path, sep="\t").rename(columns=_ANNOTATION_RENAME)
    for required in ["GO_term_name", "Systematic ID", "Name"]:
        if required not in annotation.columns:
            raise ValueError(
                f"complex annotation missing required column '{required}' (have: {list(annotation.columns)})"
            )
    return annotation[["GO_term_name", "Systematic ID", "Name"]].drop_duplicates()


def resolve_module_genes(
    module: str, config_symbols: list[str], annotation: pd.DataFrame
) -> list[str]:
    """Resolve a module to a list of Systematic IDs.

    A non-empty `config_symbols` list is treated as explicit gene SYMBOLS
    (Name) and mapped to Systematic IDs via the annotation. An empty list
    triggers keyword auto-resolution against GO_term_name using
    _MODULE_KEYWORDS.
    """
    if config_symbols:
        matched = annotation[annotation["Name"].isin(config_symbols)]
        ids = sorted(set(matched["Systematic ID"]))
        missing = set(config_symbols) - set(matched["Name"])
        if missing:
            logger.warning(f"[{module}] {len(missing)} config symbols not in annotation: {sorted(missing)[:10]}...")
        return ids

    spec = _MODULE_KEYWORDS.get(module)
    if not spec:
        logger.warning(f"[{module}] no keyword rule and empty config list; skipping.")
        return []
    names = annotation["GO_term_name"].fillna("")
    include_mask = pd.Series(False, index=annotation.index)
    for kw in spec["include"]:
        include_mask |= names.str.contains(kw, case=False, regex=False)
    for kw in spec["exclude"]:
        include_mask &= ~names.str.contains(kw, case=False, regex=False)
    ids = sorted(set(annotation[include_mask]["Systematic ID"]))
    logger.info(f"[{module}] auto-resolved {len(ids)} genes from annotation keywords")
    return ids


# =============================================================================
# PLOTTING
# =============================================================================
def plot_modules(clusters: pd.DataFrame, module_genes: dict[str, list[str]]) -> plt.Figure:
    """One feature-space subplot per module, highlighting its members among all genes."""
    modules = list(module_genes.keys())
    n = max(len(modules), 1)
    col_num = min(3, n)
    row_num = int(np.ceil(n / col_num))

    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH * col_num, AX_HEIGHT * row_num))
    axes = np.atleast_1d(axes).flatten()

    for idx, module in enumerate(modules):
        genes = module_genes[module]
        plot_given_genes_on_feature_space(
            ax=axes[idx],
            data_df=clusters,
            genes=genes,
            gene_column="Systematic ID",
            title=module.replace("_", " "),
            x_feature="DR",
            y_feature="DL",
            cmap="#9D343C",
            label=module.replace("_", " "),
            s=40,
        )

    for j in range(len(modules), len(axes)):
        fig.delaxes(axes[j])

    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: ModuleConfig) -> None:
    """Load -> resolve each module's genes -> feature-space figure + sentinel flag."""
    config.validate()

    clusters = load_final_clusters(config.final_clusters)
    annotation = load_complex_annotation(config.complex_annotation)

    # Default to the five known modules when config supplies none.
    module_specs = config.modules or {m: [] for m in _MODULE_KEYWORDS}
    module_genes = {
        module: resolve_module_genes(module, symbols, annotation)
        for module, symbols in module_specs.items()
    }
    total = sum(len(g) for g in module_genes.values())
    logger.info(f"Resolved {len(module_genes)} modules covering {total} gene-assignments")

    fig = plot_modules(clusters, module_genes)
    with PdfPages(config.output_figure) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    config.output_flag.write_text(
        "module_visualization complete\n"
        + "\n".join(f"{m}\t{len(g)}" for m, g in module_genes.items())
        + "\n"
    )
    logger.success(f"Wrote {config.output_figure} and sentinel {config.output_flag}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Visualize named complex modules in fitness space")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--complex-annotation", type=Path, required=True, help="PomBase macromolecular_complex_annotation.tsv")
    parser.add_argument("--modules", type=str, default="", help="Module -> gene-symbol list dict literal (empty lists auto-resolve)")
    parser.add_argument("--output-flag", type=Path, required=True, help="Sentinel flag written on success")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output module visualization PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = ModuleConfig(
            final_clusters=args.final_clusters,
            complex_annotation=args.complex_annotation,
            output_flag=args.output_flag,
            output_figure=args.output_figure,
            modules=parse_modules_arg(args.modules),
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
