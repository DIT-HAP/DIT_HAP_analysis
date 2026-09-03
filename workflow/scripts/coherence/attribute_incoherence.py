#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Incoherence Attribution (WHY a complex is internally dispersed)
=========================================================================

Complementary to the coherence z-score: the tightest groups are the expected
obligate machines, but the INCOHERENT ones (z > 0 = more dispersed than random)
are biologically informative too — a complex's members scatter in DR-DL fitness
space when they play different functional roles. This stage scores the candidate
causes of that dispersion for every group and labels each with its most likely
explanation, so the meaningful incoherent complexes can be found and interpreted.

Diagnostic signals (per group, from workflow/src/coherence/attribution.py)
--------------------------------------------------------------------------
- major/minor split: a 2-component GMM on the members' normalized (DR, DL/10)
  points, called a genuine split when its silhouette is high enough — a tight
  essential "core" + a looser dispensable "minority" (e.g. eIF3 core vs eIF3e).
- shared-subunit fraction: fraction of members that also belong to OTHER groups
  of the same source (cross-complex members drag the centroid apart).
- paralog fraction: fraction of members with a paralog (deletion phenotype may be
  buffered by redundancy, dampening DR and pulling the group toward WT).
The label ladder (attribute_incoherence): conditional_module (split + shared) >
major_minor_split > shared_subunits > paralog_buffered > data_limited >
intrinsic_heterogeneity (real spread, no detected cause — for manual review;
also where annotation/technical artefacts fall, which are NOT auto-labelled).

Input
-----
- --metrics: a source's coherence_metrics.tsv (source, group_id, group_name,
  term_size, covered_genes, z_score, p_fdr, ...).
- --annotation: that source's group_annotation_long.tsv (group -> member genes).
- --fitting-results: upstream fitting_results.tsv (systematic id index; DR/DL).
- --paralogs: Ensembl paralog export TSV (its "Gene stable ID" column lists genes
  with >=1 paralog).

Output
------
- --output-table: incoherence_attribution.tsv — one row per scored group with
  z_score, p_fdr, gmm_silhouette, core_size, minor_size, shared_fraction,
  paralog_fraction, attribution_label, is_incoherent, and the shared/other-group
  member detail. Sorted by z_score descending (most incoherent first).
- --output-figure: incoherence_attribution.pdf — DR-DL scatter grid for the top-N
  incoherent groups (members coloured by GMM core/minor component), each titled
  with name / z / label; plus an attribution-label frequency panel.

Usage
-----
    python attribute_incoherence.py \\
        --metrics results/coherence/{dataset}/go_macrocomplex/coherence_metrics.tsv \\
        --annotation results/coherence/{dataset}/go_macrocomplex/group_annotation_long.tsv \\
        --fitting-results .../fitting_results.tsv \\
        --paralogs resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --z-threshold 0.0 --top-n-plot 16 \\
        --output-table results/coherence/{dataset}/go_macrocomplex/incoherence_attribution.tsv \\
        --output-figure results/coherence/{dataset}/go_macrocomplex/incoherence_attribution.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from loguru import logger  # noqa: E402

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from coherence.attribution import (  # noqa: E402
    major_minor_split,
    shared_subunits,
    shared_fraction,
    paralog_fraction,
    attribute_incoherence,
)
from plotting.style import AX_HEIGHT, AX_WIDTH  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Byte-faithful to compute_coherence.py: legacy um/lam -> DR/DL, DL normalized /10.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}
_DR_NORM, _DL_NORM = 1.0, 10.0


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class AttributionConfig:
    """Inputs, outputs, and parameters for incoherence attribution."""
    metrics: Path
    annotation: Path
    fitting_results: Path
    paralogs: Path
    output_table: Path
    output_figure: Path
    z_threshold: float = 0.0
    top_n_plot: int = 16
    shared_frac_threshold: float = 0.5
    paralog_frac_threshold: float = 0.5

    def validate(self) -> None:
        """Raise ValueError on missing inputs, then make output dirs."""
        for path in [self.metrics, self.annotation, self.fitting_results, self.paralogs]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_table, self.output_figure]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
# =============================================================================
# HELPERS — data loading
# =============================================================================
def load_member_points(fitting_results: Path) -> dict[str, tuple[float, float]]:
    """Load fitting_results.tsv -> {Systematic ID: (norm_DR, norm_DL)} for DR/DL-finite genes.

    Mirrors compute_coherence.py::load_fitting_results: legacy um/lam -> DR/DL,
    +/-inf -> NaN dropped, DR against [0,1] (unchanged), DL against [0,10] (/10).
    Returns ALL fitted genes (the DR>threshold filter is applied later against the
    same member set coherence used).
    """
    fitting = pd.read_csv(fitting_results, sep="\t", index_col=0).reset_index()
    if "Systematic ID" not in fitting.columns:
        fitting = fitting.rename(columns={fitting.columns[0]: "Systematic ID"})
    rename = {o: n for o, n in _LEGACY_METRIC_RENAME.items()
              if o in fitting.columns and n not in fitting.columns}
    fitting = fitting.rename(columns=rename)
    for req in ["Systematic ID", "DR", "DL"]:
        if req not in fitting.columns:
            raise ValueError(f"fitting_results.tsv missing '{req}' (have: {list(fitting.columns)})")
    fitting = fitting.replace([np.inf, -np.inf], np.nan).dropna(subset=["DR", "DL"])
    return {
        row["Systematic ID"]: (row["DR"] / _DR_NORM, row["DL"] / _DL_NORM)
        for _, row in fitting.iterrows()
    }


def load_paralog_ids(paralogs: Path) -> set[str]:
    """Genes that actually HAVE a paralog, from the Ensembl paralog export.

    The export has one row per (gene, paralogue) pair, but genes WITHOUT a
    paralogue still appear with the paralogue columns left blank — so filtering to
    a non-empty `...paralogue gene stable ID` is essential (otherwise every gene
    counts as having a paralog and paralog_fraction is a useless 1.0 everywhere).
    """
    par = pd.read_csv(paralogs, sep="\t")
    gene_col = "Gene stable ID"
    para_col = "Schizosaccharomyces pombe paralogue gene stable ID"
    if gene_col not in par.columns or para_col not in par.columns:
        raise ValueError(f"paralog TSV missing '{gene_col}'/'{para_col}' (have: {list(par.columns)[:5]}...)")
    return set(par.loc[par[para_col].notna(), gene_col].dropna().astype(str))


# =============================================================================
# CORE LOGIC — per-group attribution
# =============================================================================
def group_member_points(
    long_table: pd.DataFrame, group_id: str, points: dict[str, tuple[float, float]]
) -> tuple[list[str], np.ndarray]:
    """The group's members that have fitness points, as (ids, (n,2) array).

    Uses the same long-table -> point-cloud join coherence used, so the GMM sees
    exactly the member set the z-score was computed on (members without a fitted
    DR/DL are dropped, matching compute_coherence's inner merge).
    """
    members = long_table.loc[long_table["group_id"] == group_id, "Systematic ID"].unique()
    ids = [m for m in members if m in points]
    X = np.array([points[m] for m in ids]) if ids else np.empty((0, 2))
    return ids, X


def attribute_all(
    metrics: pd.DataFrame,
    long_table: pd.DataFrame,
    points: dict[str, tuple[float, float]],
    paralog_ids: set[str],
    config: AttributionConfig,
) -> tuple[pd.DataFrame, dict[str, dict]]:
    """One attribution row per scored group + a {group_id: gmm_split_dict} map for plotting.

    For each group: reconstruct its member points, run the GMM major/minor split,
    compute the shared-subunit and paralog fractions, and assign the label. The
    per-group split dict (with GMM labels/core_label) is returned separately so the
    figure can colour members by component without recomputing.
    """
    rows = []
    splits: dict[str, dict] = {}
    for _, m in metrics.iterrows():
        gid = m["group_id"]
        ids, X = group_member_points(long_table, gid, points)
        split = major_minor_split(X)
        splits[gid] = {"split": split, "ids": ids, "X": X}

        shared_frac = shared_fraction(long_table, gid)
        par_frac = paralog_fraction(ids, paralog_ids)
        label = attribute_incoherence(
            split, shared_frac, par_frac,
            shared_frac_threshold=config.shared_frac_threshold,
            paralog_frac_threshold=config.paralog_frac_threshold,
        )
        sizes = split.get("component_sizes")
        core_label = split.get("core_label")
        shared_members = shared_subunits(long_table, gid)
        rows.append({
            "source": m.get("source"),
            "group_id": gid,
            "group_name": m["group_name"],
            "term_size": m["term_size"],
            "z_score": m["z_score"],
            "p_fdr": m.get("p_fdr", np.nan),
            "is_incoherent": bool(m["z_score"] > config.z_threshold),
            "gmm_silhouette": split.get("silhouette", np.nan),
            "gmm_is_split": bool(split.get("is_split")),
            "core_size": sizes[core_label] if sizes and core_label is not None else np.nan,
            "minor_size": sizes[1 - core_label] if sizes and core_label is not None else np.nan,
            "shared_fraction": shared_frac,
            "paralog_fraction": par_frac,
            "attribution_label": label,
            "n_shared_members": len(shared_members),
            "shared_members": "; ".join(shared_members["Systematic ID"]) if not shared_members.empty else "",
        })
    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.sort_values("z_score", ascending=False).reset_index(drop=True)
    return table, splits


# =============================================================================
# PLOTTING
# =============================================================================
def plot_attribution(table: pd.DataFrame, splits: dict[str, dict], top_n: int) -> plt.Figure:
    """DR-DL scatter grid for the top-N incoherent groups + a label-frequency panel.

    Each scatter shows a group's members in (norm_DR, norm_DL) space coloured by
    GMM component (core = tight blue, minor = loose red) when a 2-component split
    was fit, else a single colour; titled name / z / label. The final panel is a
    horizontal bar of attribution_label counts across the incoherent groups.
    """
    incoherent = table[table["is_incoherent"]].head(top_n)
    n_scatter = len(incoherent)
    if n_scatter == 0:
        fig, ax = plt.subplots(1, 1, figsize=(AX_WIDTH, AX_HEIGHT))
        ax.text(0.5, 0.5, "No incoherent groups (z > threshold)", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        return fig

    n_panels = n_scatter + 1  # + label-frequency panel
    ncols = 4
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(AX_WIDTH * ncols, AX_HEIGHT * nrows), squeeze=False)
    flat = axes.ravel()

    for ax, (_, row) in zip(flat, incoherent.iterrows()):
        info = splits.get(row["group_id"], {})
        X = info.get("X", np.empty((0, 2)))
        split = info.get("split", {})
        labels = split.get("labels")
        core = split.get("core_label")
        if labels is not None and core is not None:
            ax.scatter(X[labels == core, 0], X[labels == core, 1], c="#2c6fbb", s=28, label="core", alpha=0.85)
            ax.scatter(X[labels != core, 0], X[labels != core, 1], c="#c0392b", s=28, label="minor", alpha=0.85)
            ax.legend(fontsize=6, loc="best", frameon=False)
        elif X.shape[0]:
            ax.scatter(X[:, 0], X[:, 1], c="#555555", s=28, alpha=0.85)
        ax.set_xlabel("norm DR")
        ax.set_ylabel("norm DL/10")
        sil = row["gmm_silhouette"]
        sil_str = f", sil={sil:.2f}" if pd.notna(sil) else ""
        ax.set_title(f"{str(row['group_name'])[:30]}\nz={row['z_score']:.2f}{sil_str} [{row['attribution_label']}]", fontsize=7)

    # Label-frequency panel over ALL incoherent groups (not just the plotted top-N).
    ax_freq = flat[n_scatter]
    counts = table[table["is_incoherent"]]["attribution_label"].value_counts()
    ax_freq.barh(counts.index[::-1], counts.values[::-1], color="#6b99df")
    ax_freq.set_xlabel("Number of incoherent groups")
    ax_freq.set_title("Attribution label frequency", fontsize=8)

    for ax in flat[n_panels:]:
        fig.delaxes(ax)
    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: AttributionConfig) -> None:
    """Load -> per-group attribution -> TSV + figure."""
    config.validate()
    metrics = pd.read_csv(config.metrics, sep="\t")
    for req in ["group_id", "group_name", "term_size", "z_score"]:
        if req not in metrics.columns:
            raise ValueError(f"metrics TSV missing '{req}' (have: {list(metrics.columns)})")
    long_table = pd.read_csv(config.annotation, sep="\t")
    points = load_member_points(config.fitting_results)
    paralog_ids = load_paralog_ids(config.paralogs)

    if metrics.empty:
        logger.warning("metrics table is empty; writing empty attribution outputs")
        pd.DataFrame().to_csv(config.output_table, sep="\t", index=False)
        with PdfPages(config.output_figure) as pdf:
            pdf.savefig(plot_attribution(pd.DataFrame(columns=["is_incoherent"]), {}, config.top_n_plot))
        return

    table, splits = attribute_all(metrics, long_table, points, paralog_ids, config)
    table.to_csv(config.output_table, sep="\t", index=False)
    fig = plot_attribution(table, splits, config.top_n_plot)
    with PdfPages(config.output_figure) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)

    n_inc = int(table["is_incoherent"].sum())
    label_counts = table[table["is_incoherent"]]["attribution_label"].value_counts().to_dict()
    logger.success(
        f"{len(table):,} groups scored, {n_inc:,} incoherent (z>{config.z_threshold}); "
        f"labels={label_counts}; wrote {config.output_table}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Attribute the cause of coherence incoherence per group")
    parser.add_argument("--metrics", type=Path, required=True, help="A source's coherence_metrics.tsv")
    parser.add_argument("--annotation", type=Path, required=True, help="That source's group_annotation_long.tsv")
    parser.add_argument("--fitting-results", type=Path, required=True, help="Upstream fitting_results.tsv (DR/DL)")
    parser.add_argument("--paralogs", type=Path, required=True, help="Ensembl paralog export TSV")
    parser.add_argument("--z-threshold", type=float, default=0.0, help="z_score above this = incoherent")
    parser.add_argument("--top-n-plot", type=int, default=16, help="How many top-incoherent groups to scatter")
    parser.add_argument("--shared-frac-threshold", type=float, default=0.5, help="shared_fraction >= this triggers the shared-subunit label")
    parser.add_argument("--paralog-frac-threshold", type=float, default=0.5, help="paralog_fraction >= this triggers the paralog-buffered label")
    parser.add_argument("--output-table", type=Path, required=True, help="Output attribution TSV")
    parser.add_argument("--output-figure", type=Path, required=True, help="Output attribution PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run attribution, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = AttributionConfig(
            metrics=args.metrics,
            annotation=args.annotation,
            fitting_results=args.fitting_results,
            paralogs=args.paralogs,
            output_table=args.output_table,
            output_figure=args.output_figure,
            z_threshold=args.z_threshold,
            top_n_plot=args.top_n_plot,
            shared_frac_threshold=args.shared_frac_threshold,
            paralog_frac_threshold=args.paralog_frac_threshold,
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
