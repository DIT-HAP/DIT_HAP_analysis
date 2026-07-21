#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Non-coding RNA Depletion Analysis
===================================

Per-dataset: merges non-coding-gene DIT-HAP fitting results with the ncRNA
genome-region bed, GtRNAdb tRNA annotations (matched by chr+start+end, NOT by
name), and Marguerat 2012 mRNA abundance, then characterizes nuclear tRNA
depletion (tRNA copy number, amino-acid/anticodon breakdown, DR distribution).
Ported from
DIT_HAP_pipeline/workflow/notebooks/non_coding_RNA_analysis.ipynb.

This is a simplified port: the notebook's telomere/centromere Location
categorization and the per-tRNA flanking-sequence nucleotide heatmap (which
needs the genome FASTA + BioPython) are out of scope here — only the
tRNA-annotation merge, copy-number, and depletion summary the task calls for.

Input
-----
- Non-coding-gene fitting_results.tsv (Systematic ID + per-gene depletion
  stats). Legacy releases ship the pre-rename um/lam headers instead of DR/DL;
  normalized on load (same quirk as coverage.smk / verification's load_gene_level).
- ncRNA genome-region bed (`#Chr Start End ... Feature Systematic ID Type Name
  ...`) — provides the genomic coordinates + Feature type per ncRNA gene.
- GtRNAdb tRNA bed (schiPomb_972H-tRNAs.bed, headerless 12-col BED) — carries
  the GtRNAdb_Name (e.g. "tRNA-Ala-AGC-1-1") the anticodon is parsed from. Its
  chromosome names are "chrI/chrII/chrIII" and are normalized to "I/II/III"
  before the positional merge (source-notebook quirk).
- Marguerat 2012 xlsx (Table_S2, comment="#"): MM1/MM2/MN1/MN2.tot.cpc_ex
  columns indexed on Systematic.name; mean per condition
  (EMM_Proliferating = mean(MM1,MM2), EMM_Nitrogen_Starved = mean(MN1,MN2)).

Output
------
- ncrna_stats.tsv: per-nuclear-tRNA table (Systematic ID, GtRNAdb_Name,
  Amino_Acid, Anticodon, tRNA_copy_number, DR/DL, mRNA abundance means).
- ncrna_analysis.pdf: ncRNA Feature-type donut + tRNA copy-number distribution
  + DR-by-copy-number scatter.

Usage
-----
    python analyze_noncoding_rna.py \\
        --ncrna-fitting .../Non_coding_genes_Gene_level_statistics_fitted.tsv \\
        --ncrna-bed resources/external/pombase/{ver}/genome_region/non_coding_rna.bed \\
        --gtrnadb-bed resources/external/pombase/schiPomb_972H-tRNAs.bed \\
        --marguerat-excel resources/literature/margueratQuantitativeAnalysisFission2012.xlsx \\
        --output-stats results/noncoding_rna/{dataset}/ncrna_stats.tsv \\
        --output-figures results/noncoding_rna/{dataset}/ncrna_analysis.pdf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-21
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import re
import sys
from dataclasses import dataclass
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
from workflow.src.plotting.generic import donut_chart  # noqa: E402
from workflow.src.plotting.style import AX_HEIGHT, AX_WIDTH, COLORS  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Byte-faithful to the source notebook's Config: GtRNAdb bed chromosome names
# are "chrI/chrII/chrIII" but the ncRNA bed uses bare "I/II/III", so the
# GtRNAdb chr column is normalized before the positional merge.
_CHROMOSOME_NAME_MAP = {"chrI": "I", "chrII": "II", "chrIII": "III"}

# Legacy -> current metric column names, same quirk as coverage.smk /
# verification: the non-coding-gene fitting_results.tsv still ships the
# pre-rename um/lam headers instead of DR/DL.
_LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Headerless GtRNAdb BED column names (standard 12-column BED); only the first
# 6 are used downstream (GtRNAdb_Name carries the amino-acid/anticodon label).
_GTRNADB_BED_COLUMNS = [
    "#Chr", "Start", "End", "GtRNAdb_Name", "Score", "Strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts",
]

# Marguerat 2012 Table_S2 replicate columns -> condition mean, byte-faithful to
# the notebook's mRNA_abundance cell (proliferating = MM, nitrogen-starved = MN).
_MARGUERAT_CONDITIONS = {
    "EMM_Proliferating": ["MM1.tot.cpc_ex", "MM2.tot.cpc_ex"],
    "EMM_Nitrogen_Starved": ["MN1.tot.cpc_ex", "MN2.tot.cpc_ex"],
}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class NoncodingRNAConfig:
    """Inputs and outputs for the non-coding RNA depletion analysis."""
    ncrna_fitting: Path
    ncrna_bed: Path
    gtrnadb_bed: Path
    marguerat_excel: Path
    output_stats: Path
    output_figures: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dirs exist."""
        for path in [self.ncrna_fitting, self.ncrna_bed, self.gtrnadb_bed, self.marguerat_excel]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        for out in [self.output_stats, self.output_figures]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


def load_ncrna_fitting(ncrna_fitting_path: Path) -> pd.DataFrame:
    """Load non-coding-gene fitting statistics, normalizing legacy um/lam -> DR/DL columns."""
    fitting = pd.read_csv(ncrna_fitting_path, sep="\t")
    rename = {
        old: new
        for old, new in _LEGACY_METRIC_RENAME.items()
        if old in fitting.columns and new not in fitting.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        fitting = fitting.rename(columns=rename)
    return fitting


def load_gtrnadb(gtrnadb_bed_path: Path) -> pd.DataFrame:
    """Load the headerless GtRNAdb tRNA bed and normalize its chromosome names."""
    gtrnadb = pd.read_csv(
        gtrnadb_bed_path, sep="\t", header=None, names=_GTRNADB_BED_COLUMNS
    )[["#Chr", "Start", "End", "GtRNAdb_Name", "Score", "Strand"]]
    return normalize_chromosome_names(gtrnadb)


def load_marguerat_abundance(marguerat_excel_path: Path) -> pd.DataFrame:
    """Load Marguerat 2012 Table_S2 and reduce to per-condition mean mRNA abundance.

    Indexed on Systematic.name; the returned frame has one column per condition
    (EMM_Proliferating, EMM_Nitrogen_Starved) holding the mean of that
    condition's two replicate `.tot.cpc_ex` columns.
    """
    abundance = pd.read_excel(
        marguerat_excel_path, sheet_name="Table_S2", comment="#"
    ).set_index("Systematic.name")
    means = pd.DataFrame(
        {condition: abundance[cols].mean(axis=1) for condition, cols in _MARGUERAT_CONDITIONS.items()}
    )
    means.columns = [f"{c}_RNA_Abundance_mean" for c in means.columns]
    return means


# =============================================================================
# CORE LOGIC — merge + parse (unit-tested)
# =============================================================================
def normalize_chromosome_names(df: pd.DataFrame) -> pd.DataFrame:
    """Replace chrI/chrII/chrIII with I/II/III in the #Chr column (source-notebook quirk)."""
    return df.replace({"#Chr": _CHROMOSOME_NAME_MAP})


def merge_gtrnadb_by_position(ncrna: pd.DataFrame, gtrnadb: pd.DataFrame) -> pd.DataFrame:
    """Left-merge the GtRNAdb_Name onto the ncRNA table by genomic position.

    Merge key is ["#Chr", "Start", "End"] — NOT gene name — because ncRNA
    Systematic IDs and GtRNAdb names use different vocabularies but share exact
    genomic coordinates. ncRNA rows without a coordinate match get NaN
    GtRNAdb_Name (left join).
    """
    return ncrna.merge(
        gtrnadb[["#Chr", "Start", "End", "GtRNAdb_Name"]],
        on=["#Chr", "Start", "End"],
        how="left",
    )


def extract_tRNA_amino_acid_and_anticodon(row: pd.Series) -> pd.Series:
    """Parse (Amino_Acid, Anticodon) for one tRNA row.

    Amino acid from the Systematic ID via ``TRNA(\\w+)\\.`` (e.g. "SPATRNAPRO.01"
    -> "PRO"); anticodon from the GtRNAdb_Name's 3rd hyphen field (e.g.
    "tRNA-Ala-AGC-1-1" -> "AGC"). Byte-faithful to the notebook. Rows where
    GtRNAdb_Name is NaN, or the sysID lacks a TRNA<AA>. token, or the name has
    fewer than 3 hyphen fields, yield None for the missing piece rather than
    raising.
    """
    sys_id = row["Systematic ID"]
    trna_name = row["GtRNAdb_Name"]

    amino_acid = None
    anticodon = None
    if pd.notna(trna_name):
        match = re.search(r"TRNA(\w+)\.", str(sys_id))
        amino_acid = match.group(1) if match else None
        parts = str(trna_name).split("-")
        anticodon = parts[2] if len(parts) > 2 else None

    return pd.Series({"Amino_Acid": amino_acid, "Anticodon": anticodon})


def compute_tRNA_copy_number(df: pd.DataFrame) -> pd.DataFrame:
    """Add tRNA_copy_number: count of tRNAs sharing the same Amino_Acid + Anticodon."""
    result = df.copy()
    result["tRNA_copy_number"] = result.groupby(["Amino_Acid", "Anticodon"])["Systematic ID"].transform("count")
    return result


# =============================================================================
# CORE LOGIC — assembly
# =============================================================================
def build_ncrna_table(
    ncrna_fitting: pd.DataFrame,
    ncrna_bed: pd.DataFrame,
    gtrnadb: pd.DataFrame,
    marguerat_means: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble the full annotated ncRNA table (bed + GtRNAdb + fitting + abundance).

    Mirrors the notebook's Config ordering: GtRNAdb_Name is merged onto the bed
    by position first, then the per-gene fitting stats are merged on Systematic
    ID (outer join, so unfitted annotated genes and unannotated fitted genes are
    both kept), then Marguerat abundance means are left-joined on Systematic ID.
    """
    bed_cols = ["Systematic ID", "Name", "#Chr", "Start", "End", "Strand", "Feature", "Type"]
    meta = merge_gtrnadb_by_position(ncrna_bed, gtrnadb)
    meta = meta[[c for c in bed_cols if c in meta.columns] + ["GtRNAdb_Name"]]

    fitting = ncrna_fitting.drop(columns=["Name"], errors="ignore")
    combined = meta.merge(fitting, on="Systematic ID", how="outer")

    combined = combined.merge(marguerat_means, left_on="Systematic ID", right_index=True, how="left")
    return combined


def select_nuclear_tRNAs(combined: pd.DataFrame) -> pd.DataFrame:
    """Filter to nuclear tRNAs, then annotate amino acid / anticodon / copy number.

    Nuclear tRNA = ``Feature == 'tRNA' and #Chr != 'mitochondrial'`` (notebook
    quirk). Sorted by copy number then DR (descending) for a stable, readable
    stats table.
    """
    # Boolean indexing (not .query) because the "#Chr" column name breaks
    # pandas' backtick parser in some versions.
    nuclear = combined[(combined["Feature"] == "tRNA") & (combined["#Chr"] != "mitochondrial")].copy()
    nuclear[["Amino_Acid", "Anticodon"]] = nuclear.apply(extract_tRNA_amino_acid_and_anticodon, axis=1)
    nuclear = compute_tRNA_copy_number(nuclear)
    sort_cols = [c for c in ["tRNA_copy_number", "DR"] if c in nuclear.columns]
    if sort_cols:
        nuclear = nuclear.sort_values(sort_cols, ascending=[True] + [False] * (len(sort_cols) - 1))
    return nuclear.reset_index(drop=True)


# =============================================================================
# PLOTTING
# =============================================================================
def plot_feature_type_donut(combined: pd.DataFrame) -> plt.Figure:
    """Donut chart of annotated ncRNA genes per Feature type (tRNA, lncRNA, ...)."""
    counts = combined["Feature"].dropna().value_counts()
    palette = (COLORS * (len(counts) // len(COLORS) + 1))[: len(counts)]
    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    donut_chart(
        values=list(counts.values),
        labels=list(counts.index),
        colors=palette,
        center_text=f"Total\n{int(counts.sum()):,}\nncRNA genes",
        ax=ax,
    )
    ax.set_title("Non-coding RNA gene types")
    fig.tight_layout()
    return fig


def plot_trna_summary(nuclear_trnas: pd.DataFrame) -> plt.Figure:
    """tRNA copy-number distribution (bar) + DR-by-copy-number scatter."""
    fig, (ax_hist, ax_scatter) = plt.subplots(1, 2, figsize=(AX_WIDTH * 2, AX_HEIGHT))

    copy_counts = nuclear_trnas["tRNA_copy_number"].dropna().astype(int).value_counts().sort_index()
    ax_hist.bar(copy_counts.index, copy_counts.values, color=COLORS[0])
    ax_hist.set_xlabel("tRNA copy number (shared amino acid + anticodon)")
    ax_hist.set_ylabel("Number of tRNA genes")
    ax_hist.set_title("tRNA copy-number distribution")

    if "DR" in nuclear_trnas.columns:
        fitted = nuclear_trnas.dropna(subset=["DR", "tRNA_copy_number"])
        rng = np.random.default_rng(42)
        jitter = rng.uniform(-0.2, 0.2, size=len(fitted))
        ax_scatter.scatter(
            fitted["tRNA_copy_number"] + jitter, fitted["DR"],
            alpha=0.5, s=12, color=COLORS[1],
        )
        ax_scatter.set_xlabel("tRNA copy number")
        ax_scatter.set_ylabel("Depletion Rate (DR)")
        ax_scatter.set_title(f"tRNA depletion by copy number (n={len(fitted):,})")

    fig.tight_layout()
    return fig


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
@logger.catch(reraise=True)
def run(config: NoncodingRNAConfig) -> None:
    """Load -> merge -> characterize nuclear tRNAs -> save TSV + figures."""
    config.validate()

    ncrna_fitting = load_ncrna_fitting(config.ncrna_fitting)
    ncrna_bed = pd.read_csv(config.ncrna_bed, sep="\t")
    gtrnadb = load_gtrnadb(config.gtrnadb_bed)
    marguerat_means = load_marguerat_abundance(config.marguerat_excel)

    combined = build_ncrna_table(ncrna_fitting, ncrna_bed, gtrnadb, marguerat_means)
    nuclear_trnas = select_nuclear_tRNAs(combined)

    nuclear_trnas.to_csv(config.output_stats, sep="\t", index=False)

    fig_donut = plot_feature_type_donut(combined)
    fig_trna = plot_trna_summary(nuclear_trnas)
    with PdfPages(config.output_figures) as pdf:
        pdf.savefig(fig_donut, dpi=300, bbox_inches="tight")
        pdf.savefig(fig_trna, dpi=300, bbox_inches="tight")
    plt.close(fig_donut)
    plt.close(fig_trna)

    fitted = int(nuclear_trnas["DR"].notna().sum()) if "DR" in nuclear_trnas.columns else 0
    logger.success(
        f"Non-coding RNA: {len(nuclear_trnas):,} nuclear tRNAs "
        f"({fitted:,} with DR), {nuclear_trnas['Anticodon'].nunique():,} distinct anticodons"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Analyze non-coding RNA (tRNA) depletion patterns")
    parser.add_argument("--ncrna-fitting", type=Path, required=True, help="Non-coding-gene fitting_results.tsv")
    parser.add_argument("--ncrna-bed", type=Path, required=True, help="ncRNA genome-region bed")
    parser.add_argument("--gtrnadb-bed", type=Path, required=True, help="GtRNAdb tRNA bed (schiPomb_972H-tRNAs.bed)")
    parser.add_argument("--marguerat-excel", type=Path, required=True, help="Marguerat 2012 abundance xlsx")
    parser.add_argument("--output-stats", type=Path, required=True, help="Output ncRNA stats TSV")
    parser.add_argument("--output-figures", type=Path, required=True, help="Output ncRNA figures PDF")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the analysis, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = NoncodingRNAConfig(
            ncrna_fitting=args.ncrna_fitting,
            ncrna_bed=args.ncrna_bed,
            gtrnadb_bed=args.gtrnadb_bed,
            marguerat_excel=args.marguerat_excel,
            output_stats=args.output_stats,
            output_figures=args.output_figures,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
