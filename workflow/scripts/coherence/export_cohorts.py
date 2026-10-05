#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Cohort Export (threshold-filtered Excel tables) — Export Only
=======================================================================

Turns the coherence metrics table into the two cohort tables the analysis
actually reads — the COHERENT groups (significantly tighter than random) and the
INCOHERENT ones (dispersed) — as one Excel workbook with a sheet per cohort plus
the full annotated table. The cohort cuts are the same ones the figures label
with (config/analysis.yaml's `coherence:` block):

- coherent   = q_value <= coherent_q_max AND median_pairwise_distance_z < coherent_z_threshold
- incoherent = median_pairwise_distance_z > incoherent_z_threshold — a DESCRIPTIVE
  cut, not an FDR-controlled one: at n_permutations=1000 the upper-tail p floor is
  1/1001, so no group reaches q_up <= 0.05 and the incoherent side cannot carry a
  significance claim. Everything else is `other`.

Input
-----
- --metrics: a coherence metrics table. The dataset-level
  combined/coherence_metrics.parquet is the intended one (all sources, q
  re-derived over the union, which is the q the de-duplication and attribution
  stages use); a per-source coherence_metrics.parquet works too — same columns,
  but its q is the per-source BH and not cross-source comparable.
- --source (optional): keep one source's rows only, writing that source's workbook
  into its own folder. The per-source rule runs this on that source's OWN
  coherence_metrics.parquet, so the cohort is called from the per-source BH family
  and the workbook is self-contained — no `--dedup-terms`, no de-duplication
  column. Slicing the COMBINED table with this flag also works, giving the same
  cut on the pooled q instead (58 of the 322 coherent groups flip between the two,
  2026-09-29); that is the job of the combined/ workbook, not this one.
- --dedup-terms (optional): dedup/coherence_terms_deduplicated.tsv — every term
  with its cluster and representative flags. Adds `in_dedup_set` (the
  is_representative column) so a reader can see which rows of a cohort survive
  de-duplication, and carries `moonlighting_fraction` over with it. Joined on
  (source, group_id) — group_id alone is NOT unique across sources (173 collide),
  so the pair is the key everywhere. A table without is_representative (e.g. the
  representatives-only TSV) still works: presence in it IS the flag.

Output
------
- --output: coherence_cohorts.xlsx — sheets `coherent` (z ascending, tightest
  first), `incoherent` (z descending, most dispersed first), `all_terms` (every
  row plus the `cohort` label, original order) and `thresholds` (the cuts, the row
  counts and the source subset, so the workbook stays readable without the config).
  List-valued cells (`scored_member_names`) are joined with ", " on the way out:
  Excel would otherwise render the raw numpy repr, quotes and line breaks included.
  The Parquet keeps the real list — it is only the display that flattens.
  Every other column of the metrics table rides along untouched, so the descriptive
  metrics appear without this script knowing them: `paralog_fraction` is already a
  column there, and `moonlighting_fraction` arrives through --dedup-terms.

Usage
-----
    # The `combined` view: every source's rows, one q family.
    python export_cohorts.py \\
        --metrics results/3a_coherence/{dataset}/combined/coherence_metrics.parquet \\
        --dedup-terms results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv \\
        --q-max 0.05 --coherent-z -2.0 --incoherent-z 1.0 \\
        --output results/3a_coherence/{dataset}/combined/coherence_cohorts.xlsx

    # The `dedup` view, over the representative set alone.
    python export_cohorts.py \\
        --metrics results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv \\
        --output results/3a_coherence/{dataset}/dedup/coherence_cohorts.xlsx

    # A single source: its own metrics table, its own q family.
    python export_cohorts.py \\
        --metrics results/3a_coherence/{dataset}/go_cc/coherence_metrics.parquet \\
        --source go_cc \\
        --output results/3a_coherence/{dataset}/go_cc/coherence_cohorts.xlsx

Author:   Yusheng Yang (guidance) + Claude Sonnet 4.6 (implementation)
Date:     2026-09-29
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from io_table import read_file  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
class Cohort(StrEnum):
    """The three cohort labels written into the `cohort` column."""
    COHERENT = "coherent"
    INCOHERENT = "incoherent"
    OTHER = "other"


_Z_COLUMN = "median_pairwise_distance_z"
_Q_COLUMN = "q_value"
_KEY_COLUMNS = ["source", "group_id"]

# Column widths for the workbook (characters). Only the two long text columns and
# the source key are widened; everything else is numeric and fits the default.
_COLUMN_WIDTHS = {"group_name": 46, "scored_member_names": 60, "source": 16}
_DEFAULT_WIDTH = 13


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ExportConfig:
    """Inputs, cohort thresholds, source subset and output path for the cohort export."""
    metrics: Path
    output: Path
    dedup_terms: Path | None = None
    source: str | None = None
    q_max: float = 0.05
    coherent_z: float = -2.0
    incoherent_z: float = 1.0

    def validate(self) -> None:
        """Raise ValueError if an input is missing, then create the output dir."""
        if not self.metrics.exists():
            raise ValueError(f"Required input not found: {self.metrics}")
        if self.dedup_terms is not None and not self.dedup_terms.exists():
            raise ValueError(f"Required input not found: {self.dedup_terms}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
def read_metrics(config: ExportConfig) -> pd.DataFrame:
    """Read the metrics table, optionally keep one source, and check the cohort columns."""
    table = read_file(config.metrics)
    missing = [c for c in (_Z_COLUMN, _Q_COLUMN) if c not in table.columns]
    if missing:
        raise ValueError(f"{config.metrics.name} is missing required column(s): {missing}")
    logger.info(f"Read {len(table):,} groups from {config.metrics.name}")

    if config.source is not None:
        if "source" not in table.columns:
            raise ValueError(f"--source given but {config.metrics.name} has no `source` column")
        table = table[table["source"] == config.source].reset_index(drop=True)
        if table.empty:
            logger.warning(f"No rows for source {config.source!r} in {config.metrics.name}")
        else:
            logger.info(f"Kept {len(table):,} groups for source {config.source}")
    return table


def attach_dedup_columns(table: pd.DataFrame, dedup_terms: Path | None) -> pd.DataFrame:
    """Add `in_dedup_set` and the dedup columns worth carrying (moonlighting_fraction)."""
    # Guard: a table without both key columns (e.g. a hand-made subset) has nothing
    # to join on, so the columns are simply not added rather than silently mis-joined.
    if dedup_terms is None:
        return table
    if not all(key in table.columns for key in _KEY_COLUMNS):
        logger.warning(f"Skipping --dedup-terms: no {_KEY_COLUMNS} in the metrics table")
        return table

    dedup = read_file(dedup_terms)
    # `paralog_fraction` needs no carrying: it is already a column of the metrics
    # table itself. `moonlighting_fraction` exists only downstream, in the dedup
    # table, so it is joined here or not at all. `is_representative` IS the flag when
    # the all-terms table is passed; a representatives-only table has no such column,
    # and presence in it IS the flag — hence the two branches below.
    flag_column = "is_representative" if "is_representative" in dedup.columns else None
    carried = [column for column in ("moonlighting_fraction",) if column in dedup.columns]
    kept = dedup[_KEY_COLUMNS + [c for c in (flag_column, *carried) if c]].drop_duplicates(subset=_KEY_COLUMNS)
    if flag_column is None:
        kept["in_dedup_set"] = True
    else:
        kept = kept.rename(columns={flag_column: "in_dedup_set"})
    merged = table.merge(kept, on=_KEY_COLUMNS, how="left")
    merged["in_dedup_set"] = merged["in_dedup_set"].fillna(False).astype(bool)
    logger.info(f"{int(merged['in_dedup_set'].sum()):,} of {len(merged):,} groups are in the de-duplicated set")
    return merged


def label_cohorts(table: pd.DataFrame, config: ExportConfig) -> pd.DataFrame:
    """Append the `cohort` label from the configured thresholds."""
    z = table[_Z_COLUMN]
    coherent = (table[_Q_COLUMN] <= config.q_max) & (z < config.coherent_z)
    incoherent = z > config.incoherent_z    # disjoint from `coherent`: z < -2 vs z > 1

    labelled = table.copy()
    labelled["cohort"] = Cohort.OTHER.value
    labelled.loc[incoherent, "cohort"] = Cohort.INCOHERENT.value
    labelled.loc[coherent, "cohort"] = Cohort.COHERENT.value
    logger.info(
        f"Cohorts: {int(coherent.sum()):,} coherent "
        f"({_Q_COLUMN}<={config.q_max} and z<{config.coherent_z}), "
        f"{int(incoherent.sum()):,} incoherent (z>{config.incoherent_z})"
    )
    return labelled


def build_sheets(labelled: pd.DataFrame, config: ExportConfig) -> dict[str, pd.DataFrame]:
    """The four sheets, in workbook order."""
    coherent = labelled[labelled["cohort"] == Cohort.COHERENT.value].sort_values(_Z_COLUMN)
    incoherent = labelled[labelled["cohort"] == Cohort.INCOHERENT.value].sort_values(_Z_COLUMN, ascending=False)
    thresholds = pd.DataFrame(
        {
            "cut": ["coherent_q_max", "coherent_z_threshold", "incoherent_z_threshold",
                    "n_coherent", "n_incoherent", "n_total", "source_subset"],
            "value": [config.q_max, config.coherent_z, config.incoherent_z,
                      len(coherent), len(incoherent), len(labelled),
                      config.source or "all sources"],
        }
    )
    return {
        Cohort.COHERENT.value: coherent,
        Cohort.INCOHERENT.value: incoherent,
        "all_terms": labelled,
        "thresholds": thresholds,
    }


def flatten_list_cells(frame: pd.DataFrame) -> pd.DataFrame:
    """Join iterable cells (a list column such as scored_member_names) into one comma-separated string."""
    out = frame.copy()
    for column in out.columns:
        sample = out[column].dropna()
        # Duck-typed instead of isinstance(list, ndarray): a Parquet list column comes
        # back as a numpy array, and `str` is the one iterable that must NOT be split.
        if len(sample) and not isinstance(sample.iloc[0], str) and hasattr(sample.iloc[0], "__iter__"):
            out[column] = [
                ", ".join(str(gene) for gene in cell)
                if hasattr(cell, "__iter__") and not isinstance(cell, str)
                else ""
                for cell in out[column]
            ]
    return out


def write_workbook(sheets: dict[str, pd.DataFrame], output: Path) -> None:
    """Write one sheet per table, with a frozen header and readable column widths."""
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            flatten_list_cells(frame).to_excel(writer, sheet_name=name, index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            for (header,) in sheet.iter_cols(min_row=1, max_row=1):
                sheet.column_dimensions[header.column_letter].width = _COLUMN_WIDTHS.get(
                    str(header.value), _DEFAULT_WIDTH
                )
    logger.info(f"Wrote {output} ({output.stat().st_size / 1024:.0f} KB, {len(sheets)} sheets)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Export threshold-filtered coherence cohort tables to Excel.")
    parser.add_argument("--metrics", type=Path, required=True, help="Coherence metrics table (parquet/tsv)")
    parser.add_argument("--dedup-terms", type=Path, default=None,
                        help="dedup/coherence_terms_deduplicated.tsv; adds in_dedup_set + moonlighting_fraction")
    parser.add_argument("--source", default=None,
                        help="Keep one source's rows only (run this on the combined table)")
    parser.add_argument("--q-max", type=float, default=0.05, help="Coherent cohort: q_value at or below this")
    parser.add_argument("--coherent-z", type=float, default=-2.0, help="Coherent cohort: z below this")
    parser.add_argument("--incoherent-z", type=float, default=1.0, help="Incoherent cohort: z above this")
    parser.add_argument("--output", type=Path, required=True, help="Output .xlsx path")
    parser.add_argument("--verbose", action="store_true", help="Enable debug-level logging")
    return parser.parse_args()


def main() -> int:
    """Read the metrics table, label the cohorts and write the workbook."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    try:
        config = ExportConfig(
            metrics=args.metrics,
            output=args.output,
            dedup_terms=args.dedup_terms,
            source=args.source,
            q_max=args.q_max,
            coherent_z=args.coherent_z,
            incoherent_z=args.incoherent_z,
        )
        config.validate()
        table = attach_dedup_columns(read_metrics(config), config.dedup_terms)
        write_workbook(build_sheets(label_cohorts(table, config), config), config.output)
    except Exception as e:
        logger.exception(f"Cohort export failed: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
