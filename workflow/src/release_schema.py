"""
Upstream Release Schema
========================

The spellings, signs and column filters of the tables this repo reads from
DIT_HAP_snakemake's ``release/`` output and from ``resources/curated/``. They used
to be byte-identical copies in ten modules, each carrying a "same quirk as ..."
comment pointing at the others; they live here so a schema change has one place
to land.

Input
-----
- Nothing: three constants plus two frame-in/frame-out helpers.

Output
------
- ``LEGACY_METRIC_RENAME``, ``IN_GENE_FILTER``, ``GRNA_METRIC_SIGN``
- ``normalize_legacy_metrics()`` and ``read_gene_level()``, which apply the rename

Usage
-----
    from release_schema import IN_GENE_FILTER, read_gene_level
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Legacy -> current metric column names. Upstream's fit computes a depletion rate
# `um` and a lag `lam`; the columns were renamed DR/DL, but older release tables
# still carry the old headers. Renaming on load keeps every stage speaking DR/DL,
# and is a no-op once a table ships DR/DL directly.
LEGACY_METRIC_RENAME = {"um": "DR", "lam": "DL"}

# Byte-faithful to the source notebook's Config.in_gene_filter: an insertion counts
# as "in a gene" only if it is annotated as non-intergenic AND at least 5bp
# upstream of the stop codon (the >4 threshold, not >=5, is the notebook's own
# quirk — kept verbatim).
IN_GENE_FILTER = "Type != 'Intergenic region' and Distance_to_stop_codon > 4"

# The curated gRNA fitted-parameters table (resources/curated/*_gRNA_HDdata_*.tsv)
# is frozen at the pre-2026-09-17 sign convention — positive = depleted — while
# upstream flipped DIT-HAP so negative DR is now the depleted end. Flip the gRNA
# metric on the way in so both point the same way. Without this, comparison
# reports a strong ANTI-correlation: measured r = -0.92 where the two studies
# actually agree at r = +0.92 (n = 4,465). (`lam` is a lag in generations, not a
# signed magnitude, and is left alone.)
GRNA_METRIC_SIGN = -1.0


# =============================================================================
# HELPERS
# =============================================================================
def normalize_legacy_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Rename legacy um/lam metric columns to DR/DL when the new names are absent.

    The "new name absent" guard makes this idempotent: a table that already ships
    DR/DL is returned unchanged rather than gaining duplicate columns.
    """
    rename = {
        old: new
        for old, new in LEGACY_METRIC_RENAME.items()
        if old in df.columns and new not in df.columns
    }
    if rename:
        logger.info(f"Normalizing legacy metric columns: {rename}")
        df = df.rename(columns=rename)
    return df


def read_gene_level(gene_level_path: Path) -> pd.DataFrame:
    """Read a gene-level fitting_statistics/fitting_results table, um/lam normalized to DR/DL."""
    return normalize_legacy_metrics(pd.read_csv(gene_level_path, sep="\t"))
