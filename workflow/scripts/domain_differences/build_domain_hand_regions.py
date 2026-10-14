#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build Domain Hand Regions
=========================

Stage 4a, step 2: turn the hand-curated terminus spreadsheet into explicit
per-gene residue intervals.

The spreadsheet (``resources/curated/non_esential_domain_candidates.xlsx``, three
sheets) records a manual call of the form "the last X residues of this protein are
not essential" by giving the first / last affected residue. Protein length comes
from PomBase, so each row becomes a concrete interval:

    Non-essential C terminus  -> non-essential region = [affected_C, protein_len]
    Non-essential N terminus  -> non-essential region = [1, affected_N]
    both                      -> the two intervals above

The spreadsheet also records the bp distance from the affected residue to the
stop codon, which makes the reading checkable: ``(protein_len - affected_C + 1) * 3``
must reproduce it. The check is reported in the log rather than enforced — the
handful of rows that miss are annotated in the spreadsheet's Comment column as
wider / polar-effect calls.

Input
-----
- resources/curated/non_esential_domain_candidates.xlsx
- resources/external/pombase/{pombase_version}/Protein_features/peptide_stats.tsv
- release gene_level/fitting_results.tsv — gene name -> Systematic ID

Output
------
- domain_hand_regions.tsv — one row per spreadsheet row: ``Gene, sid, term,
  protein_len, regions, n_regions, affected_N, affected_C, dist_start_bp,
  dist_stop_bp, comment, bp_check``. ``regions`` is ``a-b|c-d``, ready to join
  onto the segmentation.

Usage
-----
    python workflow/scripts/domain_differences/build_domain_hand_regions.py \\
        --xlsx resources/curated/non_esential_domain_candidates.xlsx \\
        --peptide-stats resources/external/pombase/2026-06-01/Protein_features/peptide_stats.tsv \\
        --gene-level <release>/gene_level/fitting_results.tsv \\
        --output results/4a_domain_differences/HD_DIT_HAP/domain_hand_regions.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
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

# 2. Third-party Imports
from loguru import logger

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from domain_differences.core import gene_name_to_sid, read_hand_regions  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class HandRegionsConfig:
    """Runtime parameters for the hand-region build."""

    xlsx: Path
    peptide_stats: Path
    gene_level: Path
    output: Path


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: HandRegionsConfig) -> int:
    """Parse the spreadsheet into residue intervals and write them as TSV."""
    hand = read_hand_regions(config.xlsx, config.peptide_stats, gene_name_to_sid(config.gene_level))
    hand.to_csv(config.output, index=False, sep="\t")

    mapped = hand["sid"].notna()
    logger.info(f"{len(hand)} spreadsheet rows: {int(mapped.sum())} mapped to a Systematic ID, "
                f"{int(hand['n_regions'].gt(0).sum())} with a usable region")
    for term, group in hand.groupby("term", sort=False):
        logger.info(f"  {term:<5} {len(group):>3} rows, {int(group['sid'].notna().sum()):>3} mapped")
    if not mapped.all():
        logger.warning(f"{int((~mapped).sum())} rows have no Systematic ID: "
                       f"{sorted(hand.loc[~mapped, 'Gene'].tolist())}")

    checked = hand["bp_check"].ne("")
    agreed = hand["bp_check"].eq("ok")
    logger.info(f"bp check: {int(agreed.sum())}/{int(checked.sum())} rows reproduce the recorded "
                f"stop-codon distance within 3 bp")
    for _, row in hand[hand["bp_check"].str.startswith("off")].iterrows():
        logger.info(f"  off: {row['Gene']} ({row['term']}) {row['bp_check']} {row['comment']}")
    logger.info(f"wrote {config.output}")
    return 0


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Build the hand-curated non-essential regions for 4a")
    parser.add_argument("--xlsx", type=Path, required=True, help="Curated candidates spreadsheet")
    parser.add_argument("--peptide-stats", type=Path, required=True,
                        help="PomBase Protein_features/peptide_stats.tsv")
    parser.add_argument("--gene-level", type=Path, required=True,
                        help="Release gene_level/fitting_results.tsv")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the parse, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = HandRegionsConfig(
            xlsx=args.xlsx.resolve(),
            peptide_stats=args.peptide_stats.resolve(),
            gene_level=args.gene_level.resolve(),
            output=args.output.resolve(),
        )
        config.output.parent.mkdir(parents=True, exist_ok=True)
        run(config)
    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
