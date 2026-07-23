#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare PCR QC Data
====================

Stage 1 of the PCR / library-prep QC split: load the 6 raw merged-reads /
spike-in inputs, merge the technical- and biological-replicate pairs, then
write 4 parquet intermediates consumed by the figure-rendering rule:

- pbl_pbr.parquet: panel (a) PBL vs PBR reads of one library.
- tech.parquet: panel (b) technical replicate, merged on (Chr, Coordinate, Strand).
- bio.parquet: panel (c) biological replicate, merged on (Chr, Coordinate, Strand).
- spikein.parquet: panel (d) spike-in dilution table (currently a placeholder).

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-22
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.io import write_parquet  # noqa: E402
from workflow.src.pcr_qc.core import read_merged_reads  # noqa: E402


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PCRQCConfig:
    """Resolved input/output paths for the PCR QC data preparation.

    `output` doubles as the pbl_pbr.parquet output path: it and the other three
    parquet outputs (tech/bio/spikein, passed separately to run()) all land
    under the same results/pcr_qc/_work/ directory, so mkdir'ing its parent
    covers all four.
    """
    pbl_pbr: Path
    tech_rep_1: Path
    tech_rep_2: Path
    bio_rep_1: Path
    bio_rep_2: Path
    spikein: Path
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any input is missing, then ensure the output dir exists."""
        for path in [self.pbl_pbr, self.tech_rep_1, self.tech_rep_2,
                     self.bio_rep_1, self.bio_rep_2, self.spikein]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: PCRQCConfig, *, output_pbl_pbr: Path, output_tech: Path, output_bio: Path, output_spikein: Path) -> None:
    """Load -> merge -> write the four parquet intermediates."""
    config.validate()

    # Panel (a): PBL vs PBR of a single library.
    pbl_pbr = read_merged_reads(config.pbl_pbr)

    # Panel (b): technical replicate — same sample, two upstream projects.
    tech = pd.merge(
        read_merged_reads(config.tech_rep_1), read_merged_reads(config.tech_rep_2),
        left_index=True, right_index=True, suffixes=("_1", "_2"),
    )

    # Panel (c): biological replicate — two samples, one project.
    bio = pd.merge(
        read_merged_reads(config.bio_rep_1), read_merged_reads(config.bio_rep_2),
        left_index=True, right_index=True, suffixes=("_1", "_2"),
    )

    # Panel (d): spike-in linearity.
    spikein = pd.read_csv(config.spikein, sep="\t")

    write_parquet(pbl_pbr, output_pbl_pbr)
    write_parquet(tech, output_tech)
    write_parquet(bio, output_bio)
    write_parquet(spikein, output_spikein)

    logger.success(
        f"Prepared PCR QC tables: {len(pbl_pbr):,} pbl_pbr rows, "
        f"{len(tech):,} tech rows, {len(bio):,} bio rows, {len(spikein):,} spikein rows"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Prepare PCR QC parquet intermediates")
    parser.add_argument("--pbl-pbr", type=Path, required=True, help="Panel (a): merged reads TSV (PBL vs PBR)")
    parser.add_argument("--tech-rep-1", type=Path, required=True, help="Panel (b): technical replicate 1 merged reads TSV")
    parser.add_argument("--tech-rep-2", type=Path, required=True, help="Panel (b): technical replicate 2 merged reads TSV")
    parser.add_argument("--bio-rep-1", type=Path, required=True, help="Panel (c): biological replicate 1 merged reads TSV")
    parser.add_argument("--bio-rep-2", type=Path, required=True, help="Panel (c): biological replicate 2 merged reads TSV")
    parser.add_argument("--spikein", type=Path, required=True, help="Panel (d): spike-in results TSV")
    parser.add_argument("--output-pbl-pbr", type=Path, required=True, help="Output pbl_pbr.parquet")
    parser.add_argument("--output-tech", type=Path, required=True, help="Output tech.parquet")
    parser.add_argument("--output-bio", type=Path, required=True, help="Output bio.parquet")
    parser.add_argument("--output-spikein", type=Path, required=True, help="Output spikein.parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run the preparation, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = PCRQCConfig(
            pbl_pbr=args.pbl_pbr,
            tech_rep_1=args.tech_rep_1,
            tech_rep_2=args.tech_rep_2,
            bio_rep_1=args.bio_rep_1,
            bio_rep_2=args.bio_rep_2,
            spikein=args.spikein,
            output=args.output_pbl_pbr,
        )
        run(
            config,
            output_pbl_pbr=args.output_pbl_pbr,
            output_tech=args.output_tech,
            output_bio=args.output_bio,
            output_spikein=args.output_spikein,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
