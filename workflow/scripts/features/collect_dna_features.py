#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DNA-Level Feature Collection (the spine)
=========================================

Builds the gffutils DB from a PomBase GFF3, computes per-mRNA DNA-level
features, and emits the anti-codon usage matrix. This is the first stage of
the feature pipeline: its Gene_id column enumerates the coding-gene set that
every other level filters against, so the DNA table is written as a parquet
(preserving the bool `Primary_candidate` dtype) for downstream levels to read.

Input
-----
- A PomBase version directory (genome FASTA/GFF3, peptide_stats.tsv)

Output
------
- dna_features.parquet: one row per mRNA (full set, before Primary_candidate filter)
- codon_usage_matrix.tsv: gene x anti-codon count matrix

Usage
-----
    python collect_dna_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --output results/features/2025-10-01/_levels/dna_features.parquet \\
        --codon-usage-output results/features/2025-10-01/codon_usage_matrix.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
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
import gffutils
from loguru import logger

# 3. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.features.assembly import collect_dna_level_features
from workflow.src.features.genome import PombaseGenomeConfig, calculate_anticodon_usage_matrix
from workflow.src.io import write_parquet


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class DnaConfig:
    """Inputs/outputs for DNA-level feature collection."""
    pombase_dir: Path
    genome_landmarks: Path
    output_dna: Path
    output_codon_usage: Path

    def validate(self) -> None:
        """Raise ValueError if a required input is missing, then ensure output dirs exist."""
        for path in [self.pombase_dir, self.genome_landmarks]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_dna.parent.mkdir(parents=True, exist_ok=True)
        self.output_codon_usage.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: DnaConfig) -> None:
    """Build the gffutils DB, collect DNA-level features, and write the codon usage matrix."""
    logger.info(f"Building gffutils DB from {config.pombase_dir}")
    genome_cfg = PombaseGenomeConfig.from_pombase_dir(config.pombase_dir, config.genome_landmarks)
    gffutils.create_db(genome_cfg.gff3_file, genome_cfg.database_file, force=True, merge_strategy="create_unique")
    db = gffutils.FeatureDB(genome_cfg.database_file)

    logger.info("Collecting DNA-level features")
    dna_df, coding_genes = collect_dna_level_features(db, genome_cfg)
    write_parquet(dna_df, config.output_dna)
    logger.success(f"Wrote {len(dna_df)} mRNA rows ({len(coding_genes)} coding genes) to {config.output_dna}")

    logger.info("Writing codon usage matrix")
    codon_usage_matrix = calculate_anticodon_usage_matrix(db, genome_cfg)
    codon_usage_matrix.to_csv(config.output_codon_usage, sep="\t", index=True)
    logger.success(f"Wrote codon usage matrix to {config.output_codon_usage}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect DNA-level pombe features (spine stage)")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--genome-landmarks", type=Path, required=True, help="Genome-landmarks YAML (telomere/centromere coordinates)")
    parser.add_argument("--output", type=Path, required=True, dest="output_dna", help="Output DNA-level features pickle")
    parser.add_argument("--codon-usage-output", type=Path, required=True, dest="output_codon_usage", help="Output codon usage matrix path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run DNA-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DnaConfig(
            pombase_dir=args.pombase_dir,
            genome_landmarks=args.genome_landmarks,
            output_dna=args.output_dna,
            output_codon_usage=args.output_codon_usage,
        )
        config.validate()
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
