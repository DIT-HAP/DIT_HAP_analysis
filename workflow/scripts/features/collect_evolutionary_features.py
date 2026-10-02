#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Evolutionary-Level Feature Collection
=======================================

Assembles ortholog/paralog counts, evolutionary rate, and phyloP/divergence
scores per coding gene. Reads the coding-gene set from the DNA-level parquet.

`paralog_count` has two selectable sources (config.features.paralog_source):
- `ensembl`          — the BioMart paralog export, every homology type
- `deletion_library` — the Hayles-2013 `Paralogues` column (curated, within-species)

Input
-----
- A PomBase version directory (curated_orthologs, gene metadata)
- An Ensembl paralog export TSV, and/or the long deletion-library paralog parquet
  (built by the build_deletion_library_paralogs rule)
- Literature tables (Rhind 2011, Grech 2019)
- DNA-level features parquet (for the coding-gene set)

Output
------
- evolutionary_features.parquet: per-gene evolutionary feature table (indexed by gene id)

Usage
-----
    python collect_evolutionary_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --literature-dir resources/literature \\
        --paralog-source deletion_library \\
        --ensembl-paralogs-tsv resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --deletion-library-paralogs results/1b_features/2025-10-01/deletion_library_paralogs.parquet \\
        --dna-features results/1b_features/2025-10-01/_levels/dna_features.parquet \\
        --output results/1b_features/2025-10-01/_levels/evolutionary_features.parquet

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
from loguru import logger

# 3. Local Imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from io_table import read_parquet, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402
from features.assembly import (  # noqa: E402
    PARALOG_SOURCES,
    collect_evolutionary_level_features,
    load_paralogs,
    load_phyloP_and_divergence,
    read_coding_genes,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class EvolutionaryConfig:
    """Inputs/outputs for evolutionary-level feature collection."""
    pombase_dir: Path
    literature_dir: Path
    paralog_source: str
    ensembl_paralogs_tsv: Path
    deletion_library_paralogs: Path
    dna_features: Path
    output_evolutionary: Path

    @property
    def paralog_file(self) -> Path:
        """The file `paralog_source` selects; the other one is never read."""
        return (
            self.ensembl_paralogs_tsv
            if self.paralog_source == "ensembl"
            else self.deletion_library_paralogs
        )

    def validate(self) -> None:
        """Raise ValueError if the source is unknown or any required input is missing, then ensure the output dir exists."""
        if self.paralog_source not in PARALOG_SOURCES:
            raise ValueError(
                f"features.paralog_source must be one of {PARALOG_SOURCES}, "
                f"got {self.paralog_source!r}"
            )
        for path in [self.pombase_dir, self.literature_dir, self.paralog_file, self.dna_features]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_evolutionary.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used for update_sysIDs()."""
        return self.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv"


# =============================================================================
# HELPERS
# =============================================================================
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: EvolutionaryConfig) -> None:
    """Collect evolutionary-level features filtered to the DNA-level coding-gene set."""
    coding_genes = read_coding_genes(config.dna_features)
    phyloP_and_divergence = load_phyloP_and_divergence(config.literature_dir, config.gene_meta_file)
    paralog_pairs = load_paralogs(
        config.paralog_source,
        config.ensembl_paralogs_tsv,
        config.deletion_library_paralogs,
        config.gene_meta_file,
    )

    logger.info(f"Collecting evolutionary-level features (paralogs from {config.paralog_source})")
    evolutionary_df = collect_evolutionary_level_features(
        config.pombase_dir, config.literature_dir, paralog_pairs,
        config.gene_meta_file, coding_genes, phyloP_and_divergence,
    )
    write_parquet(evolutionary_df, config.output_evolutionary)
    logger.success(f"Wrote {len(evolutionary_df)} gene rows to {config.output_evolutionary}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect evolutionary-level pombe features")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--paralog-source", required=True, choices=PARALOG_SOURCES, help="Which paralog source feeds paralog_count (config.features.paralog_source)")
    parser.add_argument("--ensembl-paralogs-tsv", type=Path, required=True, help="Ensembl paralog export table (read only when --paralog-source=ensembl)")
    parser.add_argument("--deletion-library-paralogs", type=Path, required=True, help="Long (gene, paralog) parquet from build_deletion_library_paralogs (read only when --paralog-source=deletion_library)")
    parser.add_argument("--dna-features", type=Path, required=True, help="DNA-level features parquet (for coding-gene set)")
    parser.add_argument("--output", type=Path, required=True, dest="output_evolutionary", help="Output evolutionary-level features parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run evolutionary-level collection, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = EvolutionaryConfig(
            pombase_dir=args.pombase_dir, literature_dir=args.literature_dir,
            paralog_source=args.paralog_source,
            ensembl_paralogs_tsv=args.ensembl_paralogs_tsv,
            deletion_library_paralogs=args.deletion_library_paralogs,
            dna_features=args.dna_features, output_evolutionary=args.output_evolutionary,
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
