#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enrichment Finalize (aggregate)
=================================

Final stage of the split GO/FYPO/MONDO enrichment pipeline: reads the three
per-ontology frame pickles from enrich_one_ontology and writes the persisted
concat TSVs + the filtered GO target consumed by the network rule and downstream
notebooks. Reproduces run_ontology_enrichment.py's TSV outputs byte-for-byte.

Input
-----
- _work/{GO,FYPO,MONDO}_frames.pkl (from enrich_one_ontology)

Output
------
- go_enrichment_full.tsv, go_enrichment_slim.tsv, fypo_enrichment_full.tsv,
  mondo_enrichment_full.tsv (per-ontology concat tables)
- go_enrichment_full_filtered.tsv (design-doc target: pop_count<max, no MF)

Usage
-----
    python finalize_enrichment.py \\
        --work-dir results/enrichment/raw/{dataset}/{version}/_work \\
        --output-dir results/enrichment/raw/{dataset}/{version}

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

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.cluster_enrichment import ONTOLOGIES, POP_COUNT_MAX, filter_go_full


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class FinalizeConfig:
    """Work dir (frame pickles), output dir, and the GO pop-count filter."""
    work_dir: Path
    output_dir: Path
    pop_count_max: int = POP_COUNT_MAX

    def validate(self) -> None:
        """Raise ValueError if any per-ontology frame pickle is missing; ensure output dir exists."""
        for onto in ONTOLOGIES:
            p = self.work_dir / f"{onto}_frames.pkl"
            if not p.exists():
                raise ValueError(f"Required input not found: {p}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


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
def run(config: FinalizeConfig) -> None:
    """Read the per-ontology frames and write the concat TSVs + filtered GO target."""
    config.validate()
    frames = {onto: pd.read_pickle(config.work_dir / f"{onto}_frames.pkl") for onto in ONTOLOGIES}
    go, fypo, mondo = frames["GO"], frames["FYPO"], frames["MONDO"]

    # Filtered GO (design-doc deterministic target). Empty-safe.
    filtered = filter_go_full(go["full"], config.pop_count_max)
    filtered.to_csv(config.output_dir / "go_enrichment_full_filtered.tsv", sep="\t", index=False)

    # Persist per-ontology concat tables for the network rule / downstream notebooks.
    for name, df in [
        ("go_enrichment_full.tsv", go["full"]), ("go_enrichment_slim.tsv", go["slim"]),
        ("fypo_enrichment_full.tsv", fypo["full"]), ("mondo_enrichment_full.tsv", mondo["full"]),
    ]:
        df.to_csv(config.output_dir / name, sep="\t", index=False)

    logger.success(f"Finalized enrichment TSVs -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Finalize enrichment: concat TSVs + filtered GO target")
    parser.add_argument("--work-dir", type=Path, required=True, help="Work dir with per-ontology frame pickles")
    parser.add_argument("--output-dir", type=Path, required=True, help="Raw results dir for TSV outputs")
    parser.add_argument("--pop-count-max", type=int, default=POP_COUNT_MAX, help="Max pop_count for filtered GO output (default 400)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, finalize enrichment, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = FinalizeConfig(work_dir=args.work_dir, output_dir=args.output_dir, pop_count_max=args.pop_count_max)
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
