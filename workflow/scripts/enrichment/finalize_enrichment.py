#!/usr/bin/env python3

"""
Enrichment Finalize (aggregate)
===============================

Final step of the deterministic enrichment chain: read the three per-ontology
frame pickles written by the ontology jobs and persist the concat tables plus
the filtered GO target the plotting scripts read.

Input
-----
- _work/{GO,FYPO,MONDO}_frames.pkl (from enrich_one_ontology)

Output
------
- go_enrichment_full.tsv, go_enrichment_slim.tsv, fypo_enrichment_full.tsv,
  mondo_enrichment_full.tsv (per-ontology concat tables)
- go_enrichment_full_filtered.tsv (pop_count below the cap, MF dropped)

Usage
-----
    python finalize_enrichment.py \\
        --work-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01/_work \\
        --output-dir results/3c_enrichment/HD_DIT_HAP/kmeans_direct9/2026-06-01 \\
        --pop-count-max 400

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-05
Version:  2.0.0
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
# Standard library
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

from loguru import logger

# Third-party
import pandas as pd

# Local project (resolved relative to this script, never the current directory)
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from enrichment.cluster_enrichment import ONTOLOGIES, filter_go_full  # noqa: E402
from logging_setup import setup_logger  # noqa: E402

# =============================================================================
# 2. CONFIGURATION
# =============================================================================
# Output file name -> (ontology, frame key); the filtered GO table is written separately.
_CONCAT_TABLES = {
    "go_enrichment_full.tsv": ("GO", "full"),
    "go_enrichment_slim.tsv": ("GO", "slim"),
    "fypo_enrichment_full.tsv": ("FYPO", "full"),
    "mondo_enrichment_full.tsv": ("MONDO", "full"),
}


@dataclass(kw_only=True, frozen=True)
class FinalizeConfig:
    """Work dir (frame pickles), output dir, and the GO pop-count cap."""

    work_dir: Path
    output_dir: Path
    pop_count_max: int

    def validate(self) -> None:
        """Raise if a per-ontology frame pickle is missing, then create the output dir."""
        for ontology in ONTOLOGIES:
            path = self.work_dir / f"{ontology}_frames.pkl"
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 3. CORE LOGIC
# =============================================================================
def run(config: FinalizeConfig) -> None:
    """Read the per-ontology frames and write the concat tables plus the filtered GO target."""
    frames = {
        ontology: pd.read_pickle(config.work_dir / f"{ontology}_frames.pkl")
        for ontology in ONTOLOGIES
    }

    filtered = filter_go_full(frames["GO"]["full"], config.pop_count_max)
    filtered.to_csv(config.output_dir / "go_enrichment_full_filtered.tsv", sep="\t", index=False)

    for file_name, (ontology, frame_key) in _CONCAT_TABLES.items():
        frames[ontology][frame_key].to_csv(config.output_dir / file_name, sep="\t", index=False)

    logger.success(f"Finalized enrichment TSVs -> {config.output_dir}")


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Finalize enrichment: concat TSVs + filtered GO target"
    )
    parser.add_argument(
        "--work-dir", type=Path, required=True, help="Work dir with the frame pickles"
    )
    parser.add_argument("--output-dir", type=Path, required=True, help="Results dir for the TSVs")
    parser.add_argument(
        "--pop-count-max", type=int, required=True, help="Max pop_count for the filtered GO table"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    return parser.parse_args()


def main() -> int:
    """Validate the boundary once, then write the finalized tables."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = FinalizeConfig(
        work_dir=args.work_dir,
        output_dir=args.output_dir,
        pop_count_max=args.pop_count_max,
    )

    try:
        config.validate()
        run(config)
    except Exception as error:
        logger.exception(f"Enrichment finalize failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
