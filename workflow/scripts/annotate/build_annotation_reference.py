#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene Annotation Reference Construction
========================================

Parses PomBase and SGD sources once into a single per-gene annotation table, so
that annotating a user table is later just a join. Combines blocks in order:
1. PomBase gene metadata (systematic_id, name, product, feature_type, characterisation_status, etc.)
2. Deletion library categories (this study's essentiality classification)
3. Gene-level HD_DIT_HAP depletion (DR/DL)
4. gRNA-level depletion (DR/DL)
5. S. cerevisiae ortholog info (id/name/essentiality/description)
6. Human ortholog symbols
7. Functional annotation (GO-slim, complex membership)

All blocks are joined on systematic_id (pombe gene systematic name).

Building this separately is what keeps the annotation CLI fast: loading the GO
OBO/GAF and the 200k-row SGD phenotype table costs far more than the join does.

Input
-----
- A PomBase version directory (curated_orthologs, Gene_metadata, ontologies_and_associations, Protein_features)
- An SGD version directory (SGD_features.tab, phenotype_data.tab) from fetch_sgd_data.sh
- resources/curated/deletion_library_categories.xlsx
- HD_DIT_HAP dataset gene-level fitting results (gene DR/DL)
- resources/curated/*_gRNA_HDdata_fitted_parameters.tsv (gRNA-level DR/DL)

Output
------
- gene_annotation_reference.{gene_type}.parquet: one row per pombe gene, one column per annotation field
  where {gene_type} can be: all, protein, lncRNA, tRNA, rRNA, snoRNA, sncRNA, snRNA, noncoding

Usage
-----
    python build_annotation_reference.py \\
        --pombase-dir resources/external/pombase/2026-06-01 \\
        --sgd-dir resources/external/sgd/2026-08-11 \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --hd-dithap-dataset HD_DIT_HAP \\
        --grna-parameters-tsv resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv \\
        --gene-type protein \\
        --output results/1c_annotation/2026-06-01/2026-08-11/gene_annotation_reference.protein.parquet

Author:   Yusheng Yang (guidance) + Claude Opus 5 (implementation)
Date:     2026-09-03
Version:  2.1.0
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
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from annotation.core import (  # noqa: E402
    VALID_GENE_TYPES,
    add_gene_status_column,
    assemble_annotation_reference_split,
    build_complex_block,
    build_deletion_library_block,
    build_gene_level_depletion_block,
    build_go_slim_block,
    build_grna_block,
    build_hs_ortholog_block,
    build_pombase_metadata_block,
    build_sc_essentiality,
    build_sc_gene_info,
    build_sc_ortholog_block,
    build_verification_block,
    read_ortholog_file,
    read_sgd_feature_data,
    read_sgd_phenotype_data,
    report_gene_status_summary,
)
from enrichment.ontology import OntologyDataConfig, load_ontology_data  # noqa: E402
from enrichment.pipeline import get_slim_ns2assoc  # noqa: E402
from io_table import read_file, write_parquet  # noqa: E402
from logging_setup import setup_logger  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================

@dataclass(kw_only=True, slots=True, frozen=True)
class AnnotationReferenceConfig:
    """Inputs/outputs for annotation-reference construction."""

    pombase_dir: Path
    sgd_dir: Path
    deletion_library_xlsx: Path
    verification_csv: Path
    hd_dithap_dataset: str
    grna_parameters_tsv: Path
    gene_type: str
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        if self.gene_type not in VALID_GENE_TYPES:
            raise ValueError(
                f"Invalid gene_type '{self.gene_type}'. Valid options: {list(VALID_GENE_TYPES.keys())}"
            )
        for path in [
            self.pombase_dir,
            self.sgd_dir,
            self.deletion_library_xlsx,
            self.verification_csv,
            self.grna_parameters_tsv,
            self.sgd_features,
            self.sgd_phenotypes,
            self.gene_ids_parquet,
        ]:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)

    @property
    def sgd_features(self) -> Path:
        """SGD_features.tab (per-ORF standard name / qualifier / description)."""
        return self.sgd_dir / "SGD_features.tab"

    @property
    def sgd_phenotypes(self) -> Path:
        """phenotype_data.tab (per-allele phenotype records)."""
        return self.sgd_dir / "phenotype_data.tab"

    @property
    def gene_ids_parquet(self) -> Path:
        """PomBase gene_ids_and_details.parquet (replaces gene_IDs_names_products.tsv)."""
        return self.pombase_dir / "Gene_metadata" / "gene_ids_and_details.parquet"

    @property
    def orthologs_dir(self) -> Path:
        """PomBase curated_orthologs directory."""
        return self.pombase_dir / "curated_orthologs"

    @property
    def ontologies_dir(self) -> Path:
        """PomBase ontologies_and_associations directory."""
        return self.pombase_dir / "ontologies_and_associations"


# =============================================================================
# CORE LOGIC
# =============================================================================

@logger.catch(reraise=True)
def build_ortholog_blocks(config: AnnotationReferenceConfig) -> list[pd.DataFrame]:
    """Build the S. cerevisiae and human ortholog blocks from PomBase + SGD sources."""
    logger.info("Reading SGD feature table")
    sc_gene_info = build_sc_gene_info(read_sgd_feature_data(config.sgd_features))
    logger.info(f"  {len(sc_gene_info):,} S. cerevisiae ORFs")

    logger.info("Deriving S. cerevisiae null-mutant essentiality")
    sc_essentiality = build_sc_essentiality(read_sgd_phenotype_data(config.sgd_phenotypes))
    calls = sc_essentiality["essentiality"].value_counts().to_dict()
    logger.info(f"  {len(sc_essentiality):,} ORFs with a viability call: {calls}")

    logger.info("Building S. cerevisiae ortholog block")
    cerevisiae = read_ortholog_file(config.orthologs_dir / "pombe_cerevisiae_orthologs.txt")
    sc_block = build_sc_ortholog_block(cerevisiae, sc_gene_info, sc_essentiality)
    with_ortholog = int((sc_block["Sc_ortholog_count"] > 0).sum())
    logger.info(f"  {with_ortholog:,}/{len(sc_block):,} pombe genes have a cerevisiae ortholog")

    logger.info("Building human ortholog block")
    human = read_ortholog_file(config.orthologs_dir / "pombe_human_orthologs.txt")
    hs_block = build_hs_ortholog_block(human)

    return [sc_block, hs_block]


@logger.catch(reraise=True)
def build_functional_blocks(config: AnnotationReferenceConfig) -> list[pd.DataFrame]:
    """Build the GO-slim and complex-membership blocks."""
    logger.info("Loading GO ontology for slim mapping (slowest step)")
    ontology = OntologyDataConfig(
        ontology_obo=config.ontologies_dir / "go-basic.obo",
        ontology_association_gaf=config.ontologies_dir / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            config.ontologies_dir / "bp_go_slim_terms.tsv",
            config.ontologies_dir / "mf_go_slim_terms.tsv",
            config.ontologies_dir / "cc_go_slim_terms.tsv",
        ],
    ).load_data()
    dag, _, ns2assoc, _, _, slim_dag = load_ontology_data(ontology)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)
    slim_term_names = {term: node.name for term, node in slim_dag.items()}

    logger.info("Building GO-slim block")
    go_block = build_go_slim_block(ns2slim_assoc["all_ancestors"], slim_term_names)
    logger.info(f"  {len(go_block):,} genes with GO-slim annotation")

    logger.info("Building complex-membership block")
    complexes = read_file(config.ontologies_dir / "macromolecular_complex_annotation.tsv")
    complex_block = build_complex_block(complexes)
    logger.info(f"  {len(complex_block):,} genes in a macromolecular complex")

    return [go_block, complex_block]


@logger.catch(reraise=True)
def run(config: AnnotationReferenceConfig) -> None:
    """Assemble every annotation block onto the PomBase gene set and write the reference table."""
    # Block 1: PomBase metadata (defines the reference row set)
    pombe_block = build_pombase_metadata_block(config.gene_ids_parquet, config.gene_type)
    logger.info(f"  {len(pombe_block):,} pombe genes define the reference row set")

    # Block 2: Deletion library categories
    deletion_block = build_deletion_library_block(config.deletion_library_xlsx)

    # Block 2b: Verification phenotype
    verification_block = build_verification_block(config.verification_csv)

    # Block 3: Gene-level HD_DIT_HAP depletion (DR/DL)
    gene_depletion_block = build_gene_level_depletion_block(config.hd_dithap_dataset)

    # Block 4: gRNA-level depletion (DR/DL)
    logger.info("Building gRNA-level depletion block")
    grna_block = build_grna_block(read_file(config.grna_parameters_tsv))
    logger.info(f"  {len(grna_block):,} genes with gRNA-level DR/DL")

    # Blocks 5-7: Orthologs and functional annotation (use left join for these)
    ortholog_blocks = build_ortholog_blocks(config)
    functional_blocks = build_functional_blocks(config)

    # Assemble: outer join for experimental data, left join for ortholog/functional annotation
    experimental_blocks = [deletion_block, verification_block, gene_depletion_block, grna_block]
    annotation_blocks = ortholog_blocks + functional_blocks

    reference = assemble_annotation_reference_split(pombe_block, experimental_blocks, annotation_blocks)

    # Fill missing deletion_essentiality with "Not_determined" after assembly (outer join)
    if "deletion_essentiality" in reference.columns:
        reference["deletion_essentiality"] = reference["deletion_essentiality"].fillna("Not_determined")

    # Add gene_status column to track gene presence and type in current PomBase
    reference["gene_status"] = add_gene_status_column(reference.index, config.gene_ids_parquet)

    # Report genes not matching expected type
    report_gene_status_summary(reference, config.gene_type)

    write_parquet(reference, config.output)
    logger.success(
        f"Wrote {len(reference):,} genes x {reference.shape[1]} annotation columns to {config.output}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Build the pombe gene annotation reference table")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument(
        "--sgd-dir", type=Path, required=True, help="SGD version directory (see fetch_sgd_data.sh)"
    )
    parser.add_argument(
        "--deletion-library-xlsx", type=Path, required=True, help="Curated deletion-library categories xlsx"
    )
    parser.add_argument(
        "--verification-csv", type=Path, required=True, help="Curated essentiality verification CSV"
    )
    parser.add_argument(
        "--hd-dithap-dataset",
        type=str,
        required=True,
        help="HD_DIT_HAP dataset name for gene-level DR/DL (e.g., 'HD_DIT_HAP')",
    )
    parser.add_argument(
        "--grna-parameters-tsv",
        type=Path,
        required=True,
        help="Curated gRNA fitted-parameters TSV (supplies gRNA-level DR/DL)",
    )
    parser.add_argument(
        "--gene-type",
        type=str,
        default="all",
        choices=list(VALID_GENE_TYPES.keys()),
        help=f"Gene type filter: {', '.join(VALID_GENE_TYPES.keys())} (default: all)",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output annotation reference parquet")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, assemble the annotation reference, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = AnnotationReferenceConfig(
            pombase_dir=args.pombase_dir,
            sgd_dir=args.sgd_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            verification_csv=args.verification_csv,
            hd_dithap_dataset=args.hd_dithap_dataset,
            grna_parameters_tsv=args.grna_parameters_tsv,
            gene_type=args.gene_type,
            output=args.output,
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
