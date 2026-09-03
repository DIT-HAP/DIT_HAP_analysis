#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene Annotation Reference Construction
========================================

Parses PomBase and SGD sources once into a single per-gene annotation table, so
that annotating a user table is later just a join. Combines blocks in order:
1. PomBase gene metadata (systematic_id, name, product, characterisation_status, etc.)
2. Deletion library categories (this study's essentiality classification)
3. Gene-level HD_DIT_HAP depletion (DR/DL)
4. gRNA-level depletion (DR/DL)
5. S. cerevisiae ortholog info (id/name/qualifier/essentiality)
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
- gene_annotation_reference.parquet: one row per pombe gene, one column per annotation field

Usage
-----
    python build_annotation_reference.py \\
        --pombase-dir resources/external/pombase/2026-06-01 \\
        --sgd-dir resources/external/sgd/2026-08-11 \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --hd-dithap-dataset HD_DIT_HAP \\
        --grna-parameters-tsv resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv \\
        --output results/annotation/2026-06-01/2026-08-11/gene_annotation_reference.parquet

Author:   Yusheng Yang (guidance) + Claude Opus 5 (implementation)
Date:     2026-09-03
Version:  2.0.0
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
    assemble_annotation_reference,
    build_complex_block,
    build_go_slim_block,
    build_grna_block,
    build_hs_ortholog_block,
    build_sc_essentiality,
    build_sc_gene_info,
    build_sc_ortholog_block,
    read_ortholog_file,
    read_sgd_feature_data,
    read_sgd_phenotype_data,
)
from data_config import load_dataset_config  # noqa: E402
from enrichment.ontology import OntologyDataConfig, load_ontology_data  # noqa: E402
from enrichment.pipeline import get_slim_ns2assoc  # noqa: E402
from io_table import read_file, read_parquet, write_parquet  # noqa: E402
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
    hd_dithap_dataset: str
    grna_parameters_tsv: Path
    output: Path

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure the output dir exists."""
        for path in [
            self.pombase_dir,
            self.sgd_dir,
            self.deletion_library_xlsx,
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

    @property
    def domains_file(self) -> Path:
        """PomBase protein_families_and_domains.tsv."""
        return self.pombase_dir / "Protein_features" / "protein_families_and_domains.tsv"


# =============================================================================
# BLOCK BUILDERS (NEW STRUCTURE)
# =============================================================================
@logger.catch(reraise=True)
def build_pombase_metadata_block(gene_ids_parquet: Path) -> pd.DataFrame:
    """Build block 1: PomBase gene metadata from gene_ids_and_details.parquet.

    Extracts: systematic_id, name, gene_product, product, characterisation_status,
    taxonomic_distribution, deletion_viability (renamed to FYPOviability).
    """
    logger.info("Reading PomBase gene metadata from parquet")
    gene_meta = read_parquet(gene_ids_parquet)

    # Select and rename columns
    columns_map = {
        "systematic_id": "systematic_id",
        "name": "gene_name",
        "gene_product": "gene_product",
        "product": "product",
        "characterisation_status": "characterisation_status",
        "taxonomic_distribution": "taxonomic_distribution",
        "deletion_viability": "FYPOviability",
    }

    available_columns = [col for col in columns_map.keys() if col in gene_meta.columns]
    block = gene_meta[available_columns].copy()
    block = block.rename(columns={k: v for k, v in columns_map.items() if k in available_columns})

    # Set systematic_id as index
    if "systematic_id" in block.columns:
        block = block.set_index("systematic_id")

    # Fill missing gene names with systematic_id
    if "gene_name" in block.columns:
        block["gene_name"] = block["gene_name"].fillna(pd.Series(block.index, index=block.index))

    logger.info(f"  {len(block):,} pombe genes in metadata block")
    return block


@logger.catch(reraise=True)
def build_deletion_library_block(deletion_library_xlsx: Path) -> pd.DataFrame:
    """Build block 2: Deletion library categories.

    Extracts: Systematic ID, Gene dispensability. This study (renamed to deletion_essentiality),
    Category, Sub_category.
    """
    logger.info("Reading deletion library categories")
    deletion_lib = read_file(deletion_library_xlsx)

    # Select and rename columns
    columns_map = {
        "Systematic ID": "systematic_id",
        "Gene dispensability. This study": "deletion_essentiality",
        "Category": "Category",
        "Sub_category": "Sub_category",
    }

    available_columns = [col for col in columns_map.keys() if col in deletion_lib.columns]
    block = deletion_lib[available_columns].copy()
    block = block.rename(columns={k: v for k, v in columns_map.items() if k in available_columns})

    # Set systematic_id as index
    if "systematic_id" in block.columns:
        block = block.set_index("systematic_id")

    logger.info(f"  {len(block):,} genes in deletion library")
    return block


@logger.catch(reraise=True)
def build_gene_level_depletion_block(dataset_name: str) -> pd.DataFrame:
    """Build block 3: Gene-level DR/DL from HD_DIT_HAP dataset.

    Extracts DR and DL from the gene-level fitting results.
    """
    logger.info(f"Reading gene-level depletion from {dataset_name}")
    dataset_config = load_dataset_config(dataset_name)

    if dataset_config.gene_level is None:
        raise ValueError(f"Dataset {dataset_name} has no gene-level data (has_time_points=False)")

    # Read fitting results which contain DR/DL
    fitting_results = read_file(dataset_config.gene_level.fitting_results)

    # Extract gene_systematic_id, DR, DL
    required_cols = ["gene_systematic_id", "DR", "DL"]
    missing = [col for col in required_cols if col not in fitting_results.columns]
    if missing:
        raise ValueError(f"Gene-level fitting results missing columns: {missing}")

    block = fitting_results[required_cols].copy()
    block = block.set_index("gene_systematic_id")

    logger.info(f"  {len(block):,} genes with gene-level DR/DL")
    return block


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
    pombe_block = build_pombase_metadata_block(config.gene_ids_parquet)
    logger.info(f"  {len(pombe_block):,} pombe genes define the reference row set")

    # Block 2: Deletion library categories
    deletion_block = build_deletion_library_block(config.deletion_library_xlsx)

    # Block 3: Gene-level HD_DIT_HAP depletion (DR/DL)
    gene_depletion_block = build_gene_level_depletion_block(config.hd_dithap_dataset)

    # Block 4: gRNA-level depletion (DR/DL)
    logger.info("Building gRNA-level depletion block")
    grna_block = build_grna_block(read_file(config.grna_parameters_tsv))
    logger.info(f"  {len(grna_block):,} genes with gRNA-level DR/DL")

    # Blocks 5-7: Orthologs and functional annotation
    ortholog_blocks = build_ortholog_blocks(config)
    functional_blocks = build_functional_blocks(config)

    # Assemble in order: deletion_lib, gene_depletion, grna, orthologs, functional
    blocks = [deletion_block, gene_depletion_block, grna_block] + ortholog_blocks + functional_blocks
    reference = assemble_annotation_reference(pombe_block, blocks)

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
            hd_dithap_dataset=args.hd_dithap_dataset,
            grna_parameters_tsv=args.grna_parameters_tsv,
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
