#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ontology Enrichment (deterministic)
====================================

Runs GO / FYPO / MONDO over-representation analysis (full + slim, goatools
FDR-BH) for each gene cluster in the curated final_clusters.tsv. This is the
DETERMINISTIC, pure-local half of
DIT_HAP_pipeline/workflow/notebooks/comprehensive_enrichment_analysis.ipynb —
the STRING-db and REVIGO network steps live in the optional network rule
(design doc §5, Phase 2 Task 5).

For efficiency each ontology's DAG + slim association is loaded ONCE and reused
across clusters (the notebook reloaded per cluster); enrichment is deterministic
given the same DAG and gene sets, so results are identical.

Input
-----
- Curated final_clusters.tsv (Systematic ID, Name, revised_cluster, um, lam)
- PomBase ontology triples (OBO + GAF/PHAF + slim tables) for GO / FYPO / MONDO
- Gene metadata + curated deletion-library table (for id->name mapping)

Output
------
- go_enrichment_full_filtered.tsv (design-doc target: pop_count<400, no MF)
- Per-ontology full/slim concat TSVs (+ nonWT variants)
- GO / FYPO / MONDO xlsx workbooks (per-namespace sheets)
- Gene-list txt files per cluster

Usage
-----
    python run_ontology_enrichment.py \\
        --final-clusters resources/curated/final_clusters.tsv \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --output-dir results/enrichment/raw/{dataset} \\
        --intermediate-dir results/enrichment/raw/{dataset}/_gaf

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-16
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
from workflow.src.enrichment.ontology import (
    GeneMetaConfig,
    OntologyDataConfig,
    format_mondo_gaf_file,
    format_phaf_file,
    load_ontology_data,
)
from workflow.src.enrichment.pipeline import (
    format_ontology_enrichment_results,
    get_slim_ns2assoc,
    ontology_enrichment,
)
from workflow.src.io import read_file

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
WT_CLUSTER = 9              # cluster 9 is the WT/background cluster (quirk)
POP_COUNT_MAX = 400         # gene-count filter for the "filtered" GO output
FDR_THRESHOLD = 0.05        # alpha for goatools FDR-BH

GO_LOAD_KWARGS = {"relationships": {"is_a", "part_of"}, "propagate_counts": True, "load_obsolete": False, "prt": None}
SIMPLE_LOAD_KWARGS = {"propagate_counts": True, "load_obsolete": False, "prt": None}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class EnrichmentConfig:
    """Inputs, output dirs, and analysis parameters for cluster enrichment."""
    final_clusters: Path
    pombase_dir: Path
    deletion_library_xlsx: Path
    output_dir: Path
    intermediate_dir: Path
    cluster_column: str = "revised_cluster"
    wt_cluster: int = WT_CLUSTER
    pop_count_max: int = POP_COUNT_MAX
    fdr_threshold: float = FDR_THRESHOLD

    @property
    def ontology_dir(self) -> Path:
        return self.pombase_dir / "ontologies_and_associations"

    def validate(self) -> None:
        """Raise ValueError if any required input is missing."""
        for path in [self.final_clusters, self.pombase_dir, self.deletion_library_xlsx]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


@logger.catch
def write_gene_lists(cluster_genes: dict, bg_genes: list[str], output_dir: Path) -> None:
    """Write DIT_HAP_all_genes.txt, per-cluster gene lists, and the cluster matrix."""
    (output_dir / "DIT_HAP_all_genes.txt").write_text("\n".join(bg_genes) + "\n")
    for cluster, genes in cluster_genes.items():
        (output_dir / f"DIT_HAP_cluster_{cluster}_genes.txt").write_text("\n".join(genes) + "\n")
    matrix = pd.DataFrame.from_dict(cluster_genes, orient="index").transpose()
    matrix.to_csv(output_dir / "DIT_HAP_cluster_genes_matrix.txt", sep="\t", index=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def enrich_all_clusters(
    ontology_data,
    cluster_genes: dict,
    bg_genes: list[str],
    nonwt_bg_genes: list[str],
    load_kwargs: dict,
    enrichment_kwargs: dict,
    format_kwargs: dict,
    wt_cluster: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load one ontology ONCE, then run full+slim enrichment for every cluster (+ nonWT for clusters < wt_cluster).

    Returns concatenated (full, slim, nonWT_full, nonWT_slim) frames with a Cluster column.
    """
    # Load the DAG + associations a single time (the notebook reloaded per cluster).
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(ontology_data, **load_kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    full_kwargs = dict(enrichment_kwargs)
    full_kwargs["godag"] = dag
    full_kwargs["ns2assoc"] = ns2assoc

    slim_kwargs = dict(enrichment_kwargs)
    slim_kwargs["godag"] = slim_dag
    slim_kwargs["ns2assoc"] = ns2slim_assoc["all_ancestors"]
    slim_kwargs["propagate_counts"] = False

    full_by_cluster, slim_by_cluster = {}, {}
    nonwt_full_by_cluster, nonwt_slim_by_cluster = {}, {}

    for cluster, genes in cluster_genes.items():
        logger.info(f"  cluster {cluster}: {len(genes)} genes")
        _, sig_full = ontology_enrichment(genes, bg_genes, **full_kwargs)
        _, sig_slim = ontology_enrichment(genes, bg_genes, **slim_kwargs)
        full_by_cluster[cluster] = format_ontology_enrichment_results("full", sig_full, **format_kwargs)
        slim_by_cluster[cluster] = format_ontology_enrichment_results("slim", sig_slim, **format_kwargs)

        if cluster < wt_cluster:
            _, sig_nonwt_full = ontology_enrichment(genes, nonwt_bg_genes, **full_kwargs)
            _, sig_nonwt_slim = ontology_enrichment(genes, nonwt_bg_genes, **slim_kwargs)
            nonwt_full_by_cluster[cluster] = format_ontology_enrichment_results("full", sig_nonwt_full, **format_kwargs)
            nonwt_slim_by_cluster[cluster] = format_ontology_enrichment_results("slim", sig_nonwt_slim, **format_kwargs)

    return (
        _concat_by_cluster(full_by_cluster),
        _concat_by_cluster(slim_by_cluster),
        _concat_by_cluster(nonwt_full_by_cluster),
        _concat_by_cluster(nonwt_slim_by_cluster),
    )


def _concat_by_cluster(by_cluster: dict) -> pd.DataFrame:
    """Concatenate per-cluster enrichment frames into one, adding a Cluster column (empty-safe)."""
    non_empty = {k: v for k, v in by_cluster.items() if v is not None and not v.empty}
    if not non_empty:
        return pd.DataFrame()
    return pd.concat(non_empty, axis=0).droplevel(1).rename_axis("Cluster").reset_index()


@logger.catch
def write_ontology_workbook(output_path: Path, full_df: pd.DataFrame, slim_df: pd.DataFrame, full_label: str, slim_label: str) -> None:
    """Write a per-ontology xlsx: full + slim sheets, then one sheet per namespace (full and slim)."""
    with pd.ExcelWriter(output_path) as writer:
        (full_df if not full_df.empty else pd.DataFrame({"note": ["no significant terms"]})).to_excel(writer, sheet_name=full_label, index=False)
        (slim_df if not slim_df.empty else pd.DataFrame({"note": ["no significant terms"]})).to_excel(writer, sheet_name=slim_label, index=False)
        if not full_df.empty and "namespace" in full_df.columns:
            for namespace, ns_df in full_df.groupby("namespace"):
                ns_df.to_excel(writer, sheet_name=str(namespace)[:31], index=False)
        if not slim_df.empty and "namespace" in slim_df.columns:
            for namespace, ns_df in slim_df.groupby("namespace"):
                ns_df.to_excel(writer, sheet_name=("Slim " + str(namespace))[:31], index=False)


@logger.catch
def run_enrichment(config: EnrichmentConfig) -> None:
    """Orchestrate gene-list prep, GAF reformatting, GO/FYPO/MONDO enrichment, and output writing."""
    config.validate()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.intermediate_dir.mkdir(parents=True, exist_ok=True)
    od = config.ontology_dir

    # --- Gene metadata (id -> name for readable study/pop items) ---
    gene_meta = GeneMetaConfig(
        gene_IDs_names_products=config.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv",
        deletion_library_essentiality=config.deletion_library_xlsx,
    ).load_data()
    format_kwargs = {"itemid2name": gene_meta.id2name}

    # --- Cluster gene sets ---
    input_data = read_file(config.final_clusters)
    col = config.cluster_column
    bg_genes = input_data["Systematic ID"].dropna().unique().tolist()
    nonwt_bg_genes = input_data[input_data[col] < config.wt_cluster]["Systematic ID"].dropna().unique().tolist()
    cluster_genes = input_data.groupby(col)["Systematic ID"].apply(list).to_dict()
    cluster_genes = {int(k): v for k, v in sorted(cluster_genes.items())}
    write_gene_lists(cluster_genes, bg_genes, config.output_dir)
    logger.info(f"{len(cluster_genes)} clusters, {len(bg_genes)} background genes")

    enrichment_kwargs = {"alpha": config.fdr_threshold, "methods": ["fdr_bh"], "prt": None}

    # --- GO ---
    logger.info("GO enrichment")
    go_cfg = OntologyDataConfig(
        ontology_obo=od / "go-basic.obo",
        ontology_association_gaf=od / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[od / "bp_go_slim_terms.tsv", od / "mf_go_slim_terms.tsv", od / "cc_go_slim_terms.tsv"],
    ).load_data()
    go_full, go_slim, go_nonwt_full, go_nonwt_slim = enrich_all_clusters(
        go_cfg, cluster_genes, bg_genes, nonwt_bg_genes,
        {**GO_LOAD_KWARGS}, {**enrichment_kwargs, "propagate_counts": True, "relationships": {"is_a", "part_of"}}, format_kwargs, config.wt_cluster,
    )
    write_ontology_workbook(config.output_dir / "gene_ontology_enrichment_results.xlsx", go_full, go_slim, "Full GO Enrichment", "GO Slim Enrichment")
    write_ontology_workbook(config.output_dir / "nonWT_gene_ontology_enrichment_results.xlsx", go_nonwt_full, go_nonwt_slim, "Full GO Enrichment", "GO Slim Enrichment")

    # Design-doc deterministic target: pop_count < max, no MF namespace. REVIGO
    # sort columns are added later by the network rule; sort deterministically here.
    if not go_full.empty:
        filtered = go_full.query(f"pop_count < {config.pop_count_max} and namespace != 'MF'").sort_values(["Cluster", "namespace", "term_id"])
        filtered.to_csv(config.output_dir / "go_enrichment_full_filtered.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(config.output_dir / "go_enrichment_full_filtered.tsv", sep="\t", index=False)

    # --- FYPO (reformat PHAF -> GO-style GAF into intermediate dir) ---
    logger.info("FYPO enrichment")
    phaf_gaf = format_phaf_file(od / "fypo-simple-pombase.obo", od / "pombase_phenotype_annotation.phaf.tsv", config.intermediate_dir / "phaf_go_style.tsv")
    fypo_cfg = OntologyDataConfig(ontology_obo=od / "fypo-simple-pombase.obo", ontology_association_gaf=phaf_gaf, slim_terms_table=[od / "fypo_slim_ids_and_names.tsv"]).load_data()
    fypo_full, fypo_slim, fypo_nonwt_full, fypo_nonwt_slim = enrich_all_clusters(
        fypo_cfg, cluster_genes, bg_genes, nonwt_bg_genes, {**SIMPLE_LOAD_KWARGS}, {**enrichment_kwargs, "propagate_counts": True}, format_kwargs, config.wt_cluster,
    )
    write_ontology_workbook(config.output_dir / "fission_yeast_phenotype_ontology_enrichment_results.xlsx", fypo_full, fypo_slim, "Full FYPO Enrichment", "FYPO Slim Enrichment")
    write_ontology_workbook(config.output_dir / "nonWT_fission_yeast_phenotype_ontology_enrichment_results.xlsx", fypo_nonwt_full, fypo_nonwt_slim, "Full FYPO Enrichment", "FYPO Slim Enrichment")

    # --- MONDO (reformat disease association -> GO-style GAF into intermediate dir) ---
    logger.info("MONDO enrichment")
    mondo_gaf = format_mondo_gaf_file(od / "mondo-simple.obo", od / "human_disease_association.tsv", config.intermediate_dir / "mondo_go_style.tsv")
    mondo_cfg = OntologyDataConfig(ontology_obo=od / "mondo-simple.obo", ontology_association_gaf=mondo_gaf, slim_terms_table=[od / "pombe_mondo_disease_slim_terms.tsv"]).load_data()
    mondo_full, mondo_slim, mondo_nonwt_full, mondo_nonwt_slim = enrich_all_clusters(
        mondo_cfg, cluster_genes, bg_genes, nonwt_bg_genes, {**SIMPLE_LOAD_KWARGS}, {**enrichment_kwargs, "propagate_counts": True}, format_kwargs, config.wt_cluster,
    )
    write_ontology_workbook(config.output_dir / "mondo_disease_ontology_enrichment_results.xlsx", mondo_full, mondo_slim, "Full MONDO Enrichment", "MONDO Slim Enrichment")
    write_ontology_workbook(config.output_dir / "nonWT_mondo_disease_ontology_enrichment_results.xlsx", mondo_nonwt_full, mondo_nonwt_slim, "Full MONDO Enrichment", "MONDO Slim Enrichment")

    # Persist per-ontology concat tables for the network rule / downstream notebooks.
    for name, df in [
        ("go_enrichment_full.tsv", go_full), ("go_enrichment_slim.tsv", go_slim),
        ("fypo_enrichment_full.tsv", fypo_full), ("mondo_enrichment_full.tsv", mondo_full),
    ]:
        df.to_csv(config.output_dir / name, sep="\t", index=False)

    logger.success(f"Enrichment complete -> {config.output_dir}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Deterministic GO/FYPO/MONDO cluster enrichment")
    parser.add_argument("--final-clusters", type=Path, required=True, help="Curated final_clusters.tsv")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion library categories xlsx")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for enrichment results")
    parser.add_argument("--intermediate-dir", type=Path, required=True, help="Directory for reformatted GAF intermediates")
    parser.add_argument("--wt-cluster", type=int, default=WT_CLUSTER, help="WT/background cluster id (default 9)")
    parser.add_argument("--pop-count-max", type=int, default=POP_COUNT_MAX, help="Max pop_count for filtered GO output (default 400)")
    parser.add_argument("--fdr-threshold", type=float, default=FDR_THRESHOLD, help="FDR-BH alpha (default 0.05)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run enrichment, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = EnrichmentConfig(
            final_clusters=args.final_clusters,
            pombase_dir=args.pombase_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            output_dir=args.output_dir,
            intermediate_dir=args.intermediate_dir,
            wt_cluster=args.wt_cluster,
            pop_count_max=args.pop_count_max,
            fdr_threshold=args.fdr_threshold,
        )
        run_enrichment(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
