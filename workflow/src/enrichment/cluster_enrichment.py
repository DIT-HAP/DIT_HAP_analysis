"""
Cluster Enrichment — Shared Helpers
=====================================

Shared logic for the split GO/FYPO/MONDO cluster-enrichment pipeline: constants,
the gene-set loader (background / nonWT / per-cluster), the per-ontology
`enrich_all_clusters` driver (load DAG once, enrich every cluster), the workbook
writer, and per-ontology config resolution. Used by the three driver scripts
(prepare_genesets / enrich_one_ontology / finalize_enrichment) so each ontology
can be its own Snakemake job while the enrichment maths lives in one place.

The three ontologies all run the SAME enrich_all_clusters over the SAME gene
sets; they differ only in the OBO/GAF/slim inputs and load kwargs — hence one
`resolve_ontology(name, ...)` config builder rather than three near-identical
code paths (byte-faithful to run_ontology_enrichment.py).

Input
-----
- Curated final_clusters.tsv (via load_cluster_genesets)
- PomBase ontology triples (OBO + GAF/PHAF + slim tables) per ontology

Output
------
- Helpers returning DataFrames; drivers persist them.

Usage
-----
    from workflow.src.enrichment.cluster_enrichment import load_cluster_genesets, resolve_ontology, enrich_all_clusters
    genesets = load_cluster_genesets(final_clusters, cluster_column, wt_cluster)
    onto = resolve_ontology("GO", ontology_dir, intermediate_dir)
    full, slim, nonwt_full, nonwt_slim = enrich_all_clusters(onto.data, genesets, ..., wt_cluster)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from workflow.src.enrichment.ontology import (
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
CLUSTER_COLUMN = "cluster"

GO_LOAD_KWARGS = {"relationships": {"is_a", "part_of"}, "propagate_counts": True, "load_obsolete": False, "prt": None}
SIMPLE_LOAD_KWARGS = {"propagate_counts": True, "load_obsolete": False, "prt": None}

# The three ontologies, in the order run_ontology_enrichment.py processed them.
ONTOLOGIES = ["GO", "FYPO", "MONDO"]


# =============================================================================
# GENE SETS
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ClusterGeneSets:
    """Background / nonWT-background / per-cluster gene sets from final_clusters.tsv."""
    cluster_genes: dict            # {int cluster -> [systematic ids]}, sorted by cluster
    bg_genes: list                 # all background genes
    nonwt_bg_genes: list           # background genes in clusters < wt_cluster


@logger.catch(reraise=True)
def load_cluster_genesets(final_clusters: Path, cluster_column: str, wt_cluster: int) -> ClusterGeneSets:
    """Load per-cluster / background / nonWT gene sets from the curated final_clusters.tsv."""
    input_data = read_file(final_clusters)
    bg_genes = input_data["Systematic ID"].dropna().unique().tolist()
    nonwt_bg_genes = input_data[input_data[cluster_column] < wt_cluster]["Systematic ID"].dropna().unique().tolist()
    cluster_genes = input_data.groupby(cluster_column)["Systematic ID"].apply(list).to_dict()
    cluster_genes = {int(k): v for k, v in sorted(cluster_genes.items())}
    return ClusterGeneSets(cluster_genes=cluster_genes, bg_genes=bg_genes, nonwt_bg_genes=nonwt_bg_genes)


@logger.catch
def write_gene_lists(cluster_genes: dict, bg_genes: list, output_dir: Path) -> None:
    """Write DIT_HAP_all_genes.txt, per-cluster gene lists, and the cluster matrix."""
    (output_dir / "DIT_HAP_all_genes.txt").write_text("\n".join(bg_genes) + "\n")
    for cluster, genes in cluster_genes.items():
        (output_dir / f"DIT_HAP_cluster_{cluster}_genes.txt").write_text("\n".join(genes) + "\n")
    matrix = pd.DataFrame.from_dict(cluster_genes, orient="index").transpose()
    matrix.to_csv(output_dir / "DIT_HAP_cluster_genes_matrix.txt", sep="\t", index=False)


# =============================================================================
# PER-ONTOLOGY CONFIG
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class OntologyPlan:
    """Everything needed to enrich one ontology: loaded data, kwargs, and output labels/names."""
    name: str
    data: object                    # loaded OntologyData
    load_kwargs: dict
    enrichment_kwargs: dict         # full kwargs incl. alpha/methods/prt + ontology extras
    workbook_name: str
    nonwt_workbook_name: str
    full_label: str
    slim_label: str


@logger.catch(reraise=True)
def resolve_ontology(name: str, ontology_dir: Path, intermediate_dir: Path, fdr_threshold: float) -> OntologyPlan:
    """Build the OntologyPlan for GO/FYPO/MONDO, reformatting PHAF/MONDO GAFs into intermediate_dir as needed."""
    od = ontology_dir
    base_kwargs = {"alpha": fdr_threshold, "methods": ["fdr_bh"], "prt": None}

    if name == "GO":
        data = OntologyDataConfig(
            ontology_obo=od / "go-basic.obo",
            ontology_association_gaf=od / "gene_ontology_annotation.gaf.tsv",
            slim_terms_table=[od / "bp_go_slim_terms.tsv", od / "mf_go_slim_terms.tsv", od / "cc_go_slim_terms.tsv"],
        ).load_data()
        return OntologyPlan(
            name=name, data=data, load_kwargs=dict(GO_LOAD_KWARGS),
            enrichment_kwargs={**base_kwargs, "propagate_counts": True, "relationships": {"is_a", "part_of"}},
            workbook_name="gene_ontology_enrichment_results.xlsx",
            nonwt_workbook_name="nonWT_gene_ontology_enrichment_results.xlsx",
            full_label="Full GO Enrichment", slim_label="GO Slim Enrichment",
        )
    if name == "FYPO":
        phaf_gaf = format_phaf_file(
            od / "fypo-simple-pombase.obo", od / "pombase_phenotype_annotation.phaf.tsv", intermediate_dir / "phaf_go_style.tsv"
        )
        data = OntologyDataConfig(
            ontology_obo=od / "fypo-simple-pombase.obo", ontology_association_gaf=phaf_gaf,
            slim_terms_table=[od / "fypo_slim_ids_and_names.tsv"],
        ).load_data()
        return OntologyPlan(
            name=name, data=data, load_kwargs=dict(SIMPLE_LOAD_KWARGS),
            enrichment_kwargs={**base_kwargs, "propagate_counts": True},
            workbook_name="fission_yeast_phenotype_ontology_enrichment_results.xlsx",
            nonwt_workbook_name="nonWT_fission_yeast_phenotype_ontology_enrichment_results.xlsx",
            full_label="Full FYPO Enrichment", slim_label="FYPO Slim Enrichment",
        )
    if name == "MONDO":
        mondo_gaf = format_mondo_gaf_file(
            od / "mondo-simple.obo", od / "human_disease_association.tsv", intermediate_dir / "mondo_go_style.tsv"
        )
        data = OntologyDataConfig(
            ontology_obo=od / "mondo-simple.obo", ontology_association_gaf=mondo_gaf,
            slim_terms_table=[od / "pombe_mondo_disease_slim_terms.tsv"],
        ).load_data()
        return OntologyPlan(
            name=name, data=data, load_kwargs=dict(SIMPLE_LOAD_KWARGS),
            enrichment_kwargs={**base_kwargs, "propagate_counts": True},
            workbook_name="mondo_disease_ontology_enrichment_results.xlsx",
            nonwt_workbook_name="nonWT_mondo_disease_ontology_enrichment_results.xlsx",
            full_label="Full MONDO Enrichment", slim_label="MONDO Slim Enrichment",
        )
    raise ValueError(f"Unknown ontology: {name!r} (expected one of {ONTOLOGIES})")


# =============================================================================
# CORE ENRICHMENT (load DAG once, enrich every cluster)
# =============================================================================
@logger.catch(reraise=True)
def enrich_all_clusters(
    ontology_data,
    cluster_genes: dict,
    bg_genes: list,
    nonwt_bg_genes: list,
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
def filter_go_full(go_full: pd.DataFrame, pop_count_max: int) -> pd.DataFrame:
    """Design-doc deterministic GO target: pop_count < max, drop MF namespace, sort by [Cluster, namespace, term_id]."""
    if go_full.empty:
        return pd.DataFrame()
    return go_full.query(f"pop_count < {pop_count_max} and namespace != 'MF'").sort_values(["Cluster", "namespace", "term_id"])

