"""
Cluster Enrichment Core
=======================

Deterministic GO/FYPO/MONDO over-representation for gene clusters: gene-set
loading, the per-ontology plan, the per-cluster enrichment driver, the workbook
writer, and the goatools formatting helpers that used to live in ``pipeline.py``.

This module is the deterministic half of the enrichment stage. The
non-deterministic half (STRING-db and REVIGO, both external APIs) lives in
``network.py``. Splitting on determinism is what lets the over-representation
rules stay pure-local while the network rules stay optional.

Every ontology runs the SAME ``enrich_all_clusters`` over the SAME gene sets;
they differ only in the OBO/GAF/slim inputs and load kwargs, hence one
``resolve_ontology(name, ...)`` config builder rather than three near-identical
code paths (byte-faithful to the original notebook-derived driver).

Input
-----
- Curated ``final_clusters.tsv`` (via ``load_cluster_genesets``)
- PomBase ontology triples (OBO + GAF/PHAF + slim tables) per ontology

Output
------
- DataFrames and file paths returned to the driver scripts; nothing is written
  here except the GO-style GAF reformats and the workbooks the drivers ask for.

Usage
-----
    from enrichment.cluster_enrichment import load_cluster_genesets, resolve_ontology

    genesets = load_cluster_genesets(final_clusters, CLUSTER_COLUMN, wt_cluster)
    plan = resolve_ontology("GO", ontology_dir, work_dir, fdr_threshold)
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

# Third-party
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from loguru import logger

# Data processing
import pandas as pd

# Local project
from enrichment.ontology import (
    OntologyData,
    OntologyDataConfig,
    format_mondo_gaf_file,
    format_phaf_file,
    load_ontology_data,
)
from enrichment.schema import ClusterCol, EnrichmentCol

# =============================================================================
# CONSTANTS
# =============================================================================
CLUSTER_COLUMN = ClusterCol.CLUSTER

# The three ontologies, in the order the original driver processed them.
ONTOLOGIES = ["GO", "FYPO", "MONDO"]

# goatools' record field names, then the tidy names the enrichment tables carry.
# An ontology with no significant terms still returns a frame with these columns
# (zero rows), so every downstream groupby/concat sees the same schema.
GOATOOLS_COLUMNS = [
    "NS",
    "GO",
    "name",
    "level",
    "depth",
    "p_uncorrected",
    "p_fdr_bh",
    "gene_ratio",
    "term_coverage",
    "study_count",
    "study_n",
    "pop_count",
    "pop_n",
    "ratio_in_study",
    "ratio_in_pop",
    "study_items",
    "pop_items",
    "defn",
]
RENAMED_COLUMNS = {"NS": "namespace", "GO": "term_id", "name": "term", "p_fdr_bh": "p_fdr"}
FORMATTED_COLUMNS = [RENAMED_COLUMNS.get(column, column) for column in GOATOOLS_COLUMNS]

GO_LOAD_KWARGS = {
    "relationships": {"is_a", "part_of"},
    "propagate_counts": True,
    "load_obsolete": False,
    "prt": None,
}
SIMPLE_LOAD_KWARGS = {"propagate_counts": True, "load_obsolete": False, "prt": None}


class OntologyName(StrEnum):
    """The ontologies the enrichment stage supports."""

    GO = "GO"
    FYPO = "FYPO"
    MONDO = "MONDO"


# =============================================================================
# GENE SETS
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ClusterGeneSets:
    """Background / nonWT-background / per-cluster gene sets from final_clusters.tsv."""

    cluster_genes: dict[int, list[str]]
    bg_genes: list[str]
    nonwt_bg_genes: list[str]


def load_cluster_genesets(
    final_clusters: pd.DataFrame, cluster_column: str, wt_cluster: int
) -> ClusterGeneSets:
    """Split a validated final_clusters frame into per-cluster, background, and nonWT gene sets."""
    input_data = final_clusters
    bg_genes = input_data[ClusterCol.GENE_ID].dropna().unique().tolist()
    is_nonwt = input_data[cluster_column] < wt_cluster
    nonwt_bg_genes = input_data.loc[is_nonwt, ClusterCol.GENE_ID].dropna().unique().tolist()
    cluster_genes = input_data.groupby(cluster_column)[ClusterCol.GENE_ID].apply(list).to_dict()
    sorted_clusters = {int(key): value for key, value in sorted(cluster_genes.items())}
    return ClusterGeneSets(
        cluster_genes=sorted_clusters, bg_genes=bg_genes, nonwt_bg_genes=nonwt_bg_genes
    )


def write_gene_lists(
    cluster_genes: dict[int, list[str]], bg_genes: list[str], output_dir: Path
) -> None:
    """Write the all-genes list, one list per cluster, and the cluster matrix."""
    (output_dir / "DIT_HAP_all_genes.txt").write_text("\n".join(bg_genes) + "\n")
    for cluster, genes in cluster_genes.items():
        (output_dir / f"DIT_HAP_cluster_{cluster}_genes.txt").write_text("\n".join(genes) + "\n")
    matrix = pd.DataFrame.from_dict(cluster_genes, orient="index").transpose()
    matrix.to_csv(output_dir / "DIT_HAP_cluster_genes_matrix.txt", sep="\t", index=False)


# =============================================================================
# SLIM ASSOCIATION MAPPING
# =============================================================================
def mapslim(term: str, dag: GODag, slim_dag: dict) -> tuple[set[str], set[str]]:
    """Map a term to its (direct, all) slim ancestors by walking every path to the root."""
    all_ancestors: set[str] = set()
    covered_ancestors: set[str] = set()

    for path in dag.paths_to_top(term):
        # Walk bottom -> top (term to root), so reverse the top-down path first.
        path.reverse()
        got_leaf = False
        for node in path:
            if node.id in slim_dag:
                all_ancestors.add(node.id)
                if got_leaf:
                    covered_ancestors.add(node.id)
                got_leaf = True

    direct_ancestors = all_ancestors - covered_ancestors
    return direct_ancestors, all_ancestors


def get_slim_ns2assoc(ns2assoc: dict, dag: GODag, slim_dag: dict) -> dict:
    """Build namespace -> gene -> slim-term associations for direct and all ancestors."""
    term2slim = {}
    for term in dag:
        direct, all_ancestors = mapslim(term, dag, slim_dag)
        term2slim[term] = {"direct_ancestors": direct, "all_ancestors": all_ancestors}

    ns2slim_assoc: dict[str, dict] = {"direct_ancestors": {}, "all_ancestors": {}}
    for kind in ("direct_ancestors", "all_ancestors"):
        for namespace, gene2terms in ns2assoc.items():
            ns2slim_assoc[kind].setdefault(namespace, {})
            for gene, terms in gene2terms.items():
                ancestors: set[str] = set()
                for term in terms:
                    ancestors |= term2slim[term][kind]
                ns2slim_assoc[kind][namespace][gene] = ancestors

    return ns2slim_assoc


# =============================================================================
# PER-ONTOLOGY CONFIG
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class OntologyPlan:
    """Everything needed to enrich one ontology: loaded data, kwargs, and output labels."""

    name: str
    data: OntologyData
    load_kwargs: dict
    enrichment_kwargs: dict
    workbook_name: str
    nonwt_workbook_name: str
    full_label: str
    slim_label: str


def _go_plan(ontology_dir: Path, fdr_threshold: float) -> OntologyPlan:
    """Build the GO plan: relationships-aware load plus the three GO slim tables."""
    data = OntologyDataConfig(
        ontology_obo=ontology_dir / "go-basic.obo",
        ontology_association_gaf=ontology_dir / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            ontology_dir / "bp_go_slim_terms.tsv",
            ontology_dir / "mf_go_slim_terms.tsv",
            ontology_dir / "cc_go_slim_terms.tsv",
        ],
    ).load_data()
    return OntologyPlan(
        name=OntologyName.GO,
        data=data,
        load_kwargs=dict(GO_LOAD_KWARGS),
        enrichment_kwargs={
            "alpha": fdr_threshold,
            "methods": ["fdr_bh"],
            "prt": None,
            "propagate_counts": True,
            "relationships": {"is_a", "part_of"},
        },
        workbook_name="gene_ontology_enrichment_results.xlsx",
        nonwt_workbook_name="nonWT_gene_ontology_enrichment_results.xlsx",
        full_label="Full GO Enrichment",
        slim_label="GO Slim Enrichment",
    )


def _fypo_plan(ontology_dir: Path, intermediate_dir: Path, fdr_threshold: float) -> OntologyPlan:
    """Build the FYPO plan, reformatting the PomBase PHAF into a GO-style GAF first."""
    phaf_gaf = format_phaf_file(
        ontology_dir / "fypo-simple-pombase.obo",
        ontology_dir / "pombase_single_locus_haploid_phenotype_annotation.phaf.tsv",
        intermediate_dir / "phaf_go_style.tsv",
    )
    data = OntologyDataConfig(
        ontology_obo=ontology_dir / "fypo-simple-pombase.obo",
        ontology_association_gaf=phaf_gaf,
        slim_terms_table=[ontology_dir / "fypo_slim_ids_and_names.tsv"],
    ).load_data()
    return OntologyPlan(
        name=OntologyName.FYPO,
        data=data,
        load_kwargs=dict(SIMPLE_LOAD_KWARGS),
        enrichment_kwargs={
            "alpha": fdr_threshold,
            "methods": ["fdr_bh"],
            "prt": None,
            "propagate_counts": True,
        },
        workbook_name="fission_yeast_phenotype_ontology_enrichment_results.xlsx",
        nonwt_workbook_name="nonWT_fission_yeast_phenotype_ontology_enrichment_results.xlsx",
        full_label="Full FYPO Enrichment",
        slim_label="FYPO Slim Enrichment",
    )


def _mondo_plan(ontology_dir: Path, intermediate_dir: Path, fdr_threshold: float) -> OntologyPlan:
    """Build the MONDO plan, reformatting the disease association into a GO-style GAF first."""
    mondo_gaf = format_mondo_gaf_file(
        ontology_dir / "mondo-simple.obo",
        ontology_dir / "human_disease_association.tsv",
        intermediate_dir / "mondo_go_style.tsv",
    )
    data = OntologyDataConfig(
        ontology_obo=ontology_dir / "mondo-simple.obo",
        ontology_association_gaf=mondo_gaf,
        slim_terms_table=[ontology_dir / "pombe_mondo_disease_slim_terms.tsv"],
    ).load_data()
    return OntologyPlan(
        name=OntologyName.MONDO,
        data=data,
        load_kwargs=dict(SIMPLE_LOAD_KWARGS),
        enrichment_kwargs={
            "alpha": fdr_threshold,
            "methods": ["fdr_bh"],
            "prt": None,
            "propagate_counts": True,
        },
        workbook_name="mondo_disease_ontology_enrichment_results.xlsx",
        nonwt_workbook_name="nonWT_mondo_disease_ontology_enrichment_results.xlsx",
        full_label="Full MONDO Enrichment",
        slim_label="MONDO Slim Enrichment",
    )


def resolve_ontology(
    name: str, ontology_dir: Path, intermediate_dir: Path, fdr_threshold: float
) -> OntologyPlan:
    """Build the GO/FYPO/MONDO plan, writing GAF reformats into ``intermediate_dir``."""
    if name not in ONTOLOGIES:
        raise ValueError(f"Unknown ontology: {name!r} (expected one of {ONTOLOGIES})")
    if name == OntologyName.GO:
        return _go_plan(ontology_dir, fdr_threshold)
    if name == OntologyName.FYPO:
        return _fypo_plan(ontology_dir, intermediate_dir, fdr_threshold)
    return _mondo_plan(ontology_dir, intermediate_dir, fdr_threshold)


# =============================================================================
# ENRICHMENT
# =============================================================================
def create_enrichment_dataframe(
    oea_results_sig: list[GOEnrichmentRecord], itemid2name: dict[str, str] | None = None
) -> pd.DataFrame:
    """Build an enrichment frame from GOEnrichmentRecords, field by field.

    This replaces goatools' ``get_goea_nts_prt`` rather than wrapping it: that
    helper derives its field list from whichever record happens to come first,
    so the output gains or loses a ``defn`` column depending on the process's
    hash seed, and it crashes outright (list + ``dict_keys`` in its own error
    path) for a term with no ``defn`` — MONDO has plenty. Reading the fields
    directly makes the schema fixed, which is why the record parameter is
    ``Any``: goatools exposes the field set dynamically.
    """
    return pd.DataFrame([_record_to_row(record, itemid2name) for record in oea_results_sig])


def _record_to_row(record: Any, itemid2name: dict[str, str] | None) -> dict[str, Any]:
    """Flatten one GOEnrichmentRecord into a table row."""
    study_items = sorted(itemid2name[item] for item in record.study_items) if itemid2name else None
    pop_items = sorted(itemid2name[item] for item in record.pop_items) if itemid2name else None
    return {
        "GO": record.GO,
        "NS": record.NS,
        "name": record.name,
        "level": record.goterm.level,
        "depth": record.goterm.depth,
        "p_uncorrected": record.p_uncorrected,
        "p_fdr_bh": record.p_fdr_bh,
        "study_count": record.study_count,
        "study_n": record.study_n,
        "pop_count": record.pop_count,
        "pop_n": record.pop_n,
        "ratio_in_study": "/".join(map(str, record.ratio_in_study)),
        "ratio_in_pop": "/".join(map(str, record.ratio_in_pop)),
        "study_items": ", ".join(study_items)
        if study_items
        else ", ".join(sorted(record.study_items)),
        "pop_items": ", ".join(pop_items) if pop_items else ", ".join(sorted(record.pop_items)),
        "defn": getattr(record.goterm, "defn", ""),
    }


def ontology_enrichment(
    query_genes: list[str], bg_genes: list[str], **kwargs: Any
) -> tuple[GOEnrichmentStudyNS, list[GOEnrichmentRecord]]:
    """Run the study and keep only significantly ENRICHED terms (FDR-BH below alpha)."""
    study = GOEnrichmentStudyNS(bg_genes, **kwargs)
    results = study.run_study(query_genes, **kwargs)
    alpha = kwargs["alpha"]
    significant = [r for r in results if (r.p_fdr_bh < alpha) and (r.enrichment == "e")]
    return study, significant


def format_ontology_enrichment_results(
    oea_results_sig: list[GOEnrichmentRecord], **kwargs: Any
) -> pd.DataFrame:
    """Format significant enrichment records into a tidy frame with gene_ratio and term_coverage."""
    formatted = create_enrichment_dataframe(oea_results_sig, kwargs.get("itemid2name"))
    if formatted.empty:
        return pd.DataFrame(columns=FORMATTED_COLUMNS)

    formatted["gene_ratio"] = round(formatted["study_count"] / formatted["study_n"], 2)
    formatted["term_coverage"] = round(formatted["study_count"] / formatted["pop_count"], 2)
    formatted["study_items"] = (
        formatted["study_items"].astype(str).str.split(", ").map(sorted).str.join(",")
    )
    formatted["pop_items"] = (
        formatted["pop_items"].astype(str).str.split(", ").map(sorted).str.join(",")
    )

    present_columns = [column for column in GOATOOLS_COLUMNS if column in formatted.columns]
    return formatted[present_columns].copy().rename(columns=RENAMED_COLUMNS)


def enrich_all_clusters(
    ontology_data: OntologyData,
    cluster_genes: dict[int, list[str]],
    bg_genes: list[str],
    nonwt_bg_genes: list[str],
    load_kwargs: dict,
    enrichment_kwargs: dict,
    format_kwargs: dict,
    wt_cluster: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run full + slim enrichment for every cluster (+ nonWT for clusters below ``wt_cluster``).

    The DAG and associations are loaded once and reused; the original notebook
    reloaded them per cluster. Returns the four concatenated frames.
    """
    dag, _, ns2assoc, _, _, slim_dag = load_ontology_data(ontology_data, **load_kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    full_kwargs = {**enrichment_kwargs, "godag": dag, "ns2assoc": ns2assoc}
    slim_kwargs = {
        **enrichment_kwargs,
        "godag": slim_dag,
        "ns2assoc": ns2slim_assoc["all_ancestors"],
        "propagate_counts": False,
    }

    full_by_cluster: dict[int, pd.DataFrame] = {}
    slim_by_cluster: dict[int, pd.DataFrame] = {}
    nonwt_full_by_cluster: dict[int, pd.DataFrame] = {}
    nonwt_slim_by_cluster: dict[int, pd.DataFrame] = {}

    for cluster, genes in cluster_genes.items():
        logger.info(f"  cluster {cluster}: {len(genes)} genes")
        _, sig_full = ontology_enrichment(genes, bg_genes, **full_kwargs)
        _, sig_slim = ontology_enrichment(genes, bg_genes, **slim_kwargs)
        full_by_cluster[cluster] = format_ontology_enrichment_results(sig_full, **format_kwargs)
        slim_by_cluster[cluster] = format_ontology_enrichment_results(sig_slim, **format_kwargs)

        if cluster < wt_cluster:
            _, sig_nonwt_full = ontology_enrichment(genes, nonwt_bg_genes, **full_kwargs)
            _, sig_nonwt_slim = ontology_enrichment(genes, nonwt_bg_genes, **slim_kwargs)
            nonwt_full_by_cluster[cluster] = format_ontology_enrichment_results(
                sig_nonwt_full, **format_kwargs
            )
            nonwt_slim_by_cluster[cluster] = format_ontology_enrichment_results(
                sig_nonwt_slim, **format_kwargs
            )

    return (
        _concat_by_cluster(full_by_cluster),
        _concat_by_cluster(slim_by_cluster),
        _concat_by_cluster(nonwt_full_by_cluster),
        _concat_by_cluster(nonwt_slim_by_cluster),
    )


def ontology_enrichment_pipeline(
    ontology_data: OntologyData,
    query_genes: list[str],
    bg_genes: list[str],
    load_kwargs: dict | None = None,
    enrichment_kwargs: dict | None = None,
    format_kwargs: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Full + slim enrichment for one gene set, loading the ontology from scratch."""
    load_kwargs = dict(load_kwargs or {})
    enrichment_kwargs = dict(enrichment_kwargs or {})
    format_kwargs = dict(format_kwargs or {})

    dag, _, ns2assoc, _, _, slim_dag = load_ontology_data(ontology_data, **load_kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    full_kwargs = {**enrichment_kwargs, "godag": dag, "ns2assoc": ns2assoc}
    slim_kwargs = {
        **enrichment_kwargs,
        "godag": slim_dag,
        "ns2assoc": ns2slim_assoc["all_ancestors"],
        "propagate_counts": False,
    }

    _, sig_full = ontology_enrichment(query_genes, bg_genes, **full_kwargs)
    _, sig_slim = ontology_enrichment(query_genes, bg_genes, **slim_kwargs)
    full_df = format_ontology_enrichment_results(sig_full, **format_kwargs)
    slim_df = format_ontology_enrichment_results(sig_slim, **format_kwargs)
    return full_df, slim_df


def _concat_by_cluster(by_cluster: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Concatenate per-cluster enrichment frames into one, adding a Cluster column."""
    non_empty = {
        cluster: frame
        for cluster, frame in by_cluster.items()
        if frame is not None and not frame.empty
    }
    if not non_empty:
        return pd.DataFrame(columns=[EnrichmentCol.CLUSTER, *FORMATTED_COLUMNS])
    return pd.concat(non_empty, axis=0).droplevel(1).rename_axis("Cluster").reset_index()


def write_ontology_workbook(
    output_path: Path,
    full_df: pd.DataFrame,
    slim_df: pd.DataFrame,
    full_label: str,
    slim_label: str,
) -> None:
    """Write a per-ontology workbook: full + slim sheets, then one sheet per namespace."""
    with pd.ExcelWriter(output_path) as writer:
        _write_sheet(writer, full_df, full_label)
        _write_sheet(writer, slim_df, slim_label)
        for namespace, namespace_df in full_df.groupby(EnrichmentCol.NAMESPACE):
            namespace_df.to_excel(writer, sheet_name=str(namespace)[:31], index=False)
        for namespace, namespace_df in slim_df.groupby(EnrichmentCol.NAMESPACE):
            namespace_df.to_excel(writer, sheet_name=("Slim " + str(namespace))[:31], index=False)


def _write_sheet(writer: pd.ExcelWriter, frame: pd.DataFrame, sheet_name: str) -> None:
    """Write one frame, substituting a one-row note when the frame is empty."""
    if frame.empty:
        pd.DataFrame({"note": ["no significant terms"]}).to_excel(
            writer, sheet_name=sheet_name, index=False
        )
        return
    frame.to_excel(writer, sheet_name=sheet_name, index=False)


def filter_go_full(go_full: pd.DataFrame, pop_count_max: int) -> pd.DataFrame:
    """Design-doc GO target: pop_count below the cap, MF dropped, deterministically sorted."""
    if go_full.empty:
        return go_full.head(0)
    kept = go_full.query(
        f"{EnrichmentCol.POP_COUNT} < {pop_count_max} and {EnrichmentCol.NAMESPACE} != 'MF'"
    )
    return kept.sort_values([EnrichmentCol.CLUSTER, EnrichmentCol.NAMESPACE, EnrichmentCol.TERM_ID])
