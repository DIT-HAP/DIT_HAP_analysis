"""
Ontology Enrichment Pipeline
============================

goatools-based GO/FYPO/MONDO enrichment (full + slim) for gene clusters,
ported from the enrichment slice of DIT_HAP_pipeline/workflow/src/enrichment_functions.py.
Pure-local and deterministic: reads OBO/GAF files, runs GOEnrichmentStudyNS with
FDR-BH, and formats results. The slim-association mapping (mapslim /
get_slim_ns2assoc) lives here rather than in ontology.py so that the feature-
collection loader (which does not need it) stays a cheap 6-tuple.

The STRING-db and REVIGO network functions are intentionally NOT here — they hit
external APIs and live in the optional network-enrichment rule (design doc §5,
Phase 2 Task 5).

Usage
-----
    from workflow.src.enrichment.pipeline import ontology_enrichment_pipeline
    full_df, slim_df, dag, objanno = ontology_enrichment_pipeline(
        ontology_data, query_genes, bg_genes, load_kwargs, enrichment_kwargs, format_kwargs
    )
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Data Processing Imports
import pandas as pd

# 2. Third-party Imports
from goatools.anno.gaf_reader import GafReader
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt

# 3. Local Imports
from workflow.src.enrichment.ontology import OntologyData, load_ontology_data


# =============================================================================
# SLIM ASSOCIATION MAPPING
# =============================================================================
def mapslim(term: str, dag: GODag, slim_dag: dict) -> tuple[set[str], set[str]]:
    """Map a term to its (direct, all) slim ancestors by walking every path to the root."""
    all_ancestors = set()
    covered_ancestors = set()

    paths = dag.paths_to_top(term)
    for path in paths:
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
    """Build a namespace->gene->slim-terms association for both direct and all ancestors."""
    term2slim = {}
    for term in dag:
        term2slim[term] = {}
        term2slim[term]["direct_ancestors"], term2slim[term]["all_ancestors"] = mapslim(term, dag, slim_dag)

    ns2slim_assoc = {"direct_ancestors": {}, "all_ancestors": {}}
    for ns, gene2terms in ns2assoc.items():
        ns2slim_assoc["direct_ancestors"][ns] = {}
        ns2slim_assoc["all_ancestors"][ns] = {}
        for gene, terms in gene2terms.items():
            ns2slim_assoc["direct_ancestors"][ns][gene] = set()
            ns2slim_assoc["all_ancestors"][ns][gene] = set()
            for term in terms:
                ns2slim_assoc["direct_ancestors"][ns][gene].update(term2slim[term]["direct_ancestors"])
                ns2slim_assoc["all_ancestors"][ns][gene].update(term2slim[term]["all_ancestors"])

    return ns2slim_assoc


# =============================================================================
# ENRICHMENT
# =============================================================================
def create_enrichment_dataframe(oea_results_sig: list[GOEnrichmentRecord], **kwargs) -> pd.DataFrame:
    """Build an enrichment DataFrame from GOEnrichmentRecord objects (fallback when get_goea_nts_prt fails)."""
    results_list = []
    for result in oea_results_sig:
        res_dict = {
            "GO": result.GO,
            "NS": result.NS,
            "name": result.name,
            "level": result.goterm.level,
            "depth": result.goterm.depth,
            "p_uncorrected": result.p_uncorrected,
            "p_fdr_bh": result.p_fdr_bh,
            "study_count": result.study_count,
            "study_n": result.study_n,
            "pop_count": result.pop_count,
            "pop_n": result.pop_n,
            "ratio_in_study": "/".join(map(str, result.ratio_in_study)),
            "ratio_in_pop": "/".join(map(str, result.ratio_in_pop)),
            "study_items": result.study_items,
            "pop_items": result.pop_items,
        }
        try:
            res_dict["defn"] = result.goterm.defn
        except Exception:
            res_dict["defn"] = ""

        if kwargs.get("itemid2name") is not None:
            study_items = sorted([kwargs["itemid2name"][item] for item in result.study_items])
            pop_items = sorted([kwargs["itemid2name"][item] for item in result.pop_items])
            res_dict["study_items"] = ", ".join(study_items)
            res_dict["pop_items"] = ", ".join(pop_items)
        else:
            res_dict["study_items"] = ", ".join(sorted(result.study_items))
            res_dict["pop_items"] = ", ".join(sorted(result.pop_items))

        results_list.append(res_dict)

    return pd.DataFrame(results_list)


def ontology_enrichment(
    query_genes: list[str], bg_genes: list[str], **kwargs
) -> tuple[GOEnrichmentStudyNS, list[GOEnrichmentRecord]]:
    """Run a GOEnrichmentStudyNS and keep only significantly ENRICHED terms (p_fdr_bh < alpha, enrichment == 'e')."""
    oea_obj = GOEnrichmentStudyNS(bg_genes, **kwargs)
    oea_results = oea_obj.run_study(query_genes, **kwargs)
    oea_results_sig = [r for r in oea_results if (r.p_fdr_bh < kwargs["alpha"]) and (r.enrichment == "e")]
    return oea_obj, oea_results_sig


def format_ontology_enrichment_results(label: str, oea_results_sig: list[GOEnrichmentRecord], **kwargs) -> pd.DataFrame:
    """Format significant enrichment records into a tidy DataFrame with gene_ratio/term_coverage."""
    try:
        oea_results_sig_prt = pd.DataFrame(get_goea_nts_prt(oea_results_sig, **kwargs))
    except Exception:
        oea_results_sig_prt = create_enrichment_dataframe(oea_results_sig, **kwargs)

    if oea_results_sig_prt.empty:
        return pd.DataFrame()

    oea_results_sig_prt["gene_ratio"] = round(
        oea_results_sig_prt["study_count"] / oea_results_sig_prt["study_n"], 2
    )
    oea_results_sig_prt["term_coverage"] = round(
        oea_results_sig_prt["study_count"] / oea_results_sig_prt["pop_count"], 2
    )
    oea_results_sig_prt["study_items"] = oea_results_sig_prt["study_items"].apply(
        lambda x: ",".join(sorted(x.split(", ")))
    )
    oea_results_sig_prt["pop_items"] = oea_results_sig_prt["pop_items"].apply(
        lambda x: ",".join(sorted(x.split(", ")))
    )

    kept_columns = [
        "NS", "GO", "name", "level", "depth", "p_uncorrected", "p_fdr_bh",
        "gene_ratio", "term_coverage", "study_count", "study_n", "pop_count",
        "pop_n", "ratio_in_study", "ratio_in_pop", "study_items", "pop_items",
    ]
    if "defn" in oea_results_sig_prt.columns:
        kept_columns.append("defn")

    return (
        oea_results_sig_prt[kept_columns]
        .copy()
        .rename(columns={"NS": "namespace", "GO": "term_id", "name": "term", "p_fdr_bh": "p_fdr"})
    )


def ontology_enrichment_pipeline(
    ontology_data: OntologyData,
    query_genes: list[str],
    bg_genes: list[str],
    load_kwargs: dict | None = None,
    enrichment_kwargs: dict | None = None,
    format_kwargs: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, GODag, GafReader]:
    """Full + slim enrichment for one gene set: load ontology, run study, format both result tables."""
    load_kwargs = dict(load_kwargs or {})
    enrichment_kwargs = dict(enrichment_kwargs or {})
    format_kwargs = dict(format_kwargs or {})

    # 6-tuple loader (feature-collection compatible); compute slim assoc here.
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(ontology_data, **load_kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    enrichment_kwargs["godag"] = dag
    enrichment_kwargs["ns2assoc"] = ns2assoc

    slim_enrichment_kwargs = enrichment_kwargs.copy()
    slim_enrichment_kwargs["godag"] = slim_dag
    slim_enrichment_kwargs["ns2assoc"] = ns2slim_assoc["all_ancestors"]
    slim_enrichment_kwargs["propagate_counts"] = False

    _, oea_results_sig = ontology_enrichment(query_genes, bg_genes, **enrichment_kwargs)
    _, oea_results_sig_slim = ontology_enrichment(query_genes, bg_genes, **slim_enrichment_kwargs)

    full_df = format_ontology_enrichment_results("full", oea_results_sig, **format_kwargs)
    slim_df = format_ontology_enrichment_results("slim", oea_results_sig_slim, **format_kwargs)

    return full_df, slim_df, dag, objanno
