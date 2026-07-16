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
# 1. Standard Library Imports
import hashlib
import json
import time
from io import StringIO
from pathlib import Path
from typing import Literal

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import requests
from goatools.anno.gaf_reader import GafReader
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt
from requests.exceptions import RequestException

# 4. Local Imports
from workflow.src.enrichment.ontology import OntologyData, load_ontology_data

# =============================================================================
# NETWORK CONSTANTS (STRING / REVIGO)
# =============================================================================
STRING_API_URL = "https://string-db.org/api"
STRING_SPECIES_ID = "4896"
STRING_CALLER_IDENTITY = "dit-hap.analysis"
STRING_SEPARATOR = "%0d"
STRING_MAX_RETRIES = 5
STRING_RETRY_DELAY = 10

REVIGO_URL = "http://revigo.irb.hr/Revigo"
REVIGO_SPECIES_TAXON = "284812"
REVIGO_MEASURE = "SIMREL"


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


# =============================================================================
# NETWORK ENRICHMENT — CACHE HELPERS
# =============================================================================
def _cache_key(*parts: str) -> str:
    """Deterministic short hash of the request-defining parts, for cache filenames."""
    return hashlib.sha1("||".join(parts).encode()).hexdigest()[:16]


def _cache_load(cache_dir: Path | None, key: str) -> pd.DataFrame | None:
    """Return a cached DataFrame for `key`, or None on a miss / when caching is off."""
    if cache_dir is None:
        return None
    path = cache_dir / f"{key}.tsv"
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def _cache_store(cache_dir: Path | None, key: str, df: pd.DataFrame) -> None:
    """Persist a DataFrame under `key` so future runs are deterministic without the network."""
    if cache_dir is None:
        return
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_dir / f"{key}.tsv", sep="\t", index=False)


# =============================================================================
# STRING-DB ENRICHMENT (network)
# =============================================================================
def stringdb_api_functions(
    output_format: Literal["tsv", "tsv-no-header", "json", "xml"] = "xml",
    method: Literal["get_string_ids", "enrichment"] = "get_string_ids",
    params: dict | None = None,
) -> pd.DataFrame:
    """POST to the STRING API with retries and parse the response into a DataFrame."""
    params = params or {}
    request_url = "/".join([STRING_API_URL, output_format, method])
    response = None
    for attempt in range(STRING_MAX_RETRIES):
        try:
            response = requests.post(request_url, data=params)
            response.raise_for_status()
            if "Something went wrong!" in response.text:
                raise RequestException("STRING communication_error: Something went wrong!")
            break
        except RequestException as e:
            if attempt < STRING_MAX_RETRIES - 1:
                time.sleep(STRING_RETRY_DELAY)
            else:
                raise ValueError(f"STRING API failed after {STRING_MAX_RETRIES} retries: {e}")

    match output_format:
        case "xml":
            return pd.read_xml(StringIO(response.text))
        case "tsv":
            return pd.read_csv(StringIO(response.text), sep="\t")
        case "tsv-no-header":
            return pd.read_csv(StringIO(response.text), sep="\t", header=None)
        case "json":
            return pd.read_json(StringIO(response.text))
        case _:
            raise ValueError(f"Invalid output format: {output_format}")


def format_string_enrichment_results(enrichment_df: pd.DataFrame, query_genes: list[str], background_genes: list[str]) -> pd.DataFrame:
    """Rename STRING enrichment columns to the shared schema and order namespaces."""
    renamed = {
        "term": "term_id", "category": "namespace", "description": "term",
        "p_value": "p_uncorrected", "fdr": "p_fdr", "gene_ratio": "gene_ratio",
        "term_coverage": "term_coverage", "number_of_genes": "study_count", "study_n": "study_n",
        "number_of_genes_in_background": "pop_count", "pop_n": "pop_n",
        "ratio_in_study": "ratio_in_study", "ratio_in_pop": "ratio_in_pop", "preferredNames": "study_items",
    }
    enrichment_df = enrichment_df.rename(columns=renamed)
    enrichment_df["study_n"] = len(query_genes)
    enrichment_df["pop_n"] = len(background_genes)
    enrichment_df["ratio_in_study"] = enrichment_df.apply(lambda r: f"{r['study_count']}/{r['study_n']}", axis=1)
    enrichment_df["ratio_in_pop"] = enrichment_df.apply(lambda r: f"{r['pop_count']}/{r['pop_n']}", axis=1)
    enrichment_df["gene_ratio"] = round(enrichment_df["study_count"] / enrichment_df["study_n"], 2)
    enrichment_df["term_coverage"] = round(enrichment_df["study_count"] / enrichment_df["pop_count"], 2)
    enrichment_df = enrichment_df[list(renamed.values())]

    namespace_description = {
        "Process": "Biological Process (Gene Ontology)",
        "Component": "Cellular Component (Gene Ontology)",
        "Function": "Molecular Function (Gene Ontology)",
        "PMID": "Reference Publications (PubMed)",
        "NetworkNeighborAL": "Local Network Cluster (STRING)",
        "KEGG": "KEGG Pathways",
        "RCTM": "Reactome Pathways",
        "COMPARTMENTS": "Subcellular Localization (COMPARTMENTS)",
        "Keyword": "Annotated Keywords (UniProt)",
        "InterPro": "Protein Domains and Features (InterPro)",
        "SMART": "Protein Domains and Features (SMART)",
    }
    enrichment_df["namespace"] = enrichment_df["namespace"].map(namespace_description)
    namespace_order = {v: i for i, v in enumerate(namespace_description.values())}
    enrichment_df["namespace_order"] = enrichment_df["namespace"].map(namespace_order)
    enrichment_df = enrichment_df.sort_values(by="namespace_order").drop("namespace_order", axis=1)
    return enrichment_df


def stringdb_enrichment(query_genes: list[str], bg_genes: list[str], cache_dir: Path | None = None) -> pd.DataFrame:
    """STRING functional enrichment for query vs background genes; cache-first, empty frame on failure."""
    key = _cache_key("string", ",".join(sorted(query_genes)), ",".join(sorted(bg_genes)))
    cached = _cache_load(cache_dir, key)
    if cached is not None:
        return cached

    get_string_id_params = {
        "identifiers": STRING_SEPARATOR.join(bg_genes),
        "species": STRING_SPECIES_ID,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    string_ids = stringdb_api_functions(output_format="xml", method="get_string_ids", params=get_string_id_params)["stringId"].tolist()

    enrichment_params = {
        "identifiers": STRING_SEPARATOR.join(query_genes),
        "species": STRING_SPECIES_ID,
        "background_string_identifiers": STRING_SEPARATOR.join(string_ids),
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    try:
        enrichment_df = stringdb_api_functions(output_format="xml", method="enrichment", params=enrichment_params)
        enrichment_df = format_string_enrichment_results(enrichment_df, query_genes, bg_genes)
        _cache_store(cache_dir, key, enrichment_df)
        return enrichment_df
    except Exception as e:
        print(f"Error performing STRING enrichment analysis: {e}")
        return pd.DataFrame()


# =============================================================================
# REVIGO ANALYSIS (network)
# =============================================================================
def revigo_analysis(
    enrich_df: pd.DataFrame,
    cut_off: float = 0.7,
    max_retries: int = 3,
    retry_delay: int = 30,
    cache_dir: Path | None = None,
) -> pd.DataFrame:
    """Summarize GO enrichment via REVIGO (semantic-similarity reduction); cache-first."""
    GOs = enrich_df["term_id"].tolist()
    padj = enrich_df["p_fdr"].tolist()
    key = _cache_key("revigo", str(cut_off), ",".join(f"{g}:{p}" for g, p in zip(GOs, padj)))
    cached = _cache_load(cache_dir, key)
    if cached is not None:
        return cached

    data = "\n".join([f"{go}\t{p}" for go, p in zip(GOs, padj)])
    payload = {
        "cutoff": f"{cut_off}",
        "valueType": "pvalue",
        "speciesTaxon": REVIGO_SPECIES_TAXON,
        "measure": REVIGO_MEASURE,
        "goList": data,
    }
    for attempt in range(max_retries):
        try:
            r = requests.post(REVIGO_URL, data=payload, timeout=60)
            r.raise_for_status()
            revigo_res = pd.read_html(StringIO(r.text))[0]
            id2name = dict(zip(revigo_res["Term ID"], revigo_res["Name"]))
            # Representative is either NaN (term is its own representative -> use Name)
            # or a zero-padded GO integer whose name we look up.
            revigo_res["Representative"] = revigo_res.apply(
                lambda row: row["Name"]
                if np.isnan(row["Representative"])
                else id2name["GO:" + str(int(row["Representative"])).rjust(7, "0")],
                axis=1,
            )
            _cache_store(cache_dir, key, revigo_res)
            return revigo_res
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
            else:
                raise ValueError(f"Error performing REVIGO analysis after {max_retries} retries: {e}")
