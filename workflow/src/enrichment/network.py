"""
Network Enrichment Core (STRING + REVIGO)
=========================================

The non-deterministic half of the enrichment stage: STRING-db functional
enrichment and REVIGO semantic-similarity annotation, plus the response cache
that makes a second run offline and deterministic.

Both halves hit external web APIs, so the rules that call them stay out of
``rule all``; every response is cached under a cache dir, and a cache hit never
touches the network.

Input
-----
- Per-cluster gene lists and ``go_enrichment_full.tsv`` from the deterministic
  half (``cluster_enrichment.py``)

Output
------
- STRING enrichment frame (written to a workbook by the driver)
- GO enrichment frame with one REVIGO Representative/Eliminated/Dispensability
  triplet per configured cutoff

Usage
-----
    from enrichment.network import StringNetworkConfig, run_string_enrichment

    config = StringNetworkConfig(enrichment_dir=..., output_dir=..., cache_dir=...)
"""

# =============================================================================
# IMPORTS
# =============================================================================
# Standard library
from dataclasses import dataclass
import hashlib
from io import StringIO
from pathlib import Path
import time
from typing import Literal

from loguru import logger

# Data processing
import numpy as np
import pandas as pd

# Third-party
import requests
from requests.exceptions import RequestException

# Local project
from enrichment.schema import EnrichmentCol

# =============================================================================
# CONSTANTS
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
REVIGO_MAX_RETRIES = 3
REVIGO_RETRY_DELAY = 30

# The STRING namespace codes and the labels the shared schema uses.
_STRING_NAMESPACE_LABELS = {
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


# =============================================================================
# CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class StringNetworkConfig:
    """Inputs, output, and cache location for STRING enrichment."""

    enrichment_dir: Path
    output_dir: Path
    cache_dir: Path

    def validate(self) -> None:
        """Raise if the upstream deterministic-enrichment directory is missing."""
        if not self.enrichment_dir.exists():
            raise ValueError(f"Enrichment input dir not found: {self.enrichment_dir}")


@dataclass(kw_only=True, frozen=True)
class RevigoNetworkConfig:
    """Inputs, output, cache location, and cutoffs for REVIGO annotation."""

    enrichment_dir: Path
    output_dir: Path
    cache_dir: Path
    revigo_cutoffs: list[float]

    def validate(self) -> None:
        """Raise if the upstream deterministic-enrichment directory is missing."""
        if not self.enrichment_dir.exists():
            raise ValueError(f"Enrichment input dir not found: {self.enrichment_dir}")


# =============================================================================
# HELPERS
# =============================================================================
def _read_gene_list(path: Path) -> list[str]:
    """Read a newline-delimited gene list into a list of non-empty ids."""
    return [line for line in path.read_text().splitlines() if line.strip()]


def _cache_key(*parts: str) -> str:
    """Deterministic short hash of the request-defining parts, for cache filenames."""
    return hashlib.sha1("||".join(parts).encode()).hexdigest()[:16]


def _cache_load(cache_dir: Path | None, key: str) -> pd.DataFrame | None:
    """Return a cached frame for ``key``, or None on a miss or when caching is off."""
    if cache_dir is None:
        return None
    path = cache_dir / f"{key}.tsv"
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def _cache_store(cache_dir: Path | None, key: str, frame: pd.DataFrame) -> None:
    """Persist a frame under ``key`` so a later run needs no network."""
    if cache_dir is None:
        return
    cache_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(cache_dir / f"{key}.tsv", sep="\t", index=False)


# =============================================================================
# STRING-DB
# =============================================================================
def _post_string(request_url: str, params: dict) -> requests.Response:
    """POST once to STRING and raise on a transport or STRING-side error."""
    response = requests.post(request_url, data=params)
    response.raise_for_status()
    if "Something went wrong!" in response.text:
        raise RequestException("STRING communication_error: Something went wrong!")
    return response


def stringdb_api_functions(
    output_format: Literal["tsv", "tsv-no-header", "json", "xml"] = "xml",
    method: Literal["get_string_ids", "enrichment"] = "get_string_ids",
    params: dict | None = None,
) -> pd.DataFrame:
    """POST to the STRING API with retries and parse the response into a frame."""
    params = params or {}
    request_url = "/".join([STRING_API_URL, output_format, method])
    response = _post_string_with_retries(request_url, params)

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


def _post_string_with_retries(request_url: str, params: dict) -> requests.Response:
    """Retry the STRING POST until it succeeds or the attempt budget runs out."""
    for attempt in range(STRING_MAX_RETRIES):
        try:
            return _post_string(request_url, params)
        except RequestException as error:
            if attempt == STRING_MAX_RETRIES - 1:
                raise ValueError(
                    f"STRING API failed after {STRING_MAX_RETRIES} retries: {error}"
                ) from error
            time.sleep(STRING_RETRY_DELAY)
    raise ValueError(f"STRING API failed after {STRING_MAX_RETRIES} retries")


def format_string_enrichment_results(
    enrichment_df: pd.DataFrame, query_genes: list[str], background_genes: list[str]
) -> pd.DataFrame:
    """Rename STRING columns to the shared schema and order the namespaces."""
    renamed = {
        "term": "term_id",
        "category": "namespace",
        "description": "term",
        "p_value": "p_uncorrected",
        "fdr": "p_fdr",
        "gene_ratio": "gene_ratio",
        "term_coverage": "term_coverage",
        "number_of_genes": "study_count",
        "study_n": "study_n",
        "number_of_genes_in_background": "pop_count",
        "pop_n": "pop_n",
        "ratio_in_study": "ratio_in_study",
        "ratio_in_pop": "ratio_in_pop",
        "preferredNames": "study_items",
    }
    formatted = enrichment_df.rename(columns=renamed)
    formatted["study_n"] = len(query_genes)
    formatted["pop_n"] = len(background_genes)
    formatted["ratio_in_study"] = (
        formatted["study_count"].astype(str) + "/" + formatted["study_n"].astype(str)
    )
    formatted["ratio_in_pop"] = (
        formatted["pop_count"].astype(str) + "/" + formatted["pop_n"].astype(str)
    )
    formatted["gene_ratio"] = round(formatted["study_count"] / formatted["study_n"], 2)
    formatted["term_coverage"] = round(formatted["study_count"] / formatted["pop_count"], 2)
    formatted = formatted[list(renamed.values())]

    formatted[EnrichmentCol.NAMESPACE] = formatted[EnrichmentCol.NAMESPACE].map(
        _STRING_NAMESPACE_LABELS
    )
    namespace_order = {
        value: index for index, value in enumerate(_STRING_NAMESPACE_LABELS.values())
    }
    formatted["namespace_order"] = formatted[EnrichmentCol.NAMESPACE].map(namespace_order)
    return formatted.sort_values(by="namespace_order").drop(columns=["namespace_order"])


def stringdb_enrichment(
    query_genes: list[str], bg_genes: list[str], cache_dir: Path | None = None
) -> pd.DataFrame:
    """STRING enrichment for query vs background genes; cache-first, empty on failure."""
    key = _cache_key("string", ",".join(sorted(query_genes)), ",".join(sorted(bg_genes)))
    cached = _cache_load(cache_dir, key)
    if cached is not None:
        return cached

    string_ids = _resolve_string_ids(bg_genes)
    enrichment_params = {
        "identifiers": STRING_SEPARATOR.join(query_genes),
        "species": STRING_SPECIES_ID,
        "background_string_identifiers": STRING_SEPARATOR.join(string_ids),
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    try:
        enrichment_df = stringdb_api_functions(
            output_format="xml", method="enrichment", params=enrichment_params
        )
    except (RequestException, ValueError, KeyError, pd.errors.ParserError) as error:
        logger.warning(f"STRING enrichment failed: {error}")
        return pd.DataFrame()

    formatted = format_string_enrichment_results(enrichment_df, query_genes, bg_genes)
    _cache_store(cache_dir, key, formatted)
    return formatted


def _resolve_string_ids(bg_genes: list[str]) -> list[str]:
    """Map PomBase systematic ids to STRING ids for the background gene set."""
    id_params = {
        "identifiers": STRING_SEPARATOR.join(bg_genes),
        "species": STRING_SPECIES_ID,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    return stringdb_api_functions(output_format="xml", method="get_string_ids", params=id_params)[
        "stringId"
    ].tolist()


# =============================================================================
# REVIGO
# =============================================================================
def _post_revigo(payload: dict) -> pd.DataFrame:
    """POST once to REVIGO and parse the semantic-similarity table."""
    response = requests.post(REVIGO_URL, data=payload, timeout=60)
    response.raise_for_status()
    return pd.read_html(StringIO(response.text))[0]


def _post_revigo_with_retries(payload: dict) -> pd.DataFrame:
    """Retry the REVIGO POST until it succeeds or the attempt budget runs out."""
    for attempt in range(REVIGO_MAX_RETRIES):
        try:
            return _post_revigo(payload)
        except RequestException as error:
            if attempt == REVIGO_MAX_RETRIES - 1:
                raise ValueError(
                    f"Error performing REVIGO analysis after {REVIGO_MAX_RETRIES} retries: {error}"
                ) from error
            time.sleep(REVIGO_RETRY_DELAY)
    raise ValueError(f"Error performing REVIGO analysis after {REVIGO_MAX_RETRIES} retries")


def _resolve_representative_names(revigo_res: pd.DataFrame) -> pd.DataFrame:
    """Replace REVIGO's zero-padded representative ids with their term names."""
    id2name = dict(zip(revigo_res["Term ID"], revigo_res["Name"], strict=True))
    is_self = revigo_res["Representative"].isna()
    padded_ids = "GO:" + revigo_res["Representative"].fillna(0).astype("int64").astype(
        str
    ).str.rjust(7, "0")
    representative = padded_ids.map(id2name)
    revigo_res["Representative"] = np.where(is_self, revigo_res["Name"], representative)
    return revigo_res


def revigo_analysis(
    enrich_df: pd.DataFrame,
    cut_off: float = 0.7,
    cache_dir: Path | None = None,
) -> pd.DataFrame:
    """Summarize GO enrichment via REVIGO (semantic-similarity reduction); cache-first."""
    go_ids = enrich_df[EnrichmentCol.TERM_ID].tolist()
    p_values = enrich_df[EnrichmentCol.P_FDR].tolist()
    cache_parts = ",".join(f"{go}:{p}" for go, p in zip(go_ids, p_values, strict=True))
    key = _cache_key("revigo", str(cut_off), cache_parts)
    cached = _cache_load(cache_dir, key)
    if cached is not None:
        return cached

    payload = {
        "cutoff": f"{cut_off}",
        "valueType": "pvalue",
        "speciesTaxon": REVIGO_SPECIES_TAXON,
        "measure": REVIGO_MEASURE,
        "goList": "\n".join(f"{go}\t{p}" for go, p in zip(go_ids, p_values, strict=True)),
    }
    revigo_res = _resolve_representative_names(_post_revigo_with_retries(payload))
    _cache_store(cache_dir, key, revigo_res)
    return revigo_res


# =============================================================================
# DRIVERS
# =============================================================================
def run_string_enrichment(config: StringNetworkConfig) -> pd.DataFrame:
    """STRING enrichment for every cluster, concatenated with a Cluster column."""
    bg_genes = _read_gene_list(config.enrichment_dir / "DIT_HAP_all_genes.txt")
    cluster_files = sorted(config.enrichment_dir.glob("DIT_HAP_cluster_*_genes.txt"))

    per_cluster: dict[int, pd.DataFrame] = {}
    for cluster_file in cluster_files:
        cluster = int(cluster_file.stem.replace("DIT_HAP_cluster_", "").replace("_genes", ""))
        genes = _read_gene_list(cluster_file)
        logger.info(f"  STRING cluster {cluster}: {len(genes)} genes")
        per_cluster[cluster] = stringdb_enrichment(genes, bg_genes, cache_dir=config.cache_dir)

    non_empty = {cluster: frame for cluster, frame in per_cluster.items() if not frame.empty}
    if not non_empty:
        return pd.DataFrame()
    return pd.concat(non_empty, axis=0).droplevel(1).rename_axis("Cluster").reset_index()


def annotate_go_with_revigo(config: RevigoNetworkConfig) -> pd.DataFrame:
    """Add one REVIGO column triplet per cutoff to the GO table."""
    go_path = config.enrichment_dir / "go_enrichment_full.tsv"
    go_df = pd.read_csv(go_path, sep="\t")
    if go_df.empty:
        return go_df

    annotated_parts = []
    for _, namespace_df in go_df.groupby([EnrichmentCol.CLUSTER, EnrichmentCol.NAMESPACE]):
        annotated_parts.append(
            _annotate_namespace(namespace_df, config.revigo_cutoffs, config.cache_dir)
        )
    return pd.concat(annotated_parts, axis=0)


def _annotate_namespace(
    namespace_df: pd.DataFrame, cutoffs: list[float], cache_dir: Path
) -> pd.DataFrame:
    """Merge every configured cutoff's REVIGO columns onto one cluster x namespace slice."""
    annotated = namespace_df.copy()
    for threshold in cutoffs:
        revigo = revigo_analysis(annotated, cut_off=threshold, cache_dir=cache_dir)
        revigo = (
            revigo.set_index("Term ID")[["Dispensability", "Eliminated", "Representative"]]
            .add_suffix(f"_{threshold}")
            .reset_index()
        )
        annotated = annotated.merge(
            revigo, left_on=EnrichmentCol.TERM_ID, right_on="Term ID", how="left"
        ).drop(columns=["Term ID"])
    return annotated
