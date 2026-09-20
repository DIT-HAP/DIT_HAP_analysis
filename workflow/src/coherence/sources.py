"""Source adapters mapping PomBase grouping databases to a unified coherence long-table.

Every adapter returns the same contract columns (LONG_TABLE_COLUMNS), so the
downstream compute/plot stages are source-agnostic. Add a database = add one
adapter here + one entry in SOURCE_LOADERS + one line in config.coherence.sources.

Call `load_source()` rather than an adapter directly: it is what reads the shared
gene-name table once and hands it to the adapter.
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_annotated_members",
]

# PomBase macrocomplex annotation column -> canonical contract name.
_MACRO_RENAME = {
    "complex_term_id": "group_id",
    "GO_term_name": "group_name",
    "systematic_id": "Systematic ID",
    "symbol": "Name",
}

# PomBase's canonical id -> name table, relative to a version directory.
_GENE_METADATA_REL = Path("Gene_metadata") / "gene_IDs_names_products.tsv"

# --- GO GAF namespace loader (go_cc / go_bp) -------------------------------
# goatools namespace string per short code; also the config source name.
_NS_LONG = {"CC": "cellular_component", "BP": "biological_process"}
_NS_SOURCE = {"CC": "go_cc", "BP": "go_bp"}
# Match the canonical GO propagation exactly: workflow/src/enrichment/cluster_enrichment.py GO_LOAD_KWARGS.
_GO_LOAD_KWARGS = {"relationships": {"is_a", "part_of"}, "propagate_counts": True,
                   "load_obsolete": False, "prt": None}


# =============================================================================
# CORE LOGIC
# =============================================================================
def load_gene_names(pombase_dir: Path) -> dict[str, str]:
    """{systematic id: gene name} from PomBase's gene_IDs_names_products.tsv."""
    # The GAF and the complex annotation carry a systematic id but often no symbol,
    # which is why go_cc / go_bp used to end up with Name == Systematic ID. About a
    # third of PomBase genes have no name at all; those are left out here and fall
    # back to the systematic id in _finalize.
    path = Path(pombase_dir) / _GENE_METADATA_REL
    meta = pd.read_csv(path, sep="\t", usecols=["gene_systematic_id", "gene_name"])
    named = meta.dropna(subset=["gene_name"]).drop_duplicates(subset="gene_systematic_id")
    logger.info(f"gene names: {len(named):,} of {len(meta):,} genes in {path.name} have one")
    return dict(zip(named["gene_systematic_id"], named["gene_name"]))


def _finalize(df: pd.DataFrame, source: str, gene_names: Mapping[str, str]) -> pd.DataFrame:
    """Fill Name from `gene_names` then the systematic ID, add source + n_annotated_members, order columns."""
    df = df.copy()
    df["source"] = source
    # Adapter's own symbol first, then PomBase's name table, then the systematic id.
    # The column is never NaN: it feeds scored_member_names and the panel titles, and a
    # silent NaN there is what made every go_cc / go_bp Name read as a systematic id.
    df["Name"] = df["Name"].where(df["Name"].notna(), df["Systematic ID"].map(gene_names))
    df["Name"] = df["Name"].fillna(df["Systematic ID"])
    # Dedup + count on group_id (the stable GO term ID), NOT group_name. The old
    # compute_complex_coherence.py grouped on GO_term_name; keying on the ID is safer
    # (two distinct term IDs could share a name) and matches how the compute stage
    # later forms groups, keeping n_annotated_members consistent with the DR-member count.
    df = df.drop_duplicates(subset=["group_id", "Systematic ID"])
    counts = df.groupby("group_id")["Systematic ID"].transform("size")
    df["n_annotated_members"] = counts
    # go2genes values are sets, so upstream row order is hash-randomized; sort for
    # deterministic output. Benefits every source (macrocomplex tests assert by
    # value/set, so ordering is irrelevant to them).
    df = df.sort_values(["group_id", "Systematic ID"]).reset_index(drop=True)
    return df[LONG_TABLE_COLUMNS]


def load_macrocomplex(pombase_dir: Path, gene_names: Mapping[str, str]) -> pd.DataFrame:
    """Flat PomBase macromolecular_complex_annotation.tsv -> unified long-table."""
    path = Path(pombase_dir) / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    raw = pd.read_csv(path, sep="\t").rename(columns=_MACRO_RENAME)
    for required in ["group_id", "group_name", "Systematic ID"]:
        if required not in raw.columns:
            raise ValueError(f"macrocomplex annotation missing '{required}' (have: {list(raw.columns)})")
    if "Name" not in raw.columns:
        raw["Name"] = pd.NA
    return _finalize(
        raw[["group_id", "group_name", "Systematic ID", "Name"]], "go_macrocomplex", gene_names
    )


def load_gaf_namespace(
    pombase_dir: Path, namespace: str, gene_names: Mapping[str, str]
) -> pd.DataFrame:
    """GO GAF for one namespace (CC/BP), goatools-propagated, -> unified long-table."""
    # Reuses enrichment/ontology.py's OBO+GAF loading (is_a/part_of propagation,
    # propagate_counts=True), then keeps only terms in the requested namespace and
    # expands the propagated go2genes dict.
    #
    # Validate before the heavy import below, so a bad namespace fails on a cheap
    # check rather than on a missing goatools.
    if namespace not in _NS_LONG:
        raise ValueError(f"namespace must be one of {sorted(_NS_LONG)}, got {namespace!r}")

    # Imported here, not at module level: enrichment.ontology pulls in goatools,
    # which only the prepare_annotation rule's env has. coherence.io imports
    # LONG_TABLE_COLUMNS from this module, and compute_coherence runs in an env
    # without goatools — a top-level import would take that whole stage down.
    from enrichment.ontology import OntologyDataConfig, load_ontology_data

    od = Path(pombase_dir) / "ontologies_and_associations"
    data = OntologyDataConfig(
        ontology_obo=od / "go-basic.obo",
        ontology_association_gaf=od / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[],  # slim table not needed for raw term->gene expansion
    ).load_data()
    dag, _objanno, _ns2assoc, _gene2go, go2genes, _slim = load_ontology_data(data, **_GO_LOAD_KWARGS)

    ns_long = _NS_LONG[namespace]
    rows = []
    for term, genes in go2genes.items():
        rec = dag.get(term)
        if rec is None or rec.namespace != ns_long:
            continue
        for gene in genes:
            rows.append({"group_id": term, "group_name": rec.name,
                         "Systematic ID": gene, "Name": pd.NA})
    df = pd.DataFrame(rows, columns=["group_id", "group_name", "Systematic ID", "Name"])
    return _finalize(df, _NS_SOURCE[namespace], gene_names)


# A registry OF the adapters above, so it has to follow them — the section order
# ("each section depends only on what came before") is what keeps it out of
# GLOBAL CONSTANTS rather than an oversight.
SOURCE_LOADERS = {
    "go_macrocomplex": load_macrocomplex,
    "go_cc": lambda d, names: load_gaf_namespace(d, "CC", names),
    "go_bp": lambda d, names: load_gaf_namespace(d, "BP", names),
}


def load_source(source: str, pombase_dir: Path) -> pd.DataFrame:
    """Dispatch to `source`'s adapter, reading the shared gene-name table once."""
    if source not in SOURCE_LOADERS:
        raise ValueError(f"unknown source {source!r} (have: {sorted(SOURCE_LOADERS)})")
    pombase_dir = Path(pombase_dir)
    return SOURCE_LOADERS[source](pombase_dir, load_gene_names(pombase_dir))
