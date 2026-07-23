"""Source adapters mapping PomBase grouping databases to a unified coherence long-table.

Every adapter returns the same contract columns (LONG_TABLE_COLUMNS), so the
downstream compute/plot stages are source-agnostic. Add a database = add one
adapter here + one entry in SOURCE_LOADERS + one line in config.coherence.sources.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]

# PomBase macrocomplex annotation column -> canonical contract name.
_MACRO_RENAME = {
    "complex_term_id": "group_id",
    "GO_term_name": "group_name",
    "systematic_id": "Systematic ID",
    "symbol": "Name",
}


def _finalize(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Fill missing Name with Systematic ID, add source + n_group_genes, order columns."""
    df = df.copy()
    df["source"] = source
    df["Name"] = df["Name"].fillna(df["Systematic ID"])
    # Dedup + count on group_id (the stable GO term ID), NOT group_name. The old
    # compute_complex_coherence.py grouped on GO_term_name; keying on the ID is safer
    # (two distinct term IDs could share a name) and matches how the compute stage
    # later forms groups, keeping n_group_genes consistent with the DR-member count.
    df = df.drop_duplicates(subset=["group_id", "Systematic ID"])
    counts = df.groupby("group_id")["Systematic ID"].transform("size")
    df["n_group_genes"] = counts
    # go2genes values are sets, so upstream row order is hash-randomized; sort for
    # deterministic output. Benefits every source (macrocomplex tests assert by
    # value/set, so ordering is irrelevant to them).
    df = df.sort_values(["group_id", "Systematic ID"]).reset_index(drop=True)
    return df[LONG_TABLE_COLUMNS]


def load_macrocomplex(pombase_dir: Path) -> pd.DataFrame:
    """Flat PomBase macromolecular_complex_annotation.tsv -> unified long-table."""
    path = Path(pombase_dir) / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    raw = pd.read_csv(path, sep="\t").rename(columns=_MACRO_RENAME)
    for required in ["group_id", "group_name", "Systematic ID"]:
        if required not in raw.columns:
            raise ValueError(f"macrocomplex annotation missing '{required}' (have: {list(raw.columns)})")
    if "Name" not in raw.columns:
        raw["Name"] = pd.NA
    return _finalize(raw[["group_id", "group_name", "Systematic ID", "Name"]], "go_macrocomplex")


# --- GO GAF namespace loader (go_cc / go_bp) -------------------------------
# goatools namespace string per short code; also the config source name.
_NS_LONG = {"CC": "cellular_component", "BP": "biological_process"}
_NS_SOURCE = {"CC": "go_cc", "BP": "go_bp"}
# Match the canonical GO propagation exactly: workflow/src/enrichment/cluster_enrichment.py GO_LOAD_KWARGS.
_GO_LOAD_KWARGS = {"relationships": {"is_a", "part_of"}, "propagate_counts": True,
                   "load_obsolete": False, "prt": None}


def load_gaf_namespace(pombase_dir: Path, namespace: str) -> pd.DataFrame:
    """GO GAF for one namespace (CC/BP), goatools-propagated, -> unified long-table.

    Reuses enrichment/ontology.py's OBO+GAF loading (is_a/part_of propagation,
    propagate_counts=True), then keeps only terms in the requested namespace and
    expands the propagated go2genes dict.
    """
    from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data

    if namespace not in _NS_LONG:
        raise ValueError(f"namespace must be one of {sorted(_NS_LONG)}, got {namespace!r}")
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
    return _finalize(df, _NS_SOURCE[namespace])


SOURCE_LOADERS = {
    "go_macrocomplex": load_macrocomplex,
    "go_cc": lambda d: load_gaf_namespace(d, "CC"),
    "go_bp": lambda d: load_gaf_namespace(d, "BP"),
}
