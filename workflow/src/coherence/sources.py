"""Source adapters mapping PomBase grouping databases to a unified coherence long-table.

Every adapter returns the same contract columns (LONG_TABLE_COLUMNS), so the
downstream compute/plot stages are source-agnostic. Add a database = add one
adapter here + one entry in SOURCE_LOADERS + one line in config.coherence.sources.

Call `load_source()` rather than an adapter directly: it is what reads the shared
gene-name table once and hands it to the adapter. Every adapter takes the same
`(pombase_dir, gene_names, kegg_dir)` triple so the registry stays uniform;
`kegg_dir` is only read by the kegg_* adapters, which need it to find the
kegg_parser derived tables.
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

# --- KEGG derived tables (kegg_parser output under resources/external/kegg) ---
# One flat file per KEGG product; see load_kegg_brite / load_kegg_pathway.
_KEGG_BRITE_FILE, _KEGG_PATHWAY_FILE = "brite_flat.tsv", "pathway_gene_mapping.tsv"


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


def load_macrocomplex(
    pombase_dir: Path, gene_names: Mapping[str, str], kegg_dir: Path | None = None,  # noqa: ARG001 - uniform registry signature
) -> pd.DataFrame:
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
    pombase_dir: Path, namespace: str, gene_names: Mapping[str, str], kegg_dir: Path | None = None,  # noqa: ARG001 - uniform registry signature
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


def load_symbol_to_systematic(pombase_dir: Path) -> dict[str, str]:
    """KEGG `Gene_Symbol` -> PomBase systematic id.

    KEGG carries PomBase's gene name when the gene has one and the systematic id
    (KEGG's `SPOM_` prefix stripped) when it does not, so this resolves both: the
    name->id table plus every systematic id as its own key. Symbols that resolve to
    neither (genes absent from this PomBase release) are simply absent from the
    map, and `Series.map` turns those into NaN for the caller to drop.
    """
    meta = pd.read_csv(
        Path(pombase_dir) / _GENE_METADATA_REL,
        sep="\t",
        usecols=["gene_systematic_id", "gene_name"],
    )
    mapping = dict(zip(meta["gene_systematic_id"], meta["gene_systematic_id"]))
    named = meta.dropna(subset=["gene_name"]).drop_duplicates(subset="gene_name")
    mapping.update(dict(zip(named["gene_name"], named["gene_systematic_id"])))
    logger.info(f"KEGG symbol map: {len(mapping):,} keys "
                f"({len(named):,} gene names + {len(meta):,} systematic ids)")
    return mapping


def _kegg_frame(
    group_id: pd.Series, group_name: pd.Series, symbols: pd.Series, resolver: Mapping[str, str]
) -> pd.DataFrame:
    """Build the four pre-`_finalize` columns, dropping genes with no PomBase id."""
    df = pd.DataFrame({
        "group_id": group_id,
        "group_name": group_name,
        "Systematic ID": symbols.map(resolver),
        "Name": symbols,
    })
    n_raw = len(df)
    df = df.dropna(subset=["Systematic ID"])
    if len(df) < n_raw:
        logger.info(f"KEGG: kept {len(df):,} of {n_raw:,} rows "
                    f"(dropped {n_raw - len(df):,} with no PomBase systematic id)")
    return df


def load_kegg_pathway(
    kegg_dir: Path, pombase_dir: Path, gene_names: Mapping[str, str]
) -> pd.DataFrame:
    """KEGG PATHWAY gene mapping -> unified long-table, one group per pathway id."""
    if kegg_dir is None:
        raise ValueError("kegg_pathway needs --kegg-dir (the kegg_parser derived tables)")
    raw = pd.read_csv(
        Path(kegg_dir) / _KEGG_PATHWAY_FILE,
        sep="\t",
        usecols=["Gene_Symbol", "Pathway_ID", "Pathway_Name"],
    )
    df = _kegg_frame(
        raw["Pathway_ID"], raw["Pathway_Name"], raw["Gene_Symbol"],
        load_symbol_to_systematic(pombase_dir),
    )
    return _finalize(df, "kegg_pathway", gene_names)


def _brite_group_ids(raw: pd.DataFrame) -> pd.Series:
    """One id per BRITE node: its KEGG id where it has one, else its full path.

    `Level_D_ID` is never empty, but it is not always an id: kegg_parser's own
    contract is that "labels carrying no id at all fall back to the label text
    itself so ids stay non-empty" — so `Level_D_ID == Level_D` is exactly the
    no-id case. It is the common case: only the nodes carrying a `[PATH:...]`
    bracket, a leading class code or an EC number have a real one (8.6k of 13.1k
    rows here; the ribosome, kinase and most of the transporter trees have none).

    The fallback has to be the PATH, not the leaf label. `Level_D` is not unique
    within a tree — four kingdom branches all end at "Large subunit" — so keying on
    it merges unrelated nodes: 68 `Level_D_ID`s cover more than one node and 17 more
    collide across trees, while the 8.6k real ids collide neither way. The tree id
    still prefixes the path, because the same path is not contractually unique
    across trees.
    """
    levels = [f"Level_{name}" for name in "ABCD"]
    path = raw[levels].fillna("").astype(str).agg(" > ".join, axis=1)
    fallback = raw["BRITE_ID"].astype(str) + ":" + path
    node_id = raw["Level_D_ID"].fillna("").astype(str)
    has_real_id = (node_id != "") & (node_id != raw["Level_D"].fillna("").astype(str))
    return node_id.where(has_real_id, fallback)


def load_kegg_brite(
    kegg_dir: Path, pombase_dir: Path, gene_names: Mapping[str, str]
) -> pd.DataFrame:
    """KEGG BRITE leaves -> unified long-table, one group per gene's deepest node.

    Level_D is the deepest classification node kegg_parser emits: it forward-fills an
    empty level from the one above, so A-D are gap-free and Level_D already means
    "the most specific category this leaf has" (it equals Level_C for 82% of rows and
    Level_B for another 10%). The three shallower levels are deliberately NOT emitted
    as their own groups - within one tree an ancestor and its descendant share members
    by construction, so that would just hand dedup a pile of nested pairs to merge
    back. The label alone is not a unique key ("Others" names a node in 11 different
    trees), and neither is `Level_D_ID` — it is the label again wherever KEGG has no
    id for the node. `_brite_group_ids` is where the id is built.
    """
    if kegg_dir is None:
        raise ValueError("kegg_brite needs --kegg-dir (the kegg_parser derived tables)")
    raw = pd.read_csv(
        Path(kegg_dir) / _KEGG_BRITE_FILE,
        sep="\t",
        usecols=["BRITE_ID", "Level_A", "Level_B", "Level_C", "Level_D",
                 "Level_D_ID", "Gene_Symbol"],
    )
    unclassified = raw["Level_D"].isna()
    if unclassified.any():
        logger.warning(f"BRITE: dropping {int(unclassified.sum()):,} rows with no classification node")
        raw = raw[~unclassified]
    df = _kegg_frame(
        _brite_group_ids(raw), raw["Level_D"], raw["Gene_Symbol"],
        load_symbol_to_systematic(pombase_dir),
    )
    return _finalize(df, "kegg_brite", gene_names)


# A registry OF the adapters above, so it has to follow them — the section order
# ("each section depends only on what came before") is what keeps it out of
# GLOBAL CONSTANTS rather than an oversight.
SOURCE_LOADERS = {
    "go_macrocomplex": load_macrocomplex,
    "go_cc": lambda d, names, kegg=None: load_gaf_namespace(d, "CC", names),  # noqa: ARG005 - uniform registry signature
    "go_bp": lambda d, names, kegg=None: load_gaf_namespace(d, "BP", names),  # noqa: ARG005 - uniform registry signature
    "kegg_pathway": lambda d, names, kegg: load_kegg_pathway(kegg, d, names),
    "kegg_brite": lambda d, names, kegg: load_kegg_brite(kegg, d, names),
}


def load_source(source: str, pombase_dir: Path, kegg_dir: Path | None = None) -> pd.DataFrame:
    """Dispatch to `source`'s adapter, reading the shared gene-name table once."""
    if source not in SOURCE_LOADERS:
        raise ValueError(f"unknown source {source!r} (have: {sorted(SOURCE_LOADERS)})")
    pombase_dir = Path(pombase_dir)
    return SOURCE_LOADERS[source](pombase_dir, load_gene_names(pombase_dir), kegg_dir)
