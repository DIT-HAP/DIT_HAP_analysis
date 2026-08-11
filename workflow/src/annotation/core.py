"""
Pombe Gene Annotation Assembly
================================

Builds a per-gene annotation reference table (budding yeast / human orthologs,
pombe essentiality, functional annotation) and joins it onto arbitrary user
tables keyed by pombe systematic ID.

Ortholog fields in PomBase's `curated_orthologs/` files carry three distinct
kinds of structure that must not be flattened: `|` separates *independent*
ortholog genes, `+` joins *fragments of one* ortholog relation (a pombe gene
fused relative to S. cerevisiae), and `(N)`/`(C)` mark which terminus a fragment
corresponds to. Independent orthologs become separate groups so that per-ortholog
columns (name, essentiality, ORF qualifier) stay positionally aligned.

Input
-----
- PomBase version directory (curated_orthologs, Gene_metadata, ontologies, Protein_features)
- SGD `SGD_features.tab` and `phenotype_data.tab`
- Curated deletion-library categories xlsx

Output
------
- A DataFrame indexed by pombe systematic ID, one column per annotation field

Author:   Yusheng Yang (guidance) + Claude Opus 5 (implementation)
Date:     2026-08-11
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# PomBase's sentinel for "this gene has no ortholog in the target species".
_NO_ORTHOLOG = "NONE"

# Separators inside a PomBase curated_orthologs field.
_INDEPENDENT_SEP = "|"  # separates independent ortholog genes
_FRAGMENT_SEP = "+"  # joins fragments of a single (fusion) ortholog relation

# Only a null (deletion) mutant speaks to essentiality; conditional / overexpression
# / reduction-of-function rows describe something else. Likewise only the exact
# viable/inviable calls count — "viability: decreased" is a graded phenotype.
_NULL_MUTANT = "null"
_INVIABLE = "inviable"
_VIABLE = "viable"
_VIABILITY_PHENOTYPES = (_INVIABLE, _VIABLE)
_CONFLICTING = "conflicting"

# phenotype_data.tab is headerless and ragged: most rows carry these 14 fields, but
# 39 rows have a spurious 15th and 19 rows are broken by an embedded newline (a
# 13-field row plus a short continuation). Names are supplied positionally and extra
# fields are tolerated so a handful of malformed rows cannot abort the whole read.
_SGD_PHENOTYPE_COLUMNS = [
    "feature_name", "feature_type", "gene_name", "sgdid", "reference",
    "experiment_type", "mutant_type", "allele", "strain_background", "phenotype",
    "chemical", "condition", "details", "reporter",
]

# Annotation column -> source column in resources/curated/deletion_library_categories.xlsx.
_DELETION_LIBRARY_COLUMNS = {
    "Sp_deletion_essentiality": "Gene dispensability. This study",
    "Sp_deletion_phenotype": "Phenotypic classification used for analysis",
    "Sp_growth_category": "Category",
}

# protein_families_and_domains.tsv mixes PFAM with PANTHER/PROSITE/etc.
_PFAM_DATABASE = "PFAM"

# goatools namespace -> annotation column. GO slim (not full GO) because a gene
# carries a dozen full-GO terms but only a handful of readable slim labels.
_GO_SLIM_NAMESPACES = {
    "BP": "GO_slim_BP",
    "CC": "GO_slim_CC",
    "MF": "GO_slim_MF",
}

# Columns with this suffix hold counts and are coerced to nullable Int64 after joining.
_COUNT_COLUMN_SUFFIX = "_count"

# The curated gRNA fitted-parameters table names depletion rate/lag with upstream's
# legacy um/lam; accept DR/DL too in case upstream renames them. Output is prefixed
# gRNA_ to keep it distinct from the gene-level DR/DL in clustering tables.
_GRNA_GENE_COLUMN = "Systematic ID"
_GRNA_DEPLETION_COLUMNS = {
    "um": "gRNA_DR",
    "lam": "gRNA_DL",
    "DR": "gRNA_DR",
    "DL": "gRNA_DL",
}
_GRNA_TARGET_COLUMNS = ["gRNA_DR", "gRNA_DL"]

# SGD_features.tab covers many feature types (CDS, intron, ARS, ...); only ORF rows
# carry the systematic name that PomBase orthologs refer to.
_SGD_ORF_TYPE = "ORF"

# SGD_features.tab is also headerless, but uniformly 16 fields.
_SGD_FEATURE_COLUMNS = [
    "sgdid", "feature_type", "qualifier", "systematic_name", "standard_name",
    "alias", "parent_feature", "secondary_sgdid", "chromosome", "start", "stop",
    "strand", "genetic_position", "coordinate_version", "sequence_version",
    "description",
]


# =============================================================================
# CORE LOGIC
# =============================================================================
def parse_ortholog_field(field: str | float) -> list[str]:
    """Split a PomBase ortholog field into one entry per independent ortholog."""
    if not isinstance(field, str):
        return []

    field = field.strip()
    if not field or field == _NO_ORTHOLOG:
        return []

    groups = []
    for group in field.split(_INDEPENDENT_SEP):
        fragments = [_strip_fragment_marker(f) for f in group.split(_FRAGMENT_SEP)]
        fragments = [f for f in fragments if f]
        if fragments:
            groups.append(_FRAGMENT_SEP.join(fragments))
    return groups


def _strip_fragment_marker(fragment: str) -> str:
    """Drop a trailing (N)/(C) terminus marker so the bare ORF id can be looked up."""
    fragment = fragment.strip()
    if fragment.endswith(")") and "(" in fragment:
        fragment = fragment[: fragment.rindex("(")]
    return fragment.strip()


def read_sgd_feature_data(feature_file: Path) -> pd.DataFrame:
    """Read SGD's headerless SGD_features.tab into named columns."""
    return pd.read_csv(
        feature_file,
        sep="\t",
        header=None,
        names=_SGD_FEATURE_COLUMNS,
        dtype=str,
        skip_blank_lines=True,
    )


def build_sc_gene_info(sgd_features: pd.DataFrame) -> pd.DataFrame:
    """Index SGD ORF rows by systematic name, keeping standard name, qualifier and description."""
    orfs = sgd_features[sgd_features["feature_type"] == _SGD_ORF_TYPE].dropna(
        subset=["systematic_name"]
    )
    return orfs.set_index("systematic_name")[["standard_name", "qualifier", "description"]]


def read_ortholog_file(ortholog_file: Path) -> pd.DataFrame:
    """Read a PomBase curated_orthologs file (pombe id + tab + ortholog field, no header)."""
    return pd.read_csv(
        ortholog_file,
        sep="\t",
        header=None,
        names=["gene_systematic_id", "orthologs"],
        dtype=str,
    )


def read_gene_viability(viability_file: Path) -> pd.Series:
    """Read PomBase gene_viability.tsv into a systematic-id -> viability Series."""
    viability = pd.read_csv(
        viability_file,
        sep="\t",
        header=None,
        names=["gene_systematic_id", "viability"],
        dtype=str,
    )
    return viability.set_index("gene_systematic_id")["viability"]


def read_sgd_phenotype_data(phenotype_file: Path) -> pd.DataFrame:
    """Read SGD's headerless, ragged phenotype_data.tab, padding/truncating rows to the 14 known fields."""
    width = len(_SGD_PHENOTYPE_COLUMNS)
    rows = [
        (line.split("\t") + [None] * width)[:width]
        for line in phenotype_file.read_text().splitlines()
        if line.strip()
    ]
    return pd.DataFrame(rows, columns=_SGD_PHENOTYPE_COLUMNS, dtype="object")


def build_sc_essentiality(phenotype_data: pd.DataFrame) -> pd.DataFrame:
    """Derive per-gene null-mutant essentiality plus a per-label evidence count from SGD phenotypes."""
    viability = phenotype_data.loc[
        (phenotype_data["mutant_type"] == _NULL_MUTANT)
        & (phenotype_data["phenotype"].isin(_VIABILITY_PHENOTYPES)),
        ["feature_name", "phenotype"],
    ]
    if viability.empty:
        return pd.DataFrame(columns=["essentiality", "essentiality_evidence"])

    counts = (
        viability.groupby(["feature_name", "phenotype"]).size().unstack(fill_value=0)
    )
    for phenotype in _VIABILITY_PHENOTYPES:
        if phenotype not in counts.columns:
            counts[phenotype] = 0

    return pd.DataFrame(
        {
            "essentiality": counts.apply(_call_essentiality, axis=1),
            "essentiality_evidence": counts.apply(_format_evidence, axis=1),
        }
    )


def _call_essentiality(counts: pd.Series) -> str:
    """Label a gene from its inviable/viable record counts, flagging disagreement as conflicting."""
    has_inviable = counts[_INVIABLE] > 0
    has_viable = counts[_VIABLE] > 0
    if has_inviable and has_viable:
        return _CONFLICTING
    return _INVIABLE if has_inviable else _VIABLE


def _format_evidence(counts: pd.Series) -> str:
    """Render supporting record counts as 'inviable:3|viable:1', omitting zero-count labels."""
    return _INDEPENDENT_SEP.join(
        f"{phenotype}:{counts[phenotype]}"
        for phenotype in _VIABILITY_PHENOTYPES
        if counts[phenotype] > 0
    )


# =============================================================================
# ORTHOLOG ANNOTATION BLOCKS
# =============================================================================
def build_sc_ortholog_block(
    orthologs: pd.DataFrame,
    sc_gene_info: pd.DataFrame,
    sc_essentiality: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble per-ortholog S. cerevisiae id/name/qualifier/essentiality columns, positionally aligned."""
    records = {}
    for pombe_id, field in zip(orthologs["gene_systematic_id"], orthologs["orthologs"]):
        groups = parse_ortholog_field(field)
        records[pombe_id] = {
            "Sc_ortholog_id": _INDEPENDENT_SEP.join(groups),
            "Sc_ortholog_name": _join_per_group(groups, sc_gene_info, _lookup_standard_name),
            "Sc_ortholog_qualifier": _join_per_group(groups, sc_gene_info, _lookup_qualifier),
            "Sc_essentiality": _join_per_group(groups, sc_essentiality, _lookup_essentiality),
            "Sc_essentiality_evidence": _join_per_group(
                groups, sc_essentiality, _lookup_essentiality_evidence
            ),
            "Sc_description": _join_per_group(groups, sc_gene_info, _lookup_description),
            "Sc_ortholog_count": len(groups),
            "Sc_ortholog_raw": field if isinstance(field, str) else "",
        }

    block = pd.DataFrame.from_dict(records, orient="index")
    block.index.name = "gene_systematic_id"
    return block


def build_hs_ortholog_block(orthologs: pd.DataFrame) -> pd.DataFrame:
    """Assemble human ortholog symbols and count straight from PomBase's curated symbols."""
    records = {}
    for pombe_id, field in zip(orthologs["gene_systematic_id"], orthologs["orthologs"]):
        groups = parse_ortholog_field(field)
        records[pombe_id] = {
            "Hs_ortholog_symbol": _INDEPENDENT_SEP.join(groups),
            "Hs_ortholog_count": len(groups),
        }

    block = pd.DataFrame.from_dict(records, orient="index")
    block.index.name = "gene_systematic_id"
    return block


def _join_per_group(groups: list[str], lookup_table: pd.DataFrame, lookup) -> str:
    """Map each ortholog group through `lookup` and pipe-join, keeping fusion fragments in one group."""
    return _INDEPENDENT_SEP.join(
        _FRAGMENT_SEP.join(lookup(fragment, lookup_table) for fragment in group.split(_FRAGMENT_SEP))
        for group in groups
    )


def _lookup_standard_name(orf: str, sc_gene_info: pd.DataFrame) -> str:
    """Return an ORF's common name, falling back to its systematic id (~1300 ORFs are unnamed)."""
    if orf not in sc_gene_info.index:
        return orf
    name = sc_gene_info.loc[orf, "standard_name"]
    return orf if pd.isna(name) else str(name)


def _lookup_qualifier(orf: str, sc_gene_info: pd.DataFrame) -> str:
    """Return an ORF's SGD qualifier (Verified / Dubious / Uncharacterized)."""
    return _lookup_field(orf, sc_gene_info, "qualifier")


def _lookup_description(orf: str, sc_gene_info: pd.DataFrame) -> str:
    """Return an ORF's SGD functional description."""
    return _lookup_field(orf, sc_gene_info, "description")


def _lookup_essentiality(orf: str, sc_essentiality: pd.DataFrame) -> str:
    """Return an ORF's null-mutant essentiality call, or empty when SGD has no viability record."""
    return _lookup_field(orf, sc_essentiality, "essentiality")


def _lookup_essentiality_evidence(orf: str, sc_essentiality: pd.DataFrame) -> str:
    """Return the per-label record counts backing an ORF's essentiality call."""
    return _lookup_field(orf, sc_essentiality, "essentiality_evidence")


def _lookup_field(key: str, table: pd.DataFrame, column: str) -> str:
    """Look up one cell, returning an empty string for absent keys so alignment is preserved."""
    if key not in table.index:
        return ""
    value = table.loc[key, column]
    return "" if pd.isna(value) else str(value)


# =============================================================================
# POMBE-SIDE BLOCK
# =============================================================================
def build_pombe_block(
    gene_meta: pd.DataFrame,
    fypo_viability: pd.Series,
    deletion_library: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble pombe identity columns plus both essentiality sources (FYPO and deletion library).

    FYPO viability and the Hayles deletion library are kept as separate columns rather
    than merged: FYPO spans the whole genome but is ~61% "unknown", while the deletion
    library covers only ~4843 genes yet has a call for every one of them.
    """
    block = gene_meta.set_index("gene_systematic_id")[
        ["gene_name", "gene_product", "synonyms"]
    ].copy()
    block["gene_name"] = block["gene_name"].fillna(pd.Series(block.index, index=block.index))
    block["Sp_FYPO_viability"] = fypo_viability.reindex(block.index)

    for column, source in _DELETION_LIBRARY_COLUMNS.items():
        block[column] = (
            deletion_library.set_index("Systematic ID")[source].reindex(block.index)
            if source in deletion_library.columns
            else pd.NA
        )

    block.index.name = "gene_systematic_id"
    return block


# =============================================================================
# FUNCTIONAL ANNOTATION BLOCKS
# =============================================================================
def build_complex_block(complex_annotation: pd.DataFrame) -> pd.DataFrame:
    """Collapse macromolecular-complex membership into one pipe-joined column per gene."""
    return _collapse_unique(
        complex_annotation, key="systematic_id", value="GO_term_name", column="complex"
    )


def build_go_slim_block(
    ns2slim_assoc: dict[str, dict[str, set[str]]],
    slim_term_names: dict[str, str],
) -> pd.DataFrame:
    """Turn goatools' namespace->gene->slim-term-ids mapping into one readable column per namespace."""
    columns = {}
    for namespace, column in _GO_SLIM_NAMESPACES.items():
        gene2terms = ns2slim_assoc.get(namespace, {})
        columns[column] = {
            gene: _INDEPENDENT_SEP.join(
                sorted(slim_term_names[term] for term in terms if term in slim_term_names)
            )
            for gene, terms in gene2terms.items()
        }

    block = pd.DataFrame(columns).reindex(columns=list(_GO_SLIM_NAMESPACES.values()))
    block.index.name = "gene_systematic_id"
    return block


def _collapse_unique(
    table: pd.DataFrame, *, key: str, value: str, column: str
) -> pd.DataFrame:
    """Group by `key` and pipe-join unique `value`s in first-seen order."""
    collapsed = (
        table.groupby(key)[value]
        .apply(lambda values: _INDEPENDENT_SEP.join(dict.fromkeys(values.dropna())))
        .to_frame(column)
    )
    collapsed.index.name = "gene_systematic_id"
    return collapsed


# =============================================================================
# gRNA-LEVEL DEPLETION PARAMETERS
# =============================================================================
def build_grna_block(grna_parameters: pd.DataFrame) -> pd.DataFrame:
    """Extract per-gene gRNA-level depletion rate and lag from the curated fitted-parameters table.

    Columns are prefixed `gRNA_` because these are NOT the gene-level DR/DL that
    clustering tables carry: they come from a single representative gRNA fit rather
    than a gene-level aggregate fit, and the two disagree substantially (DR
    correlates ~0.92 but DL only ~0.55 across ~4.5k shared genes).
    """
    available = {
        source: target
        for source, target in _GRNA_DEPLETION_COLUMNS.items()
        if source in grna_parameters.columns
    }
    if len(available) < len(_GRNA_TARGET_COLUMNS):
        raise KeyError(
            "gRNA parameter table has no depletion columns: expected legacy um/lam "
            f"or DR/DL, found {list(grna_parameters.columns)}"
        )

    genes = grna_parameters[_GRNA_GENE_COLUMN]
    if genes.duplicated().any():
        duplicates = genes[genes.duplicated()].unique().tolist()
        raise ValueError(
            f"gRNA parameter table has duplicate gene ids, which would fan out rows: {duplicates[:10]}"
        )

    block = grna_parameters.set_index(_GRNA_GENE_COLUMN)[list(available)].rename(columns=available)
    block = block[_GRNA_TARGET_COLUMNS]
    block.index.name = "gene_systematic_id"
    return block


# =============================================================================
# ASSEMBLING THE FULL REFERENCE
# =============================================================================
def assemble_annotation_reference(
    pombe_block: pd.DataFrame, annotation_blocks: list[pd.DataFrame]
) -> pd.DataFrame:
    """Left-join annotation blocks onto the pombe gene set, which alone defines the row set."""
    reference = pombe_block.copy()
    for block in annotation_blocks:
        if block.index.has_duplicates:
            duplicates = block.index[block.index.duplicated()].unique().tolist()
            raise ValueError(
                f"Annotation block has duplicate gene ids, which would fan out rows: {duplicates[:10]}"
            )
        reference = reference.join(block, how="left")

    # Genes missing from a block introduce NaN, which promotes int count columns to
    # float and renders as "1.0" in the exported table. Nullable Int64 keeps them
    # integral while still allowing a blank.
    for column in reference.columns:
        if column.endswith(_COUNT_COLUMN_SUFFIX):
            reference[column] = reference[column].astype("Int64")

    reference.index.name = "gene_systematic_id"
    return reference


# =============================================================================
# JOINING ANNOTATION ONTO A USER TABLE
# =============================================================================
def annotate_table(
    table: pd.DataFrame,
    annotation_reference: pd.DataFrame,
    *,
    gene_column: str,
    columns: list[str] | None = None,
    drop_unmatched: bool = False,
) -> pd.DataFrame:
    """Append annotation columns to `table`, matching `gene_column` against the reference index."""
    _require_gene_column(table, gene_column)

    annotation = annotation_reference
    if columns is not None:
        missing = [column for column in columns if column not in annotation.columns]
        if missing:
            raise KeyError(
                f"Annotation reference has no column(s) {missing}. "
                f"Available columns: {list(annotation.columns)}"
            )
        annotation = annotation[columns]

    # merge (not join) so repeated gene ids annotate each row in place instead of
    # fanning out, and so a caller column of the same name is suffixed rather than
    # silently overwritten.
    annotated = table.merge(
        annotation,
        how="inner" if drop_unmatched else "left",
        left_on=gene_column,
        right_index=True,
        suffixes=("", "_annotation"),
        sort=False,
    )
    return annotated.reset_index(drop=True)


def _require_gene_column(table: pd.DataFrame, gene_column: str) -> None:
    """Raise a KeyError naming the available columns, rather than letting pandas raise a bare one."""
    if gene_column not in table.columns:
        raise KeyError(
            f"Gene column {gene_column!r} not found in the input table. "
            f"Available columns: {list(table.columns)}"
        )


def summarise_match(
    table: pd.DataFrame, annotation_reference: pd.DataFrame, *, gene_column: str
) -> tuple[int, list[str]]:
    """Return the matched row count and the sorted unique gene ids that found no annotation."""
    _require_gene_column(table, gene_column)
    genes = table[gene_column]
    matched_mask = genes.isin(annotation_reference.index)
    unmatched = sorted(genes[~matched_mask].dropna().astype(str).unique())
    return int(matched_mask.sum()), unmatched
