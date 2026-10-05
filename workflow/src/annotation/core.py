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
- Curated deletion-library categories xlsx and essentiality-verification csv
- A dataset's gene-level fitting results (DR/DL) and the curated gRNA parameters

Output
------
- A DataFrame indexed by pombe systematic ID, one column per annotation field

Every annotation source is parsed by a `build_*_block` here — one builder per
source, all returning a block indexed by systematic_id — so the drivers
(build_annotation_reference.py, annotate_pombe_genes.py, build_annotated_workbook.py)
only choose which blocks to assemble and where to write them.

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

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from data_config import load_dataset_config
from io_table import read_file, read_parquet


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# feature_type values (as spelled in PomBase's gene_ids_and_details.parquet) that
# each gene_type selects; "all" means no filter. The domain of
# build_pombase_metadata_block's gene_type argument.
VALID_GENE_TYPES = {
    "all": None,  # No filtering
    "protein": ["protein"],
    "lncRNA": ["lncRNA gene"],
    "tRNA": ["tRNA gene"],
    "rRNA": ["rRNA gene"],
    "snoRNA": ["snoRNA gene"],
    "sncRNA": ["sncRNA gene"],
    "snRNA": ["snRNA gene"],
    "noncoding": ["lncRNA gene", "tRNA gene", "rRNA gene", "snoRNA gene", "sncRNA gene", "snRNA gene"],
}

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
# legacy um/lam; accept DR/DL too in case upstream renames them. The two spellings
# carry different sign conventions — see GRNA_DR_SIGN and build_grna_block. Output is
# prefixed gRNA_ to keep it distinct from the gene-level DR/DL in clustering tables.
_GRNA_GENE_COLUMN = "Systematic ID"
_GRNA_DEPLETION_COLUMNS = {
    "um": "gRNA_DR",
    "lam": "gRNA_DL",
    "DR": "gRNA_DR",
    "DL": "gRNA_DL",
}
_GRNA_TARGET_COLUMNS = ["gRNA_DR", "gRNA_DL"]

# The curated gRNA table (resources/curated/*_gRNA_HDdata_fitted_parameters.tsv) is
# frozen at the pre-2026-09-17 sign convention — positive = depleted — while
# upstream flipped DIT-HAP so negative DR is now the depleted end. Flip gRNA_DR on
# the way in, so this workbook's gRNA_DR and DR columns point the same way; read
# side by side with opposite signs, agreeing genes look contradictory.
GRNA_DR_SIGN = -1.0

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
    """Assemble per-ortholog S. cerevisiae id/name/essentiality/description columns, positionally aligned."""
    records = {}
    for pombe_id, field in zip(orthologs["gene_systematic_id"], orthologs["orthologs"]):
        groups = parse_ortholog_field(field)
        records[pombe_id] = {
            "Sc_ortholog_id": _INDEPENDENT_SEP.join(groups),
            "Sc_ortholog_name": _join_per_group(groups, sc_gene_info, _lookup_standard_name),
            "Sc_essentiality": _join_per_group(groups, sc_essentiality, _lookup_essentiality),
            "Sc_description": _join_per_group(groups, sc_gene_info, _lookup_description),
            "Sc_ortholog_count": len(groups),
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


def _lookup_description(orf: str, sc_gene_info: pd.DataFrame) -> str:
    """Return an ORF's SGD functional description."""
    return _lookup_field(orf, sc_gene_info, "description")


def _lookup_essentiality(orf: str, sc_essentiality: pd.DataFrame) -> str:
    """Return an ORF's null-mutant essentiality call, or empty when SGD has no viability record."""
    return _lookup_field(orf, sc_essentiality, "essentiality")


def _lookup_field(key: str, table: pd.DataFrame, column: str) -> str:
    """Look up one cell, returning an empty string for absent keys so alignment is preserved."""
    if key not in table.index:
        return ""
    value = table.loc[key, column]
    return "" if pd.isna(value) else str(value)


# =============================================================================
# POMBE-SIDE & EXPERIMENTAL BLOCKS
# =============================================================================
# Each builder turns one source file into one block indexed by systematic_id; the
# driver (workflow/scripts/annotate/build_annotation_reference.py) only decides
# which blocks exist and how they join. The block number in each docstring is the
# order the driver assembles them in.
def _build_block(
    table: pd.DataFrame, columns_map: dict[str, str], *, required: bool = False
) -> pd.DataFrame:
    """Select `columns_map`'s source columns, rename them, and index on systematic_id.

    Two kinds of source, told apart by `required`: in a PomBase export a column can
    legitimately be absent, so the block carries whatever exists (the column set has
    moved between releases); in a curated table we depend on, an absent column means
    the wrong file was passed, so it raises instead.
    """
    missing = [source for source in columns_map if source not in table.columns]
    if required and missing:
        raise ValueError(f"Source table is missing required column(s): {missing}")

    available = {source: target for source, target in columns_map.items() if source in table.columns}
    return table[list(available)].rename(columns=available).set_index("systematic_id")


@logger.catch(reraise=True)
def build_pombase_metadata_block(gene_ids_parquet: Path, gene_type: str) -> pd.DataFrame:
    """Build block 1: PomBase gene metadata from gene_ids_and_details.parquet.

    Extracts: systematic_id, name, gene_product, product, feature_type,
    characterisation_status, taxonomic_distribution, deletion_viability (renamed to FYPOviability).
    """
    logger.info("Reading PomBase gene metadata from parquet")
    block = _build_block(read_parquet(gene_ids_parquet), {
        "systematic_id": "systematic_id",
        "name": "gene_name",
        "gene_product": "gene_product",
        "product": "product",
        "feature_type": "feature_type",
        "characterisation_status": "characterisation_status",
        "taxonomic_distribution": "taxonomic_distribution",
        "deletion_viability": "FYPOviability",
    })

    # Apply gene type filter
    feature_types = VALID_GENE_TYPES.get(gene_type)
    if feature_types is not None and "feature_type" in block.columns:
        before_count = len(block)
        block = block[block["feature_type"].isin(feature_types)]
        logger.info(f"  Filtered {before_count:,} genes to {len(block):,} {gene_type} genes")

    # Fill missing gene names with systematic_id
    if "gene_name" in block.columns:
        block["gene_name"] = block["gene_name"].fillna(pd.Series(block.index, index=block.index))

    logger.info(f"  {len(block):,} pombe genes in metadata block")
    return block


@logger.catch(reraise=True)
def build_deletion_library_block(deletion_library_xlsx: Path) -> pd.DataFrame:
    """Build block 2: Deletion library categories.

    Extracts: Systematic ID, Gene dispensability. This study (renamed to deletion_essentiality),
    Category, Sub_category.
    """
    logger.info("Reading deletion library categories")
    block = _build_block(read_file(deletion_library_xlsx), {
        "Systematic ID": "systematic_id",
        "Gene dispensability. This study": "deletion_essentiality",
        "Category": "Category",
        "Sub_category": "Sub_category",
    })

    if "deletion_essentiality" in block.columns:
        block["deletion_essentiality"] = block["deletion_essentiality"].fillna("Not_determined")

    logger.info(f"  {len(block):,} genes in deletion library")
    return block


@logger.catch(reraise=True)
def build_verification_block(verification_csv: Path) -> pd.DataFrame:
    """Build block 2b: Essentiality verification phenotype (systematic_id, verification_phenotype)."""
    logger.info("Reading essentiality verification data")
    block = _build_block(
        read_file(verification_csv),
        {"systematic_id": "systematic_id", "verification_phenotype": "verification_phenotype"},
        required=True,
    )

    logger.info(f"  {len(block):,} genes with verification phenotype")
    return block


@logger.catch(reraise=True)
def build_gene_level_depletion_block(dataset_name: str) -> pd.DataFrame:
    """Build block 3: gene-level DR/DL from the given dataset's gene-level fitting results."""
    logger.info(f"Reading gene-level depletion from {dataset_name}")
    dataset_config = load_dataset_config(dataset_name)

    if dataset_config.gene_level is None:
        raise ValueError(f"Dataset {dataset_name} has no gene-level data (has_time_points=False)")

    # Prefixed with the dataset name, like gRNA_DR/gRNA_DL. The reference carries depletion
    # from several measurements, and a bare DR/DL collides by name with the DR/DL of whichever
    # per-dataset table the reference is later joined onto — annotate_table skips columns the
    # table already has, so that collision would pick the wrong dataset's numbers silently.
    # The gene column is "Systematic ID" (with space), not "systematic_id".
    block = _build_block(
        read_file(dataset_config.gene_level.fitting_results),
        {"Systematic ID": "systematic_id", "DR": f"{dataset_name}_DR", "DL": f"{dataset_name}_DL"},
        required=True,
    )

    logger.info(f"  {len(block):,} genes with gene-level DR/DL")
    return block


def add_gene_status_column(gene_ids: pd.Index, gene_ids_parquet: Path) -> pd.Series:
    """Build the gene_status column: each gene's presence and type in the CURRENT PomBase.

    Uses the FULL gene_ids_and_details.parquet (not filtered), so it correctly
    identifies genes of all types. Values are the gene's feature_type, or
    "not_in_pombase" when the gene has no record at all.
    """
    full_gene_meta = read_parquet(gene_ids_parquet)
    full_gene_meta = full_gene_meta.set_index("systematic_id")

    status = []
    for gene_id in gene_ids:
        if gene_id not in full_gene_meta.index:
            status.append("not_in_pombase")
        else:
            status.append(full_gene_meta.loc[gene_id, "feature_type"])

    return pd.Series(status, index=gene_ids, name="gene_status")


def report_gene_status_summary(reference: pd.DataFrame, filter_type: str) -> None:
    """Log the gene_status distribution, warning about genes the type filter did not expect."""
    status_counts = reference["gene_status"].value_counts()

    logger.info("Gene status summary:")
    for status, count in status_counts.items():
        logger.info(f"  {status}: {count:,} genes")

    # Highlight unexpected cases when filtering for specific type
    if filter_type != "all":
        not_in_pombase = (reference["gene_status"] == "not_in_pombase").sum()
        wrong_type = (
            (reference["gene_status"] != "not_in_pombase") &
            (reference["gene_status"] != filter_type)
        ).sum()

        if not_in_pombase > 0:
            logger.warning(f"  {not_in_pombase:,} genes not in current PomBase version")
        if wrong_type > 0:
            logger.warning(
                f"  {wrong_type:,} genes have different type than expected '{filter_type}'"
            )


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
    correlates ~0.92 but DL only ~0.55 across ~4.5k shared genes). gRNA_DR is
    sign-flipped on the way in WHEN the source still uses the legacy um/lam
    spelling (see GRNA_DR_SIGN), so it matches the DIT-HAP DR convention that
    upstream flipped on 2026-09-17.
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
    # Flip only the legacy um/lam spelling: that table is frozen at the pre-2026-09-17
    # convention (positive = depleted). A table shipping DR/DL comes from the post-flip
    # upstream, whose sign already matches, and flipping it again would silently invert
    # every value. See GRNA_DR_SIGN.
    if "um" in available:
        block["gRNA_DR"] = block["gRNA_DR"] * GRNA_DR_SIGN
    block.index.name = "gene_systematic_id"
    return block


# =============================================================================
# ASSEMBLING THE FULL REFERENCE
# =============================================================================
def assemble_annotation_reference_split(
    pombe_block: pd.DataFrame,
    experimental_blocks: list[pd.DataFrame],
    annotation_blocks: list[pd.DataFrame],
) -> pd.DataFrame:
    """Assemble annotation reference with different join strategies.

    experimental_blocks: Use outer join (include genes from any experimental source)
    annotation_blocks: Use left join (only annotate genes already in the reference)
    """
    reference = pombe_block.copy()

    # Outer join for experimental data blocks (deletion library, DR/DL, gRNA, etc.)
    for block in experimental_blocks:
        if block.index.has_duplicates:
            duplicates = block.index[block.index.duplicated()].unique().tolist()
            raise ValueError(
                f"Experimental block has duplicate gene ids, which would fan out rows: {duplicates[:10]}"
            )
        reference = reference.join(block, how="outer")

    # Left join for annotation blocks (orthologs, functional annotation)
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
