"""
Ontology Data Loading
=======================

Loads a GO/FYPO/MONDO-style OBO + GAF association pair into goatools objects.
Ported from the `OntologyDataConfig`/`load_ontology_data` slice of
`enrichment_functions.py` (DIT_HAP_pipeline) — only the subset
`pombe_feature_collection.ipynb` actually calls (GO term richness via
`gene2go`). The enrichment-study, STRING, and REVIGO functions in the
original module are intentionally not ported here; see Task 7's scope note
in docs/plans/2026-07-15-DIT-HAP-analysis-phase1-implementation.md.

Input
-----
- An OBO ontology file (e.g. go-basic.obo)
- A GAF-format association file (e.g. gene_ontology_annotation.gaf.tsv)
- One or more slim-term tables (Term, Description columns, no header)

Output
------
- OntologyData: validated file handles + concatenated slim term table
- load_ontology_data(...): goatools GODag/GafReader plus gene2go/go2genes dicts

Usage
-----
    from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data
    cfg = OntologyDataConfig(ontology_obo=..., ontology_association_gaf=..., slim_terms_table=[...])
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(cfg.load_data())

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag

# 4. Local Imports
from workflow.src.io import read_file

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Fixed provenance header for the reformatted GAF files. The original used
# date.today(), which broke byte-level reproducibility of the intermediate
# files; a fixed stamp keeps the reformatted GAF deterministic across runs.
_GAF_HEADER_DATE = "2025-09-01"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class OntologyData:
    """Validated ontology file paths plus the concatenated slim-term table."""
    ontology_obo_path: Path
    ontology_association_file: Path
    slim_term_dataframe: pd.DataFrame


@dataclass(kw_only=True, frozen=True)
class OntologyDataConfig:
    """Unvalidated ontology file paths, as registered in config/datasets.yaml-adjacent code."""
    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms_table: list[Path]

    def validate_paths(self) -> None:
        """Raise if the OBO, GAF, or any slim-term file is missing."""
        for file_path in [self.ontology_obo, self.ontology_association_gaf, *self.slim_terms_table]:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene ontology file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    def load_data(self) -> OntologyData:
        """Validate paths, concatenate slim-term tables, and return an OntologyData."""
        self.validate_paths()
        slim_dfs = [
            read_file(path, header=None, names=["Term", "Description"])
            for path in self.slim_terms_table
        ]
        slim_df = pd.concat(slim_dfs, ignore_index=True)
        return OntologyData(
            ontology_obo_path=self.ontology_obo,
            ontology_association_file=self.ontology_association_gaf,
            slim_term_dataframe=slim_df,
        )


# =============================================================================
# CORE LOGIC
# =============================================================================
def load_ontology_data(
    ontology_data: OntologyData, **kwargs
) -> tuple[GODag, GafReader, dict, dict, dict, dict]:
    """Load an OBO + GAF pair into a GODag/GafReader and derive gene2go/go2genes dicts."""
    try:
        dag = GODag(str(ontology_data.ontology_obo_path), optional_attrs=["def", "relationship"], load_obsolete=False)
    except KeyError:
        dag = GODag(str(ontology_data.ontology_obo_path), optional_attrs=["def"], load_obsolete=False)

    objanno = GafReader(str(ontology_data.ontology_association_file), godag=dag)

    slim_terms = ontology_data.slim_term_dataframe["Term"].to_list()
    slim_dag = {term: dag[term] for term in slim_terms if term in dag}

    ns2assoc = objanno.get_ns2assc(**kwargs)
    gene2go = objanno.get_id2gos_nss(**kwargs)
    go2genes = objanno.get_id2gos_nss(go2geneids=True, **kwargs)

    return dag, objanno, ns2assoc, gene2go, go2genes, slim_dag


@dataclass(kw_only=True, frozen=True)
class GeneMetaData:
    """Gene metadata joined with deletion-library essentiality, plus an id->name dict."""
    gene_info_with_essentiality: pd.DataFrame
    id2name: dict


@dataclass(kw_only=True, frozen=True)
class GeneMetaConfig:
    """Paths for gene metadata + curated deletion-library essentiality."""
    gene_IDs_names_products: Path
    deletion_library_essentiality: Path

    def validate_paths(self) -> None:
        """Raise if either metadata file is missing or is not a file."""
        for file_path in [self.gene_IDs_names_products, self.deletion_library_essentiality]:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene metadata file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    def load_data(self) -> GeneMetaData:
        """Load metadata, fill missing gene_name with the systematic id, and left-join essentiality."""
        gene_IDs_names_products = read_file(self.gene_IDs_names_products)
        deletion_library_essentiality = read_file(self.deletion_library_essentiality)

        gene_IDs_names_products["gene_name"] = gene_IDs_names_products["gene_name"].fillna(
            gene_IDs_names_products["gene_systematic_id"]
        )
        id2name = dict(
            zip(
                list(gene_IDs_names_products["gene_systematic_id"]),
                list(gene_IDs_names_products["gene_name"]),
            )
        )

        gene_info_with_essentiality = gene_IDs_names_products.merge(
            deletion_library_essentiality[
                [
                    "Systematic ID",
                    "Gene dispensability. This study",
                    "Deletion mutant phenotype description",
                    "Phenotypic classification used for analysis",
                    "Category",
                ]
            ],
            how="left",
            left_on="gene_systematic_id",
            right_on="Systematic ID",
        ).drop(columns=["Systematic ID"])

        return GeneMetaData(gene_info_with_essentiality=gene_info_with_essentiality, id2name=id2name)


# =============================================================================
# GAF REFORMATTING (FYPO / MONDO -> GO-style)
# =============================================================================
def assign_term_name(term_id: str, term_dag: GODag) -> str:
    """Return a term's name from the DAG, or a placeholder when the term is absent."""
    if term_id in term_dag:
        return term_dag[term_id].name
    return f"No record for {term_id}"


# GO-style GAF 2.2 column order shared by both reformatters.
_GAF_COLUMNS = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference",
    "Evidence", "With", "Aspect", "DB_Object_Name", "Synonym", "DB_Object_Type",
    "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID",
]


def _write_gaf(reformatted: pd.DataFrame, output_path: Path, url: str) -> Path:
    """Write a GO-style GAF 2.2 file: fixed header (no date.today stamp) + headerless rows."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(
            "!gaf-version: 2.2\n"
            "!generated-by: Yusheng Yang\n"
            f"!date-generated: {_GAF_HEADER_DATE}\n"
            f"!URL: {url}\n"
            "!contact: yangyusheng@nibs.ac.cn\n"
        )
    reformatted.to_csv(output_path, sep="\t", index=False, header=False, mode="a")
    return output_path


def format_phaf_file(fypo_obo_file: Path, phaf_file: Path, output_path: Path) -> Path:
    """Reformat a PomBase PHAF into a GO-style GAF (deletion/disruption alleles, standard condition only)."""
    phaf_dag = GODag(str(fypo_obo_file))
    phaf = pd.read_csv(phaf_file, sep="\t").query(
        "(`Allele type` == 'deletion' or `Allele type` == 'disruption') and Condition.str.contains('FYECO:0000005')"
    )
    phaf["DB"] = "PomBase"
    phaf["DB_Object_ID"] = phaf["Gene systematic ID"]
    phaf["DB_Object_Symbol"] = phaf["Gene symbol"]
    phaf["Qualifier"] = ""
    phaf["GO_ID"] = phaf["FYPO ID"]
    phaf["DB:Reference"] = phaf["Reference"]
    phaf["Evidence"] = phaf["Evidence"]
    phaf["With"] = ""
    phaf["Aspect"] = "FYPO"
    phaf["DB_Object_Name"] = phaf["FYPO ID"].apply(assign_term_name, term_dag=phaf_dag)
    phaf["Synonym"] = ""
    phaf["DB_Object_Type"] = "protein"
    phaf["Taxon"] = "taxon:4896"
    phaf["Date"] = phaf["Date"].str.replace("-", "")
    phaf["Assigned_By"] = phaf["#Database name"]
    phaf["Annotation_Extension"] = phaf["Extension"]
    phaf["Gene_Product_Form_ID"] = ""
    reformat_phaf = phaf[_GAF_COLUMNS].copy()
    return _write_gaf(
        reformat_phaf,
        output_path,
        "https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz",
    )


def format_mondo_gaf_file(mondo_obo_file: Path, mondo_gaf_file: Path, output_path: Path) -> Path:
    """Reformat a PomBase human-disease (MONDO) association into a GO-style GAF."""
    mondo_dag = GODag(str(mondo_obo_file))
    mondo = pd.read_csv(mondo_gaf_file, sep="\t")
    mondo["DB"] = "Pombase"
    mondo["DB_Object_ID"] = mondo["#gene_systematic_id"]
    mondo["DB_Object_Symbol"] = mondo["gene_name"]
    mondo["Qualifier"] = ""
    mondo["GO_ID"] = mondo["mondo_id"]
    mondo["DB:Reference"] = mondo["reference"]
    mondo["Evidence"] = ""
    mondo["With"] = ""
    mondo["Aspect"] = "MONDO"
    mondo["DB_Object_Name"] = mondo["mondo_id"].apply(assign_term_name, term_dag=mondo_dag)
    mondo["Synonym"] = ""
    mondo["DB_Object_Type"] = "protein"
    mondo["Taxon"] = "taxon:4896"
    mondo["Date"] = mondo["date"].fillna(_GAF_HEADER_DATE).str.replace("-", "")
    mondo["Assigned_By"] = "PomBase"
    mondo["Annotation_Extension"] = ""
    mondo["Gene_Product_Form_ID"] = ""
    reformat_mondo = mondo[_GAF_COLUMNS].copy()
    return _write_gaf(
        reformat_mondo,
        output_path,
        "https://www.pombase.org/data/annotations/Disease_associations/human_disease_association.tsv",
    )
