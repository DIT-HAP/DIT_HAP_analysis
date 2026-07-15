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
    slim_dag = {term: dag[term] for term in slim_terms}

    ns2assoc = objanno.get_ns2assc(**kwargs)
    gene2go = objanno.get_id2gos_nss(**kwargs)
    go2genes = objanno.get_id2gos_nss(go2geneids=True, **kwargs)

    return dag, objanno, ns2assoc, gene2go, go2genes, slim_dag
