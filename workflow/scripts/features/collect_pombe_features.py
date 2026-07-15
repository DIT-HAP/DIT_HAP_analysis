#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pombe Coding Gene Feature Collection
=======================================

Assembles a per-coding-gene feature matrix from 15 heterogeneous sources —
PomBase annotations, AlphaFold structures, BioGrid, Ensembl paralogs, and
8 literature supplementary tables — reproducing
DIT_HAP_pipeline/workflow/notebooks/pombe_feature_collection.ipynb exactly.
Dataset-independent: depends only on the PomBase reference version, not on
any DIT-HAP sequencing project (design doc §8).

Input
-----
- A PomBase version directory (genome FASTA/GFF3, gene metadata, protein
  features, ontology OBO/GAF files, curated orthologs)
- An AlphaFold structure directory (.pdb.gz files)
- 8 literature supplementary tables (xlsx/xls)
- Curated deletion-library and essentiality-verification tables
- BioGrid interaction table, Ensembl paralog export

Output
------
- A tab-separated per-gene feature matrix (one row per coding gene)
- A tab-separated codon (anti-codon) usage matrix

Usage
-----
    python collect_pombe_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --alphafold-dir /path/to/AlphaFold_Dataset \\
        --literature-dir resources/literature \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --essentiality-verification-csv resources/curated/essentiality_verification.csv \\
        --biogrid-tsv resources/external/biogrid/BIOGRID-....tab3.txt \\
        --ensembl-paralogs-tsv resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --output results/features/2025-10-01/pombe_coding_gene_protein_features.tsv \\
        --codon-usage-output results/features/2025-10-01/codon_usage_matrix.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import gffutils
from loguru import logger
from tqdm import tqdm

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data
from workflow.src.features.genome import DNA_level_features, PombaseGenomeConfig, calculate_anticodon_usage_matrix
from workflow.src.features.protein import extract_protein_features_from_peptide_sequence, pLDDT_statistics_report
from workflow.src.gene_ids import update_sysIDs

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Amino-acid-composition columns kept from extract_protein_features_from_peptide_sequence
# (notebook cell 37's selected_protein_features_from_peptide list)
SELECTED_PEPTIDE_FEATURE_COLUMNS = [
    "aromaticity", "aliphatic_index", "gravy", "flexibility", "instability_index",
    "aa_percent_Ala", "aa_percent_Cys", "aa_percent_Asp", "aa_percent_Glu", "aa_percent_Phe",
    "aa_percent_Gly", "aa_percent_His", "aa_percent_Ile", "aa_percent_Lys", "aa_percent_Leu",
    "aa_percent_Met", "aa_percent_Asn", "aa_percent_Pro", "aa_percent_Gln", "aa_percent_Arg",
    "aa_percent_Ser", "aa_percent_Thr", "aa_percent_Val", "aa_percent_Trp", "aa_percent_Tyr",
]
# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths for the feature collection pipeline."""
    pombase_dir: Path
    alphafold_dir: Path
    literature_dir: Path
    deletion_library_xlsx: Path
    essentiality_verification_csv: Path
    biogrid_tsv: Path
    ensembl_paralogs_tsv: Path
    output_features: Path
    output_codon_usage: Path

    def __post_init__(self) -> None:
        """Validate all input paths exist, then ensure output directories exist."""
        required_inputs = [
            self.pombase_dir, self.alphafold_dir, self.literature_dir,
            self.deletion_library_xlsx, self.essentiality_verification_csv,
            self.biogrid_tsv, self.ensembl_paralogs_tsv,
        ]
        for path in required_inputs:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_features.parent.mkdir(parents=True, exist_ok=True)
        self.output_codon_usage.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used throughout for update_sysIDs()."""
        return self.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv"


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_gene_meta(gene_meta_file: Path) -> tuple[pd.DataFrame, dict]:
    """Load gene metadata and build a uniprot_id -> gene_systematic_id map."""
    gene_meta = pd.read_csv(gene_meta_file, sep="\t")
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])
    uniprot2id = dict(zip(gene_meta["uniprot_id"], gene_meta["gene_systematic_id"]))
    return gene_meta, uniprot2id


@logger.catch
def get_ortholog_counts(ortholog_file: Path) -> pd.Series:
    """Count pipe-separated orthologs per gene from a PomBase curated_orthologs file."""
    ortholog_df = pd.read_csv(
        ortholog_file, sep="\t", index_col=0, header=None,
        names=["gene_systematic_id", "orthologs"], na_values="NONE",
    )
    ortholog_df.index = ortholog_df.index.str.split("(").str[0]
    return ortholog_df["orthologs"].str.split("|").apply(lambda x: len(x) if isinstance(x, list) else 0)


@logger.catch
def collect_dna_level_features(db: gffutils.FeatureDB, genome_cfg: PombaseGenomeConfig) -> tuple[pd.DataFrame, list[str]]:
    """Compute DNA-level features for every mRNA, and return the list of coding gene IDs."""
    mRNAs = list(db.features_of_type("mRNA"))
    records = [
        DNA_level_features.from_gffutils_feature(mRNA, db, genome_cfg)
        for mRNA in tqdm(mRNAs, desc="DNA-level features")
    ]
    df = pd.DataFrame(records)
    coding_genes = df["Gene_id"].unique().tolist()
    return df, coding_genes
@logger.catch
def collect_rna_level_features(literature_dir: Path, gene_meta_file: Path, coding_genes: list[str]) -> pd.DataFrame:
    """Assemble mRNA abundance (Marguerat 2012) and mRNA kinetics (Harigaya 2016) features."""
    abundance = pd.read_excel(
        literature_dir / "margueratQuantitativeAnalysisFission2012.xlsx",
        sheet_name="Table_S2", comment="#",
    ).set_index("Systematic.name")
    abundance = abundance[["MM1.tot.cpc_ex", "MM2.tot.cpc_ex", "MN1.tot.cpc_ex", "MN2.tot.cpc_ex"]].copy()
    abundance.columns = pd.MultiIndex.from_tuples(
        [
            ("EMM_Proliferating_Cell_RNA_Abundance", "replicate1"),
            ("EMM_Proliferating_Cell_RNA_Abundance", "replicate2"),
            ("EMM_Nitrogen_Starved_Cell_RNA_Abundance", "replicate1"),
            ("EMM_Nitrogen_Starved_Cell_RNA_Abundance", "replicate2"),
        ],
        names=["Condition", "Replicate"],
    )
    mean_ = abundance.T.groupby(level="Condition").mean().T
    std_ = abundance.T.groupby(level="Condition").std().T
    cv_ = std_ / mean_
    abundance_stats = pd.concat([mean_, std_, cv_], axis=1, keys=["mean", "std", "cv"])
    abundance_stats.index = update_sysIDs(abundance_stats.index.tolist(), gene_meta_file)
    abundance_stats = abundance_stats[abundance_stats.index.isin(coding_genes)].copy().dropna().round(3)
    abundance_stats.columns = ["_".join(col).strip() for col in abundance_stats.columns.values]
    abundance_stats = (
        abundance_stats.rename_axis("gene_systematic_id")
        .reset_index()
        .drop_duplicates(subset=["gene_systematic_id"])
        .set_index("gene_systematic_id")
    )

    kinetics = pd.read_excel(literature_dir / "harigayaAnalysisAssociationCodon2016.xls", sheet_name="Table")
    kinetics = kinetics[["Gene ID", "tAIg", "HL - Mata (5)", "SR - Mata (5)"]].set_index("Gene ID")
    kinetics.columns = ["tAIg", "mRNA_half_life_minutes", "mRNA_synthesis_rate_per_minute"]
    kinetics.index = update_sysIDs(kinetics.index.tolist(), gene_meta_file)
    kinetics = kinetics[kinetics.index.isin(coding_genes)].copy().dropna().round(3)
    kinetics = (
        kinetics.rename_axis("gene_systematic_id")
        .reset_index()
        .drop_duplicates(subset=["gene_systematic_id"])
        .set_index("gene_systematic_id")
    )

    return pd.concat([abundance_stats, kinetics], axis=1, join="outer")
@logger.catch
def collect_protein_level_features(
    pombase_dir: Path,
    alphafold_dir: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    protein_metadata: pd.DataFrame,
    uniprot2id: dict,
    coding_genes: list[str],
) -> pd.DataFrame:
    """Assemble peptide-sequence, abundance, turnover, pLDDT, and PFAM-domain protein features."""
    peptide_features = extract_protein_features_from_peptide_sequence(
        pombase_dir / "genome_sequence_and_features" / "peptide.fa"
    )

    gene_abundance = pd.read_csv(pombase_dir / "RNA_metadata" / "quantitative_gene_expression.tsv", sep="\t")
    proliferating = gene_abundance.query(
        "reference == 'PMID:23101633' and type == 'protein' and condition == 'glucose MM,standard temperature'"
    )[["gene_systematic_id", "copies_per_cell"]].dropna().astype({"copies_per_cell": float})
    quiescent = gene_abundance.query(
        "reference == 'PMID:23101633' and type == 'protein' and condition == 'glucose MM,nitrogen absent,standard temperature'"
    )[["gene_systematic_id", "copies_per_cell"]].dropna().astype({"copies_per_cell": float})
    protein_abundance = pd.merge(
        proliferating, quiescent, on="gene_systematic_id",
        suffixes=("_EMM_Proliferating_Cell", "_EMMN_Quiescent_Cell"),
    ).set_index("gene_systematic_id")

    protein_kinetics = pd.read_excel(
        literature_dir / "christianoGlobalProteomeTurnover2014.xlsx", na_values=["n.d."]
    ).dropna(subset=["Degradation rates (min-1)", "t1/2 (min)"]).rename(
        columns={"t1/2 (min)": "protein_half_life_minutes"}
    )
    protein_kinetics["ENSG"] = protein_kinetics["ENSG"].fillna(protein_kinetics["Gene name"])
    protein_kinetics = protein_kinetics[["ENSG", "protein_half_life_minutes"]].set_index("protein_half_life_minutes")
    protein_kinetics = protein_kinetics["ENSG"].str.split(";").explode().reset_index()
    protein_kinetics["gene_systematic_id"] = update_sysIDs(protein_kinetics["ENSG"].tolist(), gene_meta_file)
    protein_kinetics = protein_kinetics.drop_duplicates(subset=["gene_systematic_id"])

    pLDDTs = pLDDT_statistics_report(alphafold_dir, structure_format="pdb.gz")
    pLDDTs["Systematic_ID"] = pLDDTs["uniprot_id"].map(uniprot2id)

    protein_domains = pd.read_csv(pombase_dir / "Protein_features" / "protein_families_and_domains.tsv", sep="\t")
    pfam_domain_counts = (
        protein_domains.query("database == 'PFAM'").groupby("systematic_id").size().rename("PFAM_domain_count")
    )

    merged = (
        protein_metadata
        .merge(peptide_features.set_index("Gene_id")[SELECTED_PEPTIDE_FEATURE_COLUMNS], left_index=True, right_index=True, how="outer")
        .merge(protein_abundance, left_index=True, right_index=True, how="outer")
        .merge(protein_kinetics.set_index("gene_systematic_id")[["protein_half_life_minutes"]], left_index=True, right_index=True, how="outer")
        .merge(pLDDTs.set_index("Systematic_ID")[["mean_pLDDT", "std_pLDDT", "cv_pLDDT"]], left_index=True, right_index=True, how="outer")
        .join(pfam_domain_counts)
    )
    merged = merged[merged.index.isin(coding_genes)].copy()
    merged["PFAM_domain_count"] = merged["PFAM_domain_count"].fillna(0).astype(int)
    return merged
@logger.catch
def collect_evolutionary_level_features(
    pombase_dir: Path,
    ensembl_paralogs_tsv: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    coding_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble ortholog/paralog counts, evolutionary rate, and phyloP/divergence scores.

    Returns (evolutionary_features_df, phyloP_and_divergence_df) — the second is reused
    by collect_phenotype_level_features for transposon-insertion-density features that
    live in the same source table (grechFitnessLandscapeFission2019.xlsx).
    """
    orthologs_dir = pombase_dir / "curated_orthologs"
    num_japonicus = get_ortholog_counts(orthologs_dir / "pombe_japonicus_orthologs.txt")
    num_cerevisiae = get_ortholog_counts(orthologs_dir / "pombe_cerevisiae_orthologs.txt")
    num_human = get_ortholog_counts(orthologs_dir / "pombe_human_orthologs.txt")

    pombe_paralogs = pd.read_csv(ensembl_paralogs_tsv, sep="\t")
    paralog_count = (
        pombe_paralogs.query("`Gene stable ID` in @coding_genes")
        .groupby(["Gene stable ID", "Gene name"])
        .apply(lambda sub_df: sub_df.shape[0], include_groups=False)
        .to_frame("paralog_count")
    )
    paralog_count["gene_systematic_id"] = update_sysIDs(
        paralog_count.index.get_level_values("Gene stable ID").tolist(), gene_meta_file
    )

    evolutionary_rate = pd.read_excel(
        literature_dir / "rhindComparativeFunctionalGenomics2011.xls", sheet_name="S30", skiprows=[0, 1]
    )
    evolutionary_rate["gene_systematic_id"] = update_sysIDs(
        evolutionary_rate["Genes"].str.split(";").str[0].tolist(), gene_meta_file
    )

    phyloP_and_divergence = pd.read_excel(
        literature_dir / "grechFitnessLandscapeFission2019.xlsx", sheet_name="Table 2", skiprows=list(range(14))
    ).drop_duplicates(subset=["gene"])
    phyloP_and_divergence["gene_systematic_id"] = update_sysIDs(phyloP_and_divergence["gene"].tolist(), gene_meta_file)

    evolutionary_df = pd.concat(
        [
            num_japonicus.rename("japonicus_ortholog_count"),
            num_cerevisiae.rename("cerevisiae_ortholog_count"),
            num_human.rename("human_ortholog_count"),
            paralog_count.set_index("gene_systematic_id")[["paralog_count"]],
            evolutionary_rate.set_index("gene_systematic_id")[["Rate"]].rename(columns={"Rate": "evolutionary_rate"}),
            phyloP_and_divergence.set_index("gene_systematic_id")[
                ["mean.phylop", "diversity.S", "diversity.Pi", "diversity.Theta", "diversity.Tajima_D"]
            ],
        ],
        join="outer", axis=1,
    )
    evolutionary_df = evolutionary_df[evolutionary_df.index.isin(coding_genes)].copy()
    evolutionary_df["paralog_count"] = evolutionary_df["paralog_count"].fillna(0).astype(int)

    return evolutionary_df, phyloP_and_divergence
@logger.catch
def collect_network_level_features(
    pombase_dir: Path,
    biogrid_tsv: Path,
    gene2go: dict,
    coding_genes: list[str],
) -> pd.DataFrame:
    """Assemble GO term richness and BioGrid PPI/GI degree.

    NOTE: PPI_degree/GI_degree only count rows grouped by "Systematic Name Interactor A" —
    a gene that appears solely as "Interactor B" in BioGrid gets degree 0 here, even if it
    has real interactions. This is the original notebook's behavior (undocumented upstream,
    not a bug introduced by this port) and is preserved rather than symmetrized, since fixing
    it would change every downstream degree value against the reference this task verifies.
    """
    go_richness = {gene: len(set(terms)) for gene, terms in gene2go.items()}
    go_richness_df = pd.DataFrame.from_dict(go_richness, orient="index", columns=["GO_term_richness"])

    biogrid_data = pd.read_csv(biogrid_tsv, sep="\t")
    PPI_and_GI = biogrid_data[
        [
            "Systematic Name Interactor A", "Systematic Name Interactor B",
            "Official Symbol Interactor A", "Official Symbol Interactor B",
            "Experimental System Type",
        ]
    ].drop_duplicates()
    PPI = PPI_and_GI.query("`Experimental System Type` == 'physical'")
    GI = PPI_and_GI.query("`Experimental System Type` == 'genetic'")
    PPI_degrees = PPI.groupby("Systematic Name Interactor A").size().rename("PPI_degree")
    GI_degrees = GI.groupby("Systematic Name Interactor A").size().rename("GI_degree")

    network_df = pd.concat([go_richness_df, PPI_degrees, GI_degrees], join="outer", axis=1)
    network_df = network_df[network_df.index.isin(coding_genes)].copy()
    network_df = network_df.fillna(0).astype({"PPI_degree": int, "GI_degree": int})
    return network_df
@logger.catch
def collect_phenotype_level_features(
    pombase_dir: Path,
    deletion_library_xlsx: Path,
    essentiality_verification_csv: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    coding_genes: list[str],
    phyloP_and_divergence: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble FYPO viability, deletion-library essentiality, bar-seq fitness, transposon
    insertion density, and CRISPRi growth phenotypes.
    """
    FYPO_viability = pd.read_csv(
        pombase_dir / "Gene_metadata" / "gene_viability.tsv", sep="\t",
        header=None, names=["gene_systematic_id", "FYPOviability"],
    ).set_index("gene_systematic_id")

    DeletionLibrary_essentiality = pd.read_excel(deletion_library_xlsx)[
        ["Updated_Systematic_ID", "Gene dispensability. This study", "Category"]
    ].set_index("Updated_Systematic_ID").rename(
        columns={"Gene dispensability. This study": "DeletionLibrary_essentiality", "Category": "DeletionLibrary_category"}
    )

    revised_essentiality_map = (
        pd.read_csv(essentiality_verification_csv)[["systematic_id", "verification_essentiality"]]
        .set_index("systematic_id")["verification_essentiality"]
        .to_dict()
    )
    # Updated_essentiality intentionally ends up with TWO relevant columns:
    # "DeletionLibrary_essentiality" (carried over from the .copy() below) and the new
    # "RevisedDeletionLibrary_essentiality". Concatenating it against DeletionLibrary_essentiality
    # (Step 9) therefore produces a duplicate "DeletionLibrary_essentiality" column in the final
    # matrix — see this task's header note; both are kept to match the reference output byte-for-byte.
    Updated_essentiality = DeletionLibrary_essentiality[["DeletionLibrary_essentiality"]].copy()
    Updated_essentiality["RevisedDeletionLibrary_essentiality"] = Updated_essentiality.apply(
        lambda row: revised_essentiality_map.get(row.name, row["DeletionLibrary_essentiality"]), axis=1
    )

    bar_seq_fitness = pd.read_excel(literature_dir / "comp_fitness_QianWenFeng_Koch-1.xlsx").rename(
        columns={"yes": "Barseq_from_dulab", "SM fitness defect from Koch et al": "Barseq_from_koch"}
    ).dropna(subset=["Barseq_from_dulab", "Barseq_from_koch"])
    bar_seq_fitness["gene_systematic_id"] = update_sysIDs(bar_seq_fitness["gene"].tolist(), gene_meta_file)

    ins_density = pd.read_excel(literature_dir / "guoIntegrationProfilingGene2013.xls", sheet_name="TableS2")
    ins_density["gene_systematic_id"] = update_sysIDs(
        ins_density["Gene name"].str.strip().apply(lambda row: sorted(row.split(" "))[0]).tolist(), gene_meta_file
    )
    ins_density = ins_density.drop_duplicates(subset=["gene_systematic_id"])

    ins_grech = phyloP_and_divergence[["gene_systematic_id", "ipkm", "uipkm", "Malecki2016.KO.colony.size"]].copy().rename(
        columns={"Malecki2016.KO.colony.size": "colony_size_Malecki2016"}
    )

    CRISPRi_data = pd.read_excel(literature_dir / "ishikawaArrayedCRISPRiLibrary2024.xlsx").iloc[1:].dropna(
        subset=["Max Growth Rate", "Colony Formation"]
    )
    CRISPRi_data["gene_systematic_id"] = update_sysIDs(CRISPRi_data["Systematic ID"].tolist(), gene_meta_file)

    phenotype_df = pd.concat(
        [
            FYPO_viability,
            DeletionLibrary_essentiality,
            Updated_essentiality,
            bar_seq_fitness.set_index("gene_systematic_id")[["Barseq_from_dulab", "Barseq_from_koch"]],
            ins_density.set_index("gene_systematic_id")[["Integration density, in-vivo (integrations/kb/million inserts)"]],
            ins_grech.set_index("gene_systematic_id")[["ipkm", "uipkm", "colony_size_Malecki2016"]],
            CRISPRi_data.set_index("gene_systematic_id")[["Max Growth Rate", "Colony Formation"]],
        ],
        join="outer", axis=1,
    )
    return phenotype_df[phenotype_df.index.isin(coding_genes)].copy()
@logger.catch
def merge_all_features(
    dna_df: pd.DataFrame,
    rna_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    evolutionary_df: pd.DataFrame,
    network_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    gene_meta: pd.DataFrame,
) -> pd.DataFrame:
    """Outer-join all six feature groups on gene_systematic_id and fill category-column NAs."""
    pombe_features = pd.concat(
        [
            dna_df.query("Primary_candidate == True").set_index("Gene_id"),
            rna_df,
            protein_df,
            evolutionary_df,
            network_df,
            phenotype_df,
        ],
        join="outer", axis=1,
    ).rename_axis("gene_systematic_id")

    pombe_features = gene_meta[["gene_systematic_id", "gene_name"]].merge(
        pombe_features.reset_index(), on="gene_systematic_id", how="right"
    )
    pombe_features[["GO_term_richness", "PPI_degree", "GI_degree"]] = (
        pombe_features[["GO_term_richness", "PPI_degree", "GI_degree"]].fillna(0).astype(int)
    )
    pombe_features["DeletionLibrary_essentiality"] = pombe_features["DeletionLibrary_essentiality"].fillna("Not_determined")
    pombe_features["DeletionLibrary_category"] = pombe_features["DeletionLibrary_category"].fillna("Not_determined")
    pombe_features["RevisedDeletionLibrary_essentiality"] = pombe_features["RevisedDeletionLibrary_essentiality"].fillna("Not_determined")
    return pombe_features
@logger.catch
def run_feature_collection(config: InputOutputConfig) -> pd.DataFrame:
    """Execute the full 6-group feature collection pipeline and write both output files."""
    logger.info(f"Building gffutils DB from {config.pombase_dir}")
    genome_dir = config.pombase_dir / "genome_sequence_and_features"
    genome_cfg = PombaseGenomeConfig.from_pombase_dir(config.pombase_dir)
    db = gffutils.create_db(genome_cfg.gff3_file, genome_cfg.database_file, force=True, merge_strategy="create_unique")
    db = gffutils.FeatureDB(genome_cfg.database_file)

    gene_meta, uniprot2id = load_gene_meta(config.gene_meta_file)

    logger.info("Collecting DNA-level features")
    dna_df, coding_genes = collect_dna_level_features(db, genome_cfg)

    logger.info("Writing codon usage matrix")
    codon_usage_matrix = calculate_anticodon_usage_matrix(db, genome_cfg)
    codon_usage_matrix.to_csv(config.output_codon_usage, sep="\t", index=True)

    logger.info("Collecting RNA-level features")
    rna_df = collect_rna_level_features(config.literature_dir, config.gene_meta_file, coding_genes)

    logger.info("Collecting protein-level features")
    protein_meta = pd.read_csv(config.pombase_dir / "Protein_features" / "peptide_stats.tsv", sep="\t", index_col=0)
    protein_df = collect_protein_level_features(
        config.pombase_dir, config.alphafold_dir, config.literature_dir,
        config.gene_meta_file, protein_meta, uniprot2id, coding_genes,
    )

    logger.info("Collecting evolutionary-level features")
    evolutionary_df, phyloP_and_divergence = collect_evolutionary_level_features(
        config.pombase_dir, config.ensembl_paralogs_tsv, config.literature_dir, config.gene_meta_file, coding_genes,
    )

    logger.info("Loading GO ontology data for network-level features")
    ontology_cfg = OntologyDataConfig(
        ontology_obo=config.pombase_dir / "ontologies_and_associations" / "go-basic.obo",
        ontology_association_gaf=config.pombase_dir / "ontologies_and_associations" / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            config.pombase_dir / "ontologies_and_associations" / "bp_go_slim_terms.tsv",
            config.pombase_dir / "ontologies_and_associations" / "mf_go_slim_terms.tsv",
            config.pombase_dir / "ontologies_and_associations" / "cc_go_slim_terms.tsv",
        ],
    )
    _, _, _, gene2go, _, _ = load_ontology_data(
        ontology_cfg.load_data(),
        relationships={"is_a", "part_of"}, propagate_counts=True, load_obsolete=False, prt=None,
    )

    logger.info("Collecting network-level features")
    network_df = collect_network_level_features(config.pombase_dir, config.biogrid_tsv, gene2go, coding_genes)

    logger.info("Collecting phenotype-level features")
    phenotype_df = collect_phenotype_level_features(
        config.pombase_dir, config.deletion_library_xlsx, config.essentiality_verification_csv,
        config.literature_dir, config.gene_meta_file, coding_genes, phyloP_and_divergence,
    )

    logger.info("Merging all feature groups")
    pombe_features = merge_all_features(dna_df, rna_df, protein_df, evolutionary_df, network_df, phenotype_df, gene_meta)

    pombe_features.to_csv(config.output_features, sep="\t", index=False)
    logger.success(f"Wrote {len(pombe_features)} gene records to {config.output_features}")

    return pombe_features
# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect pombe coding gene features from 15 sources")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--alphafold-dir", type=Path, required=True, help="AlphaFold structure directory (.pdb.gz files)")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion library categories xlsx")
    parser.add_argument("--essentiality-verification-csv", type=Path, required=True, help="Curated essentiality verification csv")
    parser.add_argument("--biogrid-tsv", type=Path, required=True, help="BioGrid interaction table")
    parser.add_argument("--ensembl-paralogs-tsv", type=Path, required=True, help="Ensembl paralog export table")
    parser.add_argument("--output", type=Path, required=True, dest="output_features", help="Output feature matrix path")
    parser.add_argument("--codon-usage-output", type=Path, required=True, dest="output_codon_usage", help="Output codon usage matrix path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate paths, run the feature collection pipeline, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            pombase_dir=args.pombase_dir,
            alphafold_dir=args.alphafold_dir,
            literature_dir=args.literature_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            essentiality_verification_csv=args.essentiality_verification_csv,
            biogrid_tsv=args.biogrid_tsv,
            ensembl_paralogs_tsv=args.ensembl_paralogs_tsv,
            output_features=args.output_features,
            output_codon_usage=args.output_codon_usage,
        )
        run_feature_collection(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
