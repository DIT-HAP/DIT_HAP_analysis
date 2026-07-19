"""
Pombe Feature Assembly (per-level)
====================================

Per-level feature-assembly functions shared by the six `collect_*_features.py`
driver scripts and the `merge_features.py` merge script. Split out of the
former monolithic `collect_pombe_features.py` so each biological level
(DNA / RNA / protein / evolutionary / network / phenotype) can be a separate
Snakemake rule while the assembly logic lives in one importable place.

The DNA level is the "spine": it enumerates every coding gene, and every other
level filters its rows to that gene set. Drivers therefore read the DNA-level
table to recover `coding_genes` before assembling their own level.

Input
-----
- A PomBase version directory + literature tables (per level; see each function)
- The DNA-level table's Gene_id column, used as `coding_genes` by other levels

Output
------
- One DataFrame per level; `merge_all_features` outer-joins the six into the
  final per-gene feature matrix.

Usage
-----
    from workflow.src.features.assembly import collect_rna_level_features
    rna_df = collect_rna_level_features(literature_dir, gene_meta_file, coding_genes)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-17
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
import gffutils
from loguru import logger
from tqdm import tqdm

# 4. Local Imports
from workflow.src.features.genome import DNA_level_features, PombaseGenomeConfig
from workflow.src.features.protein import (
    extract_protein_features_from_peptide_sequence,
    pLDDT_statistics_report,
)
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
# SHARED HELPERS
# =============================================================================
@logger.catch
def load_gene_meta(gene_meta_file: Path) -> tuple[pd.DataFrame, dict]:
    """Load gene metadata and build a uniprot_id -> gene_systematic_id map."""
    gene_meta = pd.read_csv(gene_meta_file, sep="\t")
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])
    if "uniprot_id" in gene_meta.columns:
        uniprot_col = "uniprot_id"
    elif "external_id" in gene_meta.columns:
        uniprot_col = "external_id"
    else:
        raise KeyError("Neither 'uniprot_id' nor 'external_id' column found in gene metadata")
    uniprot2id = dict(zip(gene_meta[uniprot_col], gene_meta["gene_systematic_id"]))
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
def read_coding_genes(dna_pickle: Path) -> list[str]:
    """Recover the coding-gene set (unique Gene_id) from the DNA-level pickle."""
    return pd.read_pickle(dna_pickle)["Gene_id"].unique().tolist()


@logger.catch
def load_phyloP_and_divergence(literature_dir: Path, gene_meta_file: Path) -> pd.DataFrame:
    """Load the Grech 2019 phyloP/divergence table (shared by evolutionary + phenotype levels)."""
    phyloP_and_divergence = pd.read_excel(
        literature_dir / "grechFitnessLandscapeFission2019.xlsx", sheet_name="Table 2", skiprows=list(range(14))
    ).drop_duplicates(subset=["gene"])
    phyloP_and_divergence["gene_systematic_id"] = update_sysIDs(phyloP_and_divergence["gene"].tolist(), gene_meta_file)
    return phyloP_and_divergence

# =============================================================================
# DNA LEVEL (the "spine" — enumerates coding genes)
# =============================================================================
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


# =============================================================================
# RNA LEVEL
# =============================================================================
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

# =============================================================================
# PROTEIN LEVEL
# =============================================================================
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

# =============================================================================
# EVOLUTIONARY LEVEL
# =============================================================================
@logger.catch
def collect_evolutionary_level_features(
    pombase_dir: Path,
    ensembl_paralogs_tsv: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    coding_genes: list[str],
    phyloP_and_divergence: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble ortholog/paralog counts, evolutionary rate, and phyloP/divergence scores.

    `phyloP_and_divergence` (Grech 2019) is loaded once via load_phyloP_and_divergence()
    and shared with collect_phenotype_level_features, which reads transposon-insertion-density
    columns from the same source table.
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

    return evolutionary_df

# =============================================================================
# NETWORK LEVEL
# =============================================================================
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


# =============================================================================
# PHENOTYPE LEVEL
# =============================================================================
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
        ["Systematic ID", "Gene dispensability. This study", "Category", "Sub_category", "Growth_tier"]
    ].set_index("Systematic ID").rename(
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

# =============================================================================
# MERGE
# =============================================================================
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
    pombe_features["Sub_category"] = pombe_features["Sub_category"].fillna("Not_determined")
    pombe_features["Growth_tier"] = pombe_features["Growth_tier"].fillna(0).astype(int)
    return pombe_features
