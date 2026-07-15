"""
DNA-Level Gene Features
=========================

Per-mRNA DNA-level features (telomere/centromere distance, GC content, intron
structure, codon usage) extracted from a PomBase GFF3 annotation + genome
FASTA. Byte-faithful port of the DNA-level half of
`pombe_feature_functions.py` (DIT_HAP_pipeline) — protein-level functions
moved to `workflow/src/features/protein.py`.

Input
-----
- A gffutils FeatureDB built from a PomBase GFF3
- The corresponding genome FASTA
- peptide_stats.tsv (for primary-transcript peptide-length comparison)

Output
------
- One DNA_level_features record per mRNA feature in the FeatureDB

Usage
-----
    from workflow.src.features.genome import PombaseGenomeConfig, DNA_level_features
    cfg = PombaseGenomeConfig.from_pombase_dir(pombase_dir)
    db = gffutils.FeatureDB(cfg.database_file)
    features = [DNA_level_features.from_gffutils_feature(m, db, cfg) for m in db.features_of_type("mRNA")]

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC123, gc_fraction
from codonbias.scores import EffectiveNumberOfCodons
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
CHR1_LEFT_TELOMERE_END = 29663
CHR1_RIGHT_TELOMERE_START = 5554844
CHR2_LEFT_TELOMERE_END = 39186
CHR2_RIGHT_TELOMERE_START = 4500619
CHR3_LEFT_RIBOSOMAL_DNA_END = 23130
CHR3_RIGHT_RIBOSOMAL_DNA_START = 2440994

CHROMOSOME_END = {
    "left": {"I": CHR1_LEFT_TELOMERE_END, "II": CHR2_LEFT_TELOMERE_END, "III": CHR3_LEFT_RIBOSOMAL_DNA_END},
    "right": {"I": CHR1_RIGHT_TELOMERE_START, "II": CHR2_RIGHT_TELOMERE_START, "III": CHR3_RIGHT_RIBOSOMAL_DNA_START},
}

CENTROMERE_POSITIONS = {
    "I": (3753687, 3789421),
    "II": (1602418, 1644747),
    "III": (1070904, 1137003),
}

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PombaseGenomeConfig:
    """Genome sequence/annotation paths and derived lookup tables for one PomBase version."""
    fasta_file: str
    fai_file: str
    gff3_file: str
    database_file: str
    genome_sequence_dict: dict[str, Seq]
    genome_length_dict: dict[str, int]
    primary_peptide_length: dict[str, int]

    @classmethod
    def from_pombase_dir(cls, pombase_dir: Path) -> PombaseGenomeConfig:
        """Build config from a PomBase version directory (e.g. resources/external/pombase/2025-10-01)."""
        genome_dir = pombase_dir / "genome_sequence_and_features"
        fasta_file = str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.fa")
        peptide_stats = pd.read_csv(pombase_dir / "Protein_features" / "peptide_stats.tsv", sep="\t", index_col=0)
        genome_sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

        return cls(
            fasta_file=fasta_file,
            fai_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.fa.fai"),
            gff3_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.gff3"),
            database_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.db"),
            genome_sequence_dict=genome_sequence_dict,
            genome_length_dict={k: len(v) for k, v in genome_sequence_dict.items()},
            primary_peptide_length=peptide_stats["Residues"].to_dict(),
        )


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def determine_primary_candidate(gene_id: str, mRNA_id: str, peptide_length: int, primary_peptide_length: int) -> bool:
    """Determine if the mRNA is the primary transcript, with hardcoded exceptions for known mis-annotated loci."""
    if gene_id in ["SPBC119.04", "SPBC17A3.07"]:
        return mRNA_id.endswith(".1")
    elif gene_id in ["SPAC212.11", "SPAC2E12.05", "SPAC977.01", "SPMIT.03", "SPMIT.06", "SPMIT.08"]:
        return mRNA_id.endswith(".1")
    else:
        return peptide_length == primary_peptide_length


@dataclass
class DNA_level_features:
    """One mRNA's DNA-level feature record."""
    Gene_id: str
    mRNA_id: str
    Chromosome: Literal["chr_II_telomeric_gap", "I", "II", "III", "mating_type_region", "mitochondrial"]
    Start: int
    End: int
    Strand: Literal["+", "-"]
    Abs_distance_from_telomere: float
    Relative_distance_from_telomere: float
    Abs_distance_from_centromere: float
    Relative_distance_from_centromere: float
    Gene_length: int
    GC_content_of_gene: float
    CDS_number: int
    GC_content_of_CDS: float
    Fraction_of_CDS: float
    GC3: float
    Containing_intron: bool
    Intron_number: int
    GC_content_of_intron: float
    Total_intron_length: int
    Average_intron_length: float
    Length_of_first_intron: int
    GC_contents_of_first_intron: float
    ENC: float
    Peptide_length: int
    Primary_peptide_length: int
    Primary_candidate: bool

    @classmethod
    def from_gffutils_feature(cls, mRNA: gffutils.Feature, db: gffutils.FeatureDB, cfg: PombaseGenomeConfig) -> DNA_level_features:
        """Compute one mRNA's DNA-level features from its gffutils Feature."""
        if mRNA.strand == "+":
            CDSs = list(db.children(mRNA, featuretype="CDS", order_by="start"))
            start = getattr(CDSs[0], "start", 0)
            end = getattr(CDSs[-1], "end", 0)
        else:
            CDSs = list(db.children(mRNA, featuretype="CDS", reverse=True, order_by="start"))
            start = getattr(CDSs[0], "end", 0)
            end = getattr(CDSs[-1], "start", 0)

        gene_id = mRNA.attributes.get("Parent")[0]
        mRNA_id = mRNA.id
        chrom = mRNA.chrom
        strand = mRNA.strand
        midpoint = (start + end) // 2

        abs_distance_from_telomere = min(
            abs(midpoint - CHROMOSOME_END["left"].get(chrom, np.nan)),
            abs(midpoint - CHROMOSOME_END["right"].get(chrom, np.nan)),
        )
        relative_distance_from_telomere = round(abs_distance_from_telomere / cfg.genome_length_dict[chrom], 3)
        abs_distance_from_centromere = abs(midpoint - np.mean(CENTROMERE_POSITIONS.get(chrom, (np.nan, np.nan))))
        relative_distance_from_centromere = round(abs_distance_from_centromere / cfg.genome_length_dict[chrom], 3)

        gene_length = abs(end - start) + 1
        GC_content_of_gene = round(gc_fraction(cfg.genome_sequence_dict[chrom][min(start, end):max(start, end)]), 3)

        CDS_number = len(CDSs)
        CDS_sequence = "".join(cds.sequence(cfg.fasta_file) for cds in CDSs)
        GC_content_of_CDS = round(gc_fraction(CDS_sequence), 3)
        Fraction_of_CDS = round(len(CDS_sequence) / gene_length, 3)
        GC3 = round(GC123(Seq(CDS_sequence))[-1], 3)

        introns = list(db.children(mRNA, featuretype="intron", order_by="start"))
        Containing_intron = len(introns) > 0
        Intron_number = len(introns)
        intron_sequences = [intron.sequence(cfg.fasta_file) for intron in introns]
        intron_sequence = "".join(intron_sequences)
        GC_content_of_intron = round(gc_fraction(intron_sequence), 3)
        Total_intron_length = len(intron_sequence)
        Average_intron_length = Total_intron_length / Intron_number if Intron_number > 0 else 0
        Length_of_first_intron = len(intron_sequences[0]) if Intron_number > 0 else 0
        GC_contents_of_first_intron = round(gc_fraction(intron_sequences[0]), 3) if Intron_number > 0 else 0.0

        ENC_model = EffectiveNumberOfCodons(mean="unweighted")
        ENC = np.round(ENC_model.get_score(CDS_sequence), 2)

        Peptide_length = len(Seq(CDS_sequence).translate(to_stop=True))
        primary_peptide_length = cfg.primary_peptide_length[gene_id]
        primary_candidate = determine_primary_candidate(gene_id, mRNA_id, Peptide_length, primary_peptide_length)

        return cls(
            Gene_id=gene_id, mRNA_id=mRNA_id, Chromosome=chrom, Start=start, End=end, Strand=strand,
            Abs_distance_from_telomere=abs_distance_from_telomere,
            Relative_distance_from_telomere=relative_distance_from_telomere,
            Abs_distance_from_centromere=abs_distance_from_centromere,
            Relative_distance_from_centromere=relative_distance_from_centromere,
            Gene_length=gene_length, GC_content_of_gene=GC_content_of_gene,
            CDS_number=CDS_number, GC_content_of_CDS=GC_content_of_CDS,
            Fraction_of_CDS=Fraction_of_CDS, GC3=GC3,
            Containing_intron=Containing_intron, Intron_number=Intron_number,
            GC_content_of_intron=GC_content_of_intron, Total_intron_length=Total_intron_length,
            Average_intron_length=Average_intron_length, Length_of_first_intron=Length_of_first_intron,
            GC_contents_of_first_intron=GC_contents_of_first_intron, ENC=ENC,
            Peptide_length=Peptide_length, Primary_peptide_length=primary_peptide_length,
            Primary_candidate=primary_candidate,
        )


@logger.catch
def calculate_anticodon_usage_matrix(db: gffutils.FeatureDB, cfg: PombaseGenomeConfig) -> pd.DataFrame:
    """Compute a gene x anti-codon count matrix across all coding genes' concatenated CDS sequence."""
    from collections import Counter

    bases = ["A", "T", "G", "C"]
    codons = [f"{b1}{b2}{b3}" for b1 in bases for b2 in bases for b3 in bases]
    anticodons = [str(Seq(codon).reverse_complement()) for codon in codons]
    codon_to_anticodon = dict(zip(codons, anticodons))

    records = []
    for mRNA in db.features_of_type("mRNA"):
        gene_id = mRNA.attributes.get("Parent")[0]
        CDSs = list(db.children(mRNA, featuretype="CDS", order_by="start"))
        if not CDSs:
            logger.warning(f"No CDS found for {gene_id}, skipping")
            continue

        CDS_sequence = "".join(cds.sequence(cfg.fasta_file) for cds in CDSs)
        CDS_length = len(CDS_sequence) - (len(CDS_sequence) % 3)
        CDS_sequence = CDS_sequence[:CDS_length]

        codon_counts = Counter(CDS_sequence[i:i + 3] for i in range(0, CDS_length, 3))
        anticodon_counts = {
            codon_to_anticodon.get(codon): count
            for codon, count in codon_counts.items()
            if codon in codon_to_anticodon
        }
        anticodon_counts["Gene_id"] = gene_id
        records.append(anticodon_counts)

    df = pd.DataFrame(records).fillna(0)
    df = df.set_index("Gene_id")
    for anticodon in anticodons:
        if anticodon not in df.columns:
            df[anticodon] = 0
    return df[sorted(anticodons)].astype(int).reset_index()
