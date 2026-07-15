"""
Protein-Level Gene Features
==============================

Peptide-sequence-derived features (aromaticity, aliphatic index, amino acid
composition) and AlphaFold pLDDT confidence statistics. Merges the
protein-level half of `pombe_feature_functions.py` with all of
`protein_structure_functions.py` (DIT_HAP_pipeline) — design doc §7 groups
both under `features/protein.py`.

Input
-----
- A PomBase peptide FASTA (peptide.fa)
- A directory of AlphaFold structure files (.pdb.gz / .pdb / .cif / .cif.gz)

Output
------
- extract_protein_features_from_peptide_sequence: one row per peptide record
- pLDDT_statistics_report: one row per structure file, keyed by UniProt ID

Usage
-----
    from workflow.src.features.protein import extract_protein_features_from_peptide_sequence, pLDDT_statistics_report
    protein_features = extract_protein_features_from_peptide_sequence(peptide_fasta)
    pLDDTs = pLDDT_statistics_report(alphafold_dir, structure_format="pdb.gz")

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import gzip
import re
from pathlib import Path
from typing import Literal

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils import seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from loguru import logger
from tqdm import tqdm

# =============================================================================
# CORE LOGIC — peptide-sequence features
# =============================================================================
@logger.catch
def calculate_aliphatic_index_biopython(protein_sequence: str) -> float:
    """Calculate the aliphatic index (Ikai 1980) of a protein sequence."""
    analysis = ProteinAnalysis(protein_sequence)
    aa_percent = analysis.amino_acids_percent
    X_ala = aa_percent.get("A", 0) * 100
    X_val = aa_percent.get("V", 0) * 100
    X_leu = aa_percent.get("L", 0) * 100
    X_ile = aa_percent.get("I", 0) * 100
    aliphatic_index = X_ala + 2.9 * X_val + 3.9 * (X_leu + X_ile)
    return round(aliphatic_index, 3)


@logger.catch
def extract_protein_features_from_peptide_sequence(peptide_fasta_file: Path, return_redundant_meta: bool = False) -> pd.DataFrame:
    """Extract per-gene protein features (aromaticity, aliphatic index, AA composition, ...) from a peptide FASTA."""
    records = []
    aa_content = {}
    aa_percent = {}
    for record in SeqIO.parse(peptide_fasta_file, "fasta"):
        gene_id = re.search(r"(\S+)\.\d:pep$", record.id).groups()[0]
        sequence = str(record.seq).rstrip("*")
        if "*" in sequence:
            logger.warning(f"Stop codon found in sequence of {gene_id}. Truncating at first stop codon.")
            sequence = sequence.split("*")[0]
        analysis = ProteinAnalysis(sequence)
        aa_content[gene_id] = analysis.count_amino_acids()
        aa_percent[gene_id] = analysis.amino_acids_percent
        protein_features = {
            "Gene_id": gene_id,
            "aromaticity": analysis.aromaticity(),
            "aliphatic_index": calculate_aliphatic_index_biopython(sequence),
            "gravy": analysis.gravy(),
            "flexibility": np.mean(analysis.flexibility()),
            "instability_index": analysis.instability_index(),
            "monoisotopic": analysis.monoisotopic,
        }
        protein_features.update(
            dict(zip(("molar_extinction_reduced", "molar_extinction_cystines"), analysis.molar_extinction_coefficient()))
        )
        protein_features.update(
            dict(zip(("Helix_fraction", "Turn_fraction", "Sheet_fraction"), analysis.secondary_structure_fraction()))
        )
        if return_redundant_meta:
            protein_features.update({
                "charge_at_pH": analysis.charge_at_pH(7.0),
                "isoelectric_point": analysis.isoelectric_point(),
                "length": len(sequence),
                "molecular_weight(kDa)": analysis.molecular_weight() / 1000,
            })
        records.append(protein_features)

    aa_content_df = pd.DataFrame.from_dict(aa_content, orient="index")
    aa_percent_df = pd.DataFrame.from_dict(aa_percent, orient="index")
    aa_content_df.columns = [f"aa_count_{seq3(col)}" for col in aa_content_df.columns]
    aa_percent_df.columns = [f"aa_percent_{seq3(col)}" for col in aa_percent_df.columns]

    records_df = pd.DataFrame(records).set_index("Gene_id")
    records_df = records_df.join(aa_content_df).join(aa_percent_df).reset_index()
    return records_df


# =============================================================================
# CORE LOGIC — AlphaFold pLDDT statistics
# =============================================================================
@logger.catch
def extract_pLDDT(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from a PDB or mmCIF file, compressed or not."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)

    f = gzip.open(structure_file, "rt") if structure_file.name.endswith(".gz") else open(structure_file, "r")

    stem = structure_file.name.rstrip(".gz").lower()
    if stem.endswith(".pdb"):
        parser = PDBParser()
    elif stem.endswith(".cif"):
        parser = MMCIFParser()
    else:
        raise ValueError(f"Unknown file format: {structure_file.name}")
    structure = parser.get_structure(structure_file.stem, f)

    pLDDT = [residue["CA"].bfactor for residue in structure.get_residues() if residue.has_id("CA")]
    f.close()
    return pLDDT


@logger.catch
def extract_pLDDT_pdb_gz(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from a .pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)
    pLDDT = [residue["CA"].bfactor for residue in structure.get_residues()]
    f.close()
    return pLDDT


@logger.catch
def extract_pLDDT_pdb(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from an uncompressed PDB file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, structure_file)
    return [residue["CA"].bfactor for residue in structure.get_residues()]


@logger.catch
def extract_protein_seq_pdb_gz(structure_file: Path | str) -> str:
    """Extract the residue sequence from a .pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)
    ppb = PPBuilder()
    seq = ppb.build_peptides(structure)[0].get_sequence()
    f.close()
    return seq


@logger.catch
def pLDDT_statistics_report(
    structure_dir: Path,
    structure_format: Literal["pdb", "pdb.gz", "cif", "cif.gz", "mixed"] = "pdb.gz",
) -> pd.DataFrame:
    """Summarize per-residue pLDDT into per-structure mean/std/CV/disorder-fraction, keyed by UniProt ID."""
    all_pdb_files = list(structure_dir.glob(f"*.{structure_format}"))
    pLDDT_records = []

    for pdb_file in tqdm(all_pdb_files):
        uniprot_id = pdb_file.name.split("-F1-")[0].split("AF-")[1]
        match structure_format:
            case "pdb":
                pLDDT = np.array(extract_pLDDT_pdb(pdb_file))
            case "pdb.gz":
                pLDDT = np.array(extract_pLDDT_pdb_gz(pdb_file))
            case "cif" | "cif.gz" | "mixed":
                pLDDT = np.array(extract_pLDDT(pdb_file))
            case _:
                raise ValueError(f"Unsupported structure format: {structure_format}")
        length_protein = len(pLDDT)
        mean_pLDDT = np.mean(pLDDT)
        std_pLDDT = np.std(pLDDT)
        cv_pLDDT = std_pLDDT / mean_pLDDT if mean_pLDDT != 0 else np.nan
        disorder_fraction = np.sum(pLDDT < 50) / length_protein
        pLDDT_records.append({
            "uniprot_id": uniprot_id,
            "protein_length": length_protein,
            "pLDDT": ",".join(pLDDT.astype(str)),
            "mean_pLDDT": round(mean_pLDDT, 3),
            "std_pLDDT": round(std_pLDDT, 3),
            "cv_pLDDT": round(cv_pLDDT, 3),
            "disorder_fraction": round(disorder_fraction, 3),
        })

    return pd.DataFrame(pLDDT_records)
