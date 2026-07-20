"""[S0] Build the two shared foundation tables for themes A (verification) + D (coherence).

Outputs (results/verification/HD_DIT_HAP/):
  - gene_annotation_long.tsv : gene x (group_type, group_id, group_name) long table,
    unioning GO-CC macrocomplex + GO-BP (goatools-propagated) + KEGG pathway.
    Plus n_groups_per_gene (used by D2 shared-gene detection).
  - verification_master.tsv  : 411 verified genes + DR/DL/R2/RMSE + library/verification
    essentiality + QC fields (filtering_fraction / total_colonies / colony area).

This is a deterministic build; run with the `data_analysis` conda env. Iterate here,
then port to a notebook cell + (later) workflow/scripts/verification/build_foundation.py.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data  # noqa: E402

PIPE = Path("/data/c/yangyusheng_optimized/DIT_HAP_pipeline")
POMBASE = REPO / "resources/external/pombase/2026-06-01"
ONTO = POMBASE / "ontologies_and_associations"


@dataclass(frozen=True)
class Paths:
    verification_xlsx: Path = PIPE / "results/HD_DIT_HAP_generationRAW/20_essentiality_verification/organized_verification_summary.xlsx"
    dr_dl_tsv: Path = PIPE / "results/HD_DIT_HAP_generationRAW/18_gene_level_clustering/all_coding_genes_with_DIT_HAP_clustering.tsv"
    complex_cc: Path = ONTO / "macromolecular_complex_annotation.tsv"
    gaf: Path = ONTO / "gene_ontology_annotation.gaf.tsv"
    obo: Path = ONTO / "go-basic.obo"
    bp_slim: Path = ONTO / "bp_go_slim_terms.tsv"
    cc_slim: Path = ONTO / "cc_go_slim_terms.tsv"
    kegg_combined: Path = PIPE / "resources/KEGG/combined_brite_table.tsv"
    gene_ids: Path = POMBASE / "Gene_metadata/gene_IDs_names_products.tsv"
    out_dir: Path = REPO / "results/verification/HD_DIT_HAP"


P = Paths()

def build_go_cc_complex() -> pd.DataFrame:
    """GO cellular_component macromolecular-complex annotation → long rows."""
    cx = pd.read_csv(P.complex_cc, sep="\t")
    out = cx[["systematic_id", "complex_term_id", "GO_term_name"]].dropna(subset=["systematic_id"])
    out = out.rename(columns={"systematic_id": "gene", "complex_term_id": "group_id", "GO_term_name": "group_name"})
    out.insert(0, "group_type", "GO_CC_complex")
    return out.drop_duplicates()


def build_go_bp() -> pd.DataFrame:
    """GO biological_process, goatools-propagated (part_of/is_a) via load_ontology_data → long rows."""
    cfg = OntologyDataConfig(
        ontology_obo=P.obo,
        ontology_association_gaf=P.gaf,
        slim_terms_table=[P.bp_slim, P.cc_slim],
    )
    dag, _obj, _ns2assoc, _gene2go, go2genes, _slim = load_ontology_data(
        cfg.load_data(), propagate_counts=True, relationships={"part_of", "is_a"}
    )
    # go2genes is flat {GO_id: {gene, ...}}; keep terms whose DAG namespace is BP.
    rows = []
    for go_id, genes in go2genes.items():
        term = dag.get(go_id)
        if term is None or term.namespace != "biological_process":
            continue
        for g in genes:
            rows.append(("GO_BP", go_id, term.name, g))
    out = pd.DataFrame(rows, columns=["group_type", "group_id", "group_name", "gene"])
    return out.drop_duplicates()


def build_kegg_pathway() -> pd.DataFrame:
    """KEGG pathway from combined_brite_table Level_4 [PATH:spoNNNNN]; map gene Name → systematic_id."""
    kg = pd.read_csv(P.kegg_combined, sep="\t")
    lvl4 = kg["Level_4"].astype(str)
    path = lvl4.str.extract(r"\[PATH:(spo\d+)\]")[0]
    path_name = lvl4.str.replace(r"\s*\[PATH:spo\d+\]", "", regex=True).str.strip()
    kg = kg.assign(group_id=path, group_name=path_name).dropna(subset=["group_id"])
    # name -> systematic_id (include synonyms so more KEGG gene names resolve)
    ids = pd.read_csv(P.gene_ids, sep="\t")
    name2id: dict[str, str] = {}
    for _, r in ids.iterrows():
        sid = r["gene_systematic_id"]
        if isinstance(r.get("gene_name"), str):
            name2id.setdefault(r["gene_name"], sid)
        if isinstance(r.get("synonyms"), str):
            for syn in str(r["synonyms"]).split(","):
                name2id.setdefault(syn.strip(), sid)
    kg["gene"] = kg["Name"].map(name2id)
    out = kg.dropna(subset=["gene"])[["group_id", "group_name", "gene"]].copy()
    out.insert(0, "group_type", "KEGG_pathway")
    return out.drop_duplicates()


def build_annotation_long() -> pd.DataFrame:
    parts = [build_go_cc_complex(), build_go_bp(), build_kegg_pathway()]
    long = pd.concat(parts, ignore_index=True)[["group_type", "group_id", "group_name", "gene"]]
    long = long.drop_duplicates()
    # Shared-gene signal for D2 = # distinct *physical GO-CC complexes* a gene belongs to.
    # (NOT the all-types count: propagated GO-BP ancestors would inflate it to meaninglessness.)
    cc = long[long["group_type"] == "GO_CC_complex"]
    n_cc = cc.groupby("gene")["group_id"].nunique().rename("n_cc_complexes_per_gene")
    long = long.merge(n_cc, on="gene", how="left")
    long["n_cc_complexes_per_gene"] = long["n_cc_complexes_per_gene"].fillna(0).astype(int)
    return long


def build_verification_master() -> pd.DataFrame:
    """411 verified genes joined to DR/DL cluster table (DR/DL authoritative from cluster file)."""
    ver = pd.read_excel(P.verification_xlsx)
    dr = pd.read_csv(P.dr_dl_tsv, sep="\t")[
        ["Systematic ID", "DR", "DL", "Cluster", "RevisedDeletion_essentiality"]
    ]
    master = ver.merge(dr, on="Systematic ID", how="left", suffixes=("", "_cluster"))
    return master


def main() -> None:
    P.out_dir.mkdir(parents=True, exist_ok=True)

    long = build_annotation_long()
    long_path = P.out_dir / "gene_annotation_long.tsv"
    long.to_csv(long_path, sep="\t", index=False)
    print(f"[gene_annotation_long] {long.shape} -> {long_path}")
    print(long.groupby("group_type").agg(
        n_rows=("gene", "size"), n_genes=("gene", "nunique"), n_groups=("group_id", "nunique")
    ).to_string())

    master = build_verification_master()
    master_path = P.out_dir / "verification_master.tsv"
    master.to_csv(master_path, sep="\t", index=False)
    print(f"\n[verification_master] {master.shape} -> {master_path}")
    print("DR notna:", master["DR"].notna().sum(), "/", len(master))


if __name__ == "__main__":
    main()
