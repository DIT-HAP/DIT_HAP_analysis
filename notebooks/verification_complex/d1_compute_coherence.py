"""[D1] Compute DR-DL coherence for GO-CC complexes + GO-BP + KEGG pathways.

Extends the pipeline's complex-only coherence to all three group types in the
[S0] foundation table. Faithful to complex_analysis.ipynb: DR>0.3 filter on both
members and the permutation background; term-size window 3..300; median_pairwise
z-score as the primary coherence axis (negative = coherent).

Reads:  results/verification/HD_DIT_HAP/gene_annotation_long.tsv (+ DR/DL source)
Writes: results/coherence/HD_DIT_HAP/coherence_metrics_all_namespaces.tsv

Run with the `data_analysis` conda env (numpy/scipy/goatools present).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from workflow.src.coherence.metrics import (  # noqa: E402
    ZSCORE_METHODS,
    coherence_metrics,
    compute_distance_zscore,
    normalize_dr_dl,
)

PIPE = Path("/data/c/yangyusheng_optimized/DIT_HAP_pipeline")
DR_DL_TSV = PIPE / "results/HD_DIT_HAP_generationRAW/18_gene_level_clustering/all_coding_genes_with_DIT_HAP_clustering.tsv"
LONG_TSV = REPO / "results/verification/HD_DIT_HAP/gene_annotation_long.tsv"
OUT_DIR = REPO / "results/coherence/HD_DIT_HAP"

DR_FILTER = 0.3          # non-essential genes dropped before coherence (both member + bg)
TERM_SIZE_MIN = 3        # permutation z unstable below ~4; floor kept at notebook value
TERM_SIZE_MAX = 300
PRIMARY_METRIC = "median_pairwise_distance"


def load_dr_dl() -> pd.DataFrame:
    df = pd.read_csv(DR_DL_TSV, sep="\t")[["Systematic ID", "DR", "DL"]].dropna(subset=["DR", "DL"])
    df = df.rename(columns={"Systematic ID": "gene"})
    return df.query("DR > @DR_FILTER").reset_index(drop=True)


def compute_group_type(long: pd.DataFrame, dr_dl: pd.DataFrame, group_type: str) -> pd.DataFrame:
    """Coherence for every group in one group_type (CC / BP / KEGG)."""
    ann = long[long["group_type"] == group_type][["group_id", "group_name", "gene"]].drop_duplicates()
    ann = ann.merge(dr_dl, on="gene", how="inner")  # inner = keep only DR>0.3 members

    bg = normalize_dr_dl(dr_dl["DR"].values, dr_dl["DL"].values)

    rows = []
    for (gid, gname), sub in ann.groupby(["group_id", "group_name"]):
        n = sub["gene"].nunique()
        if not (TERM_SIZE_MIN <= n <= TERM_SIZE_MAX):
            continue
        X = normalize_dr_dl(sub["DR"].values, sub["DL"].values)
        row = {
            "group_type": group_type,
            "group_id": gid,
            "group_name": gname,
            "term_size": n,
            "covered_genes": ", ".join(sorted(sub["gene"])),
            **coherence_metrics(X),
        }
        for method in ZSCORE_METHODS:
            z, p = compute_distance_zscore(X, bg, method=method)
            row[f"{method}_zscore"] = z
            row[f"{method}_pvalue"] = p
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    long = pd.read_csv(LONG_TSV, sep="\t")
    dr_dl = load_dr_dl()
    print(f"DR>{DR_FILTER} background genes: {len(dr_dl)}")

    parts = [compute_group_type(long, dr_dl, gt) for gt in ["GO_CC_complex", "GO_BP", "KEGG_pathway"]]
    out = pd.concat(parts, ignore_index=True)
    out_path = OUT_DIR / "coherence_metrics_all_namespaces.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[coherence_metrics_all_namespaces] {out.shape} -> {out_path}\n")

    zc = f"{PRIMARY_METRIC}_zscore"
    print("=== groups scored per type ===")
    print(out.groupby("group_type").size().to_string())
    print(f"\n=== most coherent (lowest {zc}) per type ===")
    for gt, sub in out.groupby("group_type"):
        top = sub.nsmallest(3, zc)[["group_name", "term_size", zc]]
        print(f"\n[{gt}]")
        print(top.to_string(index=False))
    print(f"\n=== most incoherent (highest {zc}) overall ===")
    print(out.nlargest(8, zc)[["group_type", "group_name", "term_size", zc]].to_string(index=False))


if __name__ == "__main__":
    main()
