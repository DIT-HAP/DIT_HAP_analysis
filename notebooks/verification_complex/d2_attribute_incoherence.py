"""[D2] Attribute incoherence for each incoherent complex (major/minor + shared-gene).

Reads:  results/3a_coherence/HD_DIT_HAP/coherence_metrics_all_namespaces.tsv (D1)
        results/verification/HD_DIT_HAP/gene_annotation_long.tsv (S0)
        DR/DL cluster table (member coordinates)
Writes: results/3a_coherence/HD_DIT_HAP/incoherence_attribution.tsv  (one row / complex)
        results/3a_coherence/HD_DIT_HAP/shared_subunits_long.tsv      (member x other-complex)

Focus on GO-CC complexes with median_pairwise z > threshold (the pipeline's
incoherent set). Run with the `data_analysis` conda env.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from workflow.src.coherence.attribution import (  # noqa: E402
    attribute_incoherence,
    major_minor_split,
    shared_subunit_fractions,
    shared_subunits,
)

PIPE = Path("/data/c/yangyusheng_optimized/DIT_HAP_pipeline")
DR_DL_TSV = PIPE / "results/HD_DIT_HAP_generationRAW/18_gene_level_clustering/all_coding_genes_with_DIT_HAP_clustering.tsv"
COH_TSV = REPO / "results/3a_coherence/HD_DIT_HAP/coherence_metrics_all_namespaces.tsv"
LONG_TSV = REPO / "results/verification/HD_DIT_HAP/gene_annotation_long.tsv"
OUT_DIR = REPO / "results/3a_coherence/HD_DIT_HAP"

# (DR, DL/10): plain divisors, not a fitted range (DR_NORM_MAX = 1.0,
# DL_NORM_MAX = 10.0 in workflow/src/coherence/metrics.py).
DR_NORM_MAX = 1.0
DL_NORM_MAX = 10.0

Z_INCOHERENT_MIN = 0.5   # median_pairwise z above this = incoherent candidate
DR_FILTER = 0.3          # same membership filter as D1 coherence


def main() -> None:
    coh = pd.read_csv(COH_TSV, sep="\t")
    long = pd.read_csv(LONG_TSV, sep="\t")
    long_cc = long[long["group_type"] == "GO_CC_complex"][["group_name", "gene"]].drop_duplicates()
    dr = (
        pd.read_csv(DR_DL_TSV, sep="\t")[["Systematic ID", "DR", "DL"]]
        .dropna(subset=["DR", "DL"]).rename(columns={"Systematic ID": "gene"})
        .query("DR > @DR_FILTER")
    )
    dr_map = dr.set_index("gene")[["DR", "DL"]]

    zc = "median_pairwise_distance_z"
    incoherent = coh[(coh["group_type"] == "GO_CC_complex") & (coh[zc] > Z_INCOHERENT_MIN)]
    incoherent = incoherent.sort_values(zc, ascending=False)
    print(f"incoherent CC complexes (z>{Z_INCOHERENT_MIN}): {len(incoherent)}")

    attr_rows = []
    shared_long_rows = []
    for _, r in incoherent.iterrows():
        name = r["group_name"]
        members = [g for g in long_cc.loc[long_cc["group_name"] == name, "gene"] if g in dr_map.index]
        if len(members) < 3:
            continue
        X = np.column_stack([
            dr_map.loc[members, "DR"].values / DR_NORM_MAX,
            dr_map.loc[members, "DL"].values / DL_NORM_MAX,
        ])
        split = major_minor_split(X)
        sf = shared_subunit_fractions(long_cc).get(name, np.nan)
        label = attribute_incoherence(split, sf)

        attr_rows.append({
            "complex": name,
            "n_scored_members": len(members),
            "median_pairwise_z": round(r[zc], 3),
            "gmm_silhouette": round(split["silhouette"], 3) if pd.notna(split.get("silhouette")) else None,
            "gmm_is_split": split["is_split"],
            "gmm_component_sizes": str(split.get("component_sizes")),
            "shared_fraction": round(sf, 3) if pd.notna(sf) else None,
            "attribution": label,
        })

        sh = shared_subunits(long_cc, name)
        for _, s in sh.iterrows():
            shared_long_rows.append({"complex": name, **s.to_dict()})

    attr = pd.DataFrame(attr_rows)
    attr_path = OUT_DIR / "incoherence_attribution.tsv"
    attr.to_csv(attr_path, sep="\t", index=False)
    shared_df = pd.DataFrame(shared_long_rows)
    shared_path = OUT_DIR / "shared_subunits_long.tsv"
    shared_df.to_csv(shared_path, sep="\t", index=False)

    print(f"[incoherence_attribution] {attr.shape} -> {attr_path}")
    print(f"[shared_subunits_long] {shared_df.shape} -> {shared_path}\n")
    print("=== attribution breakdown ===")
    print(attr["attribution"].value_counts().to_string())
    print("\n=== per-complex ===")
    print(attr.to_string(index=False))


if __name__ == "__main__":
    main()
