#!/usr/bin/env python
"""
[X] Closed loop: contrast library-essentiality vs DR-DL coherence vs verified-essentiality

Target groups: coherent (z<-1) AND flip-heavy (≥2 flip, from A2)

Dimensions:
- Library annotation: E vs V labels (pombase + manual curation)
- DR-DL coherence: geometric-median distance z-scores from D1
- Verified essentiality: actual DR phenotypes from current study

Key question: within coherent groups, does library-annotation heterogeneity
(flip-heavy = mixed E/V) reflect measurement noise OR genuine boundary cases?

Expected output:
- For each target group: E/V/flip breakdown, coherence_z, DR/DL distribution
- Enrichment: are flip genes at coherence boundary (higher distance)?
- Contrast: do flip_V2E (library said V, verified E) differ from flip_E2V in DR/DL?

Input:
- verified_group_summary.tsv (from A2: groups × verification outcomes)
- coherence_metrics_all_namespaces.tsv (from D1: aggregate coherence per group)
- gene_outcomes.tsv (from A1: per-gene outcome labels)
- verified_genes.tsv (DR/DL values)
- gene_annotation_long.tsv (S0: group membership)

Output:
- closed_loop_groups.tsv: target groups with enrichment stats
- closed_loop_genes_long.tsv: per-gene records for target groups
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from workflow.src.coherence.metrics import geometric_median, normalize_dr_dl

# Paths
base = Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis/.claude/worktrees/followup-analysis-expansion")
results = base / "results/verification/HD_DIT_HAP"
coherence_dir = Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis/results/coherence/HD_DIT_HAP")

# (1) Load A2 summary: groups × verification outcomes
summary = pd.read_csv(results / "verified_group_summary.tsv", sep="\t")

# (2) Filter: coherent (z<-1) AND flip-heavy (≥2 flip)
scored = summary[summary["coherence_z"].notna() & (summary["n_flip"] >= 2)]
targets = scored[scored["coherence_z"] < -1].copy()
print(f"[closed_loop_groups] {len(targets)} target groups (coherent z<-1 AND flip-heavy)")

# (3) Load D1 coherence aggregates
coh_agg = pd.read_csv(coherence_dir / "coherence_metrics_all_namespaces.tsv", sep="\t")
coh_agg = coh_agg[["group_type", "group_id", "group_name", "covered_genes", "centroid_x", "centroid_y"]]

# (4) Load A1 verification table: gene → outcome + DR/DL
verification = pd.read_csv(
    Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis/results/verification/HD_DIT_HAP/verification_outcome_table.tsv"),
    sep="\t"
)
verification = verification.rename(columns={"Systematic ID": "gene"})
verification = verification[["gene", "outcome", "DR", "DL"]]

# (6) Load S0 long table: group membership
long = pd.read_csv(
    Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis/results/verification/HD_DIT_HAP/gene_annotation_long.tsv"),
    sep="\t"
)

# (7) Build target group membership + outcomes + DR/DL
target_set = set(zip(targets["group_type"], targets["group_name"]))
long_target = long[
    long.apply(lambda r: (r["group_type"], r["group_name"]) in target_set, axis=1)
].copy()

# Join outcome + DR/DL
long_target = long_target.merge(
    verification[["gene", "outcome", "DR", "DL"]],
    on="gene",
    how="left"
)

# (8) Compute per-gene distance to group centroid
def compute_distances(group_df, centroid_x, centroid_y):
    """Add distance to centroid for each gene."""
    dr_norm = group_df["DR"].values / 1.0
    dl_norm = group_df["DL"].values / 10.0
    centroid = np.array([centroid_x, centroid_y])

    distances = []
    for dr_n, dl_n in zip(dr_norm, dl_norm):
        point = np.array([dr_n, dl_n])
        dist = np.linalg.norm(point - centroid)
        distances.append(dist)

    group_df = group_df.copy()
    group_df["distance_to_centroid"] = distances
    return group_df

# Join centroids
long_target = long_target.merge(
    coh_agg[["group_type", "group_name", "centroid_x", "centroid_y"]],
    on=["group_type", "group_name"],
    how="left"
)

# Compute distances per group
result_chunks = []
for (gtype, gname), group_df in long_target.groupby(["group_type", "group_name"]):
    cx = group_df["centroid_x"].iloc[0]
    cy = group_df["centroid_y"].iloc[0]
    if pd.notna(cx) and pd.notna(cy):
        group_df = compute_distances(group_df, cx, cy)
    result_chunks.append(group_df)

long_target = pd.concat(result_chunks, ignore_index=True)

print(f"[closed_loop_genes_long] {len(long_target)} gene records across {len(targets)} groups")

# (9) Per-group enrichment: compare flip vs non-flip distance
def enrichment_stats(group_df):
    """Compute flip enrichment at coherence boundary."""
    flip_mask = group_df["outcome"].str.startswith("flip_", na=False)

    flip_dist = group_df.loc[flip_mask, "distance_to_centroid"].dropna()
    nonflip_dist = group_df.loc[~flip_mask, "distance_to_centroid"].dropna()

    if len(flip_dist) < 2 or len(nonflip_dist) < 2:
        return pd.Series({
            "flip_mean_distance": flip_dist.mean() if len(flip_dist) > 0 else np.nan,
            "nonflip_mean_distance": nonflip_dist.mean() if len(nonflip_dist) > 0 else np.nan,
            "distance_delta": np.nan,
            "enrichment_signal": "insufficient_data"
        })

    delta = flip_dist.mean() - nonflip_dist.mean()

    # Enrichment signal: flip genes further from centroid?
    if delta > 0.05:  # arbitrary threshold, 5% of normalized space
        signal = "flip_peripheral"
    elif delta < -0.05:
        signal = "flip_central"
    else:
        signal = "no_difference"

    return pd.Series({
        "flip_mean_distance": flip_dist.mean(),
        "nonflip_mean_distance": nonflip_dist.mean(),
        "distance_delta": delta,
        "enrichment_signal": signal
    })

enrich = long_target.groupby(["group_type", "group_name"]).apply(enrichment_stats).reset_index()

# (10) Merge back to targets
targets_out = targets.merge(enrich, on=["group_type", "group_name"], how="left")

# (11) Contrast: flip_V2E vs flip_E2V DR/DL distribution
flip_v2e = long_target[long_target["outcome"] == "flip_V2E"]
flip_e2v = long_target[long_target["outcome"] == "flip_E2V"]

print(f"\n=== flip contrast ===")
print(f"flip_V2E: n={len(flip_v2e)}, DR={flip_v2e['DR'].mean():.3f}±{flip_v2e['DR'].std():.3f}, DL={flip_v2e['DL'].mean():.2f}±{flip_v2e['DL'].std():.2f}")
print(f"flip_E2V: n={len(flip_e2v)}, DR={flip_e2v['DR'].mean():.3f}±{flip_e2v['DR'].std():.3f}, DL={flip_e2v['DL'].mean():.2f}±{flip_e2v['DL'].std():.2f}")

# (12) Write outputs
targets_out.to_csv(results / "closed_loop_groups.tsv", sep="\t", index=False)
print(f"\n[closed_loop_groups] ({targets_out.shape[0]}, {targets_out.shape[1]}) -> {results / 'closed_loop_groups.tsv'}")

long_target.to_csv(results / "closed_loop_genes_long.tsv", sep="\t", index=False)
print(f"[closed_loop_genes_long] ({long_target.shape[0]}, {long_target.shape[1]}) -> {results / 'closed_loop_genes_long.tsv'}")

# (13) Summary stats
print(f"\n=== enrichment signal distribution ===")
print(targets_out["enrichment_signal"].value_counts())

print(f"\n=== top 10 flip_peripheral groups (flip genes furthest from centroid) ===")
peripheral = targets_out[targets_out["enrichment_signal"] == "flip_peripheral"].sort_values("distance_delta", ascending=False).head(10)
if len(peripheral) > 0:
    print(peripheral[["group_type", "group_name", "n_flip", "coherence_z", "distance_delta"]].to_string(index=False))
else:
    print("(none)")

print(f"\n=== top 10 flip_central groups (flip genes closest to centroid) ===")
central = targets_out[targets_out["enrichment_signal"] == "flip_central"].sort_values("distance_delta").head(10)
if len(central) > 0:
    print(central[["group_type", "group_name", "n_flip", "coherence_z", "distance_delta"]].to_string(index=False))
else:
    print("(none)")
