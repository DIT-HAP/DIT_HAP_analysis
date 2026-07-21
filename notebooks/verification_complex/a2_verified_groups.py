"""[A2] Which functional groups do the verified (esp. correctly-verified) genes sit in?

Joins A1 outcomes x S0 annotations x D1 coherence, aggregating per group. Flags
groups whose members were systematically re-called by verification (the input to
the [X] closed loop). Answers question (1): which pathways/complexes hold the
correctly-verified genes, and where verification changes the picture.

Reads:  verification_outcome_table.tsv (A1), gene_annotation_long.tsv (S0),
        coherence_metrics_all_namespaces.tsv (D1)
Writes: results/verification/HD_DIT_HAP/verified_group_summary.tsv
Run with the `data_analysis` conda env.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
VER_DIR = REPO / "results/verification/HD_DIT_HAP"
COH_TSV = REPO / "results/coherence/HD_DIT_HAP/coherence_metrics_all_namespaces.tsv"

FLIP_OUTCOMES = {"flip_V2E", "flip_E2V"}


def main() -> None:
    oc = pd.read_csv(VER_DIR / "verification_outcome_table.tsv", sep="\t")
    long = pd.read_csv(VER_DIR / "gene_annotation_long.tsv", sep="\t")
    coh = pd.read_csv(COH_TSV, sep="\t")

    # gene -> outcome (use raw concordance so a flagged flip still counts as a flip here)
    oc_slim = oc[["Systematic ID", "essentiality_concordance", "is_ambiguous"]].rename(
        columns={"Systematic ID": "gene"}
    )
    ann = long[["group_type", "group_id", "group_name", "gene"]].drop_duplicates()
    joined = ann.merge(oc_slim, on="gene", how="inner")  # verified members only

    def agg(sub: pd.DataFrame) -> pd.Series:
        conc = sub["essentiality_concordance"]
        return pd.Series({
            "n_verified": sub["gene"].nunique(),
            "n_concordant_E": (conc == "concordant_E").sum(),
            "n_concordant_V": (conc == "concordant_V").sum(),
            "n_flip_V2E": (conc == "flip_V2E").sum(),
            "n_flip_E2V": (conc == "flip_E2V").sum(),
            "n_flip": conc.isin(FLIP_OUTCOMES).sum(),
            "n_ambiguous": sub["is_ambiguous"].sum(),
            "flip_genes": ", ".join(sorted(sub.loc[conc.isin(FLIP_OUTCOMES), "gene"])),
        })

    summary = joined.groupby(["group_type", "group_id", "group_name"]).apply(agg, include_groups=False).reset_index()
    summary["flip_fraction"] = (summary["n_flip"] / summary["n_verified"]).round(3)

    # overlay D1 coherence (median-pairwise z) so a flip-heavy group can be read against
    # whether it's DR-DL coherent — the [X] closed loop's core comparison. Only groups
    # in D1's 3..300 size window get a coherence_z; broad root terms are left NaN.
    zc = "median_pairwise_distance_zscore"
    summary = summary.merge(
        coh[["group_type", "group_id", zc, "term_size"]].rename(columns={zc: "coherence_z", "term_size": "coherence_term_size"}),
        on=["group_type", "group_id"], how="left",
    )

    out_path = VER_DIR / "verified_group_summary.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"[verified_group_summary] {summary.shape} -> {out_path}\n")

    print("=== groups with >=1 verified member, per type ===")
    print(summary.groupby("group_type")["n_verified"].agg(["count", "sum"]).to_string())

    # The interpretable cut for [X]: SPECIFIC groups D1 scored (have coherence_z),
    # flip-heavy, ranked by flips. Root terms (coherence_z NaN, too big) excluded.
    scored = summary[summary["coherence_z"].notna() & (summary["n_flip"] >= 2)]
    cols = ["group_type", "group_name", "n_verified", "n_flip_V2E", "n_flip_E2V", "flip_fraction", "coherence_z"]
    print("\n=== flip-heavy SCORED groups, ranked by n_flip (excl. root terms) ===")
    print(scored.sort_values("n_flip_V2E", ascending=False)[cols].head(20).to_string(index=False))
    print("\n=== [X] target: COHERENT (z<-1) AND flip-heavy — library-mixed but DR-DL-tight ===")
    xt = scored[scored["coherence_z"] < -1].sort_values("n_flip_V2E", ascending=False)
    print(xt[cols].head(20).to_string(index=False))
    print(f"\ntotal flip-heavy scored groups: {len(scored)}; of which coherent (z<-1): {len(xt)}")


if __name__ == "__main__":
    main()
