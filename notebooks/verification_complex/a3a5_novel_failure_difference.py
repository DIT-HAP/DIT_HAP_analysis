#!/usr/bin/env python
"""
[A3-A5] Novel genes, failure attribution, and difference explanation.

Three verification deep-dives sharing verification_outcome_table.tsv:

[A3] Novel genes: uncharacterized (gene_name still systematic-ID form) ∩ high R²
     ∩ clear (non-ambiguous) outcome, prioritizing the flip_V2E set. These are
     the discovery payoff — genes the deletion library called viable/uncertain
     that DIT-HAP now confidently reclassifies as essential.

[A4] Failure attribution: for ambiguous / discordant genes, split the cause into
     data-quality (filtering_fraction / total_colonies / R²) vs biology
     (condition-dependent FYPO). A gene flagged only for data reasons is a
     fixable measurement; one that is clean but condition-dependent is a genuine
     biological boundary case.

[A5] Difference explanation: the 97 raw flips (80 V2E + 17 E2V) crossed with
     FYPOviability. Condition-dependent annotations are the expected reservoir of
     flips; unconditional viable→essential flips crossed with colony-area decline
     are the strongest reclassification calls.

Reads:  results/verification/HD_DIT_HAP/verification_outcome_table.tsv
Writes: results/verification/HD_DIT_HAP/a3_novel_genes.tsv
        results/verification/HD_DIT_HAP/a4_failure_attribution.tsv
        results/verification/HD_DIT_HAP/a5_difference_explanation.tsv
"""

from pathlib import Path
import numpy as np
import pandas as pd

MAIN = Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis")
TABLE = MAIN / "results/verification/HD_DIT_HAP/verification_outcome_table.tsv"
OUT_DIR = MAIN / "results/verification/HD_DIT_HAP"

R2_HIGH = 0.98          # top-half fit; median is 0.984
COLONY_LOW = 52         # 25th percentile of total_colonies
FILTER_LOW = 0.88       # 25th percentile of filtering_fraction


def load() -> pd.DataFrame:
    df = pd.read_csv(TABLE, sep="\t")
    df = df.rename(columns={"Systematic ID": "gene"})
    # uncharacterized = gene_name still in systematic-ID form (no trivial name)
    df["is_uncharacterized"] = df["gene_name"].astype(str).str.match(r"^SP[AB]")
    return df


# ---------------------------------------------------------------- A3
def a3_novel_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Uncharacterized + high-R² + clear outcome, flip_V2E prioritized."""
    novel = df[
        df["is_uncharacterized"]
        & (df["R2"] >= R2_HIGH)
        & (~df["is_ambiguous"])
    ].copy()

    # priority: flip_V2E (discovery) > concordant_E > everything else
    priority = {"flip_V2E": 0, "concordant_E": 1, "flip_E2V": 2}
    novel["priority"] = novel["outcome"].map(priority).fillna(9).astype(int)
    novel = novel.sort_values(["priority", "R2"], ascending=[True, False])

    cols = ["gene", "gene_name", "outcome", "essentiality_concordance",
            "FYPOviability", "DR", "DL", "R2", "verification_essentiality",
            "DeletionLibrary_essentiality"]
    out = novel[cols].reset_index(drop=True)
    out.to_csv(OUT_DIR / "a3_novel_genes.tsv", sep="\t", index=False)
    print(f"[a3_novel_genes] {out.shape} -> {OUT_DIR/'a3_novel_genes.tsv'}")
    print(f"  {out.shape[0]} novel genes; {(out['outcome']=='flip_V2E').sum()} are flip_V2E reclassifications")
    if len(out):
        print(out.head(12).to_string(index=False))
    return out


# ---------------------------------------------------------------- A4
def a4_failure_attribution(df: pd.DataFrame) -> pd.DataFrame:
    """Classify ambiguous / discordant genes: data-quality vs biology."""
    amb = df[df["is_ambiguous"]].copy()

    has_data_flag = (
        amb["qc_flag_reason"].astype(str).str.contains("poor_fit|thin_colonies|noised_gene", na=False)
        | (amb["total_colonies"] < COLONY_LOW)
        | (amb["filtering_fraction"] < FILTER_LOW)
    )
    is_cond_dep = amb["FYPOviability"].astype(str).eq("condition-dependent")
    only_comment = amb["qc_flag_reason"].astype(str).eq("borderline_comment")

    def label(row_idx):
        data = has_data_flag.loc[row_idx]
        cond = is_cond_dep.loc[row_idx]
        comment_only = only_comment.loc[row_idx]
        if data and cond:
            return "both"
        if data:
            return "data_quality"
        if cond:
            return "biology_condition_dependent"
        if comment_only:
            return "borderline_observation"
        return "unattributed"

    amb["failure_cause"] = [label(i) for i in amb.index]

    cols = ["gene", "gene_name", "outcome", "essentiality_concordance",
            "failure_cause", "qc_flag_reason", "FYPOviability",
            "R2", "total_colonies", "filtering_fraction", "comments"]
    out = amb[cols].sort_values("failure_cause").reset_index(drop=True)
    out.to_csv(OUT_DIR / "a4_failure_attribution.tsv", sep="\t", index=False)
    print(f"\n[a4_failure_attribution] {out.shape} -> {OUT_DIR/'a4_failure_attribution.tsv'}")
    print(out["failure_cause"].value_counts().to_string())
    return out


# ---------------------------------------------------------------- A5
def a5_difference_explanation(df: pd.DataFrame) -> pd.DataFrame:
    """97 raw flips crossed with FYPOviability + colony-area trajectory."""
    flips = df[df["essentiality_concordance"].isin(["flip_V2E", "flip_E2V"])].copy()

    # Colony-area (day6) is an INDEPENDENT growth axis. Essential-gene deletions
    # do not grow (area≈0); viable ones do (area>0.1). This is orthogonal to the
    # DR/DL fit, so agreement between flip direction and day6 growth is a strong
    # cross-validation of the reclassification.
    AREA_GROW = 0.1
    flips["day6_grows"] = flips["median_area_day6"] > AREA_GROW

    def flip_reason(row):
        fypo = str(row["FYPOviability"])
        if fypo == "condition-dependent":
            return "condition_dependent_annotation"
        if row["essentiality_concordance"] == "flip_V2E" and fypo == "viable":
            return "unconditional_viable_to_essential"  # strongest reclassification
        if row["essentiality_concordance"] == "flip_E2V" and fypo == "inviable":
            return "unconditional_essential_to_viable"
        return "other"

    flips["flip_reason"] = flips.apply(flip_reason, axis=1)

    # Does colony-area growth agree with the reclassified essentiality?
    # V2E should NOT grow (area≈0); E2V SHOULD grow.
    flips["area_confirms_flip"] = np.where(
        flips["essentiality_concordance"] == "flip_V2E",
        ~flips["day6_grows"],   # V2E confirmed if it does not grow
        flips["day6_grows"],    # E2V confirmed if it grows
    )

    cross = pd.crosstab(flips["essentiality_concordance"], flips["FYPOviability"])
    print(f"\n[a5] flip × FYPOviability:")
    print(cross.to_string())

    cols = ["gene", "gene_name", "essentiality_concordance", "outcome",
            "flip_reason", "FYPOviability", "DR", "DL", "R2",
            "median_area_day6", "day6_grows", "area_confirms_flip", "is_ambiguous"]
    out = flips[cols].sort_values(["flip_reason", "R2"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(OUT_DIR / "a5_difference_explanation.tsv", sep="\t", index=False)
    print(f"\n[a5_difference_explanation] {out.shape} -> {OUT_DIR/'a5_difference_explanation.tsv'}")
    print(out["flip_reason"].value_counts().to_string())

    # Colony-area cross-validation: fraction of flips whose day6 growth agrees
    print(f"\n=== colony-area cross-validation (independent growth axis) ===")
    for grp in ["flip_V2E", "flip_E2V"]:
        g = out[out["essentiality_concordance"] == grp]
        print(f"  {grp}: {g['area_confirms_flip'].sum()}/{len(g)} confirmed by day6 area "
              f"({g['area_confirms_flip'].mean()*100:.0f}%)")

    # Highest-confidence reclassifications: unconditional V2E, clean, area confirms (no growth)
    gold = out[
        (out["flip_reason"] == "unconditional_viable_to_essential")
        & (~out["is_ambiguous"])
        & (out["area_confirms_flip"])
    ].sort_values("R2", ascending=False)
    print(f"\n=== gold-standard reclassifications (unconditional V→E, non-ambiguous, area-confirmed): {len(gold)} ===")
    if len(gold):
        print(gold[["gene", "gene_name", "DR", "R2", "median_area_day6"]].head(12).to_string(index=False))
    return out


def main() -> None:
    df = load()
    a3_novel_genes(df)
    a4_failure_attribution(df)
    a5_difference_explanation(df)


if __name__ == "__main__":
    main()
