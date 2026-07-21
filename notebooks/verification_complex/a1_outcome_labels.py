"""[A1] Derive verification outcome labels — the prerequisite table for A2-A5.

Reads results/verification/HD_DIT_HAP/verification_master.tsv, writes
verification_outcome_table.tsv with:
  - essentiality_concordance : raw verified-vs-library call (kept even when ambiguous,
    so A2/A5 can still see what a flagged flip WOULD have been):
      concordant_E, concordant_V, flip_V2E, flip_E2V, no_library_ref
  - is_ambiguous / qc_flag_reason : data-quality doubt (noised list, comment keywords,
    poor fit, thin colony support)
  - outcome : final label = 'ambiguous' when is_ambiguous, else essentiality_concordance.
    Six values total (the 5 designed + no_library_ref for the 7 uncharacterised genes).

Run with the `data_analysis` conda env. Iterate here, then port to notebook +
(later) workflow/scripts/verification/.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "results/verification/HD_DIT_HAP"
MASTER = OUT_DIR / "verification_master.tsv"

# The 18 data-quality exclusions hand-curated in the source notebook
# (compare_with_deletion_library.ipynb cell 19). Matched on gene_name.
NOISED_GENES = [
    "pyp3", "mam2", "hul6", "fub1", "csn5", "rdl1", "ymr1", "meu31", "msn5",
    "SPAC4C5.03", "SPAC977.03", "mal1", "num1", "dbl8", "fmo1", "erm2", "thp1", "pmc3",
]

# Chinese comment fragments that flag a borderline / low-confidence wet-lab readout.
# Covers: faint/partial visibility, size caveats, and mixed/inconsistent growth
# ("one grows one doesn't", "very few grow", irregular shape). Comments are sparse
# (81/411) and almost always note a caveat — has_any_comment (below) keeps the looser
# signal available so A4 can choose strict-vs-loose.
AMBIGUOUS_COMMENT_KEYWORDS = [
    "轻微可见", "可见但", "很小", "小",           # faint / small colonies
    "一大一小", "一个长一个不长", "不都是", "一半",  # mixed / inconsistent growth
    "极少数", "少数", "不规则", "形状",             # sparse / irregular
    "难长",                                        # hard to grow after replica-plating (assay-technical)
    "部分", "几列", "?", "？",                      # partial / uncertain marker
]

# QC thresholds (justified against the observed distributions, see A1 exploration).
R2_MIN = 0.80              # median R2 ~0.98; below this the curve fit is poor
TOTAL_COLONIES_MIN = 50    # min observed 48, 25th pct 52; thin support
FILTERING_FRACTION_MIN = 0.5  # min 0.08, 25th pct 0.88; most colonies filtered out


def essentiality_concordance(row: pd.Series) -> str:
    """verified vs deletion-library essentiality; no_library_ref when library has no call."""
    ve = str(row["verification_essentiality"]).strip()
    dl_raw = row["DeletionLibrary_essentiality"]
    dl = str(dl_raw).strip()
    if pd.isna(dl_raw) or dl not in {"E", "V"}:
        return "no_library_ref"
    if ve == dl:
        return "concordant_E" if ve == "E" else "concordant_V"
    return "flip_V2E" if (dl == "V" and ve == "E") else "flip_E2V"


def qc_reasons(row: pd.Series) -> list[str]:
    """Collect data-quality doubt flags for one gene (empty list = clean)."""
    reasons: list[str] = []
    if str(row.get("gene_name")).strip() in NOISED_GENES:
        reasons.append("noised_gene")
    comment = row.get("comments")
    if isinstance(comment, str) and any(k in comment for k in AMBIGUOUS_COMMENT_KEYWORDS):
        reasons.append("borderline_comment")
    r2 = row.get("R2")
    if pd.notna(r2) and r2 < R2_MIN:
        reasons.append("poor_fit")
    tc = row.get("total_colonies")
    if pd.notna(tc) and tc < TOTAL_COLONIES_MIN:
        reasons.append("thin_colonies")
    ff = row.get("filtering_fraction")
    if pd.notna(ff) and ff < FILTERING_FRACTION_MIN:
        reasons.append("low_filter_fraction")
    return reasons


def main() -> None:
    m = pd.read_csv(MASTER, sep="\t")

    m["essentiality_concordance"] = m.apply(essentiality_concordance, axis=1)
    reason_lists = m.apply(qc_reasons, axis=1)
    m["qc_flag_reason"] = reason_lists.apply(lambda xs: ";".join(xs))
    m["is_ambiguous"] = reason_lists.apply(bool)
    # Looser signal kept alongside: any wet-lab comment at all (A4 can widen to this).
    m["has_any_comment"] = m["comments"].notna()
    m["outcome"] = m.apply(
        lambda r: "ambiguous" if r["is_ambiguous"] else r["essentiality_concordance"], axis=1
    )

    out_path = OUT_DIR / "verification_outcome_table.tsv"
    m.to_csv(out_path, sep="\t", index=False)
    print(f"[verification_outcome_table] {m.shape} -> {out_path}\n")

    print("=== essentiality_concordance (raw, ignores ambiguity) ===")
    print(m["essentiality_concordance"].value_counts().to_string())
    print("\n=== outcome (final, ambiguous overrides) ===")
    print(m["outcome"].value_counts().to_string())
    print(f"\nambiguous total: {m['is_ambiguous'].sum()}")
    print("=== qc_flag_reason breakdown (a gene may have several) ===")
    all_reasons = m["qc_flag_reason"].str.split(";").explode()
    print(all_reasons[all_reasons != ""].value_counts().to_string())
    print("\n=== how ambiguity is distributed across the raw concordance calls ===")
    print(pd.crosstab(m["essentiality_concordance"], m["is_ambiguous"]).to_string())


if __name__ == "__main__":
    main()
