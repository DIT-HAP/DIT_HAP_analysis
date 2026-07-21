#!/usr/bin/env python
"""
[D2b] Curation concordance: does GMM geometry recover literature subunit roles?

Cross-checks the data-driven GMM major/minor split (D2) against hand-curated
literature roles (resources/curated/complex_subunit_roles.tsv). This is the
"结合文献" validation step: GMM alone flags a geometric split, but only the
literature overlay tells us WHAT the split means (essential core vs peripheral
regulator, complex-specific vs shared subunit).

Key question per complex: does the GMM "core" (tighter component) correspond to
the literature core/complex-specific subunits, and GMM "minor" to the
peripheral/shared subunits?

Reads:
- resources/curated/complex_subunit_roles.tsv (52 curated member rows)
- results/coherence/HD_DIT_HAP/incoherence_attribution.tsv (D2 GMM output)

Writes:
- results/coherence/HD_DIT_HAP/curation_concordance.tsv (per-complex verdict)
"""

from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
CURATED = REPO / "resources/curated/complex_subunit_roles.tsv"
# D1/D2 coherence outputs live in the main tree; results/ is symlinked so either
# path resolves, but use the real main-tree path to be explicit.
MAIN = Path("/data/c/yangyusheng_optimized/DIT_HAP_analysis")
ATTR = MAIN / "results/coherence/HD_DIT_HAP/incoherence_attribution.tsv"
OUT = MAIN / "results/coherence/HD_DIT_HAP/curation_concordance.tsv"

# Literature role keywords → biological category (order matters: shared wins first)
SHARED_KEYS = ("shared with", "shared ", "TAF module", "HSA-module", "heterohexamer")
PERIPHERAL_KEYS = ("peripheral", "non-essential", "regulator", "TBP-interaction")
CORE_KEYS = ("essential", "core", "scaffold", "structural", "catalytic", "specific")


def classify_role(role: str) -> str:
    """Map a free-text literature role to {shared, peripheral, core}."""
    r = role.lower()
    if any(k.lower() in r for k in SHARED_KEYS):
        return "shared"
    if any(k.lower() in r for k in PERIPHERAL_KEYS):
        return "peripheral"
    if any(k.lower() in r for k in CORE_KEYS):
        return "core"
    return "unclear"


def main() -> None:
    cur = pd.read_csv(CURATED, sep="\t")
    cur["lit_category"] = cur["literature_role"].apply(classify_role)

    attr = pd.read_csv(ATTR, sep="\t")
    split_complexes = set(attr.loc[attr["gmm_is_split"] == True, "complex"])

    rows = []
    for cx, sub in cur.groupby("complex"):
        if cx not in split_complexes:
            continue

        gmm_minor = set(sub.loc[sub["gmm_subgroup"] == "minor", "gene_name"])
        gmm_core = set(sub.loc[sub["gmm_subgroup"] == "core", "gene_name"])

        # Literature "off-core" = shared OR peripheral (the biologically loose ones)
        lit_offcore = set(sub.loc[sub["lit_category"].isin(["shared", "peripheral"]), "gene_name"])

        minor_precision = len(gmm_minor & lit_offcore) / len(gmm_minor) if gmm_minor else float("nan")
        minor_recall = len(gmm_minor & lit_offcore) / len(lit_offcore) if lit_offcore else float("nan")

        minor_lit = sub.loc[sub["gmm_subgroup"] == "minor", "lit_category"]
        minor_dominant = minor_lit.mode().iloc[0] if len(minor_lit) else "none"

        if minor_precision != minor_precision:  # NaN
            verdict = "no_minor"
        elif minor_precision >= 0.75:
            verdict = "concordant"
        elif minor_precision >= 0.5:
            verdict = "partial"
        else:
            verdict = "discordant"

        rows.append({
            "complex": cx,
            "n_members": len(sub),
            "gmm_core": ", ".join(sorted(gmm_core)),
            "gmm_minor": ", ".join(sorted(gmm_minor)),
            "minor_dominant_role": minor_dominant,
            "minor_precision": round(minor_precision, 3) if minor_precision == minor_precision else None,
            "minor_recall": round(minor_recall, 3) if minor_recall == minor_recall else None,
            "verdict": verdict,
        })

    out = pd.DataFrame(rows).sort_values("minor_precision", ascending=False, na_position="last")
    out.to_csv(OUT, sep="\t", index=False)
    print(f"[curation_concordance] ({out.shape[0]}, {out.shape[1]}) -> {OUT}\n")
    print(out.to_string(index=False))

    print("\n=== interpretation ===")
    print(f"concordant (GMM-minor ≈ literature off-core, prec≥0.75): {(out['verdict']=='concordant').sum()}")
    print(f"partial (0.5-0.75): {(out['verdict']=='partial').sum()}")
    print(f"discordant (<0.5): {(out['verdict']=='discordant').sum()}")

    # The shared-subunit paradox: same subunits, opposite GMM label across complexes
    print("\n=== shared-module label paradox (Rvb1/Rvb2/Act1/Arp4 across Ino80 vs Swr1) ===")
    shared_mod = cur[cur["gene_name"].isin(["rvb1", "rvb2", "act1", "alp5"])]
    paradox = shared_mod[shared_mod["complex"].isin(["Ino80 complex", "Swr1 complex"])]
    print(paradox[["complex", "gene_name", "gmm_subgroup", "lit_category"]].sort_values(["gene_name", "complex"]).to_string(index=False))
    print("\nSame physical module, opposite GMM subgroup across the two complexes:")
    print("→ GMM geometry alone cannot flag shared subunits; the literature overlay is required.")


if __name__ == "__main__":
    main()
