"""Diagnose WHY a complex is incoherent in DR-DL space (theme D, task D2).

An incoherent group's members are more dispersed in normalized (DR, DL/10) space
than a random gene set of equal size (coherence z-score > 0). This module scores
the candidate biological causes of that dispersion so each incoherent group can
be labelled with its most likely explanation. The diagnostic lines:

  (a) major/minor subunit split — fit a 2-component GMM to the members' normalized
      DR-DL and test whether they separate into a tight "core" + looser "minority"
      (silhouette + component-size/spread asymmetry). E.g. eIF3's essential core
      (tif301/302, DR~1.2, DL~0) vs the dispensable regulatory eIF3e/int6 (DL~7).
  (b) shared-subunit — members that also belong to OTHER groups get pulled toward
      those groups' functional centres, inflating apparent incoherence. E.g. Swr1,
      whose members are almost all shared with NuA4 / Ino80 / HAT complexes.
  (c) paralog buffering — members with a paralog can have their deletion phenotype
      masked (lower DR), pulling the group toward the WT corner and splitting it.

Not every cause is detectable from these signals: annotation artefacts (transient
members, over-broad "complex" definitions) and technical issues (sparse insertion
coverage, curve-fit quality) are NOT auto-labelled — a group with real dispersion
but none of the above signals is left as `intrinsic_heterogeneity` for manual review.

Pure functions over arrays / the coherence long table; no IO.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

GMM_RANDOM_STATE = 42
MIN_N_FOR_GMM = 6            # need enough points for a 2-cluster split to be meaningful
SILHOUETTE_SPLIT_MIN = 0.5   # >= this = a real 2-subgroup split (major/minor)


def major_minor_split(X: np.ndarray) -> dict:
    """Fit 2-component GMM to normalized DR-DL; report whether members split core/minor.

    Returns silhouette, per-component sizes, the tighter component's index ("core"),
    and a boolean is_split (silhouette high enough to call a genuine major/minor split).
    Returns is_split=False with reason for n<MIN_N_FOR_GMM.
    """
    n = X.shape[0]
    if n < MIN_N_FOR_GMM:
        return {"is_split": False, "reason": f"n<{MIN_N_FOR_GMM}", "silhouette": np.nan,
                "labels": None, "core_label": None, "component_sizes": None}
    gmm = GaussianMixture(n_components=2, random_state=GMM_RANDOM_STATE, n_init=5)
    labels = gmm.fit_predict(X)
    if len(np.unique(labels)) < 2:
        return {"is_split": False, "reason": "degenerate_single_component", "silhouette": np.nan,
                "labels": labels, "core_label": None, "component_sizes": None}
    sil = float(silhouette_score(X, labels))
    # "core" = the tighter (lower mean-distance-to-own-centroid) component
    spreads = {}
    sizes = {}
    for lab in (0, 1):
        pts = X[labels == lab]
        sizes[lab] = int(pts.shape[0])
        spreads[lab] = float(np.mean(np.linalg.norm(pts - pts.mean(axis=0), axis=1)))
    core_label = min(spreads, key=spreads.get)
    return {
        "is_split": sil >= SILHOUETTE_SPLIT_MIN,
        "reason": "gmm_2comp",
        "silhouette": sil,
        "labels": labels,
        "core_label": core_label,
        "component_sizes": sizes,
        "core_spread": spreads[core_label],
        "minor_spread": spreads[1 - core_label],
    }


def shared_subunits(long_table: pd.DataFrame, group_id: str) -> pd.DataFrame:
    """For one group, list members that also belong to OTHER groups of the same table.

    long_table: a single source's coherence long-table (contract columns
    group_id, group_name, "Systematic ID"). Returns one row per shared member with
    the count + names of the OTHER groups it participates in, sorted by that count
    (empty df if none shared). Keyed on the stable `group_id`, matching the
    coherence long-table contract.
    """
    cols = ["Systematic ID", "n_other_groups", "other_groups"]
    members = set(long_table.loc[long_table["group_id"] == group_id, "Systematic ID"])
    if not members:
        return pd.DataFrame(columns=cols)
    sub = long_table[long_table["Systematic ID"].isin(members)].drop_duplicates(
        ["group_id", "Systematic ID"]
    )
    rows = []
    for gene, grp in sub.groupby("Systematic ID"):
        others = sorted(set(grp["group_name"]) - set(
            long_table.loc[long_table["group_id"] == group_id, "group_name"]
        ))
        if others:
            rows.append({"Systematic ID": gene, "n_other_groups": len(others),
                         "other_groups": "; ".join(others)})
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows).sort_values("n_other_groups", ascending=False)


def shared_fraction(long_table: pd.DataFrame, group_id: str) -> float:
    """Fraction of a group's members that are shared with >=1 other group."""
    members = set(long_table.loc[long_table["group_id"] == group_id, "Systematic ID"])
    if not members:
        return np.nan
    shared = set(shared_subunits(long_table, group_id)["Systematic ID"])
    return len(shared) / len(members)


def paralog_fraction(members: Iterable[str], paralog_ids: set[str]) -> float:
    """Fraction of `members` that have >=1 paralog (present in paralog_ids).

    A high paralog fraction flags a group whose members' deletion phenotypes may be
    buffered by redundant paralogs (dampened DR), a candidate incoherence cause.
    Returns NaN for an empty member set.
    """
    members = list(members)
    if not members:
        return np.nan
    return sum(1 for m in members if m in paralog_ids) / len(members)


def attribute_incoherence(
    split: dict,
    shared_frac: float,
    paralog_frac: float = np.nan,
    shared_frac_threshold: float = 0.5,
    paralog_frac_threshold: float = 0.5,
) -> str:
    """Combine the diagnostics into a single attribution label (priority ladder).

    Priority, most-specific/structural first:
      1. major/minor GMM split AND high shared fraction -> `conditional_module`
         (a distinct sub-module that is also cross-shared, e.g. CLRC's shared CRL4
         scaffold + the dispensable heterochromatin-silencing module);
      2. major/minor GMM split alone -> `major_minor_split`;
      3. high shared fraction alone -> `shared_subunits`;
      4. high paralog fraction -> `paralog_buffered`;
      5. too few members to have fit a GMM -> `data_limited`;
      6. otherwise -> `intrinsic_heterogeneity` (real spread, no detected cause;
         may also be an annotation/technical artefact — flagged for manual review).
    """
    is_split = bool(split.get("is_split"))
    high_shared = pd.notna(shared_frac) and shared_frac >= shared_frac_threshold
    high_paralog = pd.notna(paralog_frac) and paralog_frac >= paralog_frac_threshold
    if is_split and high_shared:
        return "conditional_module"
    if is_split:
        return "major_minor_split"
    if high_shared:
        return "shared_subunits"
    if high_paralog:
        return "paralog_buffered"
    if split.get("reason", "").startswith("n<"):
        return "data_limited"
    return "intrinsic_heterogeneity"
