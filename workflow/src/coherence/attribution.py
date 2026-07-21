"""Diagnose WHY a complex is incoherent in DR-DL space (theme D, task D2).

Two independent diagnostic lines:
  (a) major/minor subunit split — fit a 2-component GMM to the members' normalized
      DR-DL and test whether they separate into a tight "core" + looser "minority"
      (silhouette + component-size/spread asymmetry). Cross-checked against a curated
      literature-role table when available.
  (b) shared-subunit — members that also belong to OTHER physical complexes drag the
      centroid apart, inflating apparent incoherence.

Pure functions over arrays / the foundation long table; no IO.
"""
from __future__ import annotations

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


def shared_subunits(long_cc: pd.DataFrame, complex_name: str) -> pd.DataFrame:
    """For one CC complex, list members that also belong to other CC complexes.

    long_cc: the GO_CC_complex slice of the foundation long table
             (columns group_name, gene). Returns one row per shared member with the
             list of OTHER complexes it participates in (empty df if none shared).
    """
    members = set(long_cc.loc[long_cc["group_name"] == complex_name, "gene"])
    if not members:
        return pd.DataFrame(columns=["gene", "n_other_complexes", "other_complexes"])
    sub = long_cc[long_cc["gene"].isin(members)]
    rows = []
    for gene, grp in sub.groupby("gene"):
        others = sorted(set(grp["group_name"]) - {complex_name})
        if others:
            rows.append({"gene": gene, "n_other_complexes": len(others),
                         "other_complexes": "; ".join(others)})
    return pd.DataFrame(rows).sort_values("n_other_complexes", ascending=False) if rows else \
        pd.DataFrame(columns=["gene", "n_other_complexes", "other_complexes"])


def shared_fraction(long_cc: pd.DataFrame, complex_name: str) -> float:
    """Fraction of a complex's members that are shared with >=1 other CC complex."""
    members = set(long_cc.loc[long_cc["group_name"] == complex_name, "gene"])
    if not members:
        return np.nan
    shared = set(shared_subunits(long_cc, complex_name)["gene"])
    return len(shared) / len(members)


def attribute_incoherence(
    split: dict, shared_frac: float, shared_frac_threshold: float = 0.5
) -> str:
    """Combine the two diagnostics into a single attribution label.

    Priority: a genuine major/minor GMM split explains it structurally; else a high
    shared-subunit fraction explains it as cross-complex contamination; else the
    spread is intrinsic (real biological heterogeneity) or data-limited.
    """
    if split.get("is_split"):
        return "major_minor_split"
    if pd.notna(shared_frac) and shared_frac >= shared_frac_threshold:
        return "shared_subunits"
    if split.get("reason", "").startswith("n<"):
        return "data_limited"
    return "intrinsic_heterogeneity"
