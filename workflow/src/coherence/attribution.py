"""Diagnose WHY a complex is incoherent in DR-DL space (theme D, task D2).

An incoherent group's members are more dispersed in normalized (DR, DL/10) space
than a random gene set of equal size (coherence z-score > 0). This module scores
the candidate biological causes of that dispersion so each incoherent group can
be labelled with its most likely explanation. The diagnostic lines:

  (a) major/minor subunit split — fit a 2-component GMM to the members' normalized
      DR-DL and test whether they separate into a tight "core" + looser "minority"
      (silhouette + component-size/spread asymmetry). E.g. eIF3's essential core
      (tif301/302, DR~-1.2, DL~0) vs the dispensable regulatory eIF3e/int6 (DL~7).
  (b) shared-subunit — members that also belong to OTHER groups get pulled toward
      those groups' functional centres, inflating apparent incoherence. E.g. Swr1,
      whose members are almost all shared with NuA4 / Ino80 / HAT complexes.
  (c) paralog buffering — members with a paralog can have their deletion phenotype
      masked (DR nearer 0), pulling the group toward the WT corner and splitting it.

Not every cause is detectable from these signals: annotation artefacts (transient
members, over-broad "complex" definitions) and technical issues (sparse insertion
coverage, curve-fit quality) are NOT auto-labelled — a group with real dispersion
but none of the above signals is left as `intrinsic_heterogeneity` for manual review.

Pure functions over arrays / the coherence long table; no IO.
"""
# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations

from collections.abc import Iterable

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
GMM_RANDOM_STATE = 42
MIN_N_FOR_GMM = 6            # need enough points for a 2-cluster split to be meaningful
SILHOUETTE_SPLIT_MIN = 0.5   # >= this = a real 2-subgroup split (major/minor)


# =============================================================================
# CORE LOGIC
# =============================================================================
def major_minor_split(X: np.ndarray) -> dict:
    """Fit a 2-component GMM to normalized DR-DL and report whether they split core/minor."""
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


def member_pairs(long_table: pd.DataFrame) -> pd.DataFrame:
    """Unique (source, group_id, Systematic ID) membership rows."""
    # Sharing is per SOURCE: a gene in a go_cc term and a go_bp term belongs to two
    # different groupings, not to two alternative descriptions of one complex, so
    # cross-source membership must not read as a shared subunit. This matters on the
    # pooled table the dedup stage now uses — and it is not enough on its own there,
    # because group_id is not unique across sources either; every lookup keyed by
    # group_id must carry the source along (see `shared_subunit_fractions`).
    #
    # `source` is synthesized as "" for callers whose table lacks the column, so the
    # grouping keys are the same shape either way.
    if "source" in long_table.columns:
        return long_table[["source", "group_id", "Systematic ID"]].drop_duplicates()
    return long_table[["group_id", "Systematic ID"]].drop_duplicates().assign(source="")


def member_degrees(long_table: pd.DataFrame) -> pd.Series:
    """Per (source, gene), how many distinct group_ids contain it."""
    # The definition behind every shared-subunit number in this repo: a member is
    # "shared" exactly when its degree within its source is > 1.
    pairs = member_pairs(long_table)
    return pairs.groupby(["source", "Systematic ID"])["group_id"].nunique()


def shared_subunit_fractions(long_table: pd.DataFrame) -> dict:
    """Return {key: fraction of members that also belong to >=1 other group}.

    `key` is `(source, group_id)` when the table carries a `source` column and the
    bare `group_id` when it does not, mirroring `member_pairs`. The pair is not
    decoration on a pooled table: group_id is NOT unique across sources (173 of
    them, e.g. GO:0032040, appear in both go_cc and go_macrocomplex), so a
    group_id-keyed dict lets whichever source is iterated last silently overwrite
    the other's fraction, and every lookup then returns the wrong source's number.
    """
    # One pass for every group at once, which is what the compute stage wants: the
    # per-group loop this replaces scans the whole table per group, and go_bp's
    # 120k-row table with 3.7k groups made that cost ~40s.
    pairs = member_pairs(long_table)
    degrees = member_degrees(long_table).rename("degree")
    annotated = pairs.merge(
        degrees, left_on=["source", "Systematic ID"], right_index=True, how="left"
    )
    shared = annotated[annotated["degree"] > 1]

    n_shared = shared.groupby(["source", "group_id"])["Systematic ID"].nunique()
    n_members = annotated.groupby(["source", "group_id"])["Systematic ID"].nunique()
    ratios = n_shared.reindex(n_members.index).fillna(0) / n_members
    if "source" in long_table.columns:
        return {(str(source), group_id): float(value) for (source, group_id), value in ratios.items()}
    return {group_id: float(value) for (_source, group_id), value in ratios.items()}


def shared_subunits(
    long_table: pd.DataFrame, group_id: str, source: str | None = None
) -> pd.DataFrame:
    """List one group's members that also belong to OTHER groups of the same source."""
    # `long_table` is either a single source's coherence long-table (contract columns
    # group_id, group_name, "Systematic ID") or the pooled all-sources table.
    #
    # `source` scopes the lookups. It is optional because a single-source table
    # makes the group_id unique and the source inferable; on a POOLED table it must
    # be passed, since group_id alone is ambiguous there (see
    # `shared_subunit_fractions`) and inferring it picks an arbitrary source.
    #
    # "Other" is keyed on the stable `group_id` and scoped to the member's own
    # source (see `member_pairs`), matching how `sources.py` dedups and how the
    # compute stage forms groups. Two distinct term IDs can share a group_name, and
    # keying on the name instead would make a member of such a pair read as
    # unshared. The reported `other_groups` stays a list of NAMES (that is what a
    # reader wants to see), so `n_other_groups` counts distinct other names; only
    # the membership test is ID-keyed.
    cols = ["Systematic ID", "n_other_groups", "other_groups"]
    candidates = member_pairs(long_table)
    if "source" in long_table.columns:
        if source is None:
            group_source = long_table.loc[long_table["group_id"] == group_id, "source"]
            if group_source.empty:
                return pd.DataFrame(columns=cols)
            source = group_source.iloc[0]
        candidates = candidates[candidates["source"] == source]

    members = set(candidates.loc[candidates["group_id"] == group_id, "Systematic ID"])
    if not members:
        return pd.DataFrame(columns=cols)
    sub = (
        candidates.loc[
            candidates["Systematic ID"].isin(members), ["group_id", "Systematic ID"]
        ]
        .merge(
            long_table[["group_id", "group_name"]].drop_duplicates("group_id"),
            on="group_id", how="left",
        )
        .drop_duplicates(["group_id", "Systematic ID"])
    )
    sub = sub[sub["group_id"] != group_id]
    if sub.empty:
        return pd.DataFrame(columns=cols)
    other_names = sub.groupby("Systematic ID")["group_name"].agg(lambda names: sorted(set(names)))
    rows = [
        {"Systematic ID": gene, "n_other_groups": len(names), "other_groups": "; ".join(names)}
        for gene, names in other_names.items()
    ]
    return pd.DataFrame(rows).sort_values("n_other_groups", ascending=False)


def paralog_ids_from_features(features: pd.DataFrame) -> set[str]:
    """Genes with >=1 paralog, from the feature matrix's `paralog_count` column."""
    # Whichever source fills that column is `features.paralog_source` in
    # config/analysis.yaml (the deletion library, currently); nothing here picks a
    # source of its own, so this set and the ML feature cannot drift apart. A gene
    # absent from the matrix (not a coding gene, or outside the deletion library)
    # reads as "no paralog", which is what the column itself says.
    return set(features.loc[features["paralog_count"] > 0, "gene_systematic_id"].astype(str))


def paralog_fraction(members: Iterable[str], paralog_ids: set[str]) -> float:
    """Fraction of `members` that have >=1 paralog (present in paralog_ids)."""
    # A high paralog fraction flags a group whose members' deletion phenotypes may
    # be buffered by redundant paralogs (dampened DR), a candidate incoherence
    # cause.
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
    """Combine the diagnostics into a single attribution label (priority ladder)."""
    # Priority, most-specific/structural first:
    #   1. major/minor GMM split AND high shared fraction -> `conditional_module`
    #      (a distinct sub-module that is also cross-shared, e.g. CLRC's shared CRL4
    #      scaffold + the dispensable heterochromatin-silencing module);
    #   2. major/minor GMM split alone -> `major_minor_split`;
    #   3. high shared fraction alone -> `shared_subunits`;
    #   4. high paralog fraction -> `paralog_buffered`;
    #   5. too few members to have fit a GMM -> `data_limited`;
    #   6. otherwise -> `intrinsic_heterogeneity` (real spread, no detected cause;
    #      may also be an annotation/technical artefact — flagged for manual review).
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
