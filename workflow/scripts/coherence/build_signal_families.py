#!/usr/bin/env python3

"""
Coherence Signal Families (display layer above the Jaccard de-duplication)
=========================================================================

The Jaccard de-duplication in ``deduplicate_terms.py`` collapses terms whose
*scored member sets* are near-identical. It deliberately cannot merge a small,
specific term into a much larger parent: Jaccard = |A∩B|/|A\u222aB|, so a 41-gene
child fully contained in a 218-gene parent scores only 41/218 = 0.188. That is
the right behaviour for "same member set" but not for "same biological process
seen at a different grain" — e.g. ``small-subunit processome`` (41 genes),
``preribosome`` (108), ``ribosome biogenesis`` (207) and the two KEGG ribosome
entries (155 / 53) are one pathway family, yet the 0.65 Jaccard cut leaves them
as five separate representatives.

This stage adds a SECOND, non-replacing layer on top of the existing clusters:

- A directed containment edge child -> parent requires
  ``overlap/min >= overlap_min`` (|A∩B| / min(|A|,|B|)) and
  ``size ratio >= ratio_min`` (min(|A|,|B|) / max(|A|,|B|)).
- Each child picks at most ONE parent: the smallest qualifying parent (the most
  specific one). Sizes strictly increase along a parent edge, so the graph is a
  forest — no cycles, and no connected-component chaining. This is the reason
  the pass is directed: undirected connected components over the same edges
  chain into 114-635-cluster blobs at usable settings (measured 2026-10-22).
- A family is a tree of clusters; the family representative is the member
  cluster representative with the best evidence (smallest q_value_coherent,
  then most-negative median_pairwise_distance_z, then greater dag_depth, then
  smaller n_scored_members), mirroring deduplicate_terms.py's selection.

It does NOT alter ``redundancy_cluster``, ``is_representative``, q values or any
statistic. It only annotates each cluster/term with a ``family_id`` and the
family representative, so a table or figure can choose its display grain.

Input
-----
- --deduplicated: dedup/coherence_terms_deduplicated.tsv (all terms, with
  redundancy_cluster, is_representative, q_value_coherent, ...).
- --group-members: dedup/coherence_group_members_long.tsv (cluster -> genes).

Output
------
- --output-assignments: one row per de-duplication cluster with its parent edge
  (overlap, ratio) and family columns.
- --output-terms: all input terms + family columns, for joining a q-sorted
  coherent table.
- --output-families: one row per family (number of clusters/terms, gene union,
  representative, member cluster representatives).
- --output-preview: one row per family that contains coherent representatives
  (q_value_coherent <= --coherent-q-max AND z < --coherent-z), sorted by best z,
  with the coherent member terms listed — the "one row per signal" view.
- --output-workbook (optional): a `coherence_cohorts.xlsx`-shaped workbook whose
  rows are FAMILIES, not terms. Sheets `coherent` / `incoherent` / `all_families`
  / `thresholds`; a family appears in a cohort sheet when ANY of its dedup
  representatives falls in that cohort (so a mixed family appears in both, labelled
  `mixed` in all_families). With `--workbook-scope per_cohort`, families are built
  separately inside each cohort instead, so every family is cohort-pure (no mixed
  families) and `all_families` concatenates the two cohort partitions.
  The row's z/q are the family representative's values,
  while family_min_z / family_max_z / n_*_representatives show the spread inside
  the family. `n_scored_members` is the family gene union and
  `member_cluster_representatives` lists what was collapsed (with each cluster's z).

Usage
-----
    python build_signal_families.py \\
        --deduplicated results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv \\
        --group-members results/3a_coherence/{dataset}/dedup/coherence_group_members_long.tsv \\
        --overlap-min 0.8 --ratio-min 0.2 \\
        --output-assignments tmp/coherence_signal_families/coherence_family_assignments.tsv \\
        --output-terms tmp/coherence_signal_families/coherence_terms_with_family.tsv \\
        --output-families tmp/coherence_signal_families/coherence_family_representatives.tsv \\
        --output-preview tmp/coherence_signal_families/coherent_terms_by_family.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-10-22
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys

# 3. Third-party Imports
from loguru import logger

# 2. Data Processing Imports
import pandas as pd

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
DEFAULT_OVERLAP_MIN = 0.8
DEFAULT_RATIO_MIN = 0.2
DEFAULT_COHERENT_Q_MAX = 0.05
DEFAULT_COHERENT_Z = -2.0
_KEY = ["source", "group_id"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FamilyConfig:
    """Inputs, outputs, and the two containment thresholds."""

    deduplicated: Path
    group_members: Path
    output_assignments: Path
    output_terms: Path
    output_families: Path
    output_preview: Path
    output_workbook: Path | None = None
    overlap_min: float = DEFAULT_OVERLAP_MIN
    ratio_min: float = DEFAULT_RATIO_MIN
    coherent_q_max: float = DEFAULT_COHERENT_Q_MAX
    coherent_z: float = DEFAULT_COHERENT_Z
    incoherent_z: float = 0.5
    workbook_scope: str = "global"  # "global" | "per_cohort"

    def validate(self) -> None:
        """Raise ValueError on bad inputs/params, then make output dirs."""
        for path in [self.deduplicated, self.group_members]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if not 0.0 < self.overlap_min <= 1.0:
            raise ValueError(f"overlap_min must be in (0, 1]: {self.overlap_min}")
        if not 0.0 < self.ratio_min <= 1.0:
            raise ValueError(f"ratio_min must be in (0, 1]: {self.ratio_min}")
        if self.workbook_scope not in ("global", "per_cohort"):
            raise ValueError(
                f"workbook_scope must be 'global' or 'per_cohort': {self.workbook_scope!r}"
            )
        outputs = [
            self.output_assignments,
            self.output_terms,
            self.output_families,
            self.output_preview,
        ]
        if self.output_workbook is not None:
            outputs.append(self.output_workbook)
        for output in outputs:
            output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
    )


# =============================================================================
# CORE LOGIC — loading
# =============================================================================
def member_set(value: object) -> set[str]:
    """Parse a TSV member list (comma+space joined) into a set."""
    if value is None or isinstance(value, float):
        return set()
    return {gene.strip() for gene in str(value).split(",") if gene.strip()}


def load_deduplicated(path: Path) -> pd.DataFrame:
    """Read the deduplicated table and normalise its boolean column."""
    table = pd.read_csv(path, sep="\t")
    table["is_representative"] = table["is_representative"].astype(str).str.lower().eq("true")
    return table


def cluster_gene_sets(group_members: Path) -> dict[str, set[str]]:
    """Cluster id -> union of its member genes, from the long table."""
    table = pd.read_csv(group_members, sep="\t")
    return table.groupby("redundancy_cluster")["gene"].apply(set).to_dict()


def cluster_metadata(table: pd.DataFrame) -> pd.DataFrame:
    """One row per de-duplication cluster, with its representative's evidence."""
    representatives = table.loc[
        table["is_representative"],
        [
            "redundancy_cluster",
            *_KEY,
            "group_name",
            "q_value_coherent",
            "median_pairwise_distance_z",
            "dag_depth",
            "n_scored_members",
        ],
    ].set_index("redundancy_cluster")
    counts = table.groupby("redundancy_cluster").size().rename("n_terms")
    meta = representatives.join(counts)
    if meta.index.has_duplicates:
        raise ValueError("A redundancy_cluster has more than one representative row")
    return meta


# =============================================================================
# CORE LOGIC — containment edges
# =============================================================================
def candidate_overlaps(sets: dict[str, set[str]]) -> list[tuple[str, str, int]]:
    """Every cluster pair that shares >=1 gene, with its intersection size."""
    # The inverted index avoids the full O(n^2) sweep: only co-occurring pairs
    # can have a non-zero overlap. Same construction as deduplicate_terms.py.
    index: dict[str, list[str]] = defaultdict(list)
    for cluster, genes in sets.items():
        for gene in genes:
            index[gene].append(cluster)
    pairs: dict[tuple[str, str], int] = defaultdict(int)
    for clusters in index.values():
        for first in range(len(clusters)):
            for second in range(first + 1, len(clusters)):
                a, b = clusters[first], clusters[second]
                key = (a, b) if a < b else (b, a)
                pairs[key] += 1
    return [(a, b, inter) for (a, b), inter in pairs.items()]


def assign_parents(sets: dict[str, set[str]], overlap_min: float, ratio_min: float) -> pd.DataFrame:
    """Give every cluster at most one larger parent, with its edge statistics."""
    sizes = {cluster: len(genes) for cluster, genes in sets.items()}
    best: dict[str, tuple[str, float, float]] = {}
    for a, b, inter in candidate_overlaps(sets):
        size_a, size_b = sizes[a], sizes[b]
        if size_a == size_b:
            continue
        child, parent = (a, b) if size_a < size_b else (b, a)
        small, large = min(size_a, size_b), max(size_a, size_b)
        overlap = inter / small
        ratio = small / large
        if overlap < overlap_min or ratio < ratio_min:
            continue
        current = best.get(child)
        # The smallest qualifying parent is the most specific one. Ties go to the
        # higher overlap, then the smaller cluster id, so the result is stable.
        candidate = (parent, overlap, ratio)
        if current is None or (sizes[parent], -overlap, parent) < (
            sizes[current[0]],
            -current[1],
            current[0],
        ):
            best[child] = candidate
    rows = [
        {
            "redundancy_cluster": child,
            "parent_cluster": parent,
            "parent_overlap": overlap,
            "parent_ratio": ratio,
        }
        for child, (parent, overlap, ratio) in best.items()
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "redundancy_cluster",
            "parent_cluster",
            "parent_overlap",
            "parent_ratio",
        ],
    )


def find_root(parent: dict[str, str], node: str) -> str:
    """Union-find root with path compression."""
    while parent[node] != node:
        parent[node] = parent[parent[node]]
        node = parent[node]
    return node


def family_components(
    clusters: list[str], edges: pd.DataFrame
) -> tuple[dict[str, str], dict[str, str]]:
    """Cluster -> family root, plus the child -> parent map actually used."""
    parent = {cluster: cluster for cluster in clusters}
    child_to_parent = dict(zip(edges["redundancy_cluster"], edges["parent_cluster"], strict=False))
    for child, above in child_to_parent.items():
        child_root, parent_root = find_root(parent, child), find_root(parent, above)
        if child_root != parent_root:
            parent[child_root] = parent_root
    roots = {cluster: find_root(parent, cluster) for cluster in clusters}
    return roots, child_to_parent


# =============================================================================
# CORE LOGIC — family tables
# =============================================================================
def family_representative(members: pd.DataFrame) -> pd.Series:
    """The member cluster representative with the best evidence."""
    # Same ordering as deduplicate_terms.py::pick_representative, applied to one
    # representative row per cluster instead of to individual terms.
    ordered = members.sort_values(
        by=["q_value_coherent", "median_pairwise_distance_z", "dag_depth", "n_scored_members"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    return ordered.iloc[0]


def _family_membership(
    meta: pd.DataFrame, clusters: list[str], roots: dict[str, str], sets: dict[str, set[str]]
) -> tuple[list[dict], dict[str, str], dict[str, pd.Series]]:
    """Family rows plus the cluster -> family id and cluster -> representative maps."""
    family_members: dict[str, list[str]] = defaultdict(list)
    for cluster in clusters:
        family_members[roots[cluster]].append(cluster)
    family_rows = []
    cluster_to_family: dict[str, str] = {}
    cluster_to_family_rep: dict[str, pd.Series] = {}
    for members in family_members.values():
        member_meta = meta.loc[sorted(members)]
        representative = family_representative(member_meta)
        family_id = f"family:{representative['source']}:{representative['group_id']}"
        genes = set().union(*(sets[cluster] for cluster in members))
        n_terms = int(member_meta["n_terms"].sum())
        for cluster in members:
            cluster_to_family[cluster] = family_id
            cluster_to_family_rep[cluster] = representative
        family_rows.append(
            {
                "family_id": family_id,
                "family_representative_source": representative["source"],
                "family_representative_group_id": representative["group_id"],
                "family_representative_name": representative["group_name"],
                "n_clusters": len(members),
                "n_terms": n_terms,
                "gene_n": len(genes),
                "best_q_value_coherent": float(member_meta["q_value_coherent"].min()),
                "best_median_pairwise_distance_z": float(
                    member_meta["median_pairwise_distance_z"].min()
                ),
                "n_dedup_representatives": len(members),
                "member_cluster_representatives": " | ".join(
                    member_meta["group_name"].astype(str).tolist()
                ),
            }
        )
    return family_rows, cluster_to_family, cluster_to_family_rep


def _family_table(family_rows: list[dict]) -> pd.DataFrame:
    """The family table, largest families first."""
    return (
        pd.DataFrame(family_rows)
        .sort_values(
            by=["n_clusters", "best_q_value_coherent", "family_id"],
            ascending=[False, True, True],
        )
        .reset_index(drop=True)
    )


def _family_assignments(
    meta: pd.DataFrame,
    edges: pd.DataFrame,
    clusters: list[str],
    sets: dict[str, set[str]],
    cluster_to_family: dict[str, str],
    cluster_to_family_rep: dict[str, pd.Series],
    family_sizes: pd.DataFrame,
) -> pd.DataFrame:
    """One row per de-duplication cluster with its parent edge and family columns."""
    assignments = meta.reset_index()
    assignments = assignments.rename(
        columns={
            "source": "cluster_representative_source",
            "group_id": "cluster_representative_group_id",
            "group_name": "cluster_representative_name",
            "q_value_coherent": "cluster_best_q_value_coherent",
            "median_pairwise_distance_z": "cluster_best_median_pairwise_distance_z",
            "n_terms": "cluster_n_terms",
        }
    )
    assignments["cluster_gene_n"] = assignments["redundancy_cluster"].map(
        {cluster: len(sets[cluster]) for cluster in clusters}
    )
    assignments = assignments.merge(edges, on="redundancy_cluster", how="left")
    assignments["family_id"] = assignments["redundancy_cluster"].map(cluster_to_family)
    assignments["family_representative_source"] = assignments["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["source"] for c in clusters}
    )
    assignments["family_representative_group_id"] = assignments["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["group_id"] for c in clusters}
    )
    assignments["family_representative_name"] = assignments["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["group_name"] for c in clusters}
    )
    assignments = assignments.join(family_sizes, on="family_id")
    assignments["is_family_representative"] = (
        assignments["redundancy_cluster"].map({c: cluster_to_family_rep[c].name for c in clusters})
        == assignments["redundancy_cluster"]
    )
    return assignments.sort_values(
        by=["family_id", "cluster_best_q_value_coherent", "redundancy_cluster"]
    ).reset_index(drop=True)


def _family_terms(
    table: pd.DataFrame,
    clusters: list[str],
    cluster_to_family: dict[str, str],
    cluster_to_family_rep: dict[str, pd.Series],
    family_sizes: pd.DataFrame,
) -> pd.DataFrame:
    """All input terms with the family columns joined in (key columns + family only)."""
    terms = table[
        [
            "source",
            "group_id",
            "group_name",
            "n_scored_members",
            "median_pairwise_distance_z",
            "q_value_coherent",
            "redundancy_cluster",
            "cluster_size",
            "is_representative",
        ]
    ].copy()
    terms["family_id"] = terms["redundancy_cluster"].map(cluster_to_family)
    terms["family_representative_source"] = terms["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["source"] for c in clusters}
    )
    terms["family_representative_group_id"] = terms["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["group_id"] for c in clusters}
    )
    terms["family_representative_name"] = terms["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["group_name"] for c in clusters}
    )
    terms = terms.join(family_sizes, on="family_id")
    return terms.sort_values(
        by=["q_value_coherent", "median_pairwise_distance_z", "group_id"]
    ).reset_index(drop=True)


def _coherent_preview(
    table: pd.DataFrame,
    clusters: list[str],
    cluster_to_family: dict[str, str],
    cluster_to_family_rep: dict[str, pd.Series],
    config: FamilyConfig,
) -> pd.DataFrame:
    """One row per family with at least one coherent representative."""
    coherent = table.loc[
        (table["q_value_coherent"] <= config.coherent_q_max)
        & (table["median_pairwise_distance_z"] < config.coherent_z)
    ].copy()
    coherent["family_id"] = coherent["redundancy_cluster"].map(cluster_to_family)
    coherent["family_representative_name"] = coherent["redundancy_cluster"].map(
        {c: cluster_to_family_rep[c]["group_name"] for c in clusters}
    )
    preview_rows = []
    for family_id, sub in coherent.groupby("family_id"):
        representatives = sub.loc[sub["is_representative"]].sort_values(
            by=["median_pairwise_distance_z", "q_value_coherent"]
        )
        pieces = [
            f"{row['group_name']} (z={row['median_pairwise_distance_z']:.2f})"
            for _, row in representatives.iterrows()
        ]
        preview_rows.append(
            {
                "family_id": family_id,
                "family_representative_name": sub["family_representative_name"].iloc[0],
                "n_coherent_representatives": len(representatives),
                "n_coherent_terms": len(sub),
                "best_median_pairwise_distance_z": float(sub["median_pairwise_distance_z"].min()),
                "best_q_value_coherent": float(sub["q_value_coherent"].min()),
                "coherent_member_terms": " | ".join(pieces),
            }
        )
    if preview_rows:
        return (
            pd.DataFrame(preview_rows)
            .sort_values(by=["best_median_pairwise_distance_z", "family_id"])
            .reset_index(drop=True)
        )
    return pd.DataFrame(
        columns=[
            "family_id",
            "family_representative_name",
            "n_coherent_representatives",
            "n_coherent_terms",
            "best_median_pairwise_distance_z",
            "best_q_value_coherent",
            "coherent_member_terms",
        ]
    )


def build_family_tables(
    table: pd.DataFrame, sets: dict[str, set[str]], config: FamilyConfig
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (assignments, terms_with_family, families, coherent_preview)."""
    meta = cluster_metadata(table)
    clusters = sorted(meta.index)
    edges = assign_parents(sets, config.overlap_min, config.ratio_min)
    roots, _child_to_parent = family_components(clusters, edges)

    family_rows, cluster_to_family, cluster_to_family_rep = _family_membership(
        meta, clusters, roots, sets
    )
    families = _family_table(family_rows)
    family_sizes = families.set_index("family_id")[["n_clusters", "n_terms", "gene_n"]]

    assignments = _family_assignments(
        meta, edges, clusters, sets, cluster_to_family, cluster_to_family_rep, family_sizes
    )
    terms = _family_terms(table, clusters, cluster_to_family, cluster_to_family_rep, family_sizes)
    preview = _coherent_preview(table, clusters, cluster_to_family, cluster_to_family_rep, config)
    return assignments, terms, families, preview


# =============================================================================
# CORE LOGIC — family-level cohort workbook
# =============================================================================
_WORKBOOK_WIDTHS = {
    "family_id": 34,
    "source": 16,
    "group_name": 46,
    "scored_member_names": 60,
    "member_cluster_representatives": 90,
}
_DEFAULT_WIDTH = 13


def build_family_cohort_rows(
    table: pd.DataFrame,
    sets: dict[str, set[str]],
    assignments: pd.DataFrame,
    config: FamilyConfig,
) -> pd.DataFrame:
    """One workbook row per family, with the representative's evidence and its members."""
    # A family's cohort membership is driven by its MEMBERS, not only its
    # representative: a family with 14 coherent reps and one dispersed rep must not
    # disappear from the incoherent sheet. `cohort` therefore reads
    # coherent / incoherent / mixed (both ends present) / other, and the workbook
    # sheets select on the n_*_representatives counters. The row's own z/q remain
    # the representative's values; family_min_z / family_max_z show the spread.
    representatives = table.loc[table["is_representative"]].set_index(_KEY)
    rows = []
    for family_id, sub in assignments.groupby("family_id", sort=False):
        representative_cluster = sub.loc[sub["is_family_representative"]].iloc[0]
        key = (
            representative_cluster["cluster_representative_source"],
            representative_cluster["cluster_representative_group_id"],
        )
        representative = representatives.loc[key]
        genes = set().union(*(sets[cluster] for cluster in sub["redundancy_cluster"]))
        z = sub["cluster_best_median_pairwise_distance_z"].astype(float)
        q = sub["cluster_best_q_value_coherent"].astype(float)
        members = sub.sort_values("cluster_best_median_pairwise_distance_z")
        rows.append(
            {
                "family_id": family_id,
                "source": representative.name[0],
                "group_id": representative.name[1],
                "group_name": representative["group_name"],
                "n_clusters": len(sub),
                "n_terms": int(sub["cluster_n_terms"].sum()),
                "n_scored_members": len(genes),
                "representative_n_scored_members": int(representative["n_scored_members"]),
                "median_pairwise_distance_z": float(representative["median_pairwise_distance_z"]),
                "q_value_coherent": float(representative["q_value_coherent"]),
                "family_min_z": float(z.min()),
                "family_max_z": float(z.max()),
                "family_median_z": float(z.median()),
                "family_best_q_value_coherent": float(q.min()),
                "n_coherent_representatives": int(
                    ((q <= config.coherent_q_max) & (z < config.coherent_z)).sum()
                ),
                "n_incoherent_representatives": int((z > config.incoherent_z).sum()),
                "representative_paralog_fraction": representative.get("paralog_fraction"),
                "representative_moonlighting_fraction": representative.get("moonlighting_fraction"),
                "representative_abundance_cv": representative.get("abundance_cv"),
                "representative_conservation_cv": representative.get("conservation_cv"),
                "scored_member_names": ", ".join(sorted(genes)),
                "member_cluster_representatives": " | ".join(
                    f"{row['cluster_representative_name']} "
                    f"(z={row['cluster_best_median_pairwise_distance_z']:.2f})"
                    for _, row in members.iterrows()
                ),
            }
        )
    frame = pd.DataFrame(rows)
    coherent_any = frame["n_coherent_representatives"] > 0
    incoherent_any = frame["n_incoherent_representatives"] > 0
    frame["cohort"] = "other"
    frame.loc[coherent_any & ~incoherent_any, "cohort"] = "coherent"
    frame.loc[incoherent_any & ~coherent_any, "cohort"] = "incoherent"
    frame.loc[coherent_any & incoherent_any, "cohort"] = "mixed"
    return frame[
        [
            "family_id",
            "cohort",
            "source",
            "group_id",
            "group_name",
            "n_clusters",
            "n_terms",
            "n_scored_members",
            "representative_n_scored_members",
            "median_pairwise_distance_z",
            "q_value_coherent",
            "family_min_z",
            "family_max_z",
            "family_median_z",
            "family_best_q_value_coherent",
            "n_coherent_representatives",
            "n_incoherent_representatives",
            "representative_paralog_fraction",
            "representative_moonlighting_fraction",
            "representative_abundance_cv",
            "representative_conservation_cv",
            "scored_member_names",
            "member_cluster_representatives",
        ]
    ]


def build_cohort_family_rows(
    table: pd.DataFrame,
    sets: dict[str, set[str]],
    config: FamilyConfig,
) -> pd.DataFrame:
    """One row per family, with families built separately inside each cohort."""
    # The directed-containment rule is applied to the cohort's OWN representative
    # clusters, so every family is cohort-pure by construction: no mixed families,
    # and no edge can chain a coherent cluster to an incoherent one.
    meta = cluster_metadata(table)
    representatives = table.loc[table["is_representative"]].set_index(_KEY)
    cohorts = [
        (
            "coherent",
            (table["q_value_coherent"] <= config.coherent_q_max)
            & (table["median_pairwise_distance_z"] < config.coherent_z),
        ),
        ("incoherent", table["median_pairwise_distance_z"] > config.incoherent_z),
    ]
    rows = []
    for cohort_name, mask in cohorts:
        clusters = sorted(
            table.loc[table["is_representative"] & mask, "redundancy_cluster"].unique()
        )
        if not clusters:
            logger.warning(f"No representatives in the {cohort_name} cohort; skipping its families")
            continue
        subset_sets = {cluster: sets[cluster] for cluster in clusters}
        edges = assign_parents(subset_sets, config.overlap_min, config.ratio_min)
        roots, _ = family_components(clusters, edges)
        family_members: dict[str, list[str]] = defaultdict(list)
        for cluster in clusters:
            family_members[roots[cluster]].append(cluster)
        for members in family_members.values():
            member_meta = meta.loc[sorted(members)]
            representative_meta = family_representative(member_meta)
            key = (representative_meta["source"], representative_meta["group_id"])
            representative = representatives.loc[key]
            genes = set().union(*(sets[cluster] for cluster in members))
            z = member_meta["median_pairwise_distance_z"].astype(float)
            q = member_meta["q_value_coherent"].astype(float)
            sorted_members = member_meta.sort_values("median_pairwise_distance_z")
            rows.append(
                {
                    "family_id": (
                        f"{cohort_name}:family:"
                        f"{representative_meta['source']}:{representative_meta['group_id']}"
                    ),
                    "cohort": cohort_name,
                    "source": representative.name[0],
                    "group_id": representative.name[1],
                    "group_name": representative["group_name"],
                    "n_clusters": len(members),
                    "n_terms": int(member_meta["n_terms"].sum()),
                    "n_scored_members": len(genes),
                    "representative_n_scored_members": int(representative["n_scored_members"]),
                    "median_pairwise_distance_z": float(
                        representative["median_pairwise_distance_z"]
                    ),
                    "q_value_coherent": float(representative["q_value_coherent"]),
                    "family_min_z": float(z.min()),
                    "family_max_z": float(z.max()),
                    "family_median_z": float(z.median()),
                    "family_best_q_value_coherent": float(q.min()),
                    "representative_paralog_fraction": representative.get("paralog_fraction"),
                    "representative_moonlighting_fraction": representative.get(
                        "moonlighting_fraction"
                    ),
                    "representative_abundance_cv": representative.get("abundance_cv"),
                    "representative_conservation_cv": representative.get("conservation_cv"),
                    "scored_member_names": ", ".join(sorted(genes)),
                    "member_cluster_representatives": " | ".join(
                        f"{row['group_name']} (z={row['median_pairwise_distance_z']:.2f})"
                        for _, row in sorted_members.iterrows()
                    ),
                }
            )
    columns = [
        "family_id",
        "cohort",
        "source",
        "group_id",
        "group_name",
        "n_clusters",
        "n_terms",
        "n_scored_members",
        "representative_n_scored_members",
        "median_pairwise_distance_z",
        "q_value_coherent",
        "family_min_z",
        "family_max_z",
        "family_median_z",
        "family_best_q_value_coherent",
        "representative_paralog_fraction",
        "representative_moonlighting_fraction",
        "representative_abundance_cv",
        "representative_conservation_cv",
        "scored_member_names",
        "member_cluster_representatives",
    ]
    return pd.DataFrame(rows, columns=columns)


def build_cohort_workbook_sheets(
    rows: pd.DataFrame, config: FamilyConfig, n_clusters: int, n_terms_total: int
) -> dict[str, pd.DataFrame]:
    """The four workbook sheets when families are built inside each cohort."""
    coherent = rows.loc[rows["cohort"] == "coherent"].sort_values("family_min_z")
    incoherent = rows.loc[rows["cohort"] == "incoherent"].sort_values(
        "family_max_z", ascending=False
    )
    all_families = pd.concat([coherent, incoherent], ignore_index=True)
    thresholds = pd.DataFrame(
        {
            "cut": [
                "coherent_q_max",
                "coherent_z_threshold",
                "incoherent_z_threshold",
                "family_overlap_min",
                "family_ratio_min",
                "family_scope",
                "n_dedup_clusters",
                "n_terms",
                "n_coherent_families",
                "n_coherent_clusters",
                "n_coherent_terms",
                "n_incoherent_families",
                "n_incoherent_clusters",
                "n_incoherent_terms",
                "row_grain",
                "row_z_q",
                "n_scored_members",
                "family_*_z",
                "n_clusters",
            ],
            "value": [
                config.coherent_q_max,
                config.coherent_z,
                config.incoherent_z,
                config.overlap_min,
                config.ratio_min,
                "families built separately inside each cohort; no mixed families",
                n_clusters,
                n_terms_total,
                len(coherent),
                int(coherent["n_clusters"].sum()),
                int(coherent["n_terms"].sum()),
                len(incoherent),
                int(incoherent["n_clusters"].sum()),
                int(incoherent["n_terms"].sum()),
                "one row per signal family",
                "family representative's own values",
                "family gene union",
                "z spread across the family's dedup representatives",
                "dedup clusters folded into the family",
            ],
        }
    )
    return {
        "coherent": coherent,
        "incoherent": incoherent,
        "all_families": all_families,
        "thresholds": thresholds,
    }


def build_workbook_sheets(
    rows: pd.DataFrame, config: FamilyConfig, n_clusters: int
) -> dict[str, pd.DataFrame]:
    """The four workbook sheets over families, in workbook order."""
    # A family is listed in a cohort sheet when ANY of its dedup representatives
    # falls in that cohort; mixed families appear in both sheets.
    coherent = rows.loc[rows["n_coherent_representatives"] > 0].sort_values("family_min_z")
    incoherent = rows.loc[rows["n_incoherent_representatives"] > 0].sort_values(
        "family_max_z", ascending=False
    )
    mixed = (rows["n_coherent_representatives"] > 0) & (rows["n_incoherent_representatives"] > 0)
    thresholds = pd.DataFrame(
        {
            "cut": [
                "coherent_q_max",
                "coherent_z_threshold",
                "incoherent_z_threshold",
                "family_overlap_min",
                "family_ratio_min",
                "n_families",
                "n_dedup_clusters",
                "n_terms",
                "n_coherent",
                "n_incoherent",
                "n_mixed",
                "n_other",
                "cohort_rule",
                "row_grain",
                "row_z_q",
                "n_scored_members",
                "n_*_representatives",
            ],
            "value": [
                config.coherent_q_max,
                config.coherent_z,
                config.incoherent_z,
                config.overlap_min,
                config.ratio_min,
                len(rows),
                n_clusters,
                int(rows["n_terms"].sum()),
                len(coherent),
                len(incoherent),
                int(mixed.sum()),
                len(rows)
                - int(
                    (
                        (rows["n_coherent_representatives"] > 0)
                        | (rows["n_incoherent_representatives"] > 0)
                    ).sum()
                ),
                "family listed when ANY member representative falls in the cohort; "
                "mixed families appear in both sheets",
                "one row per signal family",
                "family representative's own values",
                "family gene union",
                "coherent / incoherent dedup representatives inside the family",
            ],
        }
    )
    all_families = rows.sort_values(["family_min_z", "n_clusters"], ascending=[True, False])
    return {
        "coherent": coherent,
        "incoherent": incoherent,
        "all_families": all_families,
        "thresholds": thresholds,
    }


def write_workbook(sheets: dict[str, pd.DataFrame], output: Path) -> None:
    """Write one sheet per table, with a frozen header and readable column widths."""
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            for (header,) in sheet.iter_cols(min_row=1, max_row=1):
                sheet.column_dimensions[header.column_letter].width = _WORKBOOK_WIDTHS.get(
                    str(header.value), _DEFAULT_WIDTH
                )
    logger.info(f"Wrote {output} ({output.stat().st_size / 1024:.0f} KB, {len(sheets)} sheets)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
@logger.catch(reraise=True)
def run(config: FamilyConfig) -> None:
    """Load -> build families -> write the four tables."""
    config.validate()
    table = load_deduplicated(config.deduplicated)
    sets = cluster_gene_sets(config.group_members)
    missing = set(table["redundancy_cluster"]) - set(sets)
    if missing:
        raise ValueError(
            f"{len(missing)} cluster(s) have no group-member rows: {sorted(missing)[:5]}"
        )
    assignments, terms, families, preview = build_family_tables(table, sets, config)
    assignments.to_csv(config.output_assignments, sep="\t", index=False)
    terms.to_csv(config.output_terms, sep="\t", index=False)
    families.to_csv(config.output_families, sep="\t", index=False)
    preview.to_csv(config.output_preview, sep="\t", index=False)
    if config.output_workbook is not None:
        if config.workbook_scope == "per_cohort":
            rows = build_cohort_family_rows(table, sets, config)
            sheets = build_cohort_workbook_sheets(rows, config, len(assignments), len(table))
        else:
            rows = build_family_cohort_rows(table, sets, assignments, config)
            sheets = build_workbook_sheets(rows, config, len(assignments))
        write_workbook(sheets, config.output_workbook)
        logger.info(
            f"Workbook scope={config.workbook_scope}: "
            f"coherent {len(sheets['coherent']):,} families, "
            f"incoherent {len(sheets['incoherent']):,} families"
        )
    logger.success(
        f"{len(table):,} terms -> {len(families):,} families from "
        f"{assignments.shape[0]:,} de-duplication clusters "
        f"(overlap_min={config.overlap_min:g}, ratio_min={config.ratio_min:g}); "
        f"{len(preview):,} families in the coherent preview; wrote "
        f"{config.output_families}"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build directed-containment signal families above the Jaccard de-duplication"
    )
    parser.add_argument(
        "--deduplicated", type=Path, required=True, help="dedup/coherence_terms_deduplicated.tsv"
    )
    parser.add_argument(
        "--group-members", type=Path, required=True, help="dedup/coherence_group_members_long.tsv"
    )
    parser.add_argument(
        "--overlap-min",
        type=float,
        default=DEFAULT_OVERLAP_MIN,
        help="Minimum |A∩B| / min(|A|,|B|) for a child -> parent edge",
    )
    parser.add_argument(
        "--ratio-min",
        type=float,
        default=DEFAULT_RATIO_MIN,
        help="Minimum min(|A|,|B|) / max(|A|,|B|) for a child -> parent edge",
    )
    parser.add_argument(
        "--coherent-q-max",
        type=float,
        default=DEFAULT_COHERENT_Q_MAX,
        help="Coherent cohort q_value_coherent cut for the preview",
    )
    parser.add_argument(
        "--coherent-z",
        type=float,
        default=DEFAULT_COHERENT_Z,
        help="Coherent cohort median_pairwise_distance_z cut for the preview",
    )
    parser.add_argument(
        "--incoherent-z",
        type=float,
        default=0.5,
        help="Incoherent cohort median_pairwise_distance_z cut for the workbook",
    )
    parser.add_argument(
        "--workbook-scope",
        choices=["global", "per_cohort"],
        default="global",
        help=(
            "Build workbook families globally (mixed families possible) "
            "or separately within each cohort (no mixed)"
        ),
    )
    parser.add_argument(
        "--output-assignments", type=Path, required=True, help="Per-cluster family assignment TSV"
    )
    parser.add_argument(
        "--output-terms", type=Path, required=True, help="All terms + family columns TSV"
    )
    parser.add_argument(
        "--output-families", type=Path, required=True, help="One row per family TSV"
    )
    parser.add_argument(
        "--output-preview", type=Path, required=True, help="One row per coherent family TSV"
    )
    parser.add_argument(
        "--output-workbook",
        type=Path,
        default=None,
        help="Optional family-level coherence_cohorts.xlsx",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        run(
            FamilyConfig(
                deduplicated=args.deduplicated,
                group_members=args.group_members,
                output_assignments=args.output_assignments,
                output_terms=args.output_terms,
                output_families=args.output_families,
                output_preview=args.output_preview,
                overlap_min=args.overlap_min,
                ratio_min=args.ratio_min,
                coherent_q_max=args.coherent_q_max,
                coherent_z=args.coherent_z,
                incoherent_z=args.incoherent_z,
                workbook_scope=args.workbook_scope,
                output_workbook=args.output_workbook,
            )
        )
    except (ValueError, OSError) as error:
        logger.error(f"Error: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
