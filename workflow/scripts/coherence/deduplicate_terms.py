#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Term Redundancy Reducer (display-layer de-duplication)
================================================================

GO terms (and macromolecular complexes) are heavily redundant: parents and
children share member genes, so the "most coherent" head of
combined/coherence_metrics.parquet is dominated by many aliases of the same signal
(e.g. ribosome biogenesis / rRNA processing / preribosome / 90S / nucleolus...).
This stage collapses that redundancy into clusters and picks one representative
per cluster — WITHOUT touching the statistics: the full-set q_value from
compute_coherence.py is carried through unchanged, and every term is retained
(flagged), so the reduction is a reproducible DISPLAY layer, not a re-test.

Redundancy axis (config-driven, "both" by design)
-------------------------------------------------
- MEMBER SIMILARITY (primary): jaccard_index(A, B) = |A ∩ B| / |A ∪ B| over each
  term's coherence member set (the DR<threshold `scored_member_names`, which is
  what the coherence z-score was computed on). Two terms with similarity
  >= dedup_jaccard_threshold are redundant.
  The denominator is the UNION. It used to be the smaller set (overlap
  coefficient), which scores 1.0 for ANY term contained in a larger one — and
  since GO propagation makes containment ubiquitous, that connected the whole
  graph into one blob: 2,097 terms collapsed to 2 clusters as soon as every
  source spelled genes with the same identifier. Jaccard still rates real
  nesting highly (a child covering most of its parent scores high) but no longer
  rewards containment on its own.
- DAG LINEAGE (optional, dedup_merge_dag_lineage — OFF by default): unite an
  ancestor/descendant pair (is_a + part_of, via GODag.get_all_upper) — but ONLY
  among term pairs that already share >=1 member, so disjoint sibling terms are
  never merged. It is off because that rule is degenerate on propagated GO:
  every term shares members with all of its ancestors, so transitivity chains the
  whole DAG together. Measured on the current table it puts 2,046 of 2,097 terms
  in one cluster, and raising the similarity threshold barely moves it
  (0.5 -> 0.9 changes the largest cluster from 2,046 to 1,966).
- DAG DEPTH (semantic tiebreak): a deeper GO term is more specific; used to break
  ties when selecting a cluster representative.

Clustering is transitive (union-find). Scope is `pooled` (default; clusters
across all sources, so the same complex appearing in go_cc AND go_macrocomplex
collapses) or `per_source`.

LINKAGE (dedup_linkage, default `complete`): how the "clears the threshold"
relation is turned into clusters. A similarity threshold alone only says which
PAIRS are redundant, and making a cluster out of that needs a linkage rule.

- `single` (connected components): transitive, so A~B and B~C fuse A and C even at
  Jaccard 0. A cluster ends up only as tight as its weakest link. Measured on the
  current 2,587-term table: 37.5% of all intra-cluster pairs share no member, the
  largest cluster is 193 terms whose union is 191 genes, and it takes the whole
  cytosolic ribosome and the mitochondrial one into a single cluster.
- `complete` (agglomerative, cut at 1 - threshold): a cluster's merge height IS its
  widest internal pair, so no cluster can contain a pair below the threshold and the
  chaining is gone by construction. Largest cluster drops to 26, and the threshold
  becomes a real knob (it controls how much gets collapsed rather than how badly a
  chain runs away).
- `average`: same machinery on the mean linkage, the middle ground.

`single` is the only one that can carry the DAG-lineage edge (a similarity cannot
express it), so `dedup_merge_dag_lineage` requires it — see DedupConfig.validate.

Representative selection
------------------------
Per cluster, the default representative is the best-evidence term: smallest
q_value, ties broken by more-negative median_pairwise_distance_z, then greater
dag_depth (more specific), then smaller n_scored_members. Any group_id listed in
dedup_force_representatives overrides this for its cluster (recorded as
representative_source="forced"); auto-picked ones are "auto". You always refine
by hand afterwards — the full cluster membership is emitted so nothing is hidden.

Moonlighting (a display-layer annotation on the same clusters)
--------------------------------------------------------------
A cluster's gene set is the UNION of its terms' members. Counting, per gene, how
many clusters it falls in gives a breadth profile: the MODE of that distribution
is the typical gene's group count, and a gene above the mode moonlights — it keeps
turning up in clusters that are not redundant with each other (the de-duplication
has already collapsed the aliases, so a high count is not GO nesting). Every row
of the deduplicated/representatives tables then carries `moonlighting_fraction`:
the share of that term's own members above the mode, i.e. how much of the term is
carried by broadly-shared genes (the "several pathways at once" signal).

Input
-----
- --combined: combined/coherence_metrics.parquet (source, group_id, group_name,
  n_scored_members, scored_member_names, median_pairwise_distance_z,
  median_pairwise_distance_p, q_value, ...).
- --obo: go-basic.obo (GO DAG for depth + is_a/part_of lineage).

Output
------
- --output-all: dedup/coherence_terms_deduplicated.tsv — every input row + columns
  redundancy_cluster, cluster_size, dag_depth, is_representative,
  representative_group_id, representative_name, representative_source,
  non_representative_terms (the rest of the cluster, newline-joined, on every row
  of that cluster), moonlighting_fraction. Sorted by (cluster's best z, then
  within-cluster z).
- --output-representatives: dedup/coherence_terms_representatives.tsv — only the
  is_representative rows (the de-duplicated view for figures/tables).
- --output-group-members: dedup/coherence_group_members_long.tsv — one row per
  (group, member gene): source / group_id / group_name identify the group by its
  representative term, full_name is every term name in the group newline-joined
  (representative first), then gene, n_groups (how many groups that gene is in)
  and is_moonlighting (n_groups > the mode).

Usage
-----
    python deduplicate_terms.py \\
        --combined results/3a_coherence/{dataset}/combined/coherence_metrics.parquet \\
        --obo resources/external/pombase/<version>/ontologies_and_associations/go-basic.obo \\
        --overlap-threshold 0.5 --merge-dag-lineage --scope pooled \\
        --force-representatives GO:0042254 GO:0005762 \\
        --output-all results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv \\
        --output-representatives results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv \\
        --output-group-members results/3a_coherence/{dataset}/dedup/coherence_group_members_long.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Defaults live here, not on the dataclass, because BOTH the dataclass field and
# argparse have to read the same value and a slots=True dataclass cannot supply
# it: `DedupConfig.jaccard_threshold` is a member_descriptor, not 0.5, so
# `default=DedupConfig.jaccard_threshold` silently hands argparse a descriptor.
DEFAULT_JACCARD_THRESHOLD = 0.5
# OFF, matching config/analysis.yaml and the docstring below. It was True, which
# was already stale against both — and now that `complete` is the default linkage
# the old value would make the bare defaults self-contradictory: validate() rejects
# the lineage rule outside single linkage, so `DedupConfig()` would not validate.
DEFAULT_MERGE_DAG_LINEAGE = False
DEFAULT_SCOPE = "pooled"
# How two terms that clear the threshold are grouped. `complete` is the default:
#   single   - connected components. Fast and sparse, but transitive: A~B and B~C
#              fuse A and C even at Jaccard 0, so a cluster is only as tight as its
#              weakest link. Measured here: a 193-term cluster whose union is 191
#              genes, a within-cluster Jaccard median of 0.13, and 37.5% of all
#              intra-cluster pairs sharing no member at all.
#   complete - agglomerative complete linkage, cut at 1 - threshold. A cluster's
#              merge height IS its widest internal pair, so no cluster can contain a
#              pair below the threshold; the chaining is gone by construction. That
#              also makes the threshold a real knob: raising it only moves how much
#              gets collapsed (1,630 -> 994 non-representatives from 0.5 to 0.9)
#              instead of racing a chain (whose largest cluster runs 193 -> 31 over
#              the same range). Biggest cluster: 26 at 0.5.
#   average  - same, on the mean linkage. The middle ground, offered because the
#              call is identical.
# `single` alone can carry the DAG-lineage edge rule (a similarity cannot express
# it), so the two are mutually exclusive - see DedupConfig.validate.
DEFAULT_LINKAGE = "complete"
_LINKAGE_METHODS = ("single", "average", "complete")

# Columns of the (group, gene) long table, in write order. `redundancy_cluster` is
# what the per-view breadth counts slice on; `n_groups` / `is_moonlighting` are per
# GENE and filled after the rows are collected.
_GROUP_MEMBER_COLUMNS = ["source", "group_id", "group_name", "full_name", "gene",
                         "redundancy_cluster", "n_groups", "is_moonlighting"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class DedupConfig:
    """Inputs, outputs, and parameters for coherence term de-duplication."""
    combined: Path
    obo: Path
    output_all: Path
    output_representatives: Path
    output_group_members: Path
    jaccard_threshold: float = DEFAULT_JACCARD_THRESHOLD
    merge_dag_lineage: bool = DEFAULT_MERGE_DAG_LINEAGE
    linkage: str = DEFAULT_LINKAGE
    scope: str = DEFAULT_SCOPE  # "pooled" | "per_source"
    force_representatives: list[str] = field(default_factory=list)

    def validate(self) -> None:
        """Raise ValueError on bad inputs/params, then make output dirs."""
        for path in [self.combined, self.obo]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if not 0.0 < self.jaccard_threshold <= 1.0:
            raise ValueError(f"jaccard_threshold must be in (0, 1]: {self.jaccard_threshold}")
        if self.scope not in ("pooled", "per_source"):
            raise ValueError(f"scope must be 'pooled' or 'per_source': {self.scope!r}")
        if self.linkage not in _LINKAGE_METHODS:
            raise ValueError(f"linkage must be one of {_LINKAGE_METHODS}: {self.linkage!r}")
        if self.linkage != "single" and self.merge_dag_lineage:
            # The lineage rule adds an edge for a member-sharing ancestor/descendant
            # pair whatever their Jaccard, which a similarity-based linkage cannot
            # express. Silently dropping it would answer a different question than
            # the one the config asked, so say so instead.
            raise ValueError(
                f"merge_dag_lineage needs linkage='single' (got {self.linkage!r}): "
                "the DAG-ancestor edge is not a similarity and cannot be cut by a linkage"
            )
        for out in [self.output_all, self.output_representatives, self.output_group_members]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from io_table import read_parquet  # noqa: E402
from coherence.fractions import moonlighting_fraction, view_breadth  # noqa: E402
# =============================================================================
# CORE LOGIC — redundancy graph (member overlap + optional DAG lineage)
# =============================================================================
def candidate_pairs(member_sets: list[set[str]]) -> set[tuple[int, int]]:
    """All (i, j) term-index pairs that share >=1 member, via a gene->terms index."""
    # Building the inverted index and only emitting co-occurring pairs is cheaper
    # than the full O(n^2) sweep, but not asymptotically so: GO terms overlap
    # heavily once propagated, so pair count is sum_g C(k_g, 2) over the per-gene
    # term degree k_g. Measured on the current 1,706-group table: 555,784 pairs
    # against 1,454,365 for a full sweep — a 2.6x saving. Returned pairs are
    # ordered i < j and de-duplicated.
    gene_to_terms: dict[str, list[int]] = defaultdict(list)
    for idx, members in enumerate(member_sets):
        for gene in members:
            gene_to_terms[gene].append(idx)
    pairs: set[tuple[int, int]] = set()
    for terms in gene_to_terms.values():
        if len(terms) < 2:
            continue
        for a in range(len(terms)):
            for b in range(a + 1, len(terms)):
                i, j = terms[a], terms[b]
                pairs.add((i, j) if i < j else (j, i))
    return pairs


def jaccard_index(a: set[str], b: set[str]) -> float:
    """|A ∩ B| / |A ∪ B|; 0.0 if either set is empty."""
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def member_set(value: object) -> set[str]:
    """A `scored_member_names` cell as a set of gene names."""
    # `scored_member_names` is stored as a real list column in the metrics Parquet, so a
    # cell is a list; a null cell (no covered genes) arrives as None/NaN and reads
    # as the empty set. A comma-joined string is rejected rather than iterated:
    # iterating it would silently yield single characters, which would corrupt the
    # overlap graph without any error.
    if value is None or isinstance(value, float):
        return set()
    if isinstance(value, str):
        raise TypeError(
            f"scored_member_names is a list column in the metrics Parquet; got a string ({value[:40]!r}...). "
            "Split it into a list at the source instead of parsing it here."
        )
    return {str(gene) for gene in value}


def condensed_index(i: int, j: int, n: int) -> int:
    """Position of the (i, j) pair in scipy's condensed distance vector (i < j)."""
    return n * i - i * (i + 1) // 2 + (j - i - 1)


def build_clusters(
    sub: pd.DataFrame,
    threshold: float,
    merge_dag_lineage: bool,
    ancestors: dict[str, set[str]],
    linkage_method: str = "single",
) -> list[int]:
    """Cluster the rows of `sub` (a single scope) -> a cluster label per row."""
    # Two terms are united when their member-set Jaccard similarity >= threshold,
    # OR (single linkage + merge_dag_lineage) one is a DAG ancestor of the other AND
    # they share >=1 member — disjoint siblings are never merged.
    # `ancestors[group_id]` is the is_a+part_of ancestor set.
    #
    # `single` returns connected components; `average`/`complete` return a
    # hierarchical cut. The returned labels are small dense integers aligned to
    # sub's row order, numbered by first appearance.
    member_sets = [member_set(cg) for cg in sub["scored_member_names"]]
    n_rows = len(sub)

    if linkage_method != "single":
        labels = hierarchical_labels(member_sets, threshold, linkage_method)
        return _relabel_by_first_appearance(labels)

    group_ids = sub["group_id"].tolist()
    edges = []
    for i, j in candidate_pairs(member_sets):
        united = jaccard_index(member_sets[i], member_sets[j]) >= threshold
        if not united and merge_dag_lineage:
            # Candidate pairs already share >=1 member, so an ancestor link here
            # is a member-sharing parent/child (safe to merge), not a disjoint pair.
            gi, gj = group_ids[i], group_ids[j]
            united = gj in ancestors.get(gi, set()) or gi in ancestors.get(gj, set())
        if united:
            edges.append((i, j))

    if not edges:
        return list(range(n_rows))

    rows, cols = zip(*edges)
    graph = csr_matrix(([1] * len(edges), (rows, cols)), shape=(n_rows, n_rows))
    _, labels = connected_components(graph, directed=False)
    return _relabel_by_first_appearance(labels)


def hierarchical_labels(
    member_sets: list[set[str]], threshold: float, linkage_method: str
) -> np.ndarray:
    """Cut an agglomerative tree at `threshold` -> one label per term.

    Unlike connected components, a cut here bounds the whole cluster and not just
    its chain: with complete linkage a cluster's merge height IS its widest internal
    pair, so no cluster can end up holding two terms below the threshold. That is
    what removes the chaining — measured on the current table, single linkage at 0.5
    produced a 193-term cluster (union: 191 genes) whose members are mostly mutually
    disjoint, which complete linkage caps at 26.

    The distance matrix is dense: a pair that shares no member has Jaccard 0 and so
    sits at the maximum distance, and the hierarchy needs every pair to build the
    tree at all. Only the sharing pairs are computed; the rest stay at 1.0.
    n(n-1)/2 doubles per doubling of the term count — 26 MB at today's 2,587 terms,
    105 MB at 5,000 — so this is not the branch for a table an order of magnitude
    larger without a sparse-linkage implementation.
    """
    n_rows = len(member_sets)
    if n_rows < 2:
        return np.arange(n_rows)

    distances = np.ones(n_rows * (n_rows - 1) // 2, dtype=np.float64)
    for i, j in candidate_pairs(member_sets):
        distances[condensed_index(i, j, n_rows)] = 1.0 - jaccard_index(member_sets[i], member_sets[j])

    tree = linkage(distances, method=linkage_method)
    # `criterion="distance"` cuts at the height itself, and the height of a complete
    # linkage is 1 - Jaccard, so t = 1 - threshold is exactly "no intra-cluster pair
    # below the threshold".
    return fcluster(tree, t=1.0 - threshold, criterion="distance")


def _relabel_by_first_appearance(labels) -> list[int]:
    """Dense 0..k-1 labels in first-appearance order.

    The caller turns labels into cluster-id strings that the output sort keys on,
    so the numbering is observable and has to be stable.
    """
    remap: dict[int, int] = {}
    return [remap.setdefault(int(label), len(remap)) for label in labels]


def pick_representative(cluster_rows: pd.DataFrame, forced: set[str]) -> tuple[int, str]:
    """Return (index_label_of_representative, source) for one cluster."""
    # A forced group_id present in the cluster wins (source="forced"); ties among
    # multiple forced members fall back to the same ordering. Otherwise the default
    # is best evidence: min q_value, then most-negative median_pairwise_distance_z,
    # then greatest dag_depth (more specific), then smallest n_scored_members.
    # The DataFrame index label
    # is returned so the caller can flag that row.
    forced_here = cluster_rows[cluster_rows["group_id"].isin(forced)]
    pool = forced_here if not forced_here.empty else cluster_rows
    source = "forced" if not forced_here.empty else "auto"
    ordered = pool.sort_values(
        by=["q_value", "median_pairwise_distance_z", "dag_depth", "n_scored_members"],
        ascending=[True, True, False, True],
    )
    return ordered.index[0], source


# =============================================================================
# CORE LOGIC — DAG depth + ancestors
# =============================================================================
def load_dag_depth_ancestors(obo: Path, group_ids: list[str]) -> tuple[dict[str, int], dict[str, set[str]]]:
    """Load the GO DAG once -> {gid: depth} and {gid: is_a+part_of ancestor set}."""
    # The `relationship` optional attr is loaded so get_all_upper() returns the
    # is_a + part_of ancestors, matching how the coherence members were propagated.
    # Terms absent from the DAG get depth 0 and no ancestors (they then only
    # cluster by member overlap).
    from goatools.obo_parser import GODag  # imported here: only this rule's env has goatools

    dag = GODag(str(obo), optional_attrs={"relationship"}, prt=None)
    depth: dict[str, int] = {}
    ancestors: dict[str, set[str]] = {}
    n_missing = 0
    for gid in set(group_ids):
        rec = dag.get(gid)
        if rec is None:
            depth[gid] = 0
            ancestors[gid] = set()
            n_missing += 1
            continue
        depth[gid] = rec.depth
        # get_all_upper() walks the `relationship` attr loaded above, i.e. is_a +
        # part_of. get_all_parents() defaults to is_a only, so a fallback to it
        # would silently drop part_of ancestors; goatools has had get_all_upper
        # since well before the pinned version, so there is nothing to fall back to.
        ancestors[gid] = rec.get_all_upper()
    if n_missing:
        logger.warning(f"{n_missing} group_id(s) not found in the GO DAG; depth=0, no lineage for those")
    return depth, ancestors


# =============================================================================
# CORE LOGIC — orchestration
# =============================================================================
def deduplicate(table: pd.DataFrame, config: DedupConfig,
                depth: dict[str, int], ancestors: dict[str, set[str]]) -> pd.DataFrame:
    """Annotate `table` with cluster + representative columns (no rows dropped)."""
    table = table.copy()
    table["dag_depth"] = table["group_id"].map(depth).fillna(0).astype(int)

    # Assign globally-unique cluster ids. For per_source scope, cluster within
    # each source and prefix the label with the source so ids stay disjoint.
    scopes = [("all", table)] if config.scope == "pooled" else list(table.groupby("source"))
    table["redundancy_cluster"] = pd.NA
    for scope_name, sub in scopes:
        labels = build_clusters(
            sub, config.jaccard_threshold, config.merge_dag_lineage, ancestors, config.linkage
        )
        cluster_ids = [f"{scope_name}:{lab}" for lab in labels]
        table.loc[sub.index, "redundancy_cluster"] = cluster_ids

    forced = set(config.force_representatives)
    table["cluster_size"] = table.groupby("redundancy_cluster")["group_id"].transform("size")
    table["is_representative"] = False
    table["representative_group_id"] = pd.NA
    table["representative_name"] = pd.NA
    table["representative_source"] = pd.NA
    table["non_representative_terms"] = ""
    for _cluster, rows in table.groupby("redundancy_cluster"):
        rep_idx, rep_source = pick_representative(rows, forced)
        table.loc[rep_idx, "is_representative"] = True
        table.loc[rows.index, "representative_group_id"] = table.loc[rep_idx, "group_id"]
        table.loc[rows.index, "representative_name"] = table.loc[rep_idx, "group_name"]
        table.loc[rows.index, "representative_source"] = rep_source
        # What this cluster folded into its representative, one term per line as
        # "source:group_name (group_id)" (the source is not redundant here — pooled
        # scope puts go_cc and go_macrocomplex copies of one term in one cluster).
        # Cluster-level like the representative_* columns: written on EVERY row of
        # the cluster, so the representatives-only table carries the collapsed list
        # next to the row it describes.
        others = rows[rows.index != rep_idx].sort_values("median_pairwise_distance_z")
        table.loc[rows.index, "non_representative_terms"] = "\n".join(
            f"{source}:{name} ({gid})"
            for source, name, gid in zip(others["source"], others["group_name"], others["group_id"])
        )

    # Sort so each cluster's best (min) z leads, members grouped, best-first within.
    table["_cluster_best_z"] = table.groupby("redundancy_cluster")["median_pairwise_distance_z"].transform("min")
    table = table.sort_values(
        by=["_cluster_best_z", "redundancy_cluster", "median_pairwise_distance_z"],
    ).drop(columns="_cluster_best_z").reset_index(drop=True)
    return table


def group_member_table(annotated: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """One row per (group, member gene) + the per-gene moonlighting counts and the mode.

    A group is one de-duplication cluster, identified by its representative term
    (source + group_id + group_name) and labelled with every term it collapsed.
    Its genes are the UNION over those terms, not the representative's own set: the
    collapsed terms are one signal, so a gene only an alias carries still belongs to
    the group, and it counts once for the group rather than once per alias.

    The breadth columns come from coherence.fractions.view_breadth over this whole
    table (every cluster), which is the same computation the per-view figures run
    over a subset of the clusters — one definition, so the table and the figures
    cannot drift. The mode is returned because the per-term fractions need the same
    cut.
    """
    rows = []
    for _cluster, sub in annotated.groupby("redundancy_cluster", sort=False):
        representative = sub[sub["is_representative"]].iloc[0]
        genes: set[str] = set()
        for cell in sub["scored_member_names"]:
            genes |= member_set(cell)
        # The representative first, then everything it collapsed: the name reads as
        # "what this group is" before the aliases it stands in for.
        names = [str(representative["group_name"])]
        names += [str(name) for name in sub.loc[~sub["is_representative"], "group_name"]]
        rows.extend(
            {
                "source": representative["source"],
                "group_id": representative["group_id"],
                "group_name": representative["group_name"],
                "full_name": "\n".join(names),
                "gene": gene,
                "redundancy_cluster": sub["redundancy_cluster"].iloc[0],
            }
            for gene in sorted(genes)
        )

    table = pd.DataFrame(rows, columns=_GROUP_MEMBER_COLUMNS)
    if table.empty:
        return table, 0
    per_gene, mode = view_breadth(table, set(table["redundancy_cluster"]))
    breadth = per_gene.set_index("gene")
    table["n_groups"] = table["gene"].map(breadth["n_groups"]).astype("int64")
    table["is_moonlighting"] = table["gene"].map(breadth["is_moonlighting"]).astype(bool)
    return table, mode


def write_dedup_table(table: pd.DataFrame, path: Path) -> None:
    """Write a human-facing dedup TSV, flattening `scored_member_names` back to a string."""
    # `scored_member_names` is a real list column in the metrics Parquet (so the list
    # survives the round-trip without a join/parse step). These two TSVs are read
    # by people, not by code, so the list is joined on the way out.
    out = table.copy()
    if "scored_member_names" in out.columns:
        out["scored_member_names"] = out["scored_member_names"].map(
            lambda genes: "" if genes is None or isinstance(genes, float) else ", ".join(map(str, genes))
        )
    out.to_csv(path, sep="\t", index=False)


@logger.catch(reraise=True)
def run(config: DedupConfig) -> None:
    """Load -> annotate clusters + representatives -> write full + representatives tables."""
    config.validate()
    table = read_parquet(config.combined)
    for required in ["source", "group_id", "group_name", "n_scored_members",
                     "scored_member_names", "median_pairwise_distance_z", "q_value"]:
        if required not in table.columns:
            raise ValueError(f"combined metrics missing required column '{required}' (have: {list(table.columns)})")

    if table.empty:
        logger.warning("combined metrics table is empty; writing empty dedup outputs")
        write_dedup_table(table, config.output_all)
        write_dedup_table(table, config.output_representatives)
        group_member_table(table)[0].to_csv(config.output_group_members, sep="\t", index=False)
        return

    depth, ancestors = load_dag_depth_ancestors(config.obo, table["group_id"].tolist())
    annotated = deduplicate(table, config, depth, ancestors)
    group_members, mode = group_member_table(annotated)
    moonlighting_genes = set(group_members.loc[group_members["is_moonlighting"], "gene"])
    annotated["moonlighting_fraction"] = annotated["scored_member_names"].map(
        lambda cell: moonlighting_fraction(member_set(cell), moonlighting_genes)
    )
    write_dedup_table(annotated, config.output_all)
    group_members.to_csv(config.output_group_members, sep="\t", index=False)

    reps = annotated[annotated["is_representative"]].reset_index(drop=True)
    write_dedup_table(reps, config.output_representatives)

    n_clusters = annotated["redundancy_cluster"].nunique()
    n_forced = int((annotated["is_representative"] & (annotated["representative_source"] == "forced")).sum())
    per_gene = group_members.drop_duplicates("gene")
    logger.success(
        f"{len(annotated):,} terms -> {n_clusters:,} non-redundant clusters "
        f"({len(annotated) - n_clusters:,} collapsed; {n_forced} forced representatives); "
        f"{int(per_gene['is_moonlighting'].sum()):,}/{len(per_gene):,} genes above the mode "
        f"({mode}) in {len(group_members):,} (group, gene) rows; wrote {config.output_representatives}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    # Defaults come from the module constants above, shared with the dataclass.
    # Do NOT read them off `DedupConfig`: the dataclass is slots=True, so its class
    # attributes are member_descriptors rather than default values.
    parser = argparse.ArgumentParser(description="De-duplicate coherence terms by member overlap + GO DAG structure")
    parser.add_argument("--combined", type=Path, required=True, help="combined/coherence_metrics.parquet")
    parser.add_argument("--obo", type=Path, required=True, help="go-basic.obo (GO DAG for depth + lineage)")
    parser.add_argument("--jaccard-threshold", type=float, default=DEFAULT_JACCARD_THRESHOLD,
                        help="Member-set Jaccard similarity at or above which two terms are redundant")
    parser.add_argument("--merge-dag-lineage", action=argparse.BooleanOptionalAction,
                        default=DEFAULT_MERGE_DAG_LINEAGE,
                        help="Also merge member-sharing ancestor/descendant pairs (--no-merge-dag-lineage to disable; single linkage only)")
    parser.add_argument("--linkage", choices=list(_LINKAGE_METHODS), default=DEFAULT_LINKAGE,
                        help="How to group terms that clear the threshold: single (connected components, chains), average or complete (hierarchical cut)")
    parser.add_argument("--scope", choices=["pooled", "per_source"], default=DEFAULT_SCOPE,
                        help="Cluster across all sources (pooled) or within each source")
    parser.add_argument("--force-representatives", nargs="*", default=[], help="group_ids forced to be their cluster's representative")
    parser.add_argument("--output-all", type=Path, required=True, help="Output annotated (all terms) TSV")
    parser.add_argument("--output-representatives", type=Path, required=True, help="Output representatives-only TSV")
    parser.add_argument("--output-group-members", type=Path, required=True,
                        help="Output (group, gene) long table TSV with the moonlighting columns")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run de-duplication, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = DedupConfig(
            combined=args.combined,
            obo=args.obo,
            output_all=args.output_all,
            output_representatives=args.output_representatives,
            output_group_members=args.output_group_members,
            jaccard_threshold=args.jaccard_threshold,
            merge_dag_lineage=args.merge_dag_lineage,
            linkage=args.linkage,
            scope=args.scope,
            force_representatives=list(args.force_representatives),
        )
        run(config)
    except (ValueError, OSError) as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
