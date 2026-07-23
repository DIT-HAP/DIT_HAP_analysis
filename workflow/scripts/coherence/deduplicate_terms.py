#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Term Redundancy Reducer (display-layer de-duplication)
================================================================

GO terms (and macromolecular complexes) are heavily redundant: parents and
children share member genes, so the "most coherent" head of
coherence_metrics_combined.tsv is dominated by many aliases of the same signal
(e.g. ribosome biogenesis / rRNA processing / preribosome / 90S / nucleolus...).
This stage collapses that redundancy into clusters and picks one representative
per cluster — WITHOUT touching the statistics: the full-set p_fdr from
compute_coherence.py is carried through unchanged, and every term is retained
(flagged), so the reduction is a reproducible DISPLAY layer, not a re-test.

Redundancy axis (config-driven, "both" by design)
-------------------------------------------------
- MEMBER OVERLAP (primary): overlap_coefficient(A, B) = |A ∩ B| / min(|A|, |B|)
  over each term's coherence member set (the DR>threshold `covered_genes`, which
  is exactly what the coherence z-score was computed on). Two terms with
  overlap >= dedup_overlap_threshold are redundant. This provably catches all
  GO parent/child nesting: propagation makes a child's genes a subset of its
  parent's, so after the same DR filter their overlap coefficient is 1.0.
- DAG LINEAGE (optional safety net, dedup_merge_dag_lineage): also unite an
  ancestor/descendant pair (is_a + part_of, via GODag.get_all_upper) — but ONLY
  among term pairs that already share >=1 member, so disjoint sibling terms
  (different coherence signals) are never merged.
- DAG DEPTH (semantic tiebreak): a deeper GO term is more specific; used to break
  ties when selecting a cluster representative.

Clustering is transitive (union-find). Scope is `pooled` (default; clusters
across all sources, so the same complex appearing in go_cc AND go_macrocomplex
collapses) or `per_source`.

Representative selection
------------------------
Per cluster, the default representative is the best-evidence term: smallest
p_fdr, ties broken by more-negative z_score, then greater dag_depth (more
specific), then smaller term_size. Any group_id listed in
dedup_force_representatives overrides this for its cluster (recorded as
representative_source="forced"); auto-picked ones are "auto". You always refine
by hand afterwards — the full cluster membership is emitted so nothing is hidden.

Input
-----
- --combined: coherence_metrics_combined.tsv (source, group_id, group_name,
  term_size, covered_genes, z_score, p_value, p_fdr, ...).
- --obo: go-basic.obo (GO DAG for depth + is_a/part_of lineage).

Output
------
- --output-all: coherence_terms_deduplicated.tsv — every input row + columns
  redundancy_cluster, cluster_size, dag_depth, is_representative,
  representative_group_id, representative_name, representative_source. Sorted by
  (cluster's best z, then within-cluster z).
- --output-representatives: coherence_terms_representatives.tsv — only the
  is_representative rows (the de-duplicated view for figures/tables).

Usage
-----
    python deduplicate_terms.py \\
        --combined results/coherence/{dataset}/coherence_metrics_combined.tsv \\
        --obo resources/external/pombase/<version>/ontologies_and_associations/go-basic.obo \\
        --overlap-threshold 0.5 --merge-dag-lineage --scope pooled \\
        --force-representatives GO:0042254 GO:0005762 \\
        --output-all results/coherence/{dataset}/coherence_terms_deduplicated.tsv \\
        --output-representatives results/coherence/{dataset}/coherence_terms_representatives.tsv

Author:   Yusheng Yang (guidance) + Claude Opus 4.8 (implementation)
Date:     2026-07-23
Version:  1.0.0
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
import pandas as pd

# 3. Third-party Imports
from loguru import logger


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class DedupConfig:
    """Inputs, outputs, and parameters for coherence term de-duplication."""
    combined: Path
    obo: Path
    output_all: Path
    output_representatives: Path
    overlap_threshold: float = 0.5
    merge_dag_lineage: bool = True
    scope: str = "pooled"  # "pooled" | "per_source"
    force_representatives: list[str] = field(default_factory=list)

    def validate(self) -> None:
        """Raise ValueError on bad inputs/params, then make output dirs."""
        for path in [self.combined, self.obo]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        if not 0.0 < self.overlap_threshold <= 1.0:
            raise ValueError(f"overlap_threshold must be in (0, 1]: {self.overlap_threshold}")
        if self.scope not in ("pooled", "per_source"):
            raise ValueError(f"scope must be 'pooled' or 'per_source': {self.scope!r}")
        for out in [self.output_all, self.output_representatives]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC — redundancy graph (member overlap + optional DAG lineage)
# =============================================================================
class UnionFind:
    """Minimal disjoint-set (path compression + union by size); no deps.

    Nodes are integer row indices into the metrics table. `union` merges two
    nodes' sets; `find` returns a set's canonical root; `groups` returns the
    final {root: [members]} partition.
    """
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # path halving
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]

    def groups(self) -> dict[int, list[int]]:
        out: dict[int, list[int]] = defaultdict(list)
        for i in range(len(self.parent)):
            out[self.find(i)].append(i)
        return out


def candidate_pairs(member_sets: list[set[str]]) -> set[tuple[int, int]]:
    """All (i, j) term-index pairs that share >=1 member, via a gene->terms index.

    Building the inverted index and only emitting co-occurring pairs avoids the
    full O(n^2) sweep (most term pairs are disjoint). Returned pairs are ordered
    i < j and de-duplicated.
    """
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


def overlap_coefficient(a: set[str], b: set[str]) -> float:
    """|A ∩ B| / min(|A|, |B|); 0.0 if either set is empty."""
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def build_clusters(
    sub: pd.DataFrame,
    threshold: float,
    merge_dag_lineage: bool,
    ancestors: dict[str, set[str]],
) -> list[int]:
    """Cluster the rows of `sub` (a single scope) -> a cluster label per row.

    Two terms are united when their member overlap coefficient >= threshold, OR
    (when merge_dag_lineage) one is a DAG ancestor of the other AND they share
    >=1 member (never merges disjoint siblings). `ancestors[group_id]` is the
    is_a+part_of ancestor set. Returns a list of small dense integer labels
    aligned to sub's row order.
    """
    member_sets = [set(str(cg).split(", ")) if pd.notna(cg) else set()
                   for cg in sub["covered_genes"]]
    group_ids = sub["group_id"].tolist()
    uf = UnionFind(len(sub))
    for i, j in candidate_pairs(member_sets):
        if overlap_coefficient(member_sets[i], member_sets[j]) >= threshold:
            uf.union(i, j)
        elif merge_dag_lineage:
            gi, gj = group_ids[i], group_ids[j]
            # Candidate pairs already share >=1 member, so an ancestor link here
            # is a member-sharing parent/child (safe to merge), not a disjoint pair.
            if gj in ancestors.get(gi, set()) or gi in ancestors.get(gj, set()):
                uf.union(i, j)
    # Densify roots -> 0..k-1 labels.
    root_to_label: dict[int, int] = {}
    labels = []
    for i in range(len(sub)):
        root = uf.find(i)
        labels.append(root_to_label.setdefault(root, len(root_to_label)))
    return labels


def pick_representative(cluster_rows: pd.DataFrame, forced: set[str]) -> tuple[int, str]:
    """Return (index_label_of_representative, source) for one cluster.

    A forced group_id present in the cluster wins (source="forced"); ties among
    multiple forced members fall back to the same ordering. Otherwise the default
    is best evidence: min p_fdr, then most-negative z_score, then greatest
    dag_depth (more specific), then smallest term_size. Returns the DataFrame
    index label so the caller can flag that row.
    """
    forced_here = cluster_rows[cluster_rows["group_id"].isin(forced)]
    pool = forced_here if not forced_here.empty else cluster_rows
    source = "forced" if not forced_here.empty else "auto"
    ordered = pool.sort_values(
        by=["p_fdr", "z_score", "dag_depth", "term_size"],
        ascending=[True, True, False, True],
    )
    return ordered.index[0], source


# =============================================================================
# CORE LOGIC — DAG depth + ancestors
# =============================================================================
def load_dag_depth_ancestors(obo: Path, group_ids: list[str]) -> tuple[dict[str, int], dict[str, set[str]]]:
    """Load the GO DAG once -> {gid: depth} and {gid: is_a+part_of ancestor set}.

    Uses goatools GODag with the `relationship` optional attr so get_all_upper()
    returns the is_a + part_of ancestors (matching how the coherence members were
    propagated). Terms absent from the DAG get depth 0 and no ancestors (they
    then only cluster by member overlap).
    """
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
        ancestors[gid] = rec.get_all_upper() if hasattr(rec, "get_all_upper") else rec.get_all_parents()
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
        labels = build_clusters(sub, config.overlap_threshold, config.merge_dag_lineage, ancestors)
        cluster_ids = [f"{scope_name}:{lab}" for lab in labels]
        table.loc[sub.index, "redundancy_cluster"] = cluster_ids

    forced = set(config.force_representatives)
    table["cluster_size"] = table.groupby("redundancy_cluster")["group_id"].transform("size")
    table["is_representative"] = False
    table["representative_group_id"] = pd.NA
    table["representative_name"] = pd.NA
    table["representative_source"] = pd.NA
    for _cluster, rows in table.groupby("redundancy_cluster"):
        rep_idx, rep_source = pick_representative(rows, forced)
        table.loc[rep_idx, "is_representative"] = True
        table.loc[rows.index, "representative_group_id"] = table.loc[rep_idx, "group_id"]
        table.loc[rows.index, "representative_name"] = table.loc[rep_idx, "group_name"]
        table.loc[rows.index, "representative_source"] = rep_source

    # Sort so each cluster's best (min) z leads, members grouped, best-first within.
    table["_cluster_best_z"] = table.groupby("redundancy_cluster")["z_score"].transform("min")
    table = table.sort_values(
        by=["_cluster_best_z", "redundancy_cluster", "z_score"],
    ).drop(columns="_cluster_best_z").reset_index(drop=True)
    return table


@logger.catch(reraise=True)
def run(config: DedupConfig) -> None:
    """Load -> annotate clusters + representatives -> write full + representatives tables."""
    config.validate()
    table = pd.read_csv(config.combined, sep="\t")
    for required in ["source", "group_id", "group_name", "term_size", "covered_genes",
                     "z_score", "p_fdr"]:
        if required not in table.columns:
            raise ValueError(f"combined metrics missing required column '{required}' (have: {list(table.columns)})")

    if table.empty:
        logger.warning("combined metrics table is empty; writing empty dedup outputs")
        table.to_csv(config.output_all, sep="\t", index=False)
        table.to_csv(config.output_representatives, sep="\t", index=False)
        return

    depth, ancestors = load_dag_depth_ancestors(config.obo, table["group_id"].tolist())
    annotated = deduplicate(table, config, depth, ancestors)
    annotated.to_csv(config.output_all, sep="\t", index=False)

    reps = annotated[annotated["is_representative"]].reset_index(drop=True)
    reps.to_csv(config.output_representatives, sep="\t", index=False)

    n_clusters = annotated["redundancy_cluster"].nunique()
    n_forced = int((annotated["is_representative"] & (annotated["representative_source"] == "forced")).sum())
    logger.success(
        f"{len(annotated):,} terms -> {n_clusters:,} non-redundant clusters "
        f"({len(annotated) - n_clusters:,} collapsed; {n_forced} forced representatives); "
        f"wrote {config.output_representatives}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="De-duplicate coherence terms by member overlap + GO DAG structure")
    parser.add_argument("--combined", type=Path, required=True, help="coherence_metrics_combined.tsv")
    parser.add_argument("--obo", type=Path, required=True, help="go-basic.obo (GO DAG for depth + lineage)")
    parser.add_argument("--overlap-threshold", type=float, default=0.5, help="Overlap-coefficient cutoff to call two terms redundant")
    parser.add_argument("--merge-dag-lineage", action="store_true", help="Also merge member-sharing ancestor/descendant pairs")
    parser.add_argument("--scope", choices=["pooled", "per_source"], default="pooled", help="Cluster across all sources (pooled) or within each source")
    parser.add_argument("--force-representatives", nargs="*", default=[], help="group_ids forced to be their cluster's representative")
    parser.add_argument("--output-all", type=Path, required=True, help="Output annotated (all terms) TSV")
    parser.add_argument("--output-representatives", type=Path, required=True, help="Output representatives-only TSV")
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
            overlap_threshold=args.overlap_threshold,
            merge_dag_lineage=args.merge_dag_lineage,
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
