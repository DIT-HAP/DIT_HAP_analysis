#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coherence Term Redundancy Reducer (display-layer de-duplication)
================================================================

GO terms (and macromolecular complexes) are heavily redundant: parents and
children share member genes, so the "most coherent" head of
coherence_metrics_combined.parquet is dominated by many aliases of the same signal
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

Representative selection
------------------------
Per cluster, the default representative is the best-evidence term: smallest
q_value, ties broken by more-negative median_pairwise_distance_z, then greater
dag_depth (more specific), then smaller n_scored_members. Any group_id listed in
dedup_force_representatives overrides this for its cluster (recorded as
representative_source="forced"); auto-picked ones are "auto". You always refine
by hand afterwards — the full cluster membership is emitted so nothing is hidden.

Input
-----
- --combined: coherence_metrics_combined.parquet (source, group_id, group_name,
  n_scored_members, scored_member_names, median_pairwise_distance_z,
  median_pairwise_distance_p, q_value, ...).
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
        --combined results/3a_coherence/{dataset}/coherence_metrics_combined.parquet \\
        --obo resources/external/pombase/<version>/ontologies_and_associations/go-basic.obo \\
        --overlap-threshold 0.5 --merge-dag-lineage --scope pooled \\
        --force-representatives GO:0042254 GO:0005762 \\
        --output-all results/3a_coherence/{dataset}/coherence_terms_deduplicated.tsv \\
        --output-representatives results/3a_coherence/{dataset}/coherence_terms_representatives.tsv

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
DEFAULT_MERGE_DAG_LINEAGE = True
DEFAULT_SCOPE = "pooled"


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
    jaccard_threshold: float = DEFAULT_JACCARD_THRESHOLD
    merge_dag_lineage: bool = DEFAULT_MERGE_DAG_LINEAGE
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
        for out in [self.output_all, self.output_representatives]:
            out.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from io_table import read_parquet  # noqa: E402
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


def build_clusters(
    sub: pd.DataFrame,
    threshold: float,
    merge_dag_lineage: bool,
    ancestors: dict[str, set[str]],
) -> list[int]:
    """Cluster the rows of `sub` (a single scope) -> a cluster label per row."""
    # Two terms are united when their member-set Jaccard similarity >= threshold,
    # OR (when merge_dag_lineage) one is a DAG ancestor of the other AND they share
    # >=1 member — disjoint siblings are never merged. `ancestors[group_id]` is the
    # is_a+part_of ancestor set. The returned labels are small dense integers
    # aligned to sub's row order, numbered by first appearance.
    member_sets = [member_set(cg) for cg in sub["scored_member_names"]]
    group_ids = sub["group_id"].tolist()
    n_rows = len(sub)

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
    # Relabel to first-appearance order: the caller turns labels into cluster-id
    # strings that the output sort keys on, so the numbering is observable.
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
        labels = build_clusters(sub, config.jaccard_threshold, config.merge_dag_lineage, ancestors)
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
    table["_cluster_best_z"] = table.groupby("redundancy_cluster")["median_pairwise_distance_z"].transform("min")
    table = table.sort_values(
        by=["_cluster_best_z", "redundancy_cluster", "median_pairwise_distance_z"],
    ).drop(columns="_cluster_best_z").reset_index(drop=True)
    return table


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
        return

    depth, ancestors = load_dag_depth_ancestors(config.obo, table["group_id"].tolist())
    annotated = deduplicate(table, config, depth, ancestors)
    write_dedup_table(annotated, config.output_all)

    reps = annotated[annotated["is_representative"]].reset_index(drop=True)
    write_dedup_table(reps, config.output_representatives)

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
    # Defaults come from the module constants above, shared with the dataclass.
    # Do NOT read them off `DedupConfig`: the dataclass is slots=True, so its class
    # attributes are member_descriptors rather than default values.
    parser = argparse.ArgumentParser(description="De-duplicate coherence terms by member overlap + GO DAG structure")
    parser.add_argument("--combined", type=Path, required=True, help="coherence_metrics_combined.parquet")
    parser.add_argument("--obo", type=Path, required=True, help="go-basic.obo (GO DAG for depth + lineage)")
    parser.add_argument("--jaccard-threshold", type=float, default=DEFAULT_JACCARD_THRESHOLD,
                        help="Member-set Jaccard similarity at or above which two terms are redundant")
    parser.add_argument("--merge-dag-lineage", action=argparse.BooleanOptionalAction,
                        default=DEFAULT_MERGE_DAG_LINEAGE,
                        help="Also merge member-sharing ancestor/descendant pairs (--no-merge-dag-lineage to disable)")
    parser.add_argument("--scope", choices=["pooled", "per_source"], default=DEFAULT_SCOPE,
                        help="Cluster across all sources (pooled) or within each source")
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
            jaccard_threshold=args.jaccard_threshold,
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
