"""Global permutation-null cache keyed by group size.

The permutation null for a group depends only on (n_scored_members, background
points, seed) — not on which group it is. Every group of the same size sees the
same 1,000 (or 100,000) random draws, so the per-group permutation loop in
`compute_coherence.py` was recomputing the same null hundreds of times. This
module computes each size once, across all sources, and stores both reducers the
metrics table needs:

    nulls[n] = {"median": array(n_permutations), "mean": array(n_permutations)}

The cache is exact, not an approximation: `null_for_n` uses the same
`default_rng(random_state)` + `rng.choice(idx, size=n, replace=False)` sequence
as `coherence.metrics._permutation_nulls`, so a cached lookup reproduces the
per-group computation bit for bit.

The fingerprint lets the compute stage reject a cache built from a different
background (a changed fitting table or DR threshold) instead of silently mixing
nulls.
"""
# =============================================================================
# IMPORTS
# =============================================================================
from __future__ import annotations

# 1. Standard Library Imports
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# 2. Data Processing Imports
import numpy as np

# 3. Third-party Imports
from scipy.spatial.distance import pdist

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
METHOD_KEYS = ("median", "mean")


# =============================================================================
# CORE LOGIC
# =============================================================================
def background_fingerprint(points: np.ndarray) -> str:
    """Stable short hash of the background point cloud (row order included)."""
    arr = np.ascontiguousarray(points, dtype=np.float64)
    return hashlib.sha256(arr.tobytes()).hexdigest()[:16]


def null_for_n(
    background_points: np.ndarray,
    n: int,
    n_permutations: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray]:
    """One group size's null: (median MPD, mean MPD) over `n_permutations` draws."""
    # Deliberately the same construction as _permutation_nulls: a fresh rng per
    # size, the same choice call, one pdist shared by both reducers. If either
    # changes, cached and live p-values stop being interchangeable.
    if n <= 1:
        return np.zeros(n_permutations), np.zeros(n_permutations)
    rng = np.random.default_rng(random_state)
    idx = np.arange(background_points.shape[0])
    median = np.empty(n_permutations)
    mean = np.empty(n_permutations)
    for i in range(n_permutations):
        pick = rng.choice(idx, size=n, replace=False)
        pairwise = pdist(background_points[pick])
        median[i] = np.median(pairwise)
        mean[i] = np.mean(pairwise)
    return median, mean


def build_null_cache(
    background_points: np.ndarray,
    n_values: list[int] | tuple[int, ...],
    n_permutations: int,
    random_state: int,
    jobs: int = 1,
) -> dict[int, dict[str, np.ndarray]]:
    """Compute the null for each distinct group size, optionally in parallel."""
    n_values = sorted({int(n) for n in n_values})
    if jobs <= 1 or len(n_values) <= 1:
        return {
            n: dict(zip(METHOD_KEYS, null_for_n(background_points, n, n_permutations, random_state)))
            for n in n_values
        }

    with ProcessPoolExecutor(max_workers=min(jobs, len(n_values))) as pool:
        futures = {
            n: pool.submit(null_for_n, background_points, n, n_permutations, random_state)
            for n in n_values
        }
        return {
            n: dict(zip(METHOD_KEYS, futures[n].result()))
            for n in n_values
        }


def save_null_cache(
    path: Path,
    nulls: dict[int, dict[str, np.ndarray]],
    meta: dict,
) -> None:
    """Write one .npz with a JSON meta string and one array per (method, n)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {"meta": np.asarray(json.dumps(meta, sort_keys=True))}
    for n, stats in nulls.items():
        for method in METHOD_KEYS:
            arrays[f"{method}_{n}"] = stats[method]
    np.savez(path, **arrays)


def load_null_cache(path: Path) -> tuple[dict[int, dict[str, np.ndarray]], dict]:
    """Read a cache written by `save_null_cache`; validate array lengths."""
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(str(data["meta"]))
        n_values = sorted(
            int(key.split("_", 1)[1]) for key in data.files if key != "meta"
        )
        expected = int(meta["n_permutations"])
        nulls: dict[int, dict[str, np.ndarray]] = {}
        for n in n_values:
            stats = {}
            for method in METHOD_KEYS:
                arr = data[f"{method}_{n}"]
                if arr.shape != (expected,):
                    raise ValueError(
                        f"null cache {path}: {method}_{n} has shape {arr.shape}, "
                        f"expected ({expected},)"
                    )
                stats[method] = arr
            nulls[n] = stats
    return nulls, meta


def validate_null_cache(
    meta: dict,
    *,
    n_permutations: int,
    random_state: int,
    background_points: np.ndarray,
) -> None:
    """Raise when the cache was built for a different run than the caller's."""
    checks = {
        "n_permutations": (meta.get("n_permutations"), n_permutations),
        "random_state": (meta.get("random_state"), random_state),
        "background_n": (meta.get("background_n"), int(background_points.shape[0])),
        "background_fingerprint": (
            meta.get("background_fingerprint"), background_fingerprint(background_points)
        ),
    }
    for name, (cached, wanted) in checks.items():
        if cached != wanted:
            raise ValueError(
                f"null cache mismatch on {name}: cache has {cached!r}, run needs {wanted!r}; "
                "rebuild it with build_null_cache.py"
            )
    missing = [m for m in METHOD_KEYS if m not in meta.get("methods", [])]
    if missing:
        raise ValueError(f"null cache is missing method(s) {missing}")
