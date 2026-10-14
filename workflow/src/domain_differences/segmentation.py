"""
Intra-gene DR Segmentation — Algorithm
======================================

Library module (no CLI, no file output) implementing TGUH: tail-greedy
unbalanced Haar, i.e. bottom-up agglomeration of adjacent regions that merges the
least-different adjacent pair first and stops at a threshold.

Fryzlewicz's original stops at a noise-scaled threshold ``sigma*sqrt(2 log n)``.
The threshold here is affine in the gene's own noise instead::

    lam = z * sigma_gene + sigma_add

which keeps the per-gene adaptivity but adds an absolute floor. That floor is what
does the work — ``sigma_gene`` has median 0.058 over the real universe, an order
of magnitude below the 0.1 additive term, so the ``z`` term only tilts the rule
against unusually noisy genes. See ``tmp/domain_scan/README.md``.

The merge statistic is the segment **mean** (the original TGUH detail coefficient
is a mean difference times a size factor; ordering by the raw mean difference
keeps the threshold in DR units).

Input
-----
One gene's DR values in CDS order, as a ``numpy`` array. Segmentation uses only
that ordering, so the coordinate system never enters here.

Output
------
Terminal segments as ``(lo, hi)`` half-open index pairs, plus the per-segment
reliability columns.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-10-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import heapq

# 2. Data Processing Imports
import numpy as np

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
MIN_SEG = 1

# MAD is biased low at small n (E[1.4826*MAD_n] = c*sigma, c=0.67 at n=3, 0.93 at
# n=10). Dividing it out is what makes the within-segment spread comparable across
# segments of different size -- without it a 3-point segment looks 50% tighter than
# a 20-point one for free. Values from 200k standard-normal draws (seed 20261012);
# above n=20 the fitted 1 - 0.83/n is within 1% out to n=64.
_MAD_C = np.array([1.0, 1.0, 0.835, 0.671, 0.7347, 0.8213, 0.8409, 0.8783,
                   0.8873, 0.9079, 0.9129, 0.9262, 0.9297, 0.9382, 0.9402,
                   0.9464, 0.9479, 0.9527, 0.9544, 0.9585, 0.9595])


# =============================================================================
# CORE LOGIC — the threshold
# =============================================================================
def noise_sigma(x: np.ndarray) -> float:
    """Robust noise scale of one gene's DR, from the differences of adjacent points."""
    d = np.diff(x)
    return 0.0 if len(d) == 0 else float(1.4826 * np.median(np.abs(d)) / np.sqrt(2))


def merge_threshold(x: np.ndarray, z: float, sigma_add: float) -> float:
    """The TGUH stop rule for one gene: ``lam = z * sigma_gene + sigma_add``."""
    return max(0.0, z * noise_sigma(x) + sigma_add)


def universe_floor(min_seg: int) -> int:
    """Fewest points a gene needs to enter the scan: two segments plus a margin.

    ``2 * min_seg`` is the structural minimum (two segments of the smallest legal
    size); the extra 2 is the absolute margin that stops two- and three-point
    genes from being segmented on nothing. At ``min_seg=1`` this is 4.
    """
    return 2 + 2 * min_seg


# =============================================================================
# CORE LOGIC — TGUH
# =============================================================================
def tguh_tuned(x: np.ndarray, min_seg: int = MIN_SEG, *, z: float = 0.0,
               sigma_add: float = 0.0) -> list[tuple[int, int]]:
    """TGUH segments of x (positions kept in original order).

    Repeatedly merges the adjacent pair with the smallest |mean difference| while
    that difference is below ``lam``, and stops at the first pair which is not.
    ``min_seg`` forces a floor on segment size: a segment below it is merged into
    whichever neighbour it differs from least, before the merge loop continues.
    Any boundary the forced merges leave below ``lam`` is then dropped, re-checking
    after each collapse because merging shifts the neighbouring means.

    ``z=0, sigma_add=0`` makes the whole rule "merge everything", so a flat series
    comes back as one segment rather than n singletons.
    """
    n = len(x)
    if n == 0:
        return []
    lam = merge_threshold(x, z, sigma_add)
    if lam <= 0:
        return [(0, n)]

    start = list(range(n))
    end = [i + 1 for i in range(n)]
    ssum = list(map(float, x))
    size = [1] * n
    alive = [True] * n
    left = [i - 1 for i in range(n)]
    right = [i + 1 for i in range(n)]
    right[n - 1] = -1

    def mean(i: int) -> float:
        return ssum[i] / size[i]

    def diff(i: int) -> float:
        return abs(mean(i) - mean(right[i]))

    def merge(i: int) -> None:              # absorbs right[i] into i
        j = right[i]
        ssum[i] += ssum[j]
        size[i] += size[j]
        end[i] = end[j]
        alive[j] = False
        k = right[j]
        right[i] = k
        if k != -1:
            left[k] = i

    heap = [(diff(i), i) for i in range(n - 1)]
    heapq.heapify(heap)
    while True:
        short = [i for i in range(n) if alive[i] and size[i] < min_seg]
        if short:
            i = min(short, key=lambda k: size[k])
            cands = [(abs(mean(i) - mean(left[i])), left[i])] if left[i] != -1 else []
            if right[i] != -1:
                cands.append((abs(mean(i) - mean(right[i])), right[i]))
            if not cands:
                break
            j = min(cands)[1]
            merge(min(i, j))                 # adjacent -> lower index absorbs
            if right[min(i, j)] != -1:
                heapq.heappush(heap, (diff(min(i, j)), min(i, j)))
            continue
        if not heap:
            break
        d, i = heapq.heappop(heap)
        if not alive[i] or right[i] == -1:
            continue
        if d != diff(i):                     # stale heap entry
            heapq.heappush(heap, (diff(i), i))
            continue
        if d >= lam:                         # least-different pair already >= lam
            break
        merge(i)
        if right[i] != -1:
            heapq.heappush(heap, (diff(i), i))
        li = left[i]
        if li != -1 and right[li] != -1:
            heapq.heappush(heap, (diff(li), li))

    segs = [(start[i], end[i]) for i in range(n) if alive[i]]
    changed = True
    while changed and len(segs) > 1:
        changed = False
        means = [x[a:b].mean() for a, b in segs]
        for k in range(1, len(segs)):
            if abs(means[k - 1] - means[k]) < lam:
                segs[k - 1] = (segs[k - 1][0], segs[k][1])
                segs.pop(k)
                changed = True
                break
    return segs


# =============================================================================
# CORE LOGIC — segment reliability
# =============================================================================
def mad_consistency(n: int) -> float:
    """E[1.4826 * MAD_n] / sigma for n normal draws (-> 1 as n grows)."""
    return float(_MAD_C[n]) if n < len(_MAD_C) else 1.0 - 0.83 / n


def segment_scale(v: np.ndarray) -> float:
    """Bias-corrected robust within-segment scale (0 if the segment is degenerate)."""
    return 1.4826 * float(np.median(np.abs(v - np.median(v)))) / mad_consistency(len(v))


def segment_metrics(drs: np.ndarray, segments: list[tuple[int, int]]) -> list[dict[str, float]]:
    """Per-segment reliability columns, in the same order as ``segments``.

    H           within-segment scale / gene noise scale (``noise_sigma``).
                Dimensionless, so comparable across segments of different level and
                size. The no-structure baseline is NOT 1: the boundary is chosen to
                maximise contrast, which also shrinks the spread *inside* the
                segments, so read H >~ 1 as heterogeneity beyond what the selection
                alone buys. H = 0 means the segment is degenerate (all values equal),
                which is the normal case for a one-point segment.
    se_seg      SE of the segment mean, scale/sqrt(n). The scale is floored at the
                gene noise: a segment cannot be tighter than the gene's
                point-to-point jitter, and a small segment's raw MAD is exactly 0
                often enough that the floor is what stops it claiming infinite
                precision. (The prototype used the median's 1.253*scale/sqrt(n) from
                its DP days; the statistic here is the mean, so its own SE is used.)
    z_boundary  |mean_L - mean_R| / sqrt(se_L^2 + se_R^2) at the segment's left
                boundary (NaN for the first segment): would the boundary survive a
                resample? ~N(0,1) under no structure. Kept as a gate at |z| >= 2.
    """
    sigma_gene = noise_sigma(drs)
    scale = [segment_scale(drs[lo:hi]) for lo, hi in segments]
    se = [max(s, sigma_gene) / np.sqrt(hi - lo) for s, (lo, hi) in zip(scale, segments)]
    mean = [float(drs[lo:hi].mean()) for lo, hi in segments]
    return [
        {
            "sigma_seg": scale[k],
            "sigma_gene": sigma_gene,
            "H": scale[k] / sigma_gene if sigma_gene > 0 else np.nan,
            "se_seg": float(se[k]),
            "z_boundary": (abs(mean[k] - mean[k - 1]) / np.hypot(se[k], se[k - 1])
                           if k else np.nan),
        }
        for k in range(len(segments))
    ]


def max_cut_contrast(x: np.ndarray, min_seg: int = MIN_SEG) -> float:
    """Max over internal cuts of |mean(L) - mean(R)| — the statistic ``perm_p`` nulls."""
    n = len(x)
    return max((abs(x[:c].mean() - x[c:].mean())
                for c in range(min_seg, n - min_seg + 1)), default=0.0)


def perm_p(x: np.ndarray, obs: float, n_perm: int, rng: np.random.Generator,
           min_seg: int = MIN_SEG) -> float:
    """Permutation p for a max-over-cuts contrast: DR values shuffled against the fixed order.

    Vectorized over permutations, looping only over cuts. The null maximises over
    the same cut set as ``max_cut_contrast``, which is what makes the gene-level
    number calibrated rather than anti-conservative.
    """
    if obs <= 0:
        return 1.0
    n = len(x)
    perms = rng.permuted(np.broadcast_to(x, (n_perm, n)), axis=1)
    edges = np.zeros(n_perm)
    for cut in range(min_seg, n - min_seg + 1):
        d = np.abs(perms[:, :cut].mean(axis=1) - perms[:, cut:].mean(axis=1))
        np.maximum(edges, d, out=edges)
    return (int((edges >= obs).sum()) + 1) / (n_perm + 1)


def segment_permutation(drs: np.ndarray, segments: list[tuple[int, int]], *, n_perm: int,
                        rng: np.random.Generator, min_seg: int = MIN_SEG) -> tuple[list[float], float]:
    """The two permutation p's: per segment against the rest of the gene, and per gene.

    ``perm_p_seg`` is ANTI-CONSERVATIVE — the segment was selected for exactly this
    contrast, so it is only comparable across segments produced by the same rule.
    ``perm_p_gene`` reuses the max-over-cuts null and is the calibrated one. Both
    are floored at 1/(n_perm+1).
    """
    n = len(drs)
    perms = rng.permuted(np.broadcast_to(drs, (n_perm, n)), axis=1)
    per_segment = []
    for lo, hi in segments:
        keep = np.ones(n, bool)
        keep[lo:hi] = False
        obs = abs(drs[lo:hi].mean() - drs[keep].mean())
        null = np.abs(perms[:, lo:hi].mean(axis=1) - perms[:, keep].mean(axis=1))
        per_segment.append(float(((null >= obs).sum() + 1) / (n_perm + 1)))
    gene = perm_p(drs, max_cut_contrast(drs, min_seg), n_perm, rng, min_seg=min_seg)
    return per_segment, gene


# =============================================================================
# SELF-CHECK
# =============================================================================
def _check() -> None:
    """Self-check: the MAD constants, the threshold's degenerate cases, the p's."""
    rng = np.random.default_rng(20261012)
    z = rng.normal(size=(60_000, 21))
    for n in (3, 5, 7, 12, 21):
        med = np.median(z[:, :n], axis=1, keepdims=True)
        sim = (1.4826 * np.median(np.abs(z[:, :n] - med), axis=1)).mean()
        assert abs(sim / mad_consistency(n) - 1) < 0.02, (n, sim, mad_consistency(n))
    assert all(mad_consistency(n) < mad_consistency(n + 1) for n in range(3, 63))
    assert mad_consistency(3) < 0.70 < mad_consistency(10) < mad_consistency(60) < 1.0

    # a flat series has no noise to scale the threshold by -> one segment, not n
    flat = np.full(12, -0.5)
    assert tguh_tuned(flat, 1, z=1.5, sigma_add=0.1) == [(0, 12)]
    assert universe_floor(1) == 4 and universe_floor(2) == 6

    # two clean levels -> the change point comes back with a large boundary z and a
    # small calibrated p; H stays low because 12 points of MAD is a noisy scale
    x = np.concatenate([np.full(12, -0.9), np.full(12, -0.1)]) + rng.normal(scale=0.05, size=24)
    segs = tguh_tuned(x, 1, z=1.5, sigma_add=0.1)
    assert segs == [(0, 12), (12, 24)], segs
    metrics = segment_metrics(x, segs)
    assert np.isnan(metrics[0]["z_boundary"]) and metrics[1]["z_boundary"] > 5
    per_seg, gene = segment_permutation(x, segs, n_perm=400, rng=rng)
    assert per_seg[1] <= 0.01 and gene <= 1 / 401, (per_seg, gene)
    assert all(m["H"] < 1.5 for m in metrics), [m["H"] for m in metrics]

    # the floor is what stops a two-point gene from ever being split: its only
    # contrast is also its only difference, so lam = 1.5*1.048*|d| + 0.1 > |d|
    d = np.array([-0.9, 0.3])
    assert len(tguh_tuned(d, 1, z=1.5, sigma_add=0.1)) == 1

    # min_seg is enforced against live means, not just by growing one point
    ragged = np.array([-1.0, -0.95, 0.4, 0.45, 0.5, 0.55])
    assert all(hi - lo >= 2 for lo, hi in tguh_tuned(ragged, 2, z=1.0, sigma_add=0.0))
    print("OK  mad_consistency / merge threshold / min_seg / reliability / permutation p")


if __name__ == "__main__":
    _check()
