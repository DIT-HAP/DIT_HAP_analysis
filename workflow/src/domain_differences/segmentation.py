"""
Intra-gene DR Segmentation — Algorithm
======================================

Library module (no CLI, no file output) implementing TGUH as published
(Fryzlewicz, Ann. Statist. 2018, 46:3390): tail-greedy unbalanced Haar, a
bottom-up agglomeration of adjacent regions that merges the smallest detail first
and stops at a noise-scaled threshold.

The merge statistic is the paper's **scaled detail**::

    d = sqrt(nL * nR / (nL + nR)) * (mean_L - mean_R)

A pair is merged while ``|d| <= lambda``. Under a constant signal with iid noise,
``Var(d) = sigma^2`` for any pair of region sizes, so the ordering is
size-comparable (the raw mean difference is not: its variance is
``sigma^2 * (1/nL + 1/nR)``). The threshold is the paper's, plus one absolute
floor on the detail scale::

    lam = C * sqrt(2 * (1 + delta) * log T) * sigma_hat + sigma_add

with the production values ``C = 1``, ``delta = 0.01`` and ``sigma_add = 0.1``.
``sigma_hat`` is the MAD of first differences over ``sqrt(2)`` — the paper says
"MAD estimate" but does not pin the formula, and this is the finest-scale TGUH
detail. See ``tmp/domain_scan/README.md`` for the measurements behind the choice.

The paper runs the full T-1-merge transform and applies *connected* thresholding
(eq. 5). Merging only the sub-lambda details gives the same segmentation because
a parent is formed after its children and never has a smaller ``|d|``; the
prototype checked this against an independent implementation on all gene
profiles and ~175k synthetic series. The paper's Stage 1/2 post-processing is
off by default there and is not implemented here.

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
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
from enum import StrEnum
import itertools

# 2. Data Processing Imports
from loguru import logger
import numpy as np

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
MIN_SEG = 1
PAPER_C = 1.0  # C in lam = C*sqrt(2(1+delta)*log T)*sigma + sigma_add
PAPER_DELTA = 0.01  # the delta of the theory's threshold
PAPER_RHO = 0.01  # tail-greedy: merge ceil(rho*alpha) pairs per scale

# MAD is biased low at small n (E[1.4826*MAD_n] = c*sigma, c=0.67 at n=3, 0.93 at
# n=10). Dividing it out is what makes the within-segment spread comparable across
# segments of different size -- without it a 3-point segment looks 50% tighter than
# a 20-point one for free. Values from 200k standard-normal draws (seed 20261012);
# above n=20 the fitted 1 - 0.83/n is within 1% out to n=64.
_MAD_C = np.array(
    [
        1.0,
        1.0,
        0.835,
        0.671,
        0.7347,
        0.8213,
        0.8409,
        0.8783,
        0.8873,
        0.9079,
        0.9129,
        0.9262,
        0.9297,
        0.9382,
        0.9402,
        0.9464,
        0.9479,
        0.9527,
        0.9544,
        0.9585,
        0.9595,
    ]
)


# =============================================================================
# CORE LOGIC — the threshold
# =============================================================================
def noise_sigma(x: np.ndarray) -> float:
    """Robust noise scale of one gene's DR, from the differences of adjacent points."""
    d = np.diff(x)
    return 0.0 if len(d) == 0 else float(1.4826 * np.median(np.abs(d)) / np.sqrt(2))


def paper_lambda(
    x: np.ndarray, *, paper_c: float = PAPER_C, delta: float = PAPER_DELTA, sigma_add: float = 0.0
) -> float:
    """The paper's stop rule for one gene, in detail units.

    ``lam = C*sqrt(2(1+delta)*lnT)*sigma + sigma_add``. The ``1e-12 * max|x|``
    floor is not part of the formula: it keeps a noiseless series
    (``sigma_hat == 0``) from reading float-rounding differences between equal
    means as change points. Without it, ten -0.9s can fragment.
    """
    n = len(x)
    if n < 2:
        return 0.0
    scale = paper_c * float(np.sqrt(2.0 * (1.0 + delta) * np.log(n)))
    lam = scale * noise_sigma(x) + sigma_add
    return max(lam, 1e-12 * float(np.max(np.abs(x))))


def detail_weight(n_left: int, n_right: int) -> float:
    """``sqrt(nL*nR/(nL+nR))``: turns a mean difference into the paper's detail units."""
    return float(np.sqrt(n_left * n_right / (n_left + n_right)))


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
class _ShortOutcome(StrEnum):
    """What the short-segment forcing step did in one TGUH pass."""

    NOT_SHORT = "not_short"
    MERGED = "merged"
    STUCK = "stuck"


@dataclass
class _TGUHState:
    """Mutable region state for the tail-greedy merge loop, one live region per point."""

    start: list[int]
    end: list[int]
    ssum: list[float]
    size: list[int]
    alive: list[bool]
    left: list[int]
    right: list[int]
    alpha: int

    def mean(self, i: int) -> float:
        """Mean of region ``i``."""
        return self.ssum[i] / self.size[i]

    def detail(self, i: int) -> float:
        """Scaled detail between region ``i`` and the live region on its right."""
        j = self.right[i]
        return detail_weight(self.size[i], self.size[j]) * (self.mean(i) - self.mean(j))

    def merge(self, i: int) -> None:
        """Absorb the region on ``i``'s right into ``i``."""
        j = self.right[i]
        self.ssum[i] += self.ssum[j]
        self.size[i] += self.size[j]
        self.end[i] = self.end[j]
        self.alive[j] = False
        k = self.right[j]
        self.right[i] = k
        if k != -1:
            self.left[k] = i

    def segments(self) -> list[tuple[int, int]]:
        """Live regions as ``(lo, hi)`` half-open index pairs, in point order."""
        return [(self.start[i], self.end[i]) for i in range(len(self.size)) if self.alive[i]]


def _tguh_initial_state(x: np.ndarray) -> _TGUHState:
    """One live region per point, each seeded with that point's value as its running sum."""
    n = len(x)
    right = [i + 1 for i in range(n)]
    right[n - 1] = -1
    return _TGUHState(
        start=list(range(n)),
        end=[i + 1 for i in range(n)],
        ssum=list(map(float, x)),
        size=[1] * n,
        alive=[True] * n,
        left=[i - 1 for i in range(n)],
        right=right,
        alpha=n,
    )


def _force_short_segment(state: _TGUHState, min_seg: int) -> _ShortOutcome:
    """Merge the smallest sub-``min_seg`` region into its closer neighbour, if one exists."""
    short = [i for i in range(len(state.size)) if state.alive[i] and state.size[i] < min_seg]
    if not short:
        return _ShortOutcome.NOT_SHORT
    i = min(short, key=lambda k: state.size[k])
    candidates = (
        [(abs(state.mean(i) - state.mean(state.left[i])), state.left[i])]
        if state.left[i] != -1
        else []
    )
    if state.right[i] != -1:
        candidates.append((abs(state.mean(i) - state.mean(state.right[i])), state.right[i]))
    if not candidates:
        return _ShortOutcome.STUCK
    j = min(candidates)[1]
    state.merge(min(i, j))
    state.alpha -= 1
    return _ShortOutcome.MERGED


def _candidate_pairs(state: _TGUHState, lam: float) -> list[tuple[float, int]]:
    """Adjacent pairs whose scaled detail is within ``lam``, as ``(detail, left index)``."""
    pairs = []
    i = 0
    while i != -1:
        j = state.right[i]
        if j != -1:
            d = abs(state.detail(i))
            if d <= lam:
                pairs.append((d, i))
        i = state.right[i]
    return pairs


def _quota_merge(state: _TGUHState, pairs: list[tuple[float, int]], rho: float) -> None:
    """Merge the smallest details, up to ``ceil(rho * alpha)`` per pass, skipping stale pairs."""
    pairs.sort()
    quota = max(1, int(np.ceil(rho * state.alpha)))
    used: set[int] = set()
    merged = 0
    for _, i in pairs:
        if merged >= quota:
            break
        j = state.right[i]
        if not state.alive[i] or j == -1 or i in used or j in used:
            continue
        state.merge(i)
        used.update((i, j))
        merged += 1
        state.alpha -= 1


def _merge_residual_short_segments(
    segments: list[tuple[int, int]], x: np.ndarray, lam: float
) -> None:
    """Merge surviving neighbours whose detail is still within ``lam``, in place."""
    changed = True
    while changed and len(segments) > 1:
        changed = False
        for k in range(1, len(segments)):
            (a, b), (c, d) = segments[k - 1], segments[k]
            detail = abs(detail_weight(b - a, d - c) * (x[a:b].mean() - x[c:d].mean()))
            if detail <= lam:
                segments[k - 1] = (a, d)
                segments.pop(k)
                changed = True
                break


def _tguh_step(state: _TGUHState, *, lam: float, rho: float, min_seg: int) -> bool:
    """Run one tail-greedy pass; return False when no further merge is possible."""
    if min_seg > 1:
        outcome = _force_short_segment(state, min_seg)
        if outcome is _ShortOutcome.STUCK:
            return False
        if outcome is _ShortOutcome.MERGED:
            return True
    pairs = _candidate_pairs(state, lam)
    if pairs:
        _quota_merge(state, pairs, rho)
    return bool(pairs)


def tguh_paper(
    x: np.ndarray,
    min_seg: int = MIN_SEG,
    *,
    paper_c: float = PAPER_C,
    delta: float = PAPER_DELTA,
    rho: float = PAPER_RHO,
    sigma_add: float = 0.0,
) -> list[tuple[int, int]]:
    """TGUH as published, with an optional minimum segment size.

    Merges the adjacent pair with the smallest scaled detail while that detail is
    at most ``paper_lambda``. Tail-greedy merges at most ``ceil(rho * alpha)``
    pairs per pass; on gene-length profiles that is one, so this is the plain
    greedy loop. ``min_seg`` is the repo's addition, not the paper's: a short
    region is forced into its more similar neighbour before the threshold is
    consulted, and surviving boundaries are re-checked afterwards.

    The paper's Stage 1/2 post-processing is disabled by default there and is not
    implemented here. Merging only the sub-lambda details gives the same result as
    its full transform plus connected thresholding (eq. 5): a parent forms after
    its children and never has a smaller detail.
    """
    n = len(x)
    if n == 0:
        return []
    if n == 1:
        return [(0, 1)]
    lam = paper_lambda(x, paper_c=paper_c, delta=delta, sigma_add=sigma_add)
    state = _tguh_initial_state(x)
    while _tguh_step(state, lam=lam, rho=rho, min_seg=min_seg):
        pass

    segments = state.segments()
    if min_seg == 1:
        return segments
    _merge_residual_short_segments(segments, x, lam)
    return segments


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
    se = [
        max(s, sigma_gene) / np.sqrt(hi - lo) for s, (lo, hi) in zip(scale, segments, strict=False)
    ]
    mean = [float(drs[lo:hi].mean()) for lo, hi in segments]
    return [
        {
            "sigma_seg": scale[k],
            "sigma_gene": sigma_gene,
            "H": scale[k] / sigma_gene if sigma_gene > 0 else np.nan,
            "se_seg": float(se[k]),
            "z_boundary": (
                abs(mean[k] - mean[k - 1]) / np.hypot(se[k], se[k - 1]) if k else np.nan
            ),
        }
        for k in range(len(segments))
    ]


def region_contrasts(means: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Each segment's mean minus the weighted mean of every OTHER segment.

    The comparison is against the rest of the same gene on purpose. DR carries a
    gene-level offset — the pooled profile sits near -0.3, not 0, and only 5% of
    segments reach 0 — so an absolute cut would call almost nothing tolerant (it
    would miss 7 of 10 of the hand-curated C termini). Being a difference of
    means, this statistic is immune to that offset: adding a constant to a whole
    profile leaves it unchanged, while the boundary contrasts the segmentation
    itself uses do not have that property.

    ``weights`` are typically the segment site counts (the statistic), so a
    one-site segment cannot outvote a large one. A single segment has nothing to
    be compared against and gets 0.
    """
    means = np.asarray(means, dtype=float)
    weights = np.asarray(weights, dtype=float)
    rest_w = weights.sum() - weights
    rest_sum = float((means * weights).sum()) - means * weights
    rest = np.divide(rest_sum, rest_w, out=np.zeros_like(means), where=rest_w > 0)
    return np.where(rest_w > 0, means - rest, 0.0)


def max_cut_contrast(x: np.ndarray, min_seg: int = MIN_SEG) -> float:
    """Max over internal cuts of |mean(L) - mean(R)| — the statistic ``perm_p`` nulls."""
    n = len(x)
    return max(
        (abs(x[:c].mean() - x[c:].mean()) for c in range(min_seg, n - min_seg + 1)), default=0.0
    )


def perm_p(
    x: np.ndarray, obs: float, n_perm: int, rng: np.random.Generator, min_seg: int = MIN_SEG
) -> float:
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


def segment_permutation(
    drs: np.ndarray,
    segments: list[tuple[int, int]],
    *,
    n_perm: int,
    rng: np.random.Generator,
    min_seg: int = MIN_SEG,
) -> tuple[list[float], float]:
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

    # the threshold is the published formula on the gene's own first-difference scale
    x = np.array([0.0, 0.1, 0.2])
    published = np.sqrt(2.0 * 1.01 * np.log(3)) * noise_sigma(x)
    assert np.isclose(paper_lambda(x), published)
    assert np.isclose(paper_lambda(x, sigma_add=0.1), published + 0.1)
    assert np.isclose(detail_weight(3, 4), np.sqrt(12 / 7))

    # a flat series has no noise to scale the threshold by -> one segment, not n
    flat = np.full(12, -0.5)
    assert tguh_paper(flat, 1, sigma_add=0.1) == [(0, 12)]
    assert universe_floor(1) == 4
    assert universe_floor(2) == 6

    # two clean levels -> the change point comes back with a large boundary z and a
    # small calibrated p; H stays low because 12 points of MAD is a noisy scale
    x = np.concatenate([np.full(12, -0.9), np.full(12, -0.1)]) + rng.normal(scale=0.05, size=24)
    segs = tguh_paper(x, 1, sigma_add=0.1)
    assert segs == [(0, 12), (12, 24)], segs
    metrics = segment_metrics(x, segs)
    assert np.isnan(metrics[0]["z_boundary"])
    assert metrics[1]["z_boundary"] > 5
    per_seg, gene = segment_permutation(x, segs, n_perm=400, rng=rng)
    assert per_seg[1] <= 0.01, (per_seg, gene)
    assert gene <= 1 / 401, (per_seg, gene)
    assert all(m["H"] < 1.5 for m in metrics), [m["H"] for m in metrics]

    # a two-point gene has only its one difference to estimate noise from, and the
    # paper threshold is larger than that difference -> never split, hence the floor
    d = np.array([-0.9, 0.3])
    assert len(tguh_paper(d, 1, sigma_add=0.1)) == 1

    # min_seg is enforced against live means, not just by growing one point
    ragged = np.array([-1.0, -0.95, 0.4, 0.45, 0.5, 0.55])
    assert all(hi - lo >= 2 for lo, hi in tguh_paper(ragged, 2, sigma_add=0.1))

    # the paper knobs: larger C and a larger floor may only coarsen the partition
    noisy = rng.normal(size=20)
    by_c = [len(tguh_paper(noisy, 1, paper_c=c)) for c in (0.5, 1.0, 1.5, 3.0)]
    assert all(a >= b for a, b in itertools.pairwise(by_c)), by_c
    by_add = [len(tguh_paper(noisy, 1, sigma_add=a)) for a in (0.0, 0.1, 0.3)]
    assert all(a >= b for a, b in itertools.pairwise(by_add)), by_add
    logger.info("OK  mad_consistency / paper threshold / min_seg / reliability / permutation p")


if __name__ == "__main__":
    _check()
