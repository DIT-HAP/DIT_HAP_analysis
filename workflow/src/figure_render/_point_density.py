"""Per-point 2D kernel density estimates for scatter colouring.

Distinct from ``workflow/src/qc/density.py``, which computes per-gene
insertion-site density: this module answers "how crowded is the scatter plot
here?" for one (x, y) point cloud.

Cost is the reason this is not a two-line call. ``gaussian_kde`` evaluation is
O(n_eval x n_train), so a 50k-point panel is ~2.5e9 kernel evaluations. Above
``DENSITY_EXACT_MAX_N`` the KDE is instead fitted on a fixed-seed subsample,
evaluated once on a coarse grid, and interpolated back to every point, which
bounds the work at roughly ``subsample_n * grid_size**2`` regardless of n.

Density is computed in the coordinate space of the arrays handed in. Callers
wanting log-space density must pass explicitly log-transformed values, exactly
as they already must for the r/P statistic in ``scatter.py``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-27
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
from loguru import logger
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import gaussian_kde

# =============================================================================
# CONSTANTS
# =============================================================================
# Above this many finite pairs, exact per-point evaluation is replaced by the
# grid path. 5000 points is ~2.5e7 kernel evaluations, which runs in well under
# a second; 50k would be ~100x that.
DENSITY_EXACT_MAX_N = 5000

# Points used to fit the KDE on the grid path. The bandwidth and shape of a 2D
# density are already well determined by a few thousand draws, so a larger
# sample buys resolution the 128x128 grid cannot represent anyway.
DENSITY_KDE_SUBSAMPLE_N = 5000

# Grid resolution for the interpolated path. 128x128 keeps grid artifacts below
# the ~3 px marker size at journal panel widths.
DENSITY_GRID_SIZE = 128

# Fixed so repeat renders of the same figure are byte-identical; the pixel
# baseline harness compares PNGs exactly.
DENSITY_SEED = 0

# gaussian_kde needs at least 3 points before its covariance estimate is
# meaningful, and raises on fewer.
_MIN_POINTS_FOR_KDE = 3


# =============================================================================
# CORE LOGIC
# =============================================================================
def _kde_on_subsample(
    xf: np.ndarray,
    yf: np.ndarray,
    *,
    subsample_n: int,
    seed: int,
) -> gaussian_kde:
    """Fit a KDE to at most subsample_n of the finite pairs, sampling without replacement."""
    if len(xf) <= subsample_n:
        return gaussian_kde(np.vstack([xf, yf]))

    rng = np.random.default_rng(seed)
    picked = rng.choice(len(xf), size=subsample_n, replace=False)

    return gaussian_kde(np.vstack([xf[picked], yf[picked]]))


def _grid_interpolated_density(
    xf: np.ndarray,
    yf: np.ndarray,
    *,
    subsample_n: int,
    grid_size: int,
    seed: int,
) -> np.ndarray:
    """Evaluate a subsampled KDE once on a regular grid, then bilinearly interpolate to every point."""
    kernel = _kde_on_subsample(xf, yf, subsample_n=subsample_n, seed=seed)

    x_edges = np.linspace(xf.min(), xf.max(), grid_size)
    y_edges = np.linspace(yf.min(), yf.max(), grid_size)
    x_mesh, y_mesh = np.meshgrid(x_edges, y_edges, indexing="ij")

    grid_density = kernel(np.vstack([x_mesh.ravel(), y_mesh.ravel()])).reshape(grid_size, grid_size)

    # fill_value=None extrapolates instead of returning NaN: floating-point
    # rounding puts min/max points a hair outside the grid they defined.
    interpolate = RegularGridInterpolator(
        (x_edges, y_edges), grid_density, method="linear", bounds_error=False, fill_value=None
    )

    return interpolate(np.column_stack([xf, yf]))


def point_density(
    x: np.ndarray,
    y: np.ndarray,
    *,
    exact_max_n: int = DENSITY_EXACT_MAX_N,
    subsample_n: int = DENSITY_KDE_SUBSAMPLE_N,
    grid_size: int = DENSITY_GRID_SIZE,
    seed: int = DENSITY_SEED,
) -> np.ndarray | None:
    """Estimate 2D point density per point, or None when the cloud is too degenerate to fit a KDE.

    The result is the same length as the inputs, carrying NaN wherever either
    coordinate was non-finite, so callers can index it straight back against
    the frame the columns came from.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same shape, got {x.shape} and {y.shape}")

    finite = np.isfinite(x) & np.isfinite(y)
    xf, yf = x[finite], y[finite]

    if len(xf) < _MIN_POINTS_FOR_KDE:
        logger.warning(f"Too few finite pairs ({len(xf)}) to estimate point density")
        return None

    # A collapsed axis makes the covariance matrix singular; gaussian_kde raises
    # rather than returning anything usable.
    if np.ptp(xf) == 0 or np.ptp(yf) == 0:
        logger.warning("Cannot estimate point density: x or y has zero range")
        return None

    density = np.full(len(x), np.nan)

    try:
        if len(xf) <= exact_max_n:
            values = np.vstack([xf, yf])
            density[finite] = gaussian_kde(values)(values)
        else:
            logger.debug(
                f"Estimating density for {len(xf)} points via a {grid_size}x{grid_size} grid "
                f"fitted on {min(len(xf), subsample_n)} points"
            )
            density[finite] = _grid_interpolated_density(
                xf, yf, subsample_n=subsample_n, grid_size=grid_size, seed=seed
            )
    except np.linalg.LinAlgError as error:
        logger.warning(f"Cannot estimate point density: singular covariance matrix ({error})")
        return None

    return density
