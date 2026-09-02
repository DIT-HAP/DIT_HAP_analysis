"""Shared input schema validation for figure rendering.

Required-column checks were duplicated verbatim across nine figure loaders.
The column *lists* are figure-specific and live in the entrypoint scripts; only
the check itself is shared.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence

import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
def require_columns(
    df: pd.DataFrame,
    required: Sequence[str],
    context: str = "input",
) -> None:
    """Raise ValueError naming every column in `required` that `df` is missing."""
    missing = [column for column in required if column not in df.columns]

    if missing:
        rendered = ", ".join(repr(column) for column in missing)
        raise ValueError(f"Missing required columns in {context}: {rendered}")
