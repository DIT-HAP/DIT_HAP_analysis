"""
Logging Configuration
=====================

Shared logging setup for all DIT-HAP analysis scripts. Configures loguru with a
consistent format and log level. All scripts import and call ``setup_logger()``
at the start of ``main()``.

This is a library module imported by all pipeline scripts — no CLI or main().

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-02
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys

from loguru import logger

# =============================================================================
# CORE LOGIC
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
