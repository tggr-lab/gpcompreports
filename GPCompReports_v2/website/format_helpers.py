"""Centralized numeric formatting for user-facing tables.

All helpers emit `MISSING` for None / NaN / unparseable input so every
table cell shows the same missing-value glyph regardless of upstream
source (pandas DataFrames, CSV strings, etc.).
"""

import math

import numpy as np

MISSING = '-'

# Sparkline binning constants. Used by the landing-page top-5 cards and
# by the per-receptor report header.
DELTA_CLIP = 3.0
SPARKLINE_BINS = 40


def _is_missing(value):
    if value is None or value == '':
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def fmt_sci(value, digits=2):
    """Format as scientific notation with `digits` significant decimals."""
    if _is_missing(value):
        return MISSING
    try:
        return f"{float(value):.{digits}e}"
    except (TypeError, ValueError):
        return MISSING


def fmt_decimal(value, digits=2):
    """Format as fixed-decimal with `digits` digits after the point."""
    if _is_missing(value):
        return MISSING
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return MISSING


def fmt_int(value, default=MISSING):
    """Format as integer; missing values fall back to `default`."""
    if _is_missing(value):
        return default
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return default


def bin_signed_delta(values, n_bins=SPARKLINE_BINS, clip=DELTA_CLIP):
    """Bin signed delta_rrcs values into a symmetric histogram around 0.

    Returns a list of integer counts. The midpoint index is where the
    positive/negative split lives; consumers render bins left of midpoint
    as inactive-favoring and right of midpoint as active-favoring.
    """
    if len(values) == 0:
        return [0] * n_bins
    v = np.clip(np.asarray(values, dtype=float), -clip, clip)
    edges = np.linspace(-clip, clip, n_bins + 1)
    counts, _ = np.histogram(v, bins=edges)
    return counts.tolist()
