"""Centralized numeric formatting for user-facing tables.

All helpers emit `MISSING` for None / NaN / unparseable input so every
table cell shows the same missing-value glyph regardless of upstream
source (pandas DataFrames, CSV strings, etc.).
"""

import math

MISSING = '-'


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
