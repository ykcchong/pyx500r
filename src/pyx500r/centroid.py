"""Centroiding dispatcher for SCIEX X500R QTOF data.

This module re-exports the monolithic JIT-compiled kernel from
:mod:`pyx500r.centroid_new`.  A pure-Python fallback is available in
:mod:`pyx500r.centroid_fallback` for reference but is no longer used at
runtime.

All public symbols (``Peak``, ``add_framing_zeros``, ``moving_average_smooth``,
``centroid_spectrum``) are always available.
"""

from __future__ import annotations

# Public symbols from the fallback module (always available for reference)
from .centroid_fallback import (
    Peak,
    add_framing_zeros,
    moving_average_smooth,
)

# Fast path: monolithic numba kernel (now required)
from .centroid_new import centroid_spectrum
