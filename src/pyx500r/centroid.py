"""Centroiding dispatcher for SCIEX WIFF2 data.

This module re-exports the fastest available implementation:

* When **numba** is installed, the monolithic JIT-compiled kernel from
  :mod:`pyx500r.centroid_new` is used (~25 ms per dense spectrum).
* Otherwise, the pure-Python + NumPy fallback from
  :mod:`pyx500r.centroid_fallback` is used (~400 ms per dense spectrum).

All public symbols (``Peak``, ``add_framing_zeros``, ``moving_average_smooth``,
``centroid_spectrum``) are always available regardless of numba status.
"""

from __future__ import annotations

# Public symbols from the fallback module (always available)
from .centroid_fallback import (
    Peak,
    add_framing_zeros,
    moving_average_smooth,
)

# Fast path: monolithic numba kernel
try:
    from .centroid_new import centroid_spectrum
except Exception:  # pragma: no cover
    from .centroid_fallback import centroid_spectrum
