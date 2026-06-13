"""Shared helpers for the pyx500r command-line tools.

These small parsing utilities are used by both :mod:`pyx500r.cli` and
:mod:`pyx500r.w2searcher` (and the parallel variant) so the transition
syntax stays identical across every entry point.
"""

from __future__ import annotations


def ppm_tolerance(mz: float, ppm: float) -> float:
    """Convert a ppm tolerance to an absolute m/z tolerance."""
    return mz * ppm * 1e-6


def parse_transitions(raw: list[str] | None) -> list[tuple[float, list[float]]]:
    """Parse ``--transition`` strings like ``"250.1587:191.0857,163.0907"``.

    Returns a list of ``(precursor_mz, [product_mz, ...])`` tuples.

    Raises
    ------
    ValueError
        If a transition is missing the ``:`` separator or has no product ions.
    """
    if not raw:
        return []
    result: list[tuple[float, list[float]]] = []
    for s in raw:
        if ":" not in s:
            raise ValueError(
                f"Invalid transition format: {s!r} "
                f"(expected 'precursor:prod1,prod2,...')"
            )
        prec_str, prods_str = s.split(":", 1)
        precursor = float(prec_str)
        products = [float(x.strip()) for x in prods_str.split(",") if x.strip()]
        if not products:
            raise ValueError(f"Transition {s!r} has no product ions")
        result.append((precursor, products))
    return result
