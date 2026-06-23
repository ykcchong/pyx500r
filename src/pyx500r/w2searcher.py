"""Search a precomputed MS2 product-ion index (built by ``x500rindex``).

Two modes:

1. **TUI** (default when no query flags given): interactive REPL::

       x500rsearch index.npz

2. **CLI** (with query flags): one-shot search, same interface as
   ``x500r transitions``::

       x500rsearch index.npz --precursor-mz 456.2 --tolerance-ppm 20
       x500rsearch index.npz -t "250.1587:191.0857,163.0907" --json

TUI commands::

    t "181.0717:124.0505,78.0338" ppm 20
    t "181.07:124.05,78.03" t "609.28:397.21,195.06" ppm 20
    p 456.2 ppm 20
    rt 4.5 5.5
    q
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._cli_common import parse_transitions as _parse_transitions
from ._cli_common import ppm_tolerance as _ppm_tolerance
from .index_builder import load_index


# ------------------------------------------------------------------ #
# core search (shared by TUI and CLI)
# ------------------------------------------------------------------ #

def search_tui_transitions(
    data: dict[str, np.ndarray],
    transitions: list[tuple[float, list[float]]],
    tolerance_ppm: float = 20.0,
    rt_min: float | None = None,
    rt_max: float | None = None,
    polarity: str | None = None,
) -> list[dict[str, Any]]:
    """Per-scan transition search for the TUI: **only** returns spectra
    where ALL product ions match (no partial hits).

    Also skips entries with no file_path and filters by *polarity*
    (``"pos"`` / ``"neg"``, derived from ``_P`` / ``_N`` filename suffix).
    """
    prec = data["precursor_mz"]
    rt = data["retention_time"]
    prod_mz = data["product_mz"]
    prod_int = data["product_intensity"]
    file_idx = data["file_index"]
    file_paths = np.array(data["file_paths"])
    index_top_n = int(data["top_n"])

    result: list[dict[str, Any]] = []

    for tidx, (trans_prec, trans_prods) in enumerate(transitions):
        tol = _ppm_tolerance(trans_prec, tolerance_ppm)

        mask = (prec >= trans_prec - tol) & (prec <= trans_prec + tol)
        if rt_min is not None:
            mask &= rt >= rt_min
        if rt_max is not None:
            mask &= rt <= rt_max

        if not mask.any():
            continue

        sel_prec = prec[mask]
        sel_rt = rt[mask]
        sel_pmz = prod_mz[mask]
        sel_pint = prod_int[mask]
        sel_fidx = file_idx[mask]

        for i in range(len(sel_prec)):
            # Check file_path
            fi = int(sel_fidx[i])
            if fi < 0 or fi >= len(file_paths) or not str(file_paths[fi]).strip():
                continue

            # Polarity filter (from filename _P / _N convention)
            if polarity is not None:
                fname = Path(str(file_paths[fi])).stem.upper()
                if polarity == "pos" and not fname.endswith("_P"):
                    continue
                if polarity == "neg" and not fname.endswith("_N"):
                    continue

            # Match ALL product ions
            product_ions: list[dict[str, Any]] = []
            all_found = True
            for prod_target in trans_prods:
                prod_tol = _ppm_tolerance(prod_target, tolerance_ppm)
                best_int = 0.0
                found = False
                for j in range(index_top_n):
                    pmz = float(sel_pmz[i, j])
                    pint = float(sel_pint[i, j])
                    if pmz <= 0:
                        continue
                    if abs(pmz - prod_target) <= prod_tol:
                        found = True
                        if pint > best_int:
                            best_int = pint

                if not found:
                    all_found = False
                    break  # skip this spectrum entirely

                product_ions.append({
                    "mz": prod_target,
                    "intensity": best_int,
                })

            if not all_found:
                continue

            result.append({
                "transition_index": tidx,
                "retention_time": float(sel_rt[i]),
                "precursor_matched": float(sel_prec[i]),
                "file_index": fi,
                "file_path": str(file_paths[fi]),
                "product_ions": product_ions,
            })

    return result


def search_tui_precursor(
    data: dict[str, np.ndarray],
    precursor_mz: float,
    tolerance_ppm: float = 20.0,
    top_n: int = 10,
    rt_min: float | None = None,
    rt_max: float | None = None,
) -> list[dict[str, Any]]:
    """Precursor-centric search: find all spectra near a given precursor."""
    prec = data["precursor_mz"]
    rt = data["retention_time"]
    prod_mz = data["product_mz"]
    prod_int = data["product_intensity"]
    file_idx = data["file_index"]
    file_paths = np.array(data["file_paths"])
    index_top_n = int(data["top_n"])

    tol = _ppm_tolerance(precursor_mz, tolerance_ppm)
    mask = (prec >= precursor_mz - tol) & (prec <= precursor_mz + tol)
    if rt_min is not None:
        mask &= rt >= rt_min
    if rt_max is not None:
        mask &= rt <= rt_max

    if not mask.any():
        return []

    sel_prec = prec[mask]
    sel_rt = rt[mask]
    sel_pmz = prod_mz[mask]
    sel_pint = prod_int[mask]
    sel_fidx = file_idx[mask]

    result: list[dict[str, Any]] = []
    use_top_n = min(top_n, index_top_n)

    for i in range(len(sel_prec)):
        fi = int(sel_fidx[i])
        if fi < 0 or fi >= len(file_paths) or not str(file_paths[fi]).strip():
            continue

        products = []
        for j in range(use_top_n):
            pmz = float(sel_pmz[i, j])
            pint = float(sel_pint[i, j])
            if pmz > 0:
                products.append({"mz": round(pmz, 6), "intensity": pint})

        result.append({
            "precursor_mz": float(sel_prec[i]),
            "retention_time": float(sel_rt[i]),
            "file_index": fi,
            "file_path": str(file_paths[fi]),
            "product_ions": products,
        })

    # Sort by RT
    result.sort(key=lambda x: x["retention_time"])
    return result


# ------------------------------------------------------------------ #
# TUI
# ------------------------------------------------------------------ #

_TUI_HELP = """
Commands:
  t "PREC:PROD1,PROD2,..." [t ...] [ppm N] [pos|neg]  search transitions
  p PREC [ppm N] [top N]                               search by precursor m/z
  rt MIN MAX                                            set RT range filter (minutes)
  rt clear                                              clear RT filter
  h, help                                               show this help
  q, quit, exit                                         exit
"""


def _parse_tui_line(line: str) -> dict[str, Any]:
    """Parse a TUI command line into a dict of parameters.

    Returns dict with keys: command (str), transitions (list), ppm (float),
    precursor_mz (float), top_n (int), rt_min/rt_max (float|None).
    """
    result: dict[str, Any] = {
        "command": None,
        "transitions": [],
        "ppm": 20.0,
        "precursor_mz": None,
        "top_n": 10,
        "rt_min": None,
        "rt_max": None,
        "polarity": None,  # None, "pos", or "neg"
    }

    # Tokenize: split on whitespace but keep quoted strings intact
    tokens: list[str] = []
    i = 0
    while i < len(line):
        if line[i].isspace():
            i += 1
            continue
        if line[i] == '"':
            j = i + 1
            while j < len(line) and line[j] != '"':
                if line[j] == '\\':
                    j += 1
                j += 1
            tokens.append(line[i+1:j])
            i = j + 1
        else:
            j = i
            while j < len(line) and not line[j].isspace():
                j += 1
            tokens.append(line[i:j])
            i = j

    if not tokens:
        result["command"] = "empty"
        return result

    cmd = tokens[0].lower()

    if cmd in ("q", "quit", "exit"):
        result["command"] = "quit"
        return result

    if cmd in ("h", "help"):
        result["command"] = "help"
        return result

    if cmd == "rt":
        result["command"] = "rt"
        return result

    if cmd == "t":
        result["command"] = "transitions"
        # Parse: t "prec:prod1,prod2" [t "prec:prod1,prod2"] [ppm N]
        ti = 0
        while ti < len(tokens):
            if tokens[ti].lower() == "t":
                ti += 1
                if ti < len(tokens) and ":" in tokens[ti]:
                    try:
                        prec_str, prods_str = tokens[ti].split(":", 1)
                        precursor = float(prec_str)
                        products = [float(x.strip()) for x in prods_str.split(",") if x.strip()]
                        if products:
                            result["transitions"].append((precursor, products))
                    except ValueError:
                        pass
                ti += 1
            elif tokens[ti].lower() == "ppm":
                ti += 1
                if ti < len(tokens):
                    try:
                        result["ppm"] = float(tokens[ti])
                    except ValueError:
                        pass
                ti += 1
            elif tokens[ti].lower() in ("pos", "positive"):
                result["polarity"] = "pos"
                ti += 1
            elif tokens[ti].lower() in ("neg", "negative"):
                result["polarity"] = "neg"
                ti += 1
            else:
                ti += 1
        return result

    if cmd == "p":
        result["command"] = "precursor"
        ti = 1
        while ti < len(tokens):
            if tokens[ti].lower() == "ppm":
                ti += 1
                if ti < len(tokens):
                    try:
                        result["ppm"] = float(tokens[ti])
                    except ValueError:
                        pass
            elif tokens[ti].lower() == "top":
                ti += 1
                if ti < len(tokens):
                    try:
                        result["top_n"] = int(tokens[ti])
                    except ValueError:
                        pass
            else:
                try:
                    result["precursor_mz"] = float(tokens[ti])
                except ValueError:
                    pass
            ti += 1
        return result

    result["command"] = "unknown"
    return result


def _print_tui_result_transitions(
    results: list[dict[str, Any]],
    transitions: list[tuple[float, list[float]]],
    _print: Any = print,
) -> None:
    """Print per-scan transition results grouped by RT windows,
    with a summary at the end."""
    if not results:
        _print("  (no matches)")
        return

    # Sort all results by RT
    results.sort(key=lambda r: r["retention_time"])

    # Build RT bins (1.0 unit)
    rt_bins: dict[int, list[dict[str, Any]]] = {}
    for r in results:
        bin_key = int(r["retention_time"])
        rt_bins.setdefault(bin_key, []).append(r)

    # Transition summary for headers
    by_t: dict[int, list[dict[str, Any]]] = {}
    for r in results:
        by_t.setdefault(r["transition_index"], []).append(r)

    for tidx, (prec_mz, prod_list) in enumerate(transitions):
        entries = by_t.get(tidx, [])
        prod_labels = ", ".join(f"{p:.4f}" for p in prod_list)
        _print(f"\n  t {prec_mz:.4f} → [{prod_labels}]  ({len(entries)} scans w/ all daughters)")

    # Collect all product m/z values across all transitions for the header
    all_prods = []
    for _t_prec, t_prods in transitions:
        all_prods.extend(t_prods)

    _print(f"\n  {'═' * 60}")

    for bin_key in sorted(rt_bins.keys()):
        entries = rt_bins[bin_key]
        rt_label = f"{bin_key}.0 – {bin_key}.999 min"
        _print(f"\n  ── RT {rt_label}  ({len(entries)} scans) ──")

        # Header
        hdr_parts = [f"{'RT':>7s}", f"{'Precursor':>12s}", f"{'File':>36s}"]
        for p in all_prods:
            hdr_parts.append(f"{p:>10.4f}")
        hdr = "  " + "  ".join(hdr_parts)
        _print(hdr)
        _print("  " + "-" * (len(hdr) - 2))

        for e in entries:
            fname = Path(e["file_path"]).name
            if len(fname) > 36:
                fname = fname[:33] + "..."

            row_parts = [
                f"{e['retention_time']:7.2f}",
                f"{e['precursor_matched']:>12.4f}",
                f"{fname:36s}",
            ]
            for pi in e["product_ions"]:
                row_parts.append(f"{pi['intensity']:>10.0f}")

            _print("  " + "  ".join(row_parts))

    # Summary
    _print(f"\n  {'═' * 60}")
    _print(f"  Summary by RT window:")
    for bin_key in sorted(rt_bins.keys()):
        entries = rt_bins[bin_key]
        by_t_in_bin: dict[int, int] = {}
        for e in entries:
            by_t_in_bin[e["transition_index"]] = by_t_in_bin.get(e["transition_index"], 0) + 1
        parts = []
        for tidx in sorted(by_t_in_bin.keys()):
            prec_mz, _ = transitions[tidx]
            parts.append(f"t{prec_mz:.4f}={by_t_in_bin[tidx]}")
        _print(f"    RT {bin_key}.0–{bin_key}.999: {len(entries)} scans  ({', '.join(parts)})")

    parts_total = []
    for tidx, (prec_mz, _) in enumerate(transitions):
        n = len(by_t.get(tidx, []))
        parts_total.append(f"t{prec_mz:.4f}={n}")
    _print(f"  Total: {len(results)} scans  ({', '.join(parts_total)})")


def _print_tui_result_precursor(
    results: list[dict[str, Any]],
    _print: Any = print,
) -> None:
    """Print precursor-centric results."""
    if not results:
        _print("  (no matches)")
        return

    # Group by file
    by_file: dict[str, list[dict[str, Any]]] = {}
    for r in results:
        fname = Path(r["file_path"]).name
        by_file.setdefault(fname, []).append(r)

    for fname, entries in sorted(by_file.items()):
        _print(f"\n  {fname}  ({len(entries)} spectra)")
        hdr = f"  {'RT':>8s}  {'Precursor':>12s}  {'Top product ions':s}"
        _print(hdr)
        _print("  " + "-" * (len(hdr) - 2))
        for e in entries:
            prods = "  ".join(
                f"{p['mz']:.4f}→{p['intensity']:.0f}"
                for p in e["product_ions"][:6]
            )
            row = f"  {e['retention_time']:8.2f}  {e['precursor_mz']:>12.4f}  {prods}"
            _print(row)


def run_tui(index_path: str | Path) -> None:
    """Interactive TUI search loop."""
    index_path = Path(index_path)
    if not index_path.exists():
        print(f"Error: index not found: {index_path}", file=sys.stderr)
        return

    print(f"Loading {index_path.name}…", end="", file=sys.stderr, flush=True)
    try:
        data = load_index(index_path)
    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        return

    n_spectra = len(data["precursor_mz"])
    n_files = len(data["file_paths"])
    print(f" {n_spectra} spectra, {n_files} files", file=sys.stderr)
    print(_TUI_HELP)

    rt_min: float | None = None
    rt_max: float | None = None

    while True:
        try:
            raw = input("w2> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break

        if not raw:
            continue

        # Check for RT command (needs state tracking)
        tokens_lower = raw.lower().split()
        if tokens_lower and tokens_lower[0] == "rt":
            parts = raw.split()
            if len(parts) >= 3 and parts[1].lower() != "clear":
                try:
                    rt_min = float(parts[1])
                    rt_max = float(parts[2])
                    print(f"  RT filter: {rt_min:.2f} – {rt_max:.2f} min")
                except ValueError:
                    print("  Invalid RT values")
            elif len(parts) >= 2 and parts[1].lower() == "clear":
                rt_min = rt_max = None
                print("  RT filter cleared")
            else:
                print(f"  Current RT filter: {rt_min} – {rt_max}")
            continue

        parsed = _parse_tui_line(raw)

        if parsed["command"] == "quit":
            break
        elif parsed["command"] == "help":
            print(_TUI_HELP)
            continue
        elif parsed["command"] == "empty":
            continue
        elif parsed["command"] == "rt":
            # handled above
            continue
        elif parsed["command"] == "unknown":
            print(f"  Unknown command. Type 'h' for help.")
            continue

        ppm = parsed["ppm"]

        if parsed["command"] == "transitions" and parsed["transitions"]:
            transitions = parsed["transitions"]
            pol = parsed.get("polarity")
            results = search_tui_transitions(
                data, transitions,
                tolerance_ppm=ppm,
                rt_min=rt_min, rt_max=rt_max,
                polarity=pol,
            )
            _print_tui_result_transitions(results, transitions)

        elif parsed["command"] == "precursor" and parsed["precursor_mz"] is not None:
            results = search_tui_precursor(
                data, parsed["precursor_mz"],
                tolerance_ppm=ppm,
                top_n=parsed["top_n"],
                rt_min=rt_min, rt_max=rt_max,
            )
            _print_tui_result_precursor(results)

        else:
            print("  Missing parameters. Type 'h' for help.")


# ------------------------------------------------------------------ #
# CLI one-shot (kept for backward compatibility)
# ------------------------------------------------------------------ #

def _print_transitions_table(output: list[dict[str, Any]], _print: Any = print) -> None:
    header = f"{'Precursor m/z':>14s}  {'RT range':>16s}  {'Scans':>6s}  {'Top product ions (m/z → intensity)'}"
    _print(header)
    _print("-" * len(header))
    for entry in output:
        rt_str = f"{entry['rt_min']:.2f} – {entry['rt_max']:.2f}"
        product_strs = []
        for p in entry["product_ions"][:10]:
            product_strs.append(f"{p['mz']:.4f}→{p['max_intensity']:.0f}")
        _print(f"{entry['precursor_mz']:14.4f}  {rt_str:>16s}  {entry['scan_count']:>5d}   {'  '.join(product_strs)}")


def _print_transition_table(
    output: list[dict[str, Any]],
    transitions: list[tuple[float, list[float]]],
    _print: Any = print,
) -> None:
    by_transition: dict[int, list[dict[str, Any]]] = {}
    for entry in output:
        by_transition.setdefault(entry["transition_index"], []).append(entry)

    for tidx, (prec_mz, prod_list) in enumerate(transitions):
        entries = by_transition.get(tidx, [])
        if not entries:
            _print(f"Transition: {prec_mz:.4f} → [...]  (0 scans)")
            _print("")
            continue
        prod_labels = ", ".join(f"{p:.4f}" for p in prod_list)
        _print(f"Transition: {prec_mz:.4f} → [{prod_labels}]  ({len(entries)} scans)")
        col_header = f"{'RT':>8s}  {'Precursor':>12s}  {'File':>40s}"
        for p in prod_list:
            col_header += f"  {p:>10.4f}"
        _print(col_header)
        _print("-" * len(col_header))
        for entry in entries:
            rt = entry["retention_time"]
            prec_matched = entry["precursor_matched"]
            fpath = Path(entry.get("file_path", "")).name[:40]
            row = f"{rt:8.2f}  {prec_matched:>12.4f}  {fpath:>40s}"
            for pi in entry["product_ions"]:
                if pi["intensity"] is not None:
                    row += f"  {pi['intensity']:>10.0f}"
                else:
                    row += f"  {'—':>10s}"
            _print(row)
        _print("")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("index", help="Path to the .npz index file")
    p.add_argument("--precursor-mz", type=float, default=None)
    p.add_argument("--tolerance-ppm", type=float, default=50.0)
    p.add_argument("--min-intensity", type=float, default=0.0)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--json", dest="json_out", action="store_true")
    p.add_argument("-t", "--transition", action="append", default=None,
                   metavar="PREC:PROD1,PROD2,...")
    p.add_argument("--rt-min", type=float, default=None)
    p.add_argument("--rt-max", type=float, default=None)
    return p


def _has_cli_query(args: argparse.Namespace) -> bool:
    return bool(args.precursor_mz or args.transition)


# ------------------------------------------------------------------ #
# CLI one-shot search functions (allow partial product-ion matches)
# ------------------------------------------------------------------ #

def _search_index_cli(
    data: dict[str, np.ndarray],
    precursor_mz: float | None = None,
    tolerance_ppm: float = 50.0,
    min_intensity: float = 0.0,
    top_n: int = 20,
    rt_min: float | None = None,
    rt_max: float | None = None,
) -> list[dict[str, Any]]:
    prec = data["precursor_mz"]
    rt = data["retention_time"]
    prod_mz = data["product_mz"]
    prod_int = data["product_intensity"]
    file_idx = data["file_index"]
    file_paths = data["file_paths"]
    index_top_n = int(data["top_n"])

    mask = np.ones(len(prec), dtype=bool)
    if precursor_mz is not None:
        tol = _ppm_tolerance(precursor_mz, tolerance_ppm)
        mask &= (prec >= precursor_mz - tol) & (prec <= precursor_mz + tol)
    if rt_min is not None:
        mask &= rt >= rt_min
    if rt_max is not None:
        mask &= rt <= rt_max
    if not mask.any():
        return []

    sel_prec = prec[mask]
    sel_rt = rt[mask]
    sel_pmz = prod_mz[mask]
    sel_pint = prod_int[mask]
    sel_fidx = file_idx[mask]

    unique_precursors: dict[float, list[int]] = {}
    for i, p in enumerate(sel_prec):
        unique_precursors.setdefault(round(float(p), 4), []).append(i)

    result: list[dict[str, Any]] = []
    use_top_n = min(top_n, index_top_n)
    for prec_key, indices in sorted(unique_precursors.items()):
        rts = sel_rt[indices]
        prod_scores: dict[float, dict[str, Any]] = {}
        for idx in indices:
            for j in range(use_top_n):
                pmz = float(sel_pmz[idx, j])
                pint = float(sel_pint[idx, j])
                if pmz <= 0 or pint < min_intensity:
                    continue
                key = round(pmz, 6)
                if key not in prod_scores:
                    prod_scores[key] = {"mz": key, "total_intensity": 0.0,
                                        "max_intensity": 0.0, "count": 0,
                                        "retention_times": []}
                s = prod_scores[key]
                s["total_intensity"] += pint
                s["max_intensity"] = max(s["max_intensity"], pint)
                s["count"] += 1
                s["retention_times"].append(float(sel_rt[idx]))
        sorted_products = sorted(prod_scores.values(),
                                 key=lambda x: x["max_intensity"],
                                 reverse=True)[:use_top_n]
        rep_fidx = int(sel_fidx[indices[0]])
        result.append({
            "precursor_mz": prec_key,
            "rt_min": float(rts.min()), "rt_max": float(rts.max()),
            "scan_count": len(indices),
            "product_ions": sorted_products,
            "file_index": rep_fidx,
            "sample_index": 0,
            "file_path": str(file_paths[rep_fidx]) if rep_fidx < len(file_paths) else "",
        })
    return result


def _search_transitions_cli(
    data: dict[str, np.ndarray],
    transitions: list[tuple[float, list[float]]],
    tolerance_ppm: float = 50.0,
    min_intensity: float = 0.0,
    rt_min: float | None = None,
    rt_max: float | None = None,
) -> list[dict[str, Any]]:
    prec = data["precursor_mz"]
    rt = data["retention_time"]
    prod_mz = data["product_mz"]
    prod_int = data["product_intensity"]
    file_idx = data["file_index"]
    file_paths = data["file_paths"]
    index_top_n = int(data["top_n"])

    result: list[dict[str, Any]] = []
    for tidx, (trans_prec, trans_prods) in enumerate(transitions):
        tol = _ppm_tolerance(trans_prec, tolerance_ppm)
        mask = (prec >= trans_prec - tol) & (prec <= trans_prec + tol)
        if rt_min is not None:
            mask &= rt >= rt_min
        if rt_max is not None:
            mask &= rt <= rt_max
        if not mask.any():
            continue
        sel_prec = prec[mask]
        sel_rt = rt[mask]
        sel_pmz = prod_mz[mask]
        sel_pint = prod_int[mask]
        sel_fidx = file_idx[mask]
        for i in range(len(sel_prec)):
            entry: dict[str, Any] = {
                "transition_index": tidx,
                "retention_time": float(sel_rt[i]),
                "cycle": 0,
                "precursor_matched": float(sel_prec[i]),
                "file_index": int(sel_fidx[i]),
                "file_path": str(file_paths[int(sel_fidx[i])]) if int(sel_fidx[i]) < len(file_paths) else "",
                "product_ions": [],
            }
            for prod_target in trans_prods:
                prod_tol = _ppm_tolerance(prod_target, tolerance_ppm)
                best_int = 0.0
                found = False
                for j in range(index_top_n):
                    pmz = float(sel_pmz[i, j])
                    pint = float(sel_pint[i, j])
                    if pmz <= 0:
                        continue
                    if abs(pmz - prod_target) <= prod_tol:
                        found = True
                        if pint > best_int:
                            best_int = pint
                entry["product_ions"].append({
                    "mz": prod_target,
                    "intensity": best_int if found and best_int >= min_intensity else None,
                })
            result.append(entry)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not _has_cli_query(args):
        # TUI mode
        run_tui(args.index)
        return 0

    # CLI one-shot mode
    index_path = Path(args.index)
    if not index_path.exists():
        print(f"Error: index file not found: {index_path}", file=sys.stderr)
        return 1

    try:
        transitions = _parse_transitions(args.transition)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Loading index: {index_path}…", file=sys.stderr)
    try:
        data = load_index(index_path)
    except Exception as exc:
        print(f"Error loading index: {exc}", file=sys.stderr)
        return 1

    print(f"  {len(data['precursor_mz'])} spectra indexed", file=sys.stderr)

    if transitions:
        results = _search_transitions_cli(data, transitions, args.tolerance_ppm,
                                          args.min_intensity, args.rt_min, args.rt_max)
        if args.json_out:
            print(json.dumps(results, indent=2))
        elif results:
            _print_transition_table(results, transitions)
        else:
            print("No matching transitions found.")
    else:
        results = _search_index_cli(data, args.precursor_mz, args.tolerance_ppm,
                                    args.min_intensity, args.top_n,
                                    args.rt_min, args.rt_max)
        if args.json_out:
            print(json.dumps(results, indent=2, default=str))
        elif results:
            _print_transitions_table(results)
        else:
            print("No matching precursor ions found.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
