"""Command-line interface for pyx500r — read WIFF2 files, centroid spectra,
and query precursor-product ion pairs.

Usage::

    python -m pyx500r.cli list acquisition.wiff2
    python -m pyx500r.cli transitions acquisition.wiff2
    python -m pyx500r.cli transitions "*.wiff2" --precursor-mz 456.2 --tolerance-ppm 20
    python -m pyx500r.cli transitions file.wiff2 -t "250.1587:191.0857,163.0907,109.0443"
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from ._cli_common import parse_transitions as _parse_transitions
from ._cli_common import ppm_tolerance as _ppm_tolerance
from .centroid import centroid_spectrum
from .reader import open_wiff2


def _resolve_paths(raw: Sequence[str]) -> list[Path]:
    """Glob-expand each path argument; return sorted unique resolved Paths."""
    seen: set[Path] = set()
    paths: list[Path] = []
    for pattern in raw:
        matches = glob.glob(pattern, recursive=True)
        if not matches:
            p = Path(pattern)
            if p.exists():
                matches = [str(p)]
        for m in matches:
            resolved = Path(m).resolve()
            if resolved not in seen:
                seen.add(resolved)
                paths.append(resolved)
    return sorted(paths)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read WIFF2 files, centroid spectra, query precursor-product ion pairs.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── list ──
    list_p = sub.add_parser("list", help="List samples and experiments in WIFF2 files")
    list_p.add_argument("paths", nargs="+", help="WIFF2 file path(s) or glob(s)")

    # ── transitions ──
    trans_p = sub.add_parser(
        "transitions",
        help="Find precursor-product ion pairs from DDA (MS2) spectra",
    )
    trans_p.add_argument("paths", nargs="+", help="WIFF2 file path(s) or glob(s)")
    trans_p.add_argument(
        "--sample", type=int, default=0, help="Sample index (default: 0)"
    )
    trans_p.add_argument(
        "--experiment",
        type=int,
        default=None,
        help="Experiment index filter (default: all MS2 experiments)",
    )
    trans_p.add_argument(
        "--precursor-mz",
        type=float,
        default=None,
        help="Filter by precursor m/z",
    )
    trans_p.add_argument(
        "--tolerance-ppm",
        type=float,
        default=50.0,
        help="Match tolerance in ppm for precursor and product m/z (default: 50)",
    )
    trans_p.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Minimum product ion intensity to report (default: 0)",
    )
    trans_p.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Report top N product ions per precursor (default: 20). "
             "Ignored when --transition is used.",
    )
    trans_p.add_argument(
        "--centroid-percentage",
        type=float,
        default=50.0,
        help="Centroid peak height percentage (default: 50)",
    )
    trans_p.add_argument(
        "--json", dest="json_out", action="store_true", help="Output as JSON"
    )
    trans_p.add_argument(
        "--no-centroid",
        dest="no_centroid",
        action="store_true",
        help="Skip centroiding (use raw profile peaks)",
    )
    trans_p.add_argument(
        "--progress",
        action="store_true",
        help="Show tqdm progress bar over files (requires tqdm)",
    )
    trans_p.add_argument(
        "-t", "--transition",
        action="append",
        default=None,
        metavar="PREC:PROD1,PROD2,...",
        help="Search for specific precursor→product transitions "
             '(e.g. "250.1587:191.0857,163.0907,109.0443"). '
             "Repeat for multiple transitions.",
    )

    return parser


def _cmd_list(args: argparse.Namespace) -> int:
    paths = _resolve_paths(args.paths)
    if not paths:
        print("No WIFF2 files found matching the given paths.", file=sys.stderr)
        return 1

    for i, wiff_path in enumerate(paths):
        if i > 0:
            print()
        try:
            with open_wiff2(wiff_path) as reader:
                samples = reader.list_samples()
                print(f"File: {wiff_path}")
                print(f"Samples: {len(samples)}")
                for s in samples:
                    print(f"  [{s.index}] {s.name or '(unnamed)'}")
                    experiments = reader.get_experiments(sample_index=s.index)
                    for e in experiments:
                        print(
                            f"    [{e.index}] {e.scan_type}  "
                            f"MS{e.ms_level}  {e.polarity}  "
                            f"{e.cycle_count} cycles"
                        )
        except Exception as exc:
            print(f"Error reading {wiff_path}: {exc}", file=sys.stderr)
    return 0


def _process_file(
    wiff_path: Path,
    args: argparse.Namespace,
    transitions: list[tuple[float, list[float]]],
) -> tuple[list[dict[str, Any]], int, int]:
    """Process one WIFF2 file; returns (output_entries, exp_count, spectra_count).

    When *transitions* is non-empty, output is per-scan transition matches.
    Otherwise output is aggregated precursor→top-N product ions.
    """
    sample_index = args.sample
    precursor_filter = args.precursor_mz
    tolerance_ppm = args.tolerance_ppm
    tol = _ppm_tolerance(precursor_filter, tolerance_ppm) if precursor_filter is not None else None
    min_intensity = args.min_intensity
    top_n = args.top_n
    centroid_pct = args.centroid_percentage
    no_centroid = args.no_centroid

    with open_wiff2(wiff_path) as reader:
        experiments = reader.get_experiments(sample_index=sample_index)

        if args.experiment is not None:
            target_exps = [args.experiment]
        else:
            target_exps = [e.index for e in experiments if e.ms_level == 2]

        if not target_exps:
            return [], 0, 0

        if transitions:
            return _process_file_transitions(
                reader, wiff_path, sample_index, target_exps,
                transitions, tolerance_ppm, min_intensity,
                centroid_pct, no_centroid,
            )

        # ── default mode: aggregate top-N product ions ──
        precursor_map: dict[float, list[dict[str, Any]]] = defaultdict(list)
        precursor_rt_map: dict[float, list[float]] = defaultdict(list)
        total_spectra = 0

        for exp_idx in target_exps:
            reader.prefetch_experiment(sample_index=sample_index, experiment_index=exp_idx)

            for spectrum in reader.iter_spectra(
                sample_index=sample_index,
                experiment_index=exp_idx,
                return_arrays=True,
            ):
                if spectrum.precursor_mz is None:
                    continue

                total_spectra += 1
                precursor = spectrum.precursor_mz
                rt = spectrum.scan_time

                if precursor_filter is not None and tol is not None:
                    if abs(precursor - precursor_filter) > tol:
                        continue

                if no_centroid:
                    mz_arr = spectrum.mz
                    int_arr = spectrum.intensities
                else:
                    mz_arr, int_arr = centroid_spectrum(
                        spectrum.mz,
                        spectrum.intensities,
                        centroid_percentage=centroid_pct,
                        return_arrays=True,
                    )

                product_ions = []
                for mz_val, int_val in zip(mz_arr, int_arr):
                    if int_val >= min_intensity:
                        product_ions.append((float(mz_val), float(int_val)))

                product_ions.sort(key=lambda x: x[1], reverse=True)
                product_ions = product_ions[:top_n]

                for prod_mz, prod_int in product_ions:
                    precursor_map[precursor].append(
                        {
                            "product_mz": round(prod_mz, 6),
                            "intensity": round(prod_int, 2),
                            "retention_time": round(rt, 4),
                            "experiment": exp_idx,
                            "cycle": spectrum.cycle_index,
                        }
                    )
                    precursor_rt_map[precursor].append(rt)

            reader.clear_prefetch()

    # Build output entries
    output: list[dict[str, Any]] = []
    for precursor_mz in sorted(precursor_map.keys()):
        entries = precursor_map[precursor_mz]
        rts = precursor_rt_map[precursor_mz]
        rt_range = (min(rts), max(rts)) if rts else (0, 0)

        prod_stats: dict[float, dict[str, Any]] = {}
        for entry in entries:
            key = entry["product_mz"]
            if key not in prod_stats:
                prod_stats[key] = {
                    "product_mz": key,
                    "total_intensity": 0.0,
                    "max_intensity": 0.0,
                    "count": 0,
                    "retention_times": [],
                }
            stats = prod_stats[key]
            stats["total_intensity"] += entry["intensity"]
            stats["max_intensity"] = max(stats["max_intensity"], entry["intensity"])
            stats["count"] += 1
            stats["retention_times"].append(entry["retention_time"])

        sorted_products = sorted(
            prod_stats.values(), key=lambda x: x["max_intensity"], reverse=True
        )

        output.append(
            {
                "precursor_mz": round(precursor_mz, 6),
                "rt_min": round(rt_range[0], 4),
                "rt_max": round(rt_range[1], 4),
                "scan_count": len(set(e["cycle"] for e in entries)),
                "product_ions": [
                    {
                        "mz": s["product_mz"],
                        "max_intensity": round(s["max_intensity"], 2),
                        "mean_intensity": round(s["total_intensity"] / s["count"], 2),
                        "occurrences": s["count"],
                    }
                    for s in sorted_products
                ],
                "file": str(wiff_path),
            }
        )

    return output, len(target_exps), total_spectra


def _process_file_transitions(
    reader: Any,
    wiff_path: Path,
    sample_index: int,
    target_exps: list[int],
    transitions: list[tuple[float, list[float]]],
    tolerance_ppm: float,
    min_intensity: float,
    centroid_pct: float,
    no_centroid: bool,
) -> tuple[list[dict[str, Any]], int, int]:
    """Transition mode: per-scan matching of requested precursor→product pairs."""
    # Build precursor → requested products lookup, with tolerance pre-computed
    prec_targets: dict[int, tuple[float, float, list[float]]] = {}
    for idx, (prec_mz, prod_list) in enumerate(transitions):
        prec_targets[idx] = (prec_mz, _ppm_tolerance(prec_mz, tolerance_ppm), prod_list)

    total_spectra = 0
    output: list[dict[str, Any]] = []

    for exp_idx in target_exps:
        reader.prefetch_experiment(sample_index=sample_index, experiment_index=exp_idx)

        for spectrum in reader.iter_spectra(
            sample_index=sample_index,
            experiment_index=exp_idx,
            return_arrays=True,
        ):
            if spectrum.precursor_mz is None:
                continue

            total_spectra += 1
            precursor = spectrum.precursor_mz
            rt = spectrum.scan_time

            # Find which (if any) requested transitions match this precursor
            matched_transitions: list[tuple[int, list[float]]] = []
            for tidx, (prec_mz, tol, prod_list) in prec_targets.items():
                if abs(precursor - prec_mz) <= tol:
                    matched_transitions.append((tidx, prod_list))

            if not matched_transitions:
                continue

            # Centroid
            if no_centroid:
                mz_arr = spectrum.mz
                int_arr = spectrum.intensities
            else:
                mz_arr, int_arr = centroid_spectrum(
                    spectrum.mz,
                    spectrum.intensities,
                    centroid_percentage=centroid_pct,
                    return_arrays=True,
                )

            # For each matched transition, search for requested product ions
            for tidx, prod_list in matched_transitions:
                prec_mz, _, _ = prec_targets[tidx]
                products_found: list[dict[str, Any]] = []

                for req_mz in prod_list:
                    prod_tol = _ppm_tolerance(req_mz, tolerance_ppm)
                    best_intensity = 0.0
                    best_mz = None
                    for mz_val, int_val in zip(mz_arr, int_arr):
                        if abs(float(mz_val) - req_mz) <= prod_tol:
                            if float(int_val) > best_intensity:
                                best_intensity = float(int_val)
                                best_mz = float(mz_val)
                    if best_intensity >= min_intensity and best_mz is not None:
                        products_found.append(
                            {
                                "requested": req_mz,
                                "matched": round(best_mz, 6),
                                "intensity": round(best_intensity, 2),
                            }
                        )
                    else:
                        products_found.append(
                            {
                                "requested": req_mz,
                                "matched": None,
                                "intensity": None,
                            }
                        )

                # Skip this scan if zero product ions were found
                if not any(p["intensity"] is not None for p in products_found):
                    continue

                output.append(
                    {
                        "transition_index": tidx,
                        "precursor_requested": prec_mz,
                        "precursor_matched": round(precursor, 6),
                        "retention_time": round(rt, 4),
                        "experiment": exp_idx,
                        "cycle": spectrum.cycle_index,
                        "product_ions": products_found,
                        "file": str(wiff_path),
                    }
                )

        reader.clear_prefetch()

    return output, len(target_exps), total_spectra


def _cmd_transitions(args: argparse.Namespace) -> int:
    paths = _resolve_paths(args.paths)
    if not paths:
        print("No WIFF2 files found matching the given paths.", file=sys.stderr)
        return 1

    try:
        transitions = _parse_transitions(args.transition)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    # ── progress bar setup ──
    tqdm = None
    file_iter: Any = paths
    if args.progress:
        try:
            from tqdm import tqdm as _tqdm
            tqdm = _tqdm
        except ImportError:
            print("Warning: tqdm not installed; install with 'pip install tqdm' or 'pip install pyx500r[cli]'.",
                  file=sys.stderr)

    if tqdm is not None:
        file_iter = tqdm(paths, desc="Processing", unit="file", file=sys.stderr)

    centroid_note = "" if args.no_centroid else ", centroid applied"
    any_results = False
    first_file = True
    errors = 0

    for wiff_path in file_iter:
        try:
            output, exp_count, spectra_count = _process_file(wiff_path, args, transitions)
        except Exception as exc:
            print(f"Error processing {wiff_path}: {exc}", file=sys.stderr)
            errors += 1
            continue

        if tqdm is None:
            print(
                f"{wiff_path}: {exp_count} experiment(s), {spectra_count} MS2 spectra{centroid_note}.",
                file=sys.stderr,
            )

        if not output:
            continue

        any_results = True

        _print = print
        if tqdm is not None:
            _print = file_iter.write  # tqdm.write prints above the bar

        if args.json_out:
            _print(json.dumps(output, indent=2))
        else:
            if not first_file:
                _print("")
            first_file = False
            if len(paths) > 1:
                _print(f"── {wiff_path} ──")
            if transitions:
                _print_transition_table(output, transitions, _print)
            else:
                _print_transitions_table(output, _print)

    if not any_results:
        print("No precursor-product ion pairs found.", file=sys.stderr)

    return 1 if errors == len(paths) else 0


def _print_transitions_table(output: list[dict[str, Any]], _print: Any = print) -> None:
    """Pretty-print aggregate precursor→top-N product ions."""
    header = f"{'Precursor m/z':>14s}  {'RT range':>16s}  {'Scans':>6s}  {'Top product ions (m/z → intensity)'}"
    _print(header)
    _print("-" * len(header))

    for entry in output:
        rt_str = f"{entry['rt_min']:.2f} – {entry['rt_max']:.2f}"
        product_strs = []
        for p in entry["product_ions"][:10]:
            product_strs.append(f"{p['mz']:.4f}→{p['max_intensity']:.0f}")
        products_line = "  ".join(product_strs)

        _print(
            f"{entry['precursor_mz']:14.4f}  {rt_str:>16s}  {entry['scan_count']:>5d}   {products_line}"
        )


def _print_transition_table(
    output: list[dict[str, Any]],
    transitions: list[tuple[float, list[float]]],
    _print: Any = print,
) -> None:
    """Pretty-print per-scan transition search results."""
    # Group by transition_index
    by_transition: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for entry in output:
        by_transition[entry["transition_index"]].append(entry)

    for tidx, (prec_mz, prod_list) in enumerate(transitions):
        entries = by_transition.get(tidx, [])
        if not entries:
            continue

        # Header
        prod_labels = ", ".join(f"{p:.4f}" for p in prod_list)
        _print(f"Transition: {prec_mz:.4f} → [{prod_labels}]  ({len(entries)} scans)")

        # Column headers
        col_header = f"{'RT':>8s}  {'Cycle':>5s}  {'Precursor':>12s}"
        for p in prod_list:
            col_header += f"  {p:>12.4f}"
        _print(col_header)
        _print("-" * len(col_header))

        for entry in entries:
            rt = entry["retention_time"]
            cycle = entry["cycle"]
            prec_matched = entry["precursor_matched"]
            row = f"{rt:8.2f}  {cycle:>5d}  {prec_matched:>12.4f}"
            for pi in entry["product_ions"]:
                if pi["intensity"] is not None:
                    row += f"  {pi['intensity']:>12.0f}"
                else:
                    row += f"  {'—':>12s}"
            _print(row)

        _print("")


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "list":
        return _cmd_list(args)
    elif args.command == "transitions":
        return _cmd_transitions(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
