"""Parallel version of ``x500r`` CLI — multiprocess-accelerated transitions.

Usage::

    x500rp list "*.wiff2"
    x500rp transitions "*.wiff2" --precursor-mz 456.2 --tolerance-ppm 20
    x500rp transitions "*.wiff2" -t "250.1587:191.0857,163.0907,109.0443" -j 8

Same interface as ``x500r`` but distributes per-file work across a
``multiprocessing.Pool`` for the ``transitions`` command.  ``list`` is
unchanged (it is I/O-bound, not CPU-bound).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path
from typing import Any

from .cli import (
    _build_parser,
    _cmd_list,
    _parse_transitions,
    _print_transition_table,
    _print_transitions_table,
    _process_file,
    _resolve_paths,
)


def _worker_transitions(
    args_tuple: tuple[Path, argparse.Namespace, list[tuple[float, list[float]]]],
) -> dict[str, Any]:
    """Worker: process one file and return serialisable result."""
    wiff_path, ns, transitions = args_tuple
    try:
        output, exp_count, spectra_count = _process_file(wiff_path, ns, transitions)
        return {
            "file": str(wiff_path),
            "output": output,
            "exp_count": exp_count,
            "spectra_count": spectra_count,
            "error": None,
        }
    except Exception as exc:
        return {
            "file": str(wiff_path),
            "output": [],
            "exp_count": 0,
            "spectra_count": 0,
            "error": str(exc),
        }


def _cmd_transitions_parallel(args: argparse.Namespace) -> int:
    """Parallel transitions — distribute files across a worker pool."""
    paths = _resolve_paths(args.paths)
    if not paths:
        print("No WIFF2 files found matching the given paths.", file=sys.stderr)
        return 1

    try:
        transitions = _parse_transitions(args.transition)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    n_workers = getattr(args, "jobs", None) or mp.cpu_count()
    n_workers = min(n_workers, len(paths))

    # ── tqdm setup ──
    tqdm = None
    if args.progress:
        try:
            from tqdm import tqdm as _tqdm
            tqdm = _tqdm
        except ImportError:
            print("Warning: tqdm not installed; install with 'pip install tqdm'",
                  file=sys.stderr)

    tasks = [(p, args, transitions) for p in paths]

    t0 = time.perf_counter()
    results: list[dict[str, Any]] = []
    errors_shown: list[str] = []

    with mp.Pool(processes=n_workers) as pool:
        if tqdm is not None:
            pbar = tqdm(
                total=len(paths), desc="Processing", unit="file",
                file=sys.stderr,
            )
            for result in pool.imap_unordered(_worker_transitions, tasks):
                results.append(result)
                if result["error"]:
                    pbar.write(f"  ✗ {Path(result['file']).name}: {result['error']}")
                pbar.update(1)
            pbar.close()
        else:
            print(file=sys.stderr)  # space before progress
            for i, result in enumerate(pool.imap_unordered(_worker_transitions, tasks)):
                results.append(result)
                if result["error"]:
                    errors_shown.append(f"  ✗ {Path(result['file']).name}: {result['error']}")
                # Single-line progress, overwritten in-place
                print(
                    f"\r  Processing [{i+1}/{len(paths)}] "
                    f"{'— ' + str(result['error']) if result['error'] else ''}",
                    end="", file=sys.stderr,
                )
            # Clear progress line and print any errors
            print("\r" + " " * 80, end="\r", file=sys.stderr)
            for err in errors_shown:
                print(err, file=sys.stderr)

    elapsed = time.perf_counter() - t0

    # Re-sort results to match original file order for output stability
    file_order = {str(p): i for i, p in enumerate(paths)}
    results.sort(key=lambda r: file_order.get(r["file"], 999999))

    # Merge output
    any_results = False
    errors = 0
    first_file = True

    for result in results:
        if result["error"]:
            errors += 1
            continue

        output = result["output"]
        if not output:
            continue

        any_results = True

        # Print per-file header
        if not args.json_out:
            if not first_file:
                print("")
            first_file = False
            if len(paths) > 1:
                print(f"── {result['file']} ──")

        if args.json_out:
            print(json.dumps(output, indent=2))
        else:
            if transitions:
                _print_transition_table(output, transitions)
            else:
                _print_transitions_table(output)

    if not any_results:
        print("No precursor-product ion pairs found.", file=sys.stderr)

    print(f"\nDone in {elapsed:.1f}s "
          f"({elapsed/len(paths):.2f}s/file, {n_workers} workers)",
          file=sys.stderr)

    return 1 if errors == len(paths) else 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    # Inject --jobs flag into the transitions subparser
    trans_parser = None
    for action in parser._subparsers._group_actions:
        if action.dest == "command":
            trans_parser = action.choices.get("transitions")
    if trans_parser is not None:
        trans_parser.add_argument(
            "-j", "--jobs", type=int, default=None,
            help="Number of parallel workers (default: cpu_count)",
        )

    args = parser.parse_args(argv)

    if args.command == "list":
        return _cmd_list(args)
    elif args.command == "transitions":
        return _cmd_transitions_parallel(args)

    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
