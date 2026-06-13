#!/usr/bin/env python3
"""CLI for library search identification.

Usage::

    # Search a centroided MGF file against the library
    python -m pyx500r.libsearch_cli search data/libview.sqlite query.mgf

    # Interactive mode: enter m/z,intensity pairs
    python -m pyx500r.libsearch_cli search data/libview.sqlite -

    # Show library statistics
    python -m pyx500r.libsearch_cli stats data/libview.sqlite

    # Export a library spectrum
    python -m pyx500r.libsearch_cli export data/libview.sqlite <spectrum_id>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .libsearch import LibrarySearcher


def parse_mgf(path: str) -> list[dict[str, Any]]:
    """Parse an MGF (Mascot Generic Format) file into a list of spectra.

    Each spectrum dict has keys: ``mz``, ``intensity``, ``precursor_mz``,
    ``polarity``, ``title``.
    """
    spectra: list[dict[str, Any]] = []
    current: Optional[dict[str, Any]] = None
    mz_vals: list[float] = []
    int_vals: list[float] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == "BEGIN IONS":
                current = {}
                mz_vals = []
                int_vals = []
            elif line == "END IONS":
                if current is not None:
                    current["mz"] = np.array(mz_vals, dtype=np.float64)
                    current["intensity"] = np.array(int_vals, dtype=np.float64)
                    spectra.append(current)
                current = None
            elif "=" in line and current is not None:
                key, value = line.split("=", 1)
                key = key.strip().upper()
                value = value.strip()
                if key == "PEPMASS":
                    parts = value.split()
                    current["precursor_mz"] = float(parts[0])
                elif key == "CHARGE":
                    ch = value.rstrip("+-")
                    sign = value[-1] if value[-1] in "+-" else "+"
                    current["polarity"] = "POS" if sign == "+" else "NEG"
                elif key == "TITLE":
                    current["title"] = value
                elif key == "RTINSECONDS":
                    current["rt"] = float(value)
            elif current is not None and line and not line.startswith("#"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz_vals.append(float(parts[0]))
                        int_vals.append(float(parts[1]))
                    except ValueError:
                        pass

    return spectra


def parse_peaks_stdin() -> tuple[np.ndarray, np.ndarray, float, str]:
    """Read peaks from stdin, one per line: ``m/z intensity``.

    Returns ``(mz, intensity, precursor_mz, polarity)``.
    """
    print("Enter precursor m/z: ", end="", flush=True)
    precursor_mz = float(sys.stdin.readline().strip())

    print("Enter polarity (POS/NEG) [POS]: ", end="", flush=True)
    pol = sys.stdin.readline().strip().upper() or "POS"

    print("Enter peaks (m/z intensity), one per line. Empty line to finish:")
    mz_vals: list[float] = []
    int_vals: list[float] = []
    while True:
        line = sys.stdin.readline()
        if not line or line.strip() == "":
            break
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                mz_vals.append(float(parts[0]))
                int_vals.append(float(parts[1]))
            except ValueError:
                print(f"  Skipping: {line.strip()}", file=sys.stderr)

    return (
        np.array(mz_vals, dtype=np.float64),
        np.array(int_vals, dtype=np.float64),
        precursor_mz,
        pol,
    )


def cmd_search(args: argparse.Namespace) -> None:
    """Search spectra against the library."""
    searcher = LibrarySearcher(args.database)

    try:
        if args.input == "-":
            mz, intensity, prec_mz, pol = parse_peaks_stdin()
            spectra = [{"mz": mz, "intensity": intensity, "precursor_mz": prec_mz, "polarity": pol, "title": "stdin"}]
        else:
            spectra = parse_mgf(args.input)

        if not spectra:
            print("No spectra found in input.", file=sys.stderr)
            sys.exit(1)

        print(f"Loaded {len(spectra)} spectra from {args.input}")
        print()

        for i, spec in enumerate(spectra):
            title = spec.get("title", f"Spectrum {i+1}")
            prec_mz = spec.get("precursor_mz", 0.0)
            pol = spec.get("polarity", "POS")

            print(f"--- {title} ---")
            print(f"Precursor: {prec_mz:.4f} | Polarity: {pol} | Peaks: {len(spec['mz'])}")

            if prec_mz <= 0:
                print("  WARNING: no precursor m/z specified, using 0 (broad search)")
                prec_mz = 0.0

            results = searcher.search(
                spec["mz"],
                spec["intensity"],
                precursor_mz=prec_mz,
                polarity=pol,
                ppm_tol=args.ppm_tol,
                dot_product_ppm=args.dot_product_ppm,
                prescreen_n=args.prescreen_n,
                top_n=args.top_n,
            )

            if not results:
                print("  No matches found.")
            else:
                print(f"  {'Rank':<5} {'Score':<8} {'Name':<40} {'Formula':<20} {'PrecMz':<10} {'CE':<6} {'Peaks':<8} {'CAS'}")
                print(f"  {'-'*5} {'-'*8} {'-'*40} {'-'*20} {'-'*10} {'-'*6} {'-'*8} {'-'*15}")
                for j, r in enumerate(results):
                    print(
                        f"  {j+1:<5} {r['score']:<8.4f} "
                        f"{r['name'][:38]:<40} {r['formula'][:18]:<20} "
                        f"{r['precursor_mz']:<10.4f} {r['collision_energy']:<6.0f} "
                        f"{r['num_peaks']:<8} {r['cas'][:13]}"
                    )
            print()
    finally:
        searcher.close()


def cmd_stats(args: argparse.Namespace) -> None:
    """Print library statistics."""
    with LibrarySearcher(args.database) as searcher:
        stats = searcher.stats
        print("Library Statistics")
        print("==================")
        print(f"  Compounds:    {stats['compounds']:,}")
        print(f"  Spectra:      {stats['spectra']:,}")
        print(f"  Positive mode: {stats['positive']:,}")
        print(f"  Negative mode: {stats['negative']:,}")

        cur = searcher._conn.execute(
            "SELECT MIN(PrecursorMass1), MAX(PrecursorMass1) FROM MassSpectrum"
        ).fetchone()
        print(f"  Precursor m/z: {cur[0]:.4f} – {cur[1]:.4f} Da")

        print(f"\n  Libraries ({len(stats['libraries'])}):")
        for lib in stats["libraries"]:
            cnt = searcher._conn.execute(
                """
                SELECT COUNT(*) FROM CompoundLibrary cl
                JOIN Library l ON cl.LibraryId = l.Id
                WHERE l.Name = ?
                """,
                (lib,),
            ).fetchone()[0]
            print(f"    {lib:<50} {cnt:>6} compounds")


def cmd_export(args: argparse.Namespace) -> None:
    """Export a library spectrum to stdout (MGF format)."""
    with LibrarySearcher(args.database) as searcher:
        spec = searcher.get_spectrum(args.spectrum_id)
        if spec is None:
            print(f"Spectrum not found: {args.spectrum_id}", file=sys.stderr)
            sys.exit(1)

        print("BEGIN IONS")
        print(f"TITLE={spec['name']}")
        print(f"PEPMASS={spec['precursor_mz']:.6f}")
        print(f"CHARGE={'1+' if spec['polarity'] == 'POS' else '1-'}")
        if spec.get("formula"):
            print(f"FORMULA={spec['formula']}")
        if spec.get("cas"):
            print(f"CAS={spec['cas']}")
        print(f"COLLISION_ENERGY={spec['collision_energy']:.1f}")

        if spec["mz"] is not None and spec["intensity"] is not None:
            for mz, inten in zip(spec["mz"], spec["intensity"]):
                print(f"{mz:.6f}\t{inten:.1f}")

        print("END IONS")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Library search identification for MS/MS spectra",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # search
    p_search = sub.add_parser("search", help="Search spectra against the library")
    p_search.add_argument("database", help="Path to libview .sqlite database")
    p_search.add_argument(
        "input",
        help="MGF file with query spectra, or '-' for interactive stdin",
    )
    p_search.add_argument(
        "--ppm-tol", type=float, default=50.0,
        help="Precursor m/z tolerance in ppm (default: 50)",
    )
    p_search.add_argument(
        "--dot-product-ppm", type=float, default=20.0,
        help="Peak matching tolerance for dot product in ppm (default: 20)",
    )
    p_search.add_argument(
        "--prescreen-n", type=int, default=200,
        help="Max candidates from pre-screening (default: 200)",
    )
    p_search.add_argument(
        "--top-n", type=int, default=10,
        help="Number of top results to display (default: 10)",
    )

    # stats
    p_stats = sub.add_parser("stats", help="Show library statistics")
    p_stats.add_argument("database", help="Path to libview .sqlite database")

    # export
    p_export = sub.add_parser("export", help="Export a library spectrum as MGF")
    p_export.add_argument("database", help="Path to libview .sqlite database")
    p_export.add_argument("spectrum_id", help="UUID of the spectrum to export")

    args = parser.parse_args(argv)

    if args.command == "search":
        cmd_search(args)
    elif args.command == "stats":
        cmd_stats(args)
    elif args.command == "export":
        cmd_export(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
