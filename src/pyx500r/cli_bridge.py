#!/usr/bin/env python3
"""Interactive qsession explorer — query compounds by index, inspect peak + XIC data.

Usage::

    qsession <qsession_path> <wiff_glob> [--match-by name|position]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from pyx500r.bridge import WiffQSessionBridge
from pyx500r.qsession import QSessionReader


def _ionisation_mode(name: str) -> str:
    upper = name.upper()
    if upper.endswith("_P"): return "P"
    if upper.endswith("_N"): return "N"
    return "?"


def _print_peak(up, idx: int, verbose: bool = False, extract_xic: bool = False,
                br=None):
    """Print one UnifiedPeak as a compact line, or full detail if verbose."""
    mode = _ionisation_mode(br.samples[up.sample_index].sample_name) if br else "?"
    area = f"{up.area:>10.1f}" if up.area else "        —"
    rt = f"{up.retention_time:>8.3f}" if up.retention_time else "      —"
    mass = f"{up.found_mass:>10.4f}" if up.found_mass else "        —"
    me = f"{up.mass_error*1e6:>6.1f}" if up.mass_error is not None and up.mass_error != 0 else ("    —" if up.mass_error is None else "   0.0")
    iso = f"{up.isotope_diff:>7.1f}" if up.isotope_diff is not None else "     —"
    rd = f"{up.rt_diff:>8.3f}" if up.rt_diff is not None else "      —"

    # Best library score (*100) from non-smart-confirmation hits
    lib_score = "    —"
    if up.library_hits:
        best_fit = max((h.fit for h in up.library_hits if not h.is_smart_confirmation), default=None)
        if best_fit is not None:
            lib_score = f"{best_fit*100:>5.1f}"

    ms2 = " ✓" if up.contains_msms else "  "
    print(f"  {idx:>5d}  {up.name[:50]:50s}  {mode:>4s}  "
          f"{area}  {rt}  {rd}  {mass}  {me}  {iso}  {lib_score}{ms2}")

    if verbose:
        xic = up.xic or {}
        for key in sorted(xic):
            if key == "__class__":
                continue
            val = xic[key]
            if isinstance(val, dict) and "__class__" in val:
                val = f"<{val['__class__']}>"
            elif isinstance(val, (list, dict)):
                val = f"({len(val)} items)" if len(val) > 0 else "(empty)"
            elif val is None:
                val = "—"
            print(f"         {key}: {val}")
        c = up._compound
        print(f"         [compound] name={up.name}  group={c.group_name if c else '?'}")
        print(f"         [compound] formula={up.formula}  adduct={c.adduct_formula if c else '?'}")
        print(f"         [compound] mz_lower={up.mz_lower}  mz_upper={up.mz_upper}")
        me = up.mass_error
        if me is not None:
            print(f"         [compound] mass_error={me:.6f}  extraction_mass={up.extraction_mass}")
        iso = up.isotope_diff
        if iso is not None:
            print(f"         [compound] isotope_diff={iso:.2f}%")
        rd = up.rt_diff
        if rd is not None:
            print(f"         [compound] rt_diff={rd:.4f} min")
        print(f"         [compound] period={up.period}  experiment={up.experiment}")
        print(f"         [compound] extraction_type={c.extraction_type if c else '?'}")
        print(f"         [compound] is_analyte={up.is_analyte}  is_reportable={up.is_reportable}")
        print(f"         [compound] is_non_targeted={c.is_non_targeted if c else '?'}")
        print(f"         [compound] isotope_index={c.isotope_index if c else '?'}")
        print(f"         [compound] expected_mw={c.expected_mw if c else '?'}")
        print(f"         [compound] internal_std={up.internal_std_name}")
        print(f"         [peak] area={up.area}  corrected={up.corrected_area}")
        print(f"         [peak] height={up.height}  corrected={up.corrected_height}")
        print(f"         [peak] retention_time={up.retention_time}")
        print(f"         [peak] apex_rt={up.apex_rt}  apex_y={up.apex_y}")
        print(f"         [peak] start_rt={up.start_rt}  end_rt={up.end_rt}")
        print(f"         [peak] noise={up.noise:.1f}  sn={up.signal_to_noise:.1f}")
        print(f"         [peak] valid_integration={up.valid_integration}")
        if up.library_hits:
            print(f"         [library] ({len(up.library_hits)} hits):")
            for h in up.library_hits:
                sc = " [SMART]" if h.is_smart_confirmation else ""
                name_info = f"  [{h.name}]" if h.name else ""
                formula_info = f" {h.formula}" if h.formula else ""
                cas_info = f" CAS={h.cas}" if h.cas else ""
                print(f"           fit={h.fit:.4f}  rev={h.reverse_fit:.4f}  purity={h.purity:.4f}{sc}{name_info}{formula_info}{cas_info}")


def _run_interactive(br, args):
    """Interactive REPL."""
    compounds = br.compounds
    n_c = len(compounds)
    unified = br.unified_results()
    merged = []
    n_s = len(unified)
    for ci in range(n_c):
        best = None
        for si in range(n_s):
            up = unified[si][ci]
            if up.area and (best is None or up.area > best.area):
                best = up
        merged.append(best)

    has_msms = sum(1 for up in merged if up and up.contains_msms)
    has_calc = sum(1 for up in merged if up and up.has_been_calculated)
    has_lib = sum(1 for up in merged if up and up.library_hits)
    print(f"Compounds: {n_c}  Samples: {n_s}  Mode: {args.match_by}")
    print(f"Calculated: {has_calc}  MS/MS: {has_msms}  Library hits: {has_lib}")
    print(f"\nType a compound index (0-{n_c-1}) to inspect, or:")
    print(f"  (l)ist MS/MS  (s)earch  (e)xport CSV  (v)erbose  (q)uit")
    print()

    _HEADER = (
        f"  {'idx':>5s}  {'name':50s}  {'mode':>4s}  "
        f"{'area':>10s}  {'rt':>8s}  {'Δrt':>8s}  {'mass':>10s}  "
        f"{'me':>6s}  {'iso':>7s}  {'lib':>5s} ms2"
    )

    verbose_list = False
    while True:
        try:
            raw_cmd = input("> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if not raw_cmd:
            continue
        cmd = raw_cmd.lower()

        if cmd in ("q", "quit", "exit"):
            break
        elif cmd in ("v", "verbose"):
            verbose_list = not verbose_list
            print(f"  verbose list = {verbose_list}")
        elif cmd in ("l", "ls", "list"):
            threshold = 5000
            shown = 0
            print(_HEADER)
            for ci in range(n_c):
                up = merged[ci]
                if up is None or not up.contains_msms:
                    continue
                if up.area and up.area < threshold:
                    continue
                if up.mass_error is not None and abs(up.mass_error * 1e6) >= 20:
                    continue
                if up.isotope_diff is not None and up.isotope_diff >= 2000:
                    continue
                _print_peak(up, ci, verbose=verbose_list, br=br)
                shown += 1
            print(f"  ({shown} compounds with MS/MS, area ≥ {threshold}, |me|<20ppm, iso<2000)")
        elif cmd in ("s", "search"):
            sub = input("  substring: ").strip().lower()
            if not sub:
                continue
            threshold = 5000
            matches = [ci for ci in range(n_c)
                       if sub in compounds[ci].name.lower()]
            shown = 0
            print(_HEADER)
            for ci in matches[:100]:
                up = merged[ci]
                if up is None:
                    continue
                if up.area and up.area < threshold:
                    continue
                if up.mass_error is not None and abs(up.mass_error * 1e6) >= 20:
                    continue
                if up.isotope_diff is not None and up.isotope_diff >= 2000:
                    continue
                _print_peak(up, ci, verbose=verbose_list, br=br)
                shown += 1
            if len(matches) > 100:
                print(f"  ... ({len(matches) - 100} more)")
        elif cmd in ("e", "export"):
            path = input("  CSV path: ").strip()
            if not path:
                continue
            with open(path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["idx", "name", "mode", "area", "rt", "found_mass",
                            "found_rt", "has_calc", "has_msms", "lib_hits"])
                for ci in range(n_c):
                    up = merged[ci]
                    if up is None:
                        continue
                    mode = "?"
                    for si in range(n_s):
                        if unified[si][ci] is up:
                            mode = _ionisation_mode(br.samples[si].sample_name)
                            break
                    w.writerow([ci, up.name, mode, up.area, up.retention_time,
                                up.found_mass, up.found_rt,
                                up.has_been_calculated, up.contains_msms,
                                len(up.library_hits)])
            print(f"  wrote {path}")
        else:
            try:
                ci = int(raw_cmd)
            except ValueError:
                print(f"  unknown command: {raw_cmd}")
                continue
            if ci < 0 or ci >= n_c:
                print(f"  index out of range (0-{n_c-1})")
                continue
            up = merged[ci]
            if up is None:
                print(f"  [{ci}] {compounds[ci].name[:50]} — no data")
                continue
            _print_peak(up, ci, verbose=True, extract_xic=args.extract_xic, br=br)
            other_si = 1 - up.sample_index
            other_up = unified[other_si][ci]
            if other_up.area and other_up is not up:
                print(f"         other mode ({_ionisation_mode(br.samples[other_si].sample_name)}): "
                      f"area={other_up.area:.1f}  rt={other_up.retention_time:.3f}  "
                      f"found_mass={other_up.found_mass}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Interactive qsession explorer",
    )
    p.add_argument("qsession", help="Path to the .qsession file")
    p.add_argument(
        "wiff_glob",
        help="Glob pattern for .wiff2 files (e.g. 'data/**/*.wiff2')",
    )
    p.add_argument(
        "-m", "--match-by", choices=["name", "position"], default="name",
        help="How to pair samples (default: name)",
    )
    p.add_argument(
        "-x", "--extract-xic", action="store_true",
        help="Extract raw XIC from wiff2 when inspecting a compound",
    )
    p.add_argument(
        "--library-db",
        default=None,
        help="Path to a libview SQLite DB for resolving library hit names "
             "(optional; if omitted, library hits show GUIDs only)",
    )
    args = p.parse_args(argv)

    if args.library_db is not None and not Path(args.library_db).exists():
        print(f"Warning: library DB not found, library names will not resolve: "
              f"{args.library_db}", file=sys.stderr)

    qs_path = Path(args.qsession).resolve()
    if not qs_path.exists():
        print(f"Error: qsession not found: {qs_path}", file=sys.stderr)
        return 1

    wiff_arg = Path(args.wiff_glob)
    if wiff_arg.is_dir():
        # User passed a directory — search recursively inside it
        all_wiffs = sorted(wiff_arg.rglob("*.wiff2"))
    else:
        all_wiffs = sorted(Path().glob(args.wiff_glob))
    if not all_wiffs:
        print(f"No .wiff2 files matched by '{args.wiff_glob}'", file=sys.stderr)
        return 1

    temp_qs = QSessionReader(qs_path)
    qs_names = {s.sample_name.lower() for s in temp_qs.list_samples()}
    temp_qs.close()

    wiff_paths = [p for p in all_wiffs if p.stem.lower() in qs_names]
    if not wiff_paths:
        print("No matching .wiff2 files found.", file=sys.stderr)
        return 1

    print(f"qsession: {qs_path}")
    print(f"wiff2 files ({len(wiff_paths)}/{len(all_wiffs)} matched):")
    for wp in wiff_paths:
        print(f"  {wp}")

    with WiffQSessionBridge(qs_path, wiff_paths, match_by=args.match_by,
                            library_db=args.library_db) as br:
        if br.library_hits_resolved is not None:
            print(f"Resolved {br.library_hits_resolved} library hit name(s).")
        # Sample routing
        print("\nSamples:")
        for info in br.match_samples():
            qs = info["qsession_sample"]
            ws = info["wiff_sample"]
            wi = info["wiff_index"]
            status = f"→ wiff[{wi}] {ws.name}" if ws else "→ NO MATCH"
            print(f"  [{info['qsession_index']}] {qs.sample_name}  {status}")

        _run_interactive(br, args)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
