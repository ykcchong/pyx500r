"""Build a precomputed MS2 product-ion index from a tree of ``.wiff2`` files.

The index is a compressed NumPy ``.npz`` archive containing centroided
MS2 spectra — for each spectrum: precursor m/z, retention time, and the
top-N product ions (m/z + intensity).  Designed for fast in-memory
boolean-mask searching.

Usage::

    python -m pyx500r.index_builder gus_data/data/2025 -o index.npz -j 8

    # Then search:
    import numpy as np
    data = np.load("index.npz", allow_pickle=True)
    mask = (data["precursor_mz"] > 400) & (data["precursor_mz"] < 410)
    hits = data["product_mz"][mask]   # shape: (n_hits, n_products)
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

from .centroid import centroid_spectrum
from .reader import open_wiff2


def _process_file(args: tuple[Path, int, int, bool]) -> dict[str, np.ndarray] | str:
    """Process one wiff2 file → compact per-spectra arrays. (worker function)"""
    wiff_path, file_idx, top_n, no_centroid = args

    try:
        precursors: list[float] = []
        rt_times: list[float] = []
        prod_mz: list[list[float]] = []
        prod_int: list[list[float]] = []
        file_indices: list[int] = []
        sample_indices: list[int] = []

        with open_wiff2(wiff_path) as reader:
            for s in reader.list_samples():
                exps = reader.get_experiments(s.index)
                for e in exps:
                    if e.ms_level != 2:
                        continue
                    reader.prefetch_experiment(s.index, e.index)
                    for spec in reader.iter_spectra(s.index, e.index, return_arrays=True):
                        if spec.precursor_mz is None:
                            continue
                        mz_arr = spec.mz
                        int_arr = spec.intensities
                        if len(mz_arr) == 0:
                            continue

                        if not no_centroid:
                            mz_arr, int_arr = centroid_spectrum(
                                mz_arr, int_arr,
                                centroid_percentage=50.0,
                                return_arrays=True,
                            )

                        # Sort by intensity descending, take top_n
                        order = np.argsort(int_arr)[::-1][:top_n]
                        pmz = mz_arr[order]
                        pint = int_arr[order]

                        precursors.append(float(spec.precursor_mz))
                        rt_times.append(float(spec.scan_time))
                        prod_mz.append(pmz.astype(np.float32).tolist())
                        prod_int.append(pint.astype(np.float32).tolist())
                        file_indices.append(file_idx)
                        sample_indices.append(s.index)

                    reader.clear_prefetch()

        n = len(precursors)
        if n == 0:
            return {
                "precursor_mz": np.array([], dtype=np.float32),
                "retention_time": np.array([], dtype=np.float32),
                "product_mz": np.empty((0, top_n), dtype=np.float32),
                "product_intensity": np.empty((0, top_n), dtype=np.float32),
                "file_index": np.array([], dtype=np.int32),
                "sample_index": np.array([], dtype=np.int8),
            }

        # Pad product arrays to uniform top_n width
        pmz_arr = np.zeros((n, top_n), dtype=np.float32)
        pint_arr = np.zeros((n, top_n), dtype=np.float32)
        for i in range(n):
            k = min(len(prod_mz[i]), top_n)
            pmz_arr[i, :k] = prod_mz[i][:k]
            pint_arr[i, :k] = prod_int[i][:k]

        return {
            "precursor_mz": np.array(precursors, dtype=np.float32),
            "retention_time": np.array(rt_times, dtype=np.float32),
            "product_mz": pmz_arr,
            "product_intensity": pint_arr,
            "file_index": np.array(file_indices, dtype=np.int32),
            "sample_index": np.array(sample_indices, dtype=np.int8),
        }

    except Exception as exc:
        return f"{wiff_path.name}: {exc}"


def build_index(
    root: str | Path,
    output: str | Path = "index.npz",
    *,
    top_n: int = 10,
    jobs: int | None = None,
    no_centroid: bool = False,
    glob_pattern: str = "**/*.wiff2",
) -> Path:
    """Build a MS2 product-ion index from all ``.wiff2`` files under *root*.

    Parameters
    ----------
    root : str or Path
        Directory tree to scan for ``.wiff2`` files.
    output : str or Path
        Output ``.npz`` path (compressed).
    top_n : int
        Number of top product ions to store per spectrum.
    jobs : int or None
        Number of parallel workers (default: cpu_count).
    no_centroid : bool
        If True, use raw profile peaks instead of centroided.
    glob_pattern : str
        Glob pattern relative to *root* (default: ``**/*.wiff2``).

    Returns
    -------
    Path
        The output file path.
    """
    root = Path(root).resolve()
    paths = sorted(root.glob(glob_pattern))
    if not paths:
        raise FileNotFoundError(f"No .wiff2 files found under {root}")

    n_workers = min(jobs or mp.cpu_count(), len(paths))
    out_path = Path(output).resolve()

    print(f"Indexing {len(paths)} .wiff2 files ({n_workers} workers), top_n={top_n}",
          file=sys.stderr)
    t0 = time.perf_counter()

    # Phase 1: process files in parallel
    tasks = [(p, i, top_n, no_centroid) for i, p in enumerate(paths)]
    all_results: list[dict[str, np.ndarray]] = []
    errors = 0
    total_spectra = 0

    with mp.Pool(processes=n_workers) as pool:
        for i, result in enumerate(pool.imap_unordered(_process_file, tasks)):
            if isinstance(result, str):
                # error string
                print(f"  ✗ {result}", file=sys.stderr)
                errors += 1
            else:
                all_results.append(result)
                total_spectra += len(result["precursor_mz"])

            # Single-line progress
            done = i + 1
            print(f"\r  [{done}/{len(paths)}] {total_spectra} spectra",
                  end="", file=sys.stderr)

    print(file=sys.stderr)  # end progress line

    if not all_results:
        raise RuntimeError("No spectra indexed (all files failed)")

    # Phase 2: merge
    print("Merging…", file=sys.stderr)
    merged: dict[str, np.ndarray] = {}
    for key in all_results[0]:
        merged[key] = np.concatenate([r[key] for r in all_results])

    merged["file_paths"] = np.array([str(p) for p in paths], dtype=object)
    merged["top_n"] = np.array([top_n], dtype=np.int32)

    # Phase 3: save
    print(f"Saving {total_spectra} spectra to {out_path}…", file=sys.stderr)
    np.savez_compressed(out_path, **merged)

    elapsed = time.perf_counter() - t0
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"Done in {elapsed:.1f}s — {total_spectra} spectra, "
          f"{size_mb:.1f} MB ({errors} errors, {n_workers} workers)",
          file=sys.stderr)

    return out_path


def load_index(path: str | Path) -> dict[str, np.ndarray]:
    """Load a previously built index into memory.

    Returns a dict with keys: ``precursor_mz``, ``retention_time``,
    ``product_mz``, ``product_intensity``, ``file_index``,
    ``sample_index``, ``file_paths``, ``top_n``.
    """
    data = dict(np.load(path, allow_pickle=True))
    data["top_n"] = int(data.get("top_n", 10))
    return data


# ------------------------------------------------------------------ #
# CLI
# ------------------------------------------------------------------ #
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Build a MS2 product-ion index from a tree of .wiff2 files",
    )
    p.add_argument("root", help="Directory tree to scan for .wiff2 files")
    p.add_argument("-o", "--output", default="index.npz",
                   help="Output .npz path (default: index.npz)")
    p.add_argument("-n", "--top-n", type=int, default=10,
                   help="Top N product ions per spectrum (default: 10)")
    p.add_argument("-j", "--jobs", type=int, default=None,
                   help="Parallel workers (default: cpu_count)")
    p.add_argument("--no-centroid", action="store_true",
                   help="Use raw profile peaks (skip centroiding)")
    p.add_argument("--glob", default="**/*.wiff2",
                   help="Glob pattern relative to root (default: **/*.wiff2)")

    args = p.parse_args(argv)

    try:
        build_index(
            args.root,
            args.output,
            top_n=args.top_n,
            jobs=args.jobs,
            no_centroid=args.no_centroid,
            glob_pattern=args.glob,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
