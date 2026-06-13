"""Parser for XicManagerXic blobs embedded in the RTParts gap.

Each blob is a .NET BinaryFormatter serialised object.  Parsing is
delegated to the pure-Python ``nrbf`` library (``pip install nrbf``),
which implements the full MS-NRBF specification.
"""

from __future__ import annotations

from typing import Any

import nrbf


def parse_xic_blobs(raw: bytes | bytearray) -> list[dict[str, Any]]:
    """Parse all ``XicManagerXic`` blobs from the raw RTParts stream.

    Returns one dict per blob (Sample × Compound).  Field names use
    the original .NET conventions (e.g. ``_compoundName``,
    ``_foundAtMass``, ``_librarySearchResults``).
    """
    if isinstance(raw, bytearray):
        raw = bytes(raw)

    pattern = b"XicManagerXic"
    blobs: list[dict[str, Any]] = []
    scan_pos = 0

    while True:
        idx = raw.find(pattern, scan_pos)
        if idx == -1:
            break
        scan_pos = idx + 1

        # Backtrack to the BinaryFormatter Header (0x00).
        # First find the ClassWithMembersAndTypes (0x05) record.
        rec_start = idx
        for back in range(1, 200):
            p = idx - back
            if p < 0:
                break
            if raw[p] == 0x05:
                rec_start = p
                break
        # Then find the Header.
        header_start = rec_start
        for back in range(1, 500):
            p = rec_start - back
            if p < 0:
                break
            if raw[p] == 0x00 and p + 17 <= rec_start and raw[p + 1] == 1 and raw[p + 2] == 0:
                header_start = p
                break

        # nrbf.loads() consumes exactly one top-level object, so we
        # don't need to know the blob size in advance.
        chunk = raw[header_start : header_start + 65536]
        obj = nrbf.loads(chunk)
        if isinstance(obj, dict):
            blobs.append(obj)

    return blobs


def build_xic_index(
    blobs: list[dict[str, Any]],
    num_samples: int = 2,
) -> dict[tuple[int, int], dict[str, Any]]:
    """Build a ``(sample_index, compound_index)`` lookup table.

    The blobs are ordered by sample first, then by compound: indices
    0 … *N*−1 belong to sample 0, *N* … 2*N*−1 to sample 1, etc.
    """
    n = len(blobs) // num_samples
    index: dict[tuple[int, int], dict[str, Any]] = {}
    for s in range(num_samples):
        offset = s * n
        for j in range(n):
            blob = blobs[offset + j]
            ci = blob.get("_index", -1)
            if 0 <= ci < n:
                index[(s, ci)] = blob
    return index
