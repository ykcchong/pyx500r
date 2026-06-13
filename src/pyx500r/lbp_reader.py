"""Reader for SCIEX LibraryView .lbp snapshot files.

The .lbp format wraps a SqlCe database inside a proprietary container:

    [32-byte header] [~82 KB preamble (BinaryFormatter)] [chunks …]

The preamble contains a serialised ``LibrarySnapshot`` object tree holding
compound metadata (formulas, CAS, molecular weight, …).  Subsequent chunks
hold BinaryFormatter-serialised ``TransactXYData`` reference MS/MS spectra.

Schema reference (from ``LibViewMapper.xml``):

    Compound      — Id, Formula, MolecularWeight, MonoIsotopicMass, CAS, …
    CompoundName  — Id, CompoundId, Name, IsDefault, RegionId
    MassSpectrum  — Id, CompoundId, CentroidedXYData (byte[]), …
    MassSpectrumSignature — Mass0…Mass15 spectral signatures

Usage::

    from pyx500r.lbp_reader import LbpFile

    lbp = LbpFile("snapshot.lbp")
    for c in lbp.compounds:
        print(c.name, c.formula, c.cas)
        for spec in c.spectra:
            print(f"  {len(spec.mz)} peaks")
"""
from __future__ import annotations

import bisect
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ── binary constants ──────────────────────────────────────────────
LBP_MAGIC = 0x4F32AA43
HEADER_SIZE = 32
HEADER_FMT = struct.Struct("<IIIIIIII")

BF_MARKER = b"\x00\x01\x00\x00\x00\xff\xff\xff\xff"
ANCHOR = b"\xcc\x34\x43\x99\x63\x41\xab\x41\xa3\xf2\x23\xa3\x1c\xc2\x31\x80"

# ── data containers ────────────────────────────────────────────────

@dataclass
class LbpHeader:
    """Raw 32-byte header."""
    magic: int
    version: int
    table_count: int
    reserved1: int
    index_size: int
    format_version: int
    record_size: int
    reserved2: int


@dataclass
class ReferenceSpectrum:
    """A reference MS/MS spectrum from the library."""
    mz: list[float] = field(default_factory=list)
    intensity: list[float] = field(default_factory=list)

    @property
    def num_peaks(self) -> int:
        return len(self.mz)


@dataclass
class LibraryCompound:
    """A single compound entry from the library snapshot."""
    name: str = ""
    formula: str = ""
    cas: str = ""
    molecular_weight: float = 0.0
    monoisotopic_mass: float = 0.0
    identifier: str = ""
    comment: str = ""
    compound_guid: str = ""
    num_spectra: int = 0
    spectrum_guids: list[str] = field(default_factory=list)
    spectra: list[ReferenceSpectrum] = field(default_factory=list)


# ── preamble BinaryFormatter tree ─────────────────────────────────

def _read_compact_ui32(view: memoryview, pos: int) -> tuple[int, int]:
    """Decode a .NET BinaryReader ``Write7BitEncodedInt`` value."""
    val = 0; shift = 0
    while True:
        b = view[pos]; pos += 1
        val |= (b & 0x7F) << shift
        if (b & 0x80) == 0:
            break
        shift += 7
    return val, pos


# ── compound extraction (existing GUID-anchor method, schema-enriched) ──

SKIP_NAMES = frozenset({
    "Clearcore2", "Sciex", "System", "mscorlib", "Version", "Culture",
    "PublicKey", "http://", " OFX", "AdditionalData", "Background",
    "CADGas", "Collision", "CompoundName", "Intensity", "Spectrum",
    "centroided", "entroided", "North America", "Europe", "Japan",
    "Australia, New Zealand", "k__Backing",
})


def _extract_compound_records(data: bytes) -> list[dict[str, Any]]:
    """Extract compound records by scanning for GUID anchors.

    Each compound record is anchored by a 16-byte GUID preceded by a
    ���-byte marker.  The compound name follows 26 bytes after the anchor,
    terminated by a null byte.  Spectrum GUIDs follow after the name.
    """
    seen: set[str] = set()
    records: list[dict[str, Any]] = []
    pos = 0

    while pos < len(data) - 200:
        p = data.find(ANCHOR, pos)
        if p < 0:
            break

        # Name after anchor + 26 bytes, null-terminated
        ns = p + 26
        nz = data.find(b"\x00", ns)
        if nz < 0 or nz - ns > 200:
            pos = p + 1
            continue

        name = data[ns:nz].decode("ascii", errors="replace")
        # Filter noise
        if not any(c.islower() for c in name) or len(name) < 3:
            pos = p + 1
            continue
        if any(name.startswith(s) for s in SKIP_NAMES):
            pos = p + 1
            continue
        if name in seen:
            pos = p + 1
            continue
        seen.add(name)

        # Compound GUID at scan position
        scan = nz
        while scan < len(data) and data[scan] == 0:
            scan += 1
        if scan + 8 > len(data):
            pos = p + 1
            continue

        sc = struct.unpack_from("<I", data, scan)[0]
        gp = scan + 8
        if gp + 16 > len(data):
            pos = p + 1
            continue

        compound_guid = data[gp:gp + 16].hex()
        spec_guids: list[str] = []
        for s in range(min(sc, 50)):
            off = gp + 16 + s * 16
            if off + 16 > len(data):
                break
            spec_guids.append(data[off:off + 16].hex())

        records.append({
            "name": name,
            "compound_guid": compound_guid,
            "num_spectra": sc,
            "spectrum_guids": spec_guids,
            "pos": p,
        })
        pos = p + 50

    return records


def _extract_spectra(data: bytes) -> dict[str, ReferenceSpectrum]:
    """Extract reference MS/MS spectra from BinaryFormatter chunks.

    Looks for arrays of float64 values (m/z and intensity) in BF data.
    Chunks with matching array lengths are paired as m/z + intensity.
    """
    spectra: dict[str, ReferenceSpectrum] = {}
    pos = 0

    while pos < len(data) - 100:
        p = data.find(BF_MARKER, pos)
        if p < 0:
            break

        chunk = data[p:p + 200_000]

        # Find first double array (m/z)
        a1 = chunk.find(b"\x0f", 50)
        if a1 < 0:
            pos = p + 100
            continue
        try:
            al1 = struct.unpack_from("<I", chunk, a1 + 5)[0]
            if chunk[a1 + 9] != 6:  # double type tag
                pos = p + 100
                continue
            if a1 + 10 + al1 * 8 > len(chunk):
                pos = p + 100
                continue
            mz = list(struct.unpack(f"<{al1}d", chunk[a1 + 10:a1 + 10 + al1 * 8]))

            # Find second double array (intensity)
            a2 = chunk.find(b"\x0f", a1 + 10 + al1 * 8)
            if a2 < 0:
                pos = p + 100
                continue
            al2 = struct.unpack_from("<I", chunk, a2 + 5)[0]
            if chunk[a2 + 9] != 6:
                pos = p + 100
                continue
            if a2 + 10 + al2 * 8 > len(chunk):
                pos = p + 100
                continue
            intensity = list(struct.unpack(f"<{al2}d", chunk[a2 + 10:a2 + 10 + al2 * 8]))

            if len(mz) > 0 and len(mz) == len(intensity):
                spec = ReferenceSpectrum(mz=mz, intensity=intensity)
                key = f"spec_{p:08x}"
                spectra[key] = spec

            pos = p + max(al1 * 8, 1000)
        except Exception:
            pos = p + 100

    return spectra


def _link_metadata(
    compounds: list[dict[str, Any]],
    preamble_data: bytes,
) -> list[dict[str, Any]]:
    """Enrich compounds with formula/CAS/MW by scanning for field-name-tagged
    values in the preamble BinaryFormatter tree.

    Strategy: locate known field-name strings (``Formula``, ``MolecularWeight``,
    ``MonoIsotopicMass``, ``CAS``…), then read the .NET serialised value that
    follows according to BinaryFormatter conventions.
    """
    view = memoryview(preamble_data)
    # Build a map: position → (field_name, value)
    field_values: dict[int, tuple[str, Any]] = {}

    # ── float64 value patterns ──
    float_fields = {
        b"MolecularWeight",
        b"MonoIsotopicMass",
        b"Mass0", b"Mass1", b"Mass2",
    }
    # ── string value patterns ──
    string_fields = {
        b"Formula",
        b"CAS",
        b"Name",
        b"Identifier",
        b"Comment",
    }

    for field_name in list(float_fields | string_fields):
        search_start = 0
        while True:
            idx = preamble_data.find(field_name, search_start)
            if idx < 0:
                break
            search_start = idx + 1

            # The value typically follows within 20–100 bytes after the field name.
            # Look for the pattern: value tag next, then the data.
            probe_start = idx + len(field_name)
            probe_end = min(probe_start + 100, len(preamble_data))

            if field_name in float_fields:
                # Look for a float64 (BF tag 0x0F) within the probe window
                for off in range(probe_start, probe_end):
                    if preamble_data[off] == 0x0F and off + 9 <= len(preamble_data):
                        val = struct.unpack_from("<d", preamble_data, off + 1)[0]
                        if 0.0 < val < 1e6:
                            field_values[off] = (field_name.decode(), val)
                            break
            elif field_name in string_fields:
                # Look for a string (BF tag 0x03 or 0x01 prefix) within the probe window
                for off in range(probe_start, probe_end):
                    if preamble_data[off] == 0x03:
                        slen, spos = _read_compact_ui32(view, off + 1)
                        if spos + slen <= len(preamble_data):
                            text = bytes(view[spos:spos + slen]).decode("utf-8", errors="replace")
                            if text and not text.startswith("\x00"):
                                field_values[off] = (field_name.decode(), text)
                                break
                    elif preamble_data[off] == 0x01:
                        # Sometimes strings are preceded by 0x01 record ref
                        for off2 in range(off + 1, min(off + 20, probe_end)):
                            if preamble_data[off2] == 0x03:
                                slen, spos = _read_compact_ui32(view, off2 + 1)
                                if spos + slen <= len(preamble_data):
                                    text = bytes(view[spos:spos + slen]).decode("utf-8", errors="replace")
                                    if text and not text.startswith("\x00"):
                                        field_values[off2] = (field_name.decode(), text)
                                        break
                        break

    # Assign values to the nearest compound (by position)
    compound_positions = [(c["pos"], c) for c in compounds]
    compound_positions.sort()
    cpos_list = [cp for cp, _ in compound_positions]

    for fpos, (fname, fval) in field_values.items():
        # Find nearest compound
        idx = bisect.bisect_left(cpos_list, fpos)
        best_c = None
        best_dist = 999999
        for j in range(max(0, idx - 2), min(len(compound_positions), idx + 2)):
            cpos, c = compound_positions[j]
            dist = abs(fpos - cpos)
            if dist < best_dist and dist < 100_000:
                best_dist = dist
                best_c = c
        if best_c is not None:
            key = fname[0].lower() + fname[1:] if fname else fname
            if key == "CAS":
                best_c["cas"] = str(fval)
            elif key == "Formula":
                best_c["formula"] = str(fval)
            elif key == "MolecularWeight":
                best_c["molecular_weight"] = float(fval)
            elif key == "MonoIsotopicMass":
                best_c["monoisotopic_mass"] = float(fval)
            elif key == "Name":
                pass  # compound already has name
            else:
                best_c.setdefault("_extra", {})[fname] = fval

    return compounds


# ── public API ────────────────────────────────────────────────────

class LbpFile:
    """Reader for SCIEX LibraryView .lbp snapshot files.

    Parameters
    ----------
    path:
        Path to the .lbp file.
    lazy:
        If True, delay parsing until first property access.
    """

    def __init__(self, path: str | Path, lazy: bool = True) -> None:
        self._path = Path(path)
        self._raw: bytes | None = None
        self._header: LbpHeader | None = None
        self._compounds: list[LibraryCompound] | None = None
        self._spectra: dict[str, ReferenceSpectrum] | None = None
        if not lazy:
            self._load()

    # ── properties ──────────────────────────────────────────────────

    @property
    def header(self) -> LbpHeader:
        if self._header is None:
            self._parse_header()
        return self._header  # type: ignore[return-value]

    @property
    def compounds(self) -> list[LibraryCompound]:
        if self._compounds is None:
            self._load()
        return self._compounds or []

    @property
    def spectra(self) -> dict[str, ReferenceSpectrum]:
        if self._spectra is None:
            self._load()
        return self._spectra or {}

    def find_by_name(self, substring: str) -> list[LibraryCompound]:
        nl = substring.lower()
        return [c for c in self.compounds if nl in c.name.lower()]

    def find_by_formula(self, formula: str) -> list[LibraryCompound]:
        fu = formula.upper()
        return [c for c in self.compounds if c.formula.upper() == fu]

    # ── internals ───────────────────────────────────────────────────

    def _parse_header(self) -> None:
        with open(self._path, "rb") as f:
            raw = f.read(HEADER_SIZE)
        fields = HEADER_FMT.unpack(raw)
        if fields[0] != LBP_MAGIC:
            raise ValueError(f"Not an LBP file (magic=0x{fields[0]:08X})")
        self._header = LbpHeader(*fields)

    def _get_data(self) -> bytes:
        if self._raw is not None:
            return self._raw
        with open(self._path, "rb") as f:
            f.seek(HEADER_SIZE)
            # Read preamble to find BF start
            pre = f.read(100_000)
            bf_start = pre.find(BF_MARKER)
            if bf_start < 0:
                raise ValueError("BinaryFormatter data not found in preamble")
            # Read from BF marker onwards
            f.seek(HEADER_SIZE + bf_start)
            # Estimate remaining data size
            remaining = self._path.stat().st_size - HEADER_SIZE - bf_start
            self._raw = f.read(min(remaining, 250_000_000))
        return self._raw

    def _load(self) -> None:
        self._parse_header()
        data = self._get_data()

        # 1. Extract compound records from GUID anchors
        raw_compounds = _extract_compound_records(data)
        if not raw_compounds:
            raise ValueError("No compounds found in LBP data")

        # 2. Link metadata from preamble strings
        raw_compounds = _link_metadata(raw_compounds, data)

        # 3. Extract spectra
        raw_spectra = _extract_spectra(data)

        # 4. Build LibraryCompound objects
        clist: list[LibraryCompound] = []
        seen: set[str] = set()
        for rc in raw_compounds:
            name = rc["name"]
            if name in seen:
                continue
            seen.add(name)

            # Link spectra by GUID
            matched_spectra: list[ReferenceSpectrum] = []
            for sg in rc.get("spectrum_guids", []):
                for sk, sp in raw_spectra.items():
                    if sg in sk:
                        matched_spectra.append(sp)
                        break

            clist.append(LibraryCompound(
                name=name,
                formula=rc.get("formula", ""),
                cas=rc.get("cas", ""),
                molecular_weight=rc.get("molecular_weight", 0.0),
                compound_guid=rc.get("compound_guid", ""),
                num_spectra=rc.get("num_spectra", 0),
                spectrum_guids=rc.get("spectrum_guids", []),
                spectra=matched_spectra,
            ))

        self._compounds = clist
        self._spectra = raw_spectra

    def __repr__(self) -> str:
        return (
            f"LbpFile({self._path.name!r}, "
            f"{len(self.compounds)} compounds, "
            f"{len(self.spectra)} spectra)"
        )


def open_lbp(path: str | Path) -> LbpFile:
    """Open a LibraryView .lbp file."""
    return LbpFile(path)


# ── CLI ────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="Read SCIEX LibraryView .lbp files")
    p.add_argument("lbp_file")
    p.add_argument("-s", "--search", help="Search compounds by name")
    p.add_argument("-f", "--formula", help="Filter by formula")
    p.add_argument("-l", "--list", action="store_true", help="List all compounds")
    p.add_argument("--spectra", action="store_true", help="Show spectra")
    p.add_argument("-n", "--limit", type=int, default=50)
    args = p.parse_args(argv)

    lbp = LbpFile(args.lbp_file)
    print(lbp, file=sys.stderr)

    if args.spectra:
        print(f"Spectra: {len(lbp.spectra)}", file=sys.stderr)

    if args.search:
        for c in lbp.find_by_name(args.search)[:args.limit]:
            print(f"\n{c.name}")
            if c.formula:
                print(f"  formula: {c.formula}")
            if c.cas:
                print(f"  CAS: {c.cas}")
            if c.molecular_weight:
                print(f"  MW: {c.molecular_weight:.4f}")
            if c.num_spectra:
                print(f"  spectra: {c.num_spectra} ({len(c.spectra)} decoded)")
            if c.spectra and args.spectra:
                for i, sp in enumerate(c.spectra[:3]):
                    mz_str = ", ".join(f"{m:.4f}" for m in sp.mz[:10])
                    int_str = ", ".join(f"{i:.0f}" for i in sp.intensity[:10])
                    print(f"  spectrum {i}: {sp.num_peaks} peaks")
                    print(f"    m/z: {mz_str}{'…' if sp.num_peaks > 10 else ''}")
                    print(f"    int: {int_str}{'…' if sp.num_peaks > 10 else ''}")
    elif args.formula:
        for c in lbp.find_by_formula(args.formula)[:args.limit]:
            print(f"  {c.name}  CAS={c.cas}")
    elif args.list:
        for c in lbp.compounds[:args.limit]:
            fm = c.formula or ""
            cas = c.cas or ""
            mw = f"MW={c.molecular_weight:.1f}" if c.molecular_weight else ""
            print(f"  {c.name[:60]:60s}  {fm:20s}  {cas:15s}  {mw}")
    else:
        print(f"\n{len(lbp.compounds)} compounds, {len(lbp.spectra)} spectra.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
