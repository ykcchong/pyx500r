# SCIEX LibraryView `LibViewSnapshot_*.lbp` Snapshot Container Format

> Reverse-engineered from binary analysis of
> `gus_data/LibViewSnapshot_20260606_141157.lbp`, cross-checked against
> `monodis` disassembly of:
> - `Clearcore2.LibraryView.Packager.dll` (DB schema / `CreateStructureCeExport`)
> - `Clearcore2.LibraryView.Common.dll` (`TransactXYData`, RC2 crypto)
> - `Clearcore2.MMF.dll` (memory-mapped file helper)
> - `Clearcore2.StructuredStorage.dll`
> - `LibViewMapper.xml` (table/column mapping schema)

> **Important — two different `.lbp` layouts exist.** See
> `docs/LBP_FORMAT.md` for the *exported package* layout (an OPC/ZIP archive
> produced by `CliquidPackage.Pack`). **This document describes the
> `LibViewSnapshot_*.lbp` layout**, which is a SCIEX *paged memory-mapped
> container* wrapping a SQL Server Compact (SqlCE 3.5) database. The two are
> not interchangeable: a snapshot file is **not** a ZIP and a package file is
> **not** a paged container.

---

## 1. Identifying a snapshot file

Read the first 4 bytes:

| Offset | Bytes (LE) | uint32 value | Meaning |
|--------|------------|--------------|---------|
| 0x00   | `43 AA 32 4F` | `0x4F32AA43` | Snapshot container magic |

If the first 4 bytes are `50 4B 03 04` (`PK\x03\x04`) the file is an OPC/ZIP
*package*, not a snapshot — use `docs/LBP_FORMAT.md` instead.

---

## 2. Container header (page 0, first 32 bytes)

The container begins with a 32-byte header (this is also the start of page 0):

| Offset | Size | Type     | Field            | Example     |
|--------|------|----------|------------------|-------------|
| 0x00   | 4    | uint32   | Magic            | `0x4F32AA43`|
| 0x04   | 4    | uint32   | Reserved         | `0`         |
| 0x08   | 4    | uint32   | `table_count`    | `27`        |
| 0x0C   | 4    | uint32   | Reserved         | `0`         |
| 0x10   | 4    | uint32   | `index_size`     | `3505053`   |
| 0x14   | 4    | uint32   | `format_version` | `1`         |
| 0x18   | 4    | uint32   | `record_size`    | `578`       |
| 0x1C   | 4    | uint32   | Reserved         | `0`         |

The remainder of page 0 (offsets `0x20`+) holds the allocation-index root node
of the paged store (a B-tree mapping logical sequence numbers to physical
pages). It is not required for image extraction and is not fully documented
here.

---

## 3. Paged storage

The file is a fixed-size **paged memory-mapped store**.

| Property | Value |
|----------|-------|
| Page size | 4096 bytes (`0x1000`) |
| Per-page frame header | 16 bytes (`0x10`) |
| Per-page payload | 4080 bytes (`0xFF0`) |

Every page (including page 0) begins with a 16-byte **frame header**:

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0x00 | 4 | `checksum` | Payload hash / integrity value (per-page) |
| 0x04 | 4 | `counter` | Allocation / transaction counter — monotonic but sparse (values up to ~9.7M; **not** a page index) |
| 0x08 | 4 | `sequence` | **Logical sequence number** — orders the payloads into the logical image |
| 0x0C | 4 | `reserved` | Always `0` |

The 4080 payload bytes follow the frame header.

### Reconstructing the logical image

1. Iterate every physical page `p` (`0 … file_size/4096 - 1`).
2. Read its `sequence` (frame-header field at offset 8).
3. Build a map `sequence -> physical_page`. If two pages share a sequence
   number (a page was rewritten), **the later physical page wins**.
4. Sort by `sequence`, and concatenate each page's 4080-byte **payload**
   (skipping the 16-byte frame header).

The result is one contiguous logical byte image (~1.0 GB for the reference
file). This image holds both the SqlCE relational database and the spectrum
blobs, laid out linearly.

> **Observed values** (reference file): 257,461 physical pages,
> 257,116 unique sequence numbers (≈345 rewritten pages), sequence range
> `0 … 1,701,328`. Sequence numbering is sparse (the store pre-allocates ids).

---

## 4. What the logical image contains

### 4.1 SqlCE 3.5 relational database

The low-sequence region of the image is a **SQL Server Compact 3.5** database.
Its catalog table `__SysObjects` and all table/index/default-constraint names
are present as readable strings. The schema matches
`Clearcore2.LibraryView.Packager`'s `CreateStructureCeExport` DDL and
`LibViewMapper.xml`:

| Table | Key columns (abridged) |
|-------|------------------------|
| `Compound` | `Id`, `Identifier`, `CAS`, `Formula`, `MolecularWeight`, `MonoIsotopicMass`, `MolecularStructureSource`, `Comment`, `Active` |
| `CompoundName` | `Id`, `CompoundId`, `Name`, `RegionId`, `IsDefault` |
| `CompoundClass` | `ClassId`, `CompoundId` |
| `CompoundLibrary` | `CompoundId`, `LibraryId` |
| `CompoundFavorite` | `Username`, `CompoundId` |
| `Class` | `Id`, `Name` |
| `Region` | `Id`, `Name` |
| `Library` | `Id`, `Name`, `ParentFolderId`, `IsFolder` |
| `ScanType` | `Id`, `ScanTypeKey`, `ScanTypeName` |
| `InstrumentType` | `ModelName`, `InstrumentType`, `InstrumentKey`, `Id` |
| `Instrument` | `Id`, `InstrumentTypeId`, `Description` |
| `MassSpectrum` | `Id`, `CompoundId`, `InstrumentId`, `ScanTypeId`, `RawXYData`, `CentroidedXYData`, `PrecursorMass1`, `CollisionEnergy`, … |
| `MassSpectrumSignature` | `MassSpectrumID`, `Mass0..Mass15` |
| `UVSpectrum` | `Id`, `CompoundId`, `RawXYData`, `ProcessedXYData`, … |
| `Transition` | `Id`, `CompoundId`, `InstrumentId`, `Q1`, `Q3`, `DP`, `EP`, `FP`, `CE`, `CXP`, `Polarity`, `TransitionOrder` |
| `RetentionTime` | `Id`, `CompoundId`, `InstrumentId`, `Value`, `IsDefault` |
| `CustomAttributeDefinition` / `CustomAttributeValue` | custom attribute schema |
| `NumberAttributeValue` / `TextAttributeValue` / `BinaryAttributeValue` | typed attribute values |

Indexes seen: `IX_Compound_Active`, `IX_Compound_CAS`, `IX_Compound_Formula`,
`IX_Compound_MolecularWeight`, `IX_CompoundName_Name_RegionId`,
`IX_MassSpectrum_PrecursorMass1`, `IX_Class_Name`, `IX_RetentionTime_Value`,
plus primary keys `PK_*` and default constraints `DF__*`.

**Compound rows** are SqlCE table rows with length-prefixed columns. Example
row decoded from the reference file:

```
Identifier : "50-33-9C19H20N2O2"
Comment    : "SM041_20200228\\S1420. Verified with Sciex library."
MolecularWeight   : 362.5155
MonoIsotopicMass  : 362.2358
```

CompoundName rows carry the human-readable names (e.g. `Histidine`,
`Tripelennamine`, `Glycopyrrolate`) as SqlCE `NVarChar` strings.

> Parsing the SqlCE B-tree to fully materialise these tables is a separate
> task. This document and the extractor only recover the contiguous SqlCE
> image; downstream tooling can then parse the SqlCE pages or scan for records.

### 4.2 Spectrum blobs (`TransactXYData`)

The `MassSpectrum.CentroidedXYData` / `MassSpectrum.RawXYData` (and the UV
equivalents) are stored as **.NET BinaryFormatter (MS-NRBF)** serialized
`Clearcore2.LibraryView.Common.TransactXYData` objects. Each blob begins with
the MS-NRBF header magic `00 01 00 00 00 FF FF FF FF`.

`TransactXYData` (from disassembly):

```csharp
[Serializable]
class TransactXYData {
    double[] _xValues;   // m/z (or wavelength for UV)
    double[] _yValues;   // intensity
    string   _longName;  // e.g. "Spectrum from SM-045.wiff2 (sample 1) ... Precursor: 453.3 Da, CE: 35"
    string   _shortName;
    string   _xName;     // e.g. "Mass/Charge"
    string   _xUnits;    // e.g. "Da"
    string   _yName;     // e.g. "Intensity"
    string   _yUnits;    // e.g. "cps"
}
```

The two `double[]` arrays serialize as `ArraySinglePrimitive` records
(record type 15, primitive type 6 = Double) and are referenced from the
object's member list by `MemberReference` records. The reference file contains
~16,075 such blobs.

> A single blob frequently spans many pages, and those pages are **not**
> physically contiguous — they are scattered across the file and only linearise
> once payloads are ordered by `sequence`. This is why a naïve byte-offset scan
> of the raw file fails to recover the larger spectra.

---

## 5. Encryption note

Snapshot databases may apply SqlCE password encryption and/or column-level RC2
encryption (`Clearcore2.LibraryView.Common.Crypto`, key/IV derived from the
license string, see `docs/LBP_FORMAT.md` §6). The reference file is
unencrypted (`Encryption` columns default to `''`). Encrypted snapshots require
the SqlCE password / license key before the SqlCE image can be opened.

---

## 6. Tooling

| Script | Purpose |
|--------|---------|
| `scripts/lbp_extract_sqlce.py` | Extract the contiguous SqlCE logical image from a snapshot container |
| `scripts/lbp_reader.py` | Decode the `TransactXYData` spectrum blobs directly from a snapshot |

### Extracting the SqlCE image

```bash
# Inspect the container without writing anything
python3 scripts/lbp_extract_sqlce.py LibViewSnapshot_xxx.lbp --info

# Write the full logical image (SqlCE tables + spectrum blobs)
python3 scripts/lbp_extract_sqlce.py LibViewSnapshot_xxx.lbp -o image.bin

# Emit only the spectrum-blob region
python3 scripts/lbp_extract_sqlce.py LibViewSnapshot_xxx.lbp --blobs-only -o blobs.bin

# Physical-order reconstruction (diagnostics / cross-check)
python3 scripts/lbp_extract_sqlce.py LibViewSnapshot_xxx.lbp --physical -o phys.bin
```

The default (logical-sequence) image is the canonical one: it places the SqlCE
schema region at the start and linearises every spectrum blob.

---

## 7. Verification (reference file)

`gus_data/LibViewSnapshot_20260606_141157.lbp`:

| Check | Result |
|-------|--------|
| Magic | `0x4F32AA43` ✓ |
| File size | 1,054,560,256 bytes = 257,461 × 4096 (exact) ✓ |
| Not a ZIP (`PK\x03\x04` absent) | ✓ |
| Unique sequences | 257,116 |
| All 21 expected tables present in extracted image | ✓ |
| Compound CAS+formula records | ~3,080 |
| `TransactXYData` spectrum blobs | 16,075 (16,057 fully parse) |
| Logical image size | 1,049,033,280 bytes |
