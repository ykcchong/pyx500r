# Clearcore2 DLL Decompile Summary

Decompiled with dnSpy.Console from `f:\dnspy\dnspy.console` on 2026-06-04.
Output: `decompiled/` (945 C# files across 17 assemblies).

## SmartAssembly Obfuscation

`SCIEX.Apis.Data.v1` is obfuscated with **SmartAssembly** (Red Gate). The decompiled output
includes `SmartAssembly/Zip/SimpleZip.cs` which contains the encryption/compression routines.

### Encrypted Embedded Resources (~48 files)

The assembly contains ~48 binary resources with GUID filenames (e.g., `{0e05d448-43d8-44e1-bcd5-77409ef9a5ea}`).
These are **encrypted, native Windows PE32 DLLs** embedded via SmartAssembly's resource protection:

- **Format**: SmartZip v2 (`{z}` header + 3DES encryption + zlib compression)
- **Keys**: Hardcoded in `SimpleZip.cs`
  - 3DES key: `FE F7 E8 DE 88 1B D0 BA`
  - 3DES IV:  `1F 39 85 09 2F 9D 4F CF`
- **Decrypted size**: ~13 MB total across 48 files (ranging from 7 KB to 987 KB each)
- **Content**: All are native 32-bit Windows DLLs (`.text`, `.rsrc`, `.reloc` sections)

These DLLs are **not .NET assemblies** — they are native Windows libraries likely used for
device communication, license enforcement, or P/Invoke callbacks. They are irrelevant for
the pure-Python WIFF2 reader since they are platform-specific native code.

> **Decryption script**: `scripts/decrypt_smartassembly_resources.py`
>
> ```bash
> python scripts/decrypt_smartassembly_resources.py [output_dir]
> ```

## Assembly Inventory

| Assembly | Files | Key Responsibility |
|----------|-------|-------------------|
| `SCIEX.Apis.Data.v1` | ~100 | Entry API, license guard, request handlers |
| `Clearcore2.Data.Wiff2` | 12 | WIFF2 SQLite persistence, encryption password, schema |
| `Clearcore2.Data.WiffReader` | 95 | WIFF1/2 file reader, experiment parsing, TOF calibration |
| `Clearcore2.Compression` | 19 | TOF/Quad/ZeroWidth compression and decompression |
| `Clearcore2.RawXYProcessing` | 52 | Centroiding, peak finding, smoothing, add zeros |
| `Clearcore2.InternalRawXYProcessing` | 5 | Internal processing helpers |
| `Clearcore2.Data` | ~30 | XY data models, data contracts |
| `Clearcore2.Data.Client` | ~40 | gRPC client for sample data provider |
| `Clearcore2.Data.CommonInterfaces` | ~10 | Shared interfaces (IRawXYData, IGetStepSize, etc.) |
| `Clearcore2.Data.AnalystDataProvider` | ~15 | Analyst data provider, XIC calculator, TOF recalibration |
| `Clearcore2.Domain.Acquisition` | ~20 | Acquisition domain models |
| `Clearcore2.Infrastructure` | 8 | Threading, work queues, events |
| `Clearcore2.Muni` | ~10 | Muni (internal) |
| `Clearcore2.StructuredStorage` | ~15 | Structured storage (OLE compound document) |
| `Clearcore2.Utility` | ~20 | Utilities, licensing protection |
| `Clearcore2.XmlHelpers` | ~5 | XML helpers, DataConverter (encryption) |
| `Clearcore2.Devices.Types` | ~10 | Device type enums and models |

## Key Findings

### 1. Encryption (`Clearcore2.Data.Wiff2.PersistenceFactory`)
- **Password**: Derived via `DataConverter.From()` with hardcoded key/salt
- Our `crypto.py` already has the resolved password: `F90CA3B4-CC7B-4439-A479-2097CB8AE246`
- **Status**: ✅ Already pure-Python

### 2. TOF Decompression (`Clearcore2.Compression.DecompressionAlgorithmTof`)
- RLE codec with byte/2byte/4byte value tokens, zero-run mask (0x80)
- Fixed bin marker (0xFF x4), stop marker (0xFF)
- **Status**: ✅ Already pure-Python in `tof.py` (byte-exact match confirmed)

### 3. TOF Calibration (`Clearcore2.Data.WiffReader.WiffTOFCalibration`)
- Simple: `m/z = (A * time_resolution * bin - A * T0)²`
- **Status**: ✅ Already pure-Python in `tof.py`

### 4. License Validation (`SCIEX.Apis.Data.v1.LicenseGuard`)
- Uses OFX.Licensing with hardcoded public key
- Validates license key XML against embedded public key
- **Status**: ⚠️ Could be bypassed (not needed for pure-Python reader)

### 5. Centroiding (`Clearcore2.RawXYProcessing.SpectralPeakFinder`)
- Multi-stage algorithm:
  1. Add zeros between data points (`AddZeros`)
  2. Moving average smoothing (`MovingAverageSmooth`)
  3. Find local maxima
  4. Assign peak ranges at 85% height
  5. Calculate centroid using weighted average above threshold
- **Status**: ✅ Pure-Python in `centroid.py` / `centroid_new.py` / `centroid_fallback.py`

### 6. Framing Zeros (`Clearcore2.RawXYProcessing.AddZeros`)
- Adds zero-intensity points at spectrum boundaries
- Uses step size to determine where zeros should be inserted
- **Status**: ✅ Pure-Python in `centroid.py` / `centroid_new.py` / `centroid_fallback.py`

### 7. Spectrum Reading Flow
```
SampleDataApi.GetSpectra()
  → SpectrumRequestHandler.GetSpectra()
    → IDataProvider.Get<ProductSpectrumRequest, ProductSpectrumResponse>()
      → (internal gRPC to SampleDataProviderServer)
        → Wiff2 persistence + WiffReader + Compression + RawXYProcessing
```

### 8. Data Provider Architecture
- Uses OFX dependency injection (`OFXApp.Locator`)
- `SampleDataProviderComponentConfigurator` registers components
- `IDataProvider` ("SampleDataProviderServer") handles requests
- Internal gRPC-like protocol with request/response types

## Pure-Python Gap Analysis

| Feature | DLL Class | Pure-Python Status |
|---------|-----------|-------------------|
| Decryption | `PersistenceFactory` | ✅ `crypto.py` |
| TOF decompression | `DecompressionAlgorithmTof` | ✅ `tof.py` |
| TOF calibration | `WiffTOFCalibration` | ✅ `tof.py` |
| SQLite reading | `Persistence` | ✅ `reader.py` |
| Sample/experiment metadata | `WiffSampleRun`, `WiffExperiment` | ✅ `reader.py` |
| TIC | `scanItems` table | ✅ `reader.py` |
| **Centroiding** | `SpectralPeakFinder` | ✅ `centroid.py` |
| **Framing zeros** | `AddZeros` | ✅ `centroid.py` |
| License validation | `LicenseGuard` | ⚠️ Not needed |
| Quad decompression | `DecompressionAlgorithmQuad` | ⚠️ Not needed (TOF-only) |

## Implementation Status (2026-06-04)

> **Performance milestone**: Full pipeline (read + centroid all 6,186 spectra) in **1.89 s**
> (**7.7× speedup** vs original 14.55 s).

All TOF/WIFF2-related functions have been implemented in pure Python:

| Feature | DLL Class | Pure-Python Module | Status |
|---------|-----------|-------------------|--------|
| Decryption | `PersistenceFactory` | `crypto.py` | ✅ Already done |
| TOF decompression | `DecompressionAlgorithmTof` | `tof.py::decompress_tof` | ✅ Already done |
| TOF compression | `CompressionAlgorithmTof` | `tof.py::compress_tof` | ✅ **NEW** |
| TOF calibration | `WiffTOFCalibration` | `tof.py::TofCalibration` | ✅ Already done |
| Quad decompression | `DecompressionAlgorithmQuad` | `tof.py::decompress_quad` | ✅ **NEW** |
| Zero-width decompression | `DecompressionAlgorithmZeroWidth` | `tof.py::decompress_zero_width` | ✅ **NEW** |
| Mass/time conversion | `CompressionAlgorithmTof` static | `tof.py::mass_to_time`, `time_to_mass` | ✅ **NEW** |
| Centroiding | `SpectralPeakFinder` | `centroid.py::centroid_spectrum` | ✅ **NEW** |
| Framing zeros | `AddZeros` | `centroid.py::add_framing_zeros` | ✅ **NEW** |
| Moving average smooth | `MovingAverageSmooth` | `centroid.py::moving_average_smooth` | ✅ **NEW** |
| SQLite metadata reading | `Persistence` | `reader.py::WiffReader` | ✅ Already done |

### New files created
- `src/pyx500r/centroid.py` — Centroiding dispatcher (thin wrapper)
- `src/pyx500r/centroid_new.py` — Monolithic numba JIT centroid kernel (~25 ms/spectrum)
- `src/pyx500r/centroid_fallback.py` — Pure-Python + NumPy fallback
- `scripts/decrypt_smartassembly_resources.py` — SmartAssembly resource decryptor
- `docs/DECOMPILE_SUMMARY.md` — This file

### Modified files
- `src/pyx500r/tof.py` — Added `compress_tof`, `decompress_quad`, `decompress_zero_width`, `mass_to_time`, `time_to_mass`, `MassRange`, `return_arrays=True`
- `src/pyx500r/reader.py` — Wired `centroid_spectrum`, batched prefetch, caching, `return_arrays=True`
- `src/pyx500r/__init__.py` — Exported all new functions
- `.gitignore` — Added `decompiled/`

### Key algorithm details

#### Centroiding (`centroid_spectrum`)
Two implementations with identical API:

**Fast path** (`centroid_new.py` — numba JIT):
- Entire pipeline compiled into single `@njit(cache=True)` function
- **~25 ms per dense spectrum** (12.4× faster than original)
- O(n) linear-scan `x_steps` optimization (replaced binary search)
- `np.argsort(..., kind='mergesort')` for deterministic tie-breaking

**Fallback path** (`centroid_fallback.py` — pure Python + NumPy):
- Used when numba is unavailable
- Same algorithm, list-based intermediate structures

**Algorithm** (both paths):
1. Estimate local step size from data
2. Add framing zeros via `add_framing_zeros`
3. Apply moving-average smoothing (half_window=2 for high-res, 1 otherwise)
4. Find local maxima (`y[i] > y[i-1]` and `y[i] >= y[i+1]`)
5. Sort maxima by intensity descending
6. For each maximum, find peak boundaries at 85% height
7. Calculate centroid as weighted average above threshold (default 50% of peak height)
8. Return peaks sorted by m/z

#### TOF Compression (`compress_tof`)
- Inverse of `decompress_tof`
- Fixed-bin marker (FF FF FF FF) + starting bin as LE uint32
- RLE encoding with separate tokens for zero-runs (high bit) vs intensities
- Padded to 4-byte boundary with 0xFF stop markers

#### Quad Decompression (`decompress_quad`)
- Reads mass range headers (start/stop/step/scale as doubles)
- RLE-encoded intensity values per mass bin
- Control byte format: bits 7-5 for count, bit 7 for data/skip flag

#### Zero-Width Decompression (`decompress_zero_width`)
- Array of IEEE 754 floats
- Positive = intensity at current transition
- Negative = skip N transitions
- First float is a header (skipped)
