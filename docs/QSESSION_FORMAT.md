# QSession File Format - Complete Analysis

## Summary

A `.qsession` file is an **encrypted SQLite database** used by SCIEX MultiQuant/Analytics software to store processed quantitation results, XIC (extracted ion chromatogram) data, and audit trail information.

---

## Encryption

| Property | Value |
|----------|-------|
| **Encryption engine** | **SEE (SQLite Encryption Extension)** — the official commercial encryption from the SQLite team |
| **Cipher** | AES-256-CBC |
| **Password** | `PQS1 is not Sirius` |
| **Page size** | 1024 bytes |
| **SQLite version** | System.Data.SQLite 1.0.98.0 (~2014) |
| **Native DLL** | `SQLite.Interop.dll` (compiled with `INTEROP_INCLUDE_SEE`) |
| **Compatibility** | SEE is the basis for **wxSQLite3**'s encryption — they share the same approach |

### Important: Not Compatible with Open-Source SQLCipher

The encryption uses **SEE (SQLite Encryption Extension)**, which is fundamentally different from SQLCipher:

| Feature | SEE / wxSQLite3 | SQLCipher |
|---------|-----------------|-----------|
| Salt | None | 16 bytes at offset 0 |
| HMAC | None | Appended to each page |
| Key derivation | RC4-based / SHA | PBKDF2-HMAC-SHA1 |
| IV | Page-number derived | Page-number derived |
| Cipher | AES-256-CBC | AES-256-CBC |

The compile-time constant `INTEROP_INCLUDE_SEE` in `SQLiteDefineConstants.cs` confirms this. All attempts with Python's `sqlcipher3` fail because the file format is completely different (SEE has no HMAC, no salt at the beginning).

**To open a qsession file**, you MUST use the vendor `System.Data.SQLite.dll` and `SQLite.Interop.dll`.

### Opening with Python (pythonnet)

```python
import clr
clr.AddReference("System.Data.SQLite")
from System.Data.SQLite import SQLiteConnection, SQLiteConnectionStringBuilder

builder = SQLiteConnectionStringBuilder()
builder.DataSource = "path/to/file.qsession"
builder.Password = "PQS1 is not Sirius"
builder.ReadOnly = True

conn = SQLiteConnection(builder.ConnectionString)
conn.Open()
# ... use conn ...
conn.Close()
```

---

## Database Schema

### 1. `RTParts` — Main Results Table Data

Contains the serialized `MultiData` objects (the core results data).

| Column | Type | Description |
|--------|------|-------------|
| `ATRecordTimeStamp` | datetime | Version timestamp for result snapshot |
| `PartId` | integer | Part ordering index (0-based) |
| `PartContent` | blob | Serialized data chunk |
| `Compressed` | bool | Whether PartContent is GZip-compressed |

- Parts are assembled in `PartId` order and deserialized via `MultiData.DeserializeIteratively()`
- GZip compression is detected via magic bytes `1F 8B` and uses ICSharpCode.SharpZipLib
- The serialization format is a custom binary format (not protobuf, not .NET BinaryFormatter)
- The most recent `ATRecordTimeStamp` contains the latest results

### 2. `XicRawTable` — XIC Chromatogram Data

Contains individual extracted ion chromatograms.

| Column | Type | Description |
|--------|------|-------------|
| `ID` | text | Unique XIC identifier including compound info |
| `SampleKey` | text | Links to sample in wiff2 file |
| `Xdata` | blob | Retention time values (double[], little-endian) |
| `Ydata` | blob | Intensity values (double[], little-endian) |
| `status` | int | Status flag |

- Xdata and Ydata are raw arrays of IEEE 754 double-precision floats
- Each array is the same length (e.g., 270 doubles = 2160 bytes)
- Xdata values are in minutes, monotonically increasing

### 3. `ColumnSettings` — Column Configuration

| Column | Type | Description |
|--------|------|-------------|
| `ColumnSettingsXML` | text | XML configuration for results table columns |
| `Layout` | blob | Serialized layout information |

### 4. `VersionInformation` — File Version

| Column | Type | Description |
|--------|------|-------------|
| `version` | varchar | Version string (e.g., "MultiQuant MD") |

### 5. `QMapVersionInformation` — Quant Map Version

| Column | Type | Description |
|--------|------|-------------|
| `version` | varchar | Version string (e.g., "MultiQuant 2.0") |

### 6. `LockInformation` — File Locking

| Column | Type | Description |
|--------|------|-------------|
| `Locked` | bool | Whether file is locked by another user |

---

## Audit Trail Tables

The qsession file includes a complete audit trail subsystem:

### 7. `AE_AuditEventEntries` — Audit Events (7 rows)
Contains all recorded audit events with full details including timestamps, user info, old/new values, and digital signatures.

### 8. `AE_ValueChangeDetails` — Value Changes (2 rows)
Detailed before/after values for audited changes.

### 9. `AE_AuditMapEntries` — Audit Map Definitions (32 rows)
Defines what events are audited.

### 10. `AE_AuditMapGroups` — Audit Map Groups (1 row)
Groups of audit map entries ("Silent" group).

### 11. `AE_EntryPredefinedReasons` — Entry Predefined Reasons (0 rows)

### 12. `AE_GroupPredefinedReasons` — Group Predefined Reasons (0 rows)

### 13. `AuditEventsVersionTable` — Audit Events Schema Version

### 14. `AuditMapVersionTable` — Audit Map Schema Version

### 15. `AuditTrailRecords` — Legacy Audit Records (0 rows)

### 16. `AuditTrailRecordDescription` — Legacy Audit Descriptions (0 rows)

### 17. `AllSamplesAtTime` — Sample Snapshots (0 rows)

### 18-22. Legacy Audit Map Tables
`AuditMaps`, `ActiveAuditMap`, `ReasonsForMap`, `AuditMapEvents`, `ReasonsPerEvent` — all empty in this file.

---

## Key Source Code References

| File | Purpose |
|------|---------|
| `Sciex.MultiQuant.Data/Utility/DBHelper.cs` | Connection management, password, table creation |
| `Sciex.MultiQuant.Data/SessionPersistance/ResultsTablePersistance.cs` | RTParts read/write, GZip decompression |
| `Sciex.MultiQuant.Data/SessionPersistance/ISessionPersistance.cs` | Session persistence interface |
| `System.Data.SQLite/SQLiteConnection.cs` | Encryption via `sqlite3_key` |
| `System.Data.SQLite/UnsafeNativeMethods.cs` | Native interop to `SQLite.Interop.dll` |

---

## Exporting to Plaintext SQLite

The `scripts/export_qsession_plaintext.py` script exports a qsession file to a standard unencrypted SQLite database using the vendor DLLs:

```bash
python scripts/export_qsession_plaintext.py
```

This creates `data/qsession/<name>_plaintext.db` which can be read with any SQLite tool (Python sqlite3, DB Browser, etc.).

---

## Utility Scripts

| Script | Description |
|--------|-------------|
| `scripts/open_qsession_dotnet.py` | Open and explore qsession using vendor DLL |
| `scripts/export_qsession_plaintext.py` | Export to unencrypted SQLite |
| `scripts/explore_qsession_data.py` | Analyze RTParts and XicRawTable data |
| `scripts/analyze_plaintext.py` | Analyze exported plaintext database |

---

## Example File Statistics

---

## Example File Statistics

**File**: `example_session.qsession`
- **Size**: 121,836,544 bytes (~116 MB)
- **Tables**: 22
- **XIC records**: 3,960 (each with 270 data points)
- **RTParts snapshots**: 2 distinct timestamps, 484 total parts
- **Latest RTParts snapshot**: 20 parts, ~3.6 MB total
- **Audit events**: 7 events (results table creation, sample ID changes, report creation, file save)
- **Samples**: example_sample_N and example_sample_P
