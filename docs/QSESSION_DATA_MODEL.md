# QSession File Format — Data Model Specification

## 1. File Identity

A `.qsession` file is an **AES-128-OFB encrypted SQLite database** (SEE / SQLite
Encryption Extension, same cipher as SCIEX `.wiff2`).

| Property | Value |
|----------|-------|
| Cipher | AES-128-OFB |
| Page size | 1024 bytes |
| Reserved bytes | 12 (per-page random nonce, in the clear) |
| IV | `struct.pack("<I", page_num) + nonce[:12]` |
| Key | `b"PQS1 is not Siri"` (first 16 bytes of UTF-8 password) |
| Password | `"PQS1 is not Sirius"` |

---

## 2. Database Schema (22 tables)

### 2.1 Primary data tables

#### `RTParts` — Serialised Results Table

The main results container. Stores the complete analysis state as a sequence of
binary chunks that, when concatenated in `PartId` order, form a custom
iterative-serialisation stream (uses `IterativeSerializer` / `IterativeDeserializer`,
**not** the standard .NET `BinaryFormatter`).

| Column | Type | Description |
|--------|------|-------------|
| `ATRecordTimeStamp` | `datetime` | Session-save identifier (`yyyy-MM-dd HH:mm:ss.fff`). Groups all parts from one save operation. |
| `PartId` | `integer` | 0-based ordering index within a timestamp group. |
| `PartContent` | `blob` | Serialised object-graph chunk (custom text-tag protocol). |
| `Compressed` | `bool` | Whether PartContent is GZip-compressed (magic `1F 8B`). |

**Primary key**: `(ATRecordTimeStamp, PartId)` — but no explicit SQL constraint;
enforced by application logic.

**How parts are assembled**:

1. Query most-recent timestamp:
   ```sql
   SELECT ATRecordTimeStamp FROM RTParts ORDER BY ATRecordTimeStamp DESC LIMIT 1
   ```
2. Fetch all parts for that timestamp in order:
   ```sql
   SELECT PartContent FROM RTParts WHERE ATRecordTimeStamp = ? ORDER BY PartId
   ```
3. GZip-decompress individual parts (detected by magic bytes `1F 8B`).
4. Feed the byte-stream into `MultiData.DeserializeIteratively()`.

**Important detail from `ResultsTablePersistance.ReadFromStream()`**:
- Parts whose `Length == 2` are returned as-is (these are version headers).
- GZip detection checks **only the first 2 bytes** (`0x1F 0x8B`).
- Decompression uses `ICSharpCode.SharpZipLib.GZip.GZipInputStream`.
- If the `Compressed` column is `False`, parts are stored raw (as in the sample file).

Older timestamps are retained indefinitely — each save appends a new set of
parts; the loader always picks the latest.

#### `XicRawTable` — Extracted Ion Chromatogram Cache

Caches XIC (chromatogram) data extracted from raw `.wiff` / `.wiff2` files so
that subsequent session loads do not re-extract.

| Column | Type | Description |
|--------|------|-------------|
| `ID` | `text` (PK) | Composite key — see **XicKey format** below. |
| `SampleKey` | `text` | Sample signature (UUID + numeric hash). Groups entries belonging to the same sample. |
| `Xdata` | `blob` | Retention-time axis as raw `double[]` (IEEE 754 LE, 8 bytes/value). |
| `Ydata` | `blob` | Intensity axis as raw `double[]` (same length as Xdata). |
| `status` | `int` | `0` = NoChange, `1` = ToBeAdded, `2` = ToBeRemoved. |

**Table creation SQL** (from `XicFileCache.CreateTable()`):
```sql
CREATE TABLE IF NOT EXISTS XicRawTable (
    ID TEXT, SampleKey TEXT, Xdata BLOB, Ydata BLOB,
    status INT, PRIMARY KEY (ID)
)
```

**XIC ID format** (from `XicKey.ComposeStringKeyInternal()` in
`Clearcore2.QuantLibrary.DataProvider.Cache.XicKey`):

```
[{period}:{expt}]_X {chromType} {sampleSignature} {startMass}-{endMass}
```

Example:
```
[0:0]_X XicScan 402e716a-0374-46fa-8e3f-66e5d30f29d5_5250687681961586044 260.018067747111-260.038067747111
```

Components:
| Part | Meaning |
|------|---------|
| `[0:0]` | `[period:expt]` — period and experiment index |
| `_X` | XIC marker |
| `XicScan` | `ChromType` = `XicScan` (enum value 3) |
| `402e716a-..._5250687681961586044` | `sampleSignature` — UUID + hash |
| `260.018...-260.038...` | m/z extraction window (start–end) |

Other `ChromType` values:
| Enum | Value | ID suffix |
|------|-------|-----------|
| `TicForSample` | 0 | `TIC` |
| `TicForExpt` | 1 | `TIC` |
| `XicIndexed` | 2 | `{extractionIndex}` |
| `XicScan` | 3 | `{startMass}-{endMass}` |
| `XicScanTimeRestricted` | 4 | `{startMass}-{endMass} cycles: {cs}-{ce}` |
| `DadTic` | 5 | — |
| `DadXwc` | 6 | `{startWL}-{endWL}` |
| `AdcXic` | 7 | `{extractionIndex}` |

**Xdata/Ydata encoding** (from `XicFileCache.AddUpdateEntryToDatabase()`):
```csharp
byte[] xBytes = new byte[xyData.GetActualXValues().Length * 8];
Buffer.BlockCopy(xyData.GetActualXValues(), 0, xBytes, 0, xBytes.Length);
byte[] yBytes = new byte[xyData.GetActualYValues().Length * 8];
Buffer.BlockCopy(xyData.GetActualYValues(), 0, yBytes, 0, yBytes.Length);
```
Pure-Python equivalent: `struct.unpack(f'<{len(blob)//8}d', blob)`.

**Status lifecycle** (from `XicFileCache`):
- `NoChange = 0` — existing entry
- `ToBeAdded = 1` — new entry since last save
- `ToBeRemoved = 2` — marked for deletion on next save
- On save finalisation: `DELETE WHERE status=2`; `UPDATE status=0 WHERE status=1`

### 2.2 SQLite Triggers (Referential Integrity)

`DBHelper.CreateTriggersHelper()` creates 7 triggers on the legacy audit tables:

| Trigger | Event | Action |
|---------|-------|--------|
| `DeleteMap_ReasonsForMap` | `DELETE ON AuditMaps` | Cascades delete to `ReasonsForMap` |
| `DeleteMap_AuditMapEvents` | `DELETE ON AuditMaps` | Cascades delete to `AuditMapEvents` |
| `DeleteMap_ReasonsPerEvent` | `DELETE ON AuditMapEvents` | Cascades delete to `ReasonsPerEvent` |
| `DeleteReason_ReasonsPerEvent` | `DELETE ON ReasonsForMap` | Cascades delete to `ReasonsPerEvent` |
| `InsertReason_ReasonsPerEvent` | `BEFORE INSERT ON ReasonsPerEvent` | `RAISE(FAIL)` if reason not in `ReasonsForMap` |
| `UpdateMap_ActiveAuditMap` | `BEFORE UPDATE ON ActiveAuditMap` | `RAISE(FAIL)` if map not in `AuditMaps` |
| `DeleteMap_ActiveAuditMap` | `BEFORE DELETE ON AuditMaps` | `RAISE(FAIL)` if map is currently active |

These triggers enforce referential integrity on the legacy audit-map tables
(`AuditMaps`, `ActiveAuditMap`, `ReasonsForMap`, `AuditMapEvents`, `ReasonsPerEvent`).

---

### 2.3 Configuration tables

#### `ColumnSettings`

| Column | Type | Description |
|--------|------|-------------|
| `ColumnSettingsXML` | `text` (PK) | XML defining column names, widths, visibility, sort order, number format. |
| `Layout` | `blob` | Serialised grid-control layout state. |

This is purely UI state — it preserves the user's results-table view
configuration between sessions. Not part of the analytical data.

#### `VersionInformation`

| Column | Type | Description |
|--------|------|-------------|
| `version` | `varchar` | Version string, e.g. `"MultiQuant MD"`. |

Guards against loading a qsession created by an incompatible software version.

#### `QMapVersionInformation`

| Column | Type | Description |
|--------|------|-------------|
| `version` | `varchar` | Quantitation-map version, e.g. `"MultiQuant 2.0"`. |

Tracks the version of the quantitation method definition.

#### `LockInformation`

| Column | Type | Description |
|--------|------|-------------|
| `Locked` | `bool` | Whether another user/process has the file open for writing. |

Used for network-share multi-user coordination.

**Note**: The table may contain **multiple rows** (the sample file has 2 rows,
both `False`). This occurs during database-merge operations where a `Tmp`
database is attached and its `LockInformation` rows are copied into the target.

### 2.4 Audit trail tables

The audit trail subsystem records every user action for regulated (GLP/GMP)
environments.

| Table | Purpose | Rows in sample |
|-------|---------|----------------|
| `AE_AuditEventEntries` | Recorded audit events (file save, sample rename, report generation, etc.) | 7 |
| `AE_ValueChangeDetails` | Before/after values for audited field changes | 2 |
| `AE_AuditMapEntries` | Which event types must be audited | 32 |
| `AE_AuditMapGroups` | Groups of audit map entries (e.g. "Silent" group) | 1 |
| `AE_EntryPredefinedReasons` | Predefined reason codes for audit entries | 0 |
| `AE_GroupPredefinedReasons` | Predefined reason codes for audit groups | 0 |
| `AuditEventsVersionTable` | Audit events schema version (`201711180002`) | 1 |
| `AuditMapVersionTable` | Audit map schema version (`201803190003`) | 1 |

**Legacy audit tables** (empty in current format):
`AuditTrailRecords`, `AuditTrailRecordDescription`, `AllSamplesAtTime`,
`AuditMaps`, `ActiveAuditMap`, `ReasonsForMap`, `AuditMapEvents`,
`ReasonsPerEvent`.

---

## 3. The MultiData Object Graph

`MultiData` is the root object deserialised from `RTParts.PartContent`. It is
a `[Serializable]` .NET class serialised with `BinaryFormatter`.

### 3.1 Object hierarchy

```
MultiData                              ← root container
├── QuantMethod                        ← compound definitions & method settings
│   ├── List<QuantCompound>            ← all compounds (analytes + internal standards)
│   │   ├── QuantAnalyte               ← analyte subclass
│   │   └── QuantInternalStd           ← internal-standard subclass
│   ├── List<FormulaInfo>              ← user-defined formulas (calculated columns)
│   └── IntegrationParameters          ← default peak-finding algorithm
│
├── List<MultiSample>                  ← samples, indexed by position
│   ├── SampleLocator                  ← link back to raw data file
│   ├── SampleName, SampleID, Barcode  ← user identifiers
│   ├── DateTime                       ← acquisition timestamp
│   ├── DilutionFactor, InjectionVolume
│   ├── InstrumentName, SerialNumber
│   ├── SampleSignature                ← unique hash → links to XicRawTable.SampleKey
│   └── MultiPeak[]                    ← one peak per compound
│       ├── CompoundIndex (PeakIndex)  ← index into QuantMethod.Compounds
│       ├── RetentionTime, Area, Height
│       ├── ApexRT, ApexY
│       ├── CorrectedArea, CorrectedHeight
│       ├── Noise, SignalToNoise
│       ├── ActualConcentration
│       ├── Use (bool)                 ← included in calibration?
│       ├── Profile (float[] + compression_type)
│       ├── ValidIntegration (bool)
│       └── Reportable (bool)
│
├── CalibrationData                    ← external calibration curves (optional)
├── IsotopePatterns                    ← isotope correction data
├── CustomFields                       ← user-defined column metadata
├── RenamedColumns                     ← column-rename tracking
└── Comment (string)                   ← session-level annotation
```

### 3.2 Key relationships

```
MultiSample[idx].MultiPeak[compound_idx]
    │                    │
    │                    └── QuantMethod.Compounds[compound_idx]
    │                         ├── Name, ExtractionValues1/2, Period, Experiment
    │                         └── AssociatedCompounds (if summed)
    │
    └── SampleSignature ──── XicRawTable.SampleKey
                                  │
                                  └── XicRawTable.ID contains compound index
```

The results table displayed to the user is a grid of:
- **Rows** = `MultiSample[]` (each sample)
- **Columns** = `QuantCompound[]` (each compound's metrics) + user formulas + fixed columns
- **Cells** = `MultiPeak` values (area, height, concentration, etc.)

### 3.3 Versioning

`MultiData` uses version byte **19** (current). The deserialiser reads this
first, then conditionally loads fields based on version:

| Version | Added |
|---------|-------|
| 1–9 | Basic samples, custom fields |
| 10+ | Quantitation flags, formulas |
| 14+ | External calibration, isotope patterns |
| 16+ | Custom field formulas |
| 18+ | Combined flagging rules |
| 19 | All current features |

`QuantMethod` has its own independent version (32).

### 3.4 SampleKey structure

**`SampleKey`** (from `Clearcore2.QuantLibrary.DataProvider.Cache.SampleKey`):

```csharp
public SampleKey(string path, int sampleIndex, string implTag) {
    _sampleSignature = string.Format("{0}:{1}", path, sampleIndex);
    _key = string.Format("{0} {1}", _sampleSignature, implTag);
}
```

| Component | Description |
|-----------|-------------|
| `path` | Raw data file path (e.g., `C:\Data\file.wiff2`) |
| `sampleIndex` | Sample index within the acquisition |
| `implTag` | Factory tag: `"PF"` (PersistenceFactory), `"RF"` (WiffReader), `"RL"` (RFLight) |
| `_sampleSignature` | `path:index` — base signature |
| `_key` | `path:index implTag` — full key including factory |

In the qsession file, `XicRawTable.SampleKey` stores the *resolved* sample
signature, which is a UUID + numeric hash derived from the base signature
(e.g., `402e716a-0374-46fa-8e3f-66e5d30f29d5_5250687681961586044`). This
hashing is done by a `SampleSignatureProvider` (interface implementation).

---

## 4. Serialisation Format

### 4.1 PartContent encoding

Each `PartContent` blob is a byte-stream produced by .NET `BinaryFormatter`
with a custom `SerializationBinder`. The stream is optionally GZip-compressed
(detected by magic bytes `1F 8B`, decompressed with `ICSharpCode.SharpZipLib`).

The custom `IterativeSerializer` writes objects with **named fields** using a
text-tag protocol (not the standard .NET BinaryFormatter binary format):

```
<field_name_as_utf8_string><type_tag_byte><value_bytes>
```

Type tags observed in the protocol:
| Tag | C# method | Type |
|-----|-----------|------|
| `i` | `SerializeInt` / `DeserializeInt` | `int32` |
| `I` | `SerializeNullableInt` / `DeserializeNullableInt` | `int32?` |
| `b` | `SerializeBool` / `DeserializeBool` | `bool` |
| `d` | `SerializeDouble` / `DeserializeDouble` | `double` |
| `s` | `SerializeString` / `DeserializeString` | `string` |
| `o` | `SerializeOneObject` / `DeserializeOneObject` | complex object |
| `O` | `SerializeCollection` / `DeserializeCollection` | `IEnumerable<T>` |
| `L` | `SerializeList` / `DeserializeList` | `List<T>` |
| `D` | `SerializeDateTime` / `DeserializeDateTime` | `DateTime` |

**`IterativeDeserializer` mechanics**:
- Maintains a `MemoryStream` (initial capacity 203 KB = `203000` bytes) fed
  from the part enumerator.
- Reuses a `BinaryFormatter`-backed `ObjectManager` for object-graph
  reconstruction (handles cross-references between `MultiPeak` and
  `MultiSample`).
- `DeserializeOneObject(name)` reads a tag, then the object header, then
  recurses into the object's fields.
- `DeserializeEnumerables<T>()` reads a count prefix, then N items.
- `DeserializeAllAndHookUp()` finalises the object graph after all parts
  are loaded — this restores:
  - `MultiPeak.fMultiSample` references
  - `MultiSample.fMultiData` reference
  - `MultiPeak.fQuantMethod` and `QuantMethod` → compound associations

### 4.2 Serialisation order (save path)

```
MultiData.SerializeIteratively()
  → writes version byte (19)
  → QuantMethod.SerializeIteratively()
      → version (32)
      → List<QuantCompound> (each serialised individually)
      → method parameters, formulas, plugins
  → List<MultiSample> (each serialised individually)
      → sample metadata
      → List<MultiPeak> (each serialised individually)
          → peak metrics, integration params, profile
  → custom fields
  → calibration data
  → isotope patterns
  → comment
```

During save, the byte-stream is split into chunks (`RTPart` objects). Each
chunk may be GZip-compressed in parallel (`Environment.ProcessorCount * 2`
tasks). Chunks are written to `RTParts` table with sequential `PartId`.

### 4.3 Deserialisation order (load path)

Reverse of save: reads version → QuantMethod → MultiSamples → MultiPeaks →
custom fields → calibration → isotope patterns → comment. Cross-references
(MultiPeak → QuantCompound, MultiSample → MultiData) are restored after all
objects are deserialised.

---

## 5. Data Flow: Opening a QSession

```
User opens file.qsession
        │
        ▼
DBHelper.CreateAndOpenConnection(path, journalModeOff=true, readOnly=true)
        │  password = "PQS1 is not Sirius"
        │  → System.Data.SQLite handles SEE decryption
        ▼
ResultsTablePersistance.LoadResultsTablePersistance(path)
        │
        ├── SELECT latest ATRecordTimeStamp FROM RTParts
        ├── SELECT PartContent FROM RTParts WHERE ATRecordTimeStamp = ? ORDER BY PartId
        ├── GZip-decompress individual parts (if magic 1F 8B)
        ▼
MultiData.DeserializeIteratively(byteEnumerable)
        │  Uses IterativeDeserializer + BinaryFormatter
        │  Reconstructs full object graph
        ▼
MultiData object returned to application
        │
        ├── QuantMethod → compound list for column headers
        ├── MultiSample[] → sample list for row headers
        ├── MultiPeak[] → cell values (area, conc, etc.)
        │
        └── (lazy) For each MultiPeak:
                construct XicRawTable ID key
                → SELECT Xdata, Ydata FROM XicRawTable WHERE ID = ?
                → deserialise double[] from blob
                → display chromatogram
```

---

## 6. Data Flow: Saving a QSession

```
User clicks Save (or auto-save triggers)
        │
        ▼
ResultsTablePersistance.SaveResultsTablePersistance(path, multiData, time, conn, auto)
        │
        ├── GenerateSessionId (timestamp string)
        ├── multiData.SerializeIteratively()
        │       → yields IEnumerable<byte[]> (one chunk per RTPart)
        │
        ├── For each chunk:
        │       optionally GZip-compress (parallel, N = CPU*2 threads)
        │       INSERT INTO RTParts (ATRecordTimeStamp, PartId, PartContent, Compressed)
        │
        ├── DBHelper.CreateVersionTable (if not exists)
        ├── DBHelper.OptimizePartsIfPossible (remove old auto-saved parts)
        │
        └── (separately, by caller)
                XicFileCache writes/updates XicRawTable entries
                ColumnSettings XML/Layout saved
                AuditEventEntries appended
```

---

## 7. Reading QSession Without Vendor DLLs

Since the encryption uses the same SEE/AES-128-OFB scheme as wiff2 files,
`pyx500r.crypto` can decrypt qsession files directly:

```python
from pyx500r.crypto import decrypt_database, QSESSION_PASSWORD
import sqlite3

plaintext = decrypt_database(
    "file.qsession",
    password=QSESSION_PASSWORD,   # "PQS1 is not Sirius"
    page_size=1024,               # wiff2 uses 4096
)
db = sqlite3.connect(":memory:")
db.deserialize(plaintext)
```

## 8. Binary Serialisation — Notes for Pure-Python Deserialisation

The `PartContent` format is .NET `BinaryFormatter`. To deserialise without .NET:

- **Required**: A parser for the custom text-tag protocol used by
  `IterativeSerializer` / `IterativeDeserializer`. This is a simplified
  .NET serialization that uses field-name tags rather than the full
  `BinaryFormatter` binary wire format.
- **Alternative**: Use Python.NET (`pythonnet`) to load the vendor DLLs and call
  `MultiData.DeserializeIteratively()` directly — this is the approach used in
  `scripts/open_qsession_dotnet.py`.
- **Pragmatic**: Extract XIC data directly from `XicRawTable` (raw `double[]`
  blobs — trivially parsed). The `PartContent` blobs can be inspected for
  specific fields without a full protocol parser (e.g., searching for UTF-8
  field-name strings like `"fQuantMethod"`, `"fSamples"`, etc.).

The `XicRawTable` Xdata/Ydata format is trivial:
```python
import struct
num_points = len(blob) // 8
xdata = struct.unpack(f'<{num_points}d', blob)
```
