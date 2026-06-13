# SCIEX LibraryView `.lbp` Package File Format

> Reverse-engineered from `monodis` disassembly of:
> - `Clearcore2.LibraryView.Packager.dll` (primary packager logic)
> - `Clearcore2.LibraryView.Common.dll` (crypto, data types)
> - `Clearcore2.Library.Core.dll` (domain model, DB schema)
> - `CliquidPackage.dll` (low-level OPC container engine)
> - `PackageUpgrader.Sciex.LibraryView.PackageUpgrader.exe` (upgrade tool)
> - `LibViewMapper.xml` (table/column mapping schema)

> **⚠ Two distinct `.lbp` layouts exist — check the magic bytes first.**
> - **OPC/ZIP package** (first bytes `50 4B 03 04` = `PK\x03\x04`): produced by
>   `CliquidPackage.Pack`. This is the format described **in this document**.
> - **Snapshot container** (first bytes `43 AA 32 4F` = magic `0x4F32AA43`):
>   a SCIEX paged memory-mapped store wrapping a SqlCE 3.5 database, used by
>   `LibViewSnapshot_*.lbp`. See **`docs/LBP_SNAPSHOT_FORMAT.md`** and
>   `scripts/lbp_extract_sqlce.py` / `scripts/lbp_reader.py`.
>
> Both layouts ultimately carry the same SqlCE schema and `TransactXYData`
> spectrum blobs; they differ only in how that database is wrapped.

---

## 1. Container Format: OPC (Open Packaging Convention)

A `.lbp` file is an **OPC package** — the same format used by `.docx`, `.xlsx`, etc. It is a ZIP archive with a specific internal structure. The underlying library used is `System.IO.Packaging.Package` from `WindowsBase.dll`.

The `CliquidPackage` library (namespace `Cliquid.Packaging`) wraps `System.IO.Packaging` and handles all pack/unpack operations:

```
CliquidPackage.Pack(files, baseDirectory, outputFile, header)
  → System.IO.Packaging.Package.Open(outputFile, FileMode.Create, FileAccess.ReadWrite)
  → For each file: Package.CreatePart(uri, "application/octet-stream")
  → PackagePart.GetStream() ← FileStream (binary copy)
  → Package.PackageProperties ← metadata header fields
  → Package.Dispose()
```

**File extension:** `.lbp` (LibraryView Package). The exporter checks the extension and routes to `ExportCompoundsToLibView` only when the extension is `lbp`. Other supported formats are `mdb` (Analyst) and `xls` (Excel).

---

## 2. OPC Package Metadata (PackageProperties)

The header metadata is stored in OPC `PackageProperties` (Dublin Core properties, stored in `/_rels/.rels` and `/docProps/core.xml` inside the ZIP). The `PackageHeader` class maps to these fields:

| `PackageHeader` Field | OPC `PackageProperties` Property | Description |
|---|---|---|
| `Title` | `Title` | Package display name |
| `ContentTable` | `Subject` | Semicolon-delimited list of content items (e.g. `"Directory,<path>;"`) |
| `Description` | `Description` | Free-text description |
| `ContentSize` | `ContentType` | Total size of packaged content in bytes (stored as string) |
| `Keywords` | `Keywords` | Free-text keywords |
| `Category` | `Category` | Instrument metadata string, e.g. `"2;MASS_SPECTROMETER=ModelName;PUMP=;AUTOSAMPLER=Agilent 1100 G1367;"` |
| `Language` | `Language` | Language tag, default `"en"` |
| `Revision` | `Revision` | Revision string |
| `ContentStatus` | `ContentStatus` | Status string (e.g. `"Final"`) |
| `LastModifiedBy` | `LastModifiedBy` | `"DOMAIN/username"` of the packing user |
| `Creator` | `Creator` | Machine/hardware ID (from `Utility.get_CreatorID()`, derived from MAC address) |
| `Version` | `Version` | Assembly version of packing DLL |
| `Identifier` | `Identifier` | New `Guid.NewGuid().ToString()` generated at pack time |
| `Created` | `Created` | UTC datetime of packaging |
| `Modified` | `Modified` | Set equal to `Created` at pack time |

### Category Field Format

The `Category` field encodes instrument information as a semicolon-delimited key=value string:

```
"<count>;MASS_SPECTROMETER=<model>;PUMP=<pump>;AUTOSAMPLER=<sampler>;"
```

On import, the packager parses `MASS_SPECTROMETER=` from this string to identify the instrument model.

### ContentTable Field Format

The `ContentTable` field (stored in `Subject`) is a semicolon-delimited list describing the package contents:

```
"Directory,<base_path>;Test,Test name;Report Style,Name;"
```

Each entry is `"<type>,<value>;"`.

---

## 3. Internal Package Parts (OPC Parts)

Each file packed into the LBP becomes an OPC Part with:
- **URI**: relative path derived from the file path with base directory stripped and spaces URL-encoded (`%20`)
- **Content-Type**: `application/octet-stream`

### Per-file Processing Pipeline

**On packing (export):**
1. If `DoCompression = true` (default): file is GZip-compressed → saved as `<filename>.Compressed`
2. The compressed (or plain) file is renamed to `<filename>.Part` (suffix appended, not replaced)
3. File is added as OPC Part at URI derived from path relative to `baseDirectory`
4. Original `.Part` file is deleted from disk

**On unpacking (import):**
1. Each OPC Part is extracted to `outputDirectory` using URI path (leading `/` stripped, `%20` → space)
2. If filename ends with `.Part`, the suffix is stripped
3. If the resulting name ends with `.Compressed`, GZip-decompress it and delete the `.Compressed` file
4. Final extracted file has the original name

**Settings defaults** (from `Cliquid.Packaging.Settings`):

| Setting | Default | Description |
|---|---|---|
| `PackagePartID` | `".Part"` | Suffix appended to files before adding to package |
| `CompressedID` | `".Compressed"` | Suffix indicating GZip-compressed content |
| `DoCompression` | `true` | Compress files before packing |
| `CleanUp` | `false` | Delete source files/dirs after packing |

---

## 4. Primary Content: SQL Server Compact (SqlCE) Database

The core content of an LBP is a **SQL Server Compact Edition (`.sdf`) database file**. This is what gets packed (and optionally GZip-compressed) into the OPC container.

### SqlCE Connection String Format

```
Data Source='{filename}';{optional_auth}LCID={lcid};Password={password};Encrypt={true|false};Max Database Size=4000; Persist Security Info=False;
```

- Unencrypted: `Encrypt=false`, `Password=` (empty)
- Encrypted: `Encrypt=true`, `Password=<user_password>`
- Max database size: 4000 MB

### Database Schema (complete)

The export creates a fresh SqlCE database with the following schema:

#### Core Tables

```sql
CREATE TABLE Compound (
    [Id]                    UniqueIdentifier NOT NULL DEFAULT (newid()),
    [Identifier]            NVarChar(100) NOT NULL,
    [CAS]                   NVarChar(100) NOT NULL,
    [Formula]               NVarChar(200) NOT NULL,
    [MolecularWeight]       Float NOT NULL,
    [MonoIsotopicMass]      Float,
    [MolecularStructureSource] NText NOT NULL,
    [PurityThreshold]       Float,
    [RedFlagThreshold]      Float,
    [YellowFlagThreshold]   Float,
    [Comment]               NVarChar(2000) NOT NULL,
    [CreatedDate]           DateTime NOT NULL,
    [LastUpdated]           DateTime,
    [Active]                Bit NOT NULL
)

CREATE TABLE CompoundName (
    [Id]          UniqueIdentifier NOT NULL DEFAULT (newid()),
    [CompoundId]  UniqueIdentifier NOT NULL,
    [Name]        NVarChar(200) NOT NULL,
    [RegionId]    UniqueIdentifier NOT NULL,
    [LastUpdated] DateTime,
    [IsDefault]   Bit NOT NULL
)

CREATE TABLE Region (
    [Id]   UniqueIdentifier NOT NULL,
    [Name] NVarChar(500) NOT NULL
)

CREATE TABLE Class (
    [Id]          UniqueIdentifier NOT NULL,
    [Name]        NVarChar(2000) NOT NULL,
    [LastUpdated] DateTime
)

CREATE TABLE Library (
    [Id]             UniqueIdentifier NOT NULL,
    [Name]           NVarChar(450) NOT NULL,
    [ParentFolderId] UniqueIdentifier,
    [LastUpdated]    DateTime,
    [IsFolder]       Bit NOT NULL DEFAULT ((1))
)
```

#### Junction Tables

```sql
CREATE TABLE CompoundClass (
    [ClassId]    UniqueIdentifier NOT NULL,
    [CompoundId] UniqueIdentifier NOT NULL
)

CREATE TABLE CompoundLibrary (
    [CompoundId] UniqueIdentifier NOT NULL,
    [LibraryId]  UniqueIdentifier NOT NULL
)

CREATE TABLE CompoundFavorite (
    [Username]   NVarChar(100) NOT NULL,
    [CompoundId] UniqueIdentifier NOT NULL
)
```

#### Instrument / Scan Type Tables

```sql
CREATE TABLE InstrumentType (
    [ModelName]     NVarChar(100) NOT NULL,
    [InstrumentType] NChar(10),
    [InstrumentKey] Int,
    [Id]            UniqueIdentifier NOT NULL
)

CREATE TABLE Instrument (
    [Id]               UniqueIdentifier NOT NULL,
    [InstrumentTypeId] UniqueIdentifier NOT NULL,
    [Description]      NVarChar(100) NOT NULL
)

CREATE TABLE ScanType (
    [Id]          UniqueIdentifier NOT NULL,
    [ScanTypeKey] Int NOT NULL,
    [ScanTypeName] NVarChar(100) NOT NULL
)
```

#### Spectral Data Tables

```sql
CREATE TABLE MassSpectrum (
    [Id]                    UniqueIdentifier NOT NULL,
    [CompoundId]            UniqueIdentifier NOT NULL,
    [InstrumentId]          UniqueIdentifier,
    [ScanTypeId]            UniqueIdentifier,
    [RawXYData]             Image,               -- GZip-compressed raw spectrum (nullable)
    [CentroidedXYData]      Image NOT NULL,      -- GZip-compressed centroided spectrum
    [WiffFile]              NVarChar(255) NOT NULL,
    [SampleIndex]           Int NOT NULL,
    [PeriodIndex]           Int NOT NULL,
    [ExperimentIndex]       Int NOT NULL,
    [StartRT]               Float NOT NULL,
    [EndRT]                 Float NOT NULL,
    [BackgroundSubtraction] Bit NOT NULL,
    [CreatedDate]           DateTime NOT NULL,
    [LastUpdated]           DateTime,
    [PrecursorMass1]        Float NOT NULL,
    [PrecursorMass2]        Float,
    [PrecursorChargeState1] Float,
    [PrecursorChargeState2] Float,
    [PositivePolarity]      Bit NOT NULL,
    [CollisionEnergy]       Float NOT NULL,
    [CollisionEnergy2]      Float,
    [CollisionEnergySpread] Float,
    [CollisionEnergySpread2] Float,
    [Type]                  NVarChar(30) NOT NULL,
    [Operator]              NVarChar(50),
    [CompoundPurityThreshold] Float,
    [IsolationWindow1]      Numeric(3,1),
    [SpectrumResolution]    Float,
    [SpectrumPeakWidth]     Float,
    [RawDataStepSize]       Float,
    [RawDataACalibration]   Float,
    [RawDataT0Calibration]  Float,
    [AdditionalData]        NText,
    [CADGasValue]           Float,
    [CADGasType]            NVarChar(50),
    [IonSource]             NVarChar(50),
    [Encryption]            NVarChar(15) DEFAULT (''),
    [Hash]                  NVarChar(32)
)

CREATE TABLE MassSpectrumSignature (
    [MassSpectrumID] UniqueIdentifier NOT NULL,
    [Mass0]  Float NOT NULL,  [Mass1]  Float NOT NULL,
    [Mass2]  Float NOT NULL,  [Mass3]  Float NOT NULL,
    [Mass4]  Float NOT NULL,  [Mass5]  Float NOT NULL,
    [Mass6]  Float NOT NULL,  [Mass7]  Float NOT NULL,
    [Mass8]  Float NOT NULL,  [Mass9]  Float NOT NULL,
    [Mass10] Float NOT NULL,  [Mass11] Float NOT NULL,
    [Mass12] Float NOT NULL,  [Mass13] Float NOT NULL,
    [Mass14] Float NOT NULL,  [Mass15] Float NOT NULL
)

CREATE TABLE UVSpectrum (
    [Id]              UniqueIdentifier NOT NULL,
    [CompoundId]      UniqueIdentifier NOT NULL,
    [InstrumentId]    UniqueIdentifier,
    [RawXYData]       Image NOT NULL,
    [ProcessedXYData] Image,
    [WiffFile]        NVarChar(255) NOT NULL,
    [SampleIndex]     Int NOT NULL,
    [StartRT]         Float NOT NULL,
    [EndRT]           Float NOT NULL,
    [DetectorType]    NVarChar(100) NOT NULL,
    [CreatedDate]     DateTime NOT NULL,
    [LastUpdated]     DateTime
)
```

#### Transition / Retention Time Tables

```sql
CREATE TABLE Transition (
    [Id]             UniqueIdentifier NOT NULL,
    [CompoundId]     UniqueIdentifier NOT NULL,
    [InstrumentId]   UniqueIdentifier NOT NULL,
    [Q1]             Float NOT NULL,
    [Q3]             Float NOT NULL,
    [DP]             Float,
    [EP]             Float,
    [FP]             Float,
    [CE]             Float,
    [CXP]            Float,
    [Polarity]       Bit NOT NULL,
    [TransitionOrder] Int NOT NULL,
    [LastUpdated]    DateTime,
    [Encryption]     NVarChar(15) DEFAULT (''),
    [Hash]           NVarChar(32)
)

CREATE TABLE RetentionTime (
    [Id]          UniqueIdentifier NOT NULL,
    [CompoundId]  UniqueIdentifier NOT NULL,
    [InstrumentId] UniqueIdentifier NOT NULL,
    [Value]       Float NOT NULL,
    [LastUpdated] DateTime,
    [Encryption]  NVarChar(15) DEFAULT (''),
    [Hash]        NVarChar(32),
    [IsDefault]   Bit NOT NULL DEFAULT ((0))
)
```

#### Custom Attribute Tables

```sql
CREATE TABLE CustomAttributeDefinition (
    [Id]           UniqueIdentifier NOT NULL DEFAULT (newid()),
    [Name]         NVarChar(50) NOT NULL,
    [IsRequired]   Bit NOT NULL,
    [MinimumValue] Float,
    [MaximumValue] Float,
    [MaxLength]    Int,
    [LastUpdated]  DateTime,
    [Type]         Int NOT NULL DEFAULT ((1))   -- 1=Number, 2=Formula/Text
)

CREATE TABLE CustomAttributeValue (
    [Id]           UniqueIdentifier NOT NULL DEFAULT (newid()),
    [CompoundId]   UniqueIdentifier NOT NULL,
    [DefinitionId] UniqueIdentifier NOT NULL
)

CREATE TABLE NumberAttributeValue (
    [Id]    UniqueIdentifier NOT NULL,
    [Value] Float
)

CREATE TABLE TextAttributeValue (
    [Id]    UniqueIdentifier NOT NULL,
    [Value] NVarChar(500)
)

CREATE TABLE BinaryAttributeValue (
    [Id]    UniqueIdentifier NOT NULL DEFAULT (newid()),
    [Value] Image
)
```

### Indexes

```sql
CREATE INDEX IX_Class_Name ON Class (Name)
CREATE INDEX IX_Compound_Active ON Compound (Active)
CREATE INDEX IX_Compound_CAS ON Compound (CAS)
CREATE INDEX IX_Compound_Formula ON Compound (Formula)
CREATE INDEX IX_Compound_MolecularWeight ON Compound (MolecularWeight)
CREATE INDEX IX_CompoundName_Name_RegionId ON CompoundName (Name, RegionId)
CREATE INDEX IX_MassSpectrum_PrecursorMass1 ON MassSpectrum
    (PrecursorMass1, Id, PrecursorMass2, CollisionEnergy, CollisionEnergySpread, Encryption)
CREATE INDEX ddd ON MassSpectrumSignature
    (Mass0, Mass1, Mass2, Mass3, Mass4, Mass5, Mass6, Mass7, Mass8,
     Mass9, Mass10, Mass11, Mass12, Mass13, Mass14, Mass15)
CREATE INDEX IX_RetentionTime_Value ON RetentionTime (Value, CompoundId)
```

---

## 5. Spectral XY Data Encoding

The `RawXYData` and `CentroidedXYData` columns in `MassSpectrum`, and `RawXYData`/`ProcessedXYData` in `UVSpectrum`, store spectral data as **GZip-compressed binary blobs** (`Image` type = `byte[]`).

The binary format inside the GZip stream corresponds to a .NET BinaryFormatter-serialized `TransactXYData` object (from `Clearcore2.LibraryView.Common`), containing two `double[]` arrays:
- Array 1: m/z values (or wavelength for UV)
- Array 2: intensity values

Both arrays have the same length. Typical MS/MS spectrum: 50–500 points.

---

## 6. Encryption

### Database-level Encryption (SqlCE password)

The `.sdf` database can be password-protected using SqlCE's built-in encryption. The connection string uses `Encrypt=true` and `Password=<user_password>`.

### Column-level Encryption (RC2)

Individual data columns in `MassSpectrum`, `Transition`, and `RetentionTime` can be encrypted column-by-column. The `Encryption` column (NVarChar(15)) records the encryption state; `Hash` (NVarChar(32)) stores an integrity hash.

The `Clearcore2.LibraryView.Common.Crypto` class implements this using **RC2** symmetric encryption (`RC2CryptoServiceProvider`):

| Parameter | Value |
|---|---|
| Algorithm | RC2 |
| Key | First 16 bytes of `Encoding.ASCII.GetBytes(license)`, cyclically repeated/truncated |
| IV | First 8 bytes of `Encoding.ASCII.GetBytes(license)`, cyclically repeated/truncated |
| Key derivation | Raw ASCII bytes of the license string, no PBKDF, no salt |
| Block mode | Default (CBC) |

The default hardcoded `CodeString` for `EncryptString`/`DecryptString` is `"LibView"`.

For column encryption via `EncryptPackage`, the license key passed is the user-supplied package password (or the system license string).

---

## 7. Additional Package Files

Besides the `.sdf` database, the LBP may also contain mapper XML files. These are the same files used to configure import/export column mappings. They are copied from the `Mappers/` directory relative to the packager DLL and packed alongside the database.

**Key mapper files:**
- `LibViewMapper.xml` — maps between LibraryView DB schema and the LBP SqlCE schema
- `LibViewExportMapper.xml` — used specifically during export
- `CliquidMapper.xml` — used for Cliquid `.clq` format source files

The mapper XML schema (`Mapper` element) describes:
- `<Connection>` — source connection string template
- `<Table>` — one per DB table, with `TargetTableName`, `SourceTableName`, optional `SourceCommandText`, `SourceJoinClause`
- `<Column>` — per column: `TargetColumnName`, `SourceColumnName`, `Type`, `IsPrimaryKey`, `IsFixed`, `Value`, `IsRequired`, `IsUsedForHash`, `ImportEnabled`, `DuplicationCheckColumnName`, `SourceForeignKeyColumnName`

---

## 8. Export / Import API

### ExportParameters

```csharp
class ExportParameters {
    string OutputFileName;        // path ending in .lbp, .mdb, or .xls
    List<CompoundIdentification> CompoundIds;  // compounds to export (null = all)
    string Password;              // optional SqlCE encryption password
    bool   IsEncrypted;           // enable SqlCE encryption
    string InstrumentName;        // filter by instrument
    DataTable ExportTable;        // optional custom column mapping table
    bool   IsSnapshot;            // snapshot mode (export all)
}
```

### ImportParameters

```csharp
class ImportParameters {
    string InputFileName;   // path to .lbp file
    string Password;        // decryption password (if encrypted)
    // ... additional filter/target fields
}
```

### Export Flow

1. `ExportCompounds(ExportParameters)` checks file extension → routes to `ExportCompoundsToLibView`
2. `Engine.CreateStructureCe(connectionString, fileName)` → creates empty `.sdf` with full schema
3. All selected compound data is exported from the source (SQL Server or MongoDB) into the `.sdf`
4. `DatabaseUtility.EncryptPackage(...)` encrypts columns if requested
5. `CliquidPackage.Instance.Pack(files, baseDir, outputFile, header)`:
   - Opens OPC package at `outputFile`
   - GZip-compresses the `.sdf` → adds as OPC Part
   - Copies mapper XML files → adds as OPC Parts
   - Writes `PackageHeader` into `PackageProperties`
   - Closes/disposes the OPC package

### Import Flow

1. `Engine.UnpackItemsToTemp(compoundId)`:
   - `CliquidPackage.Instance.Unpack(inputFile, filter, tempFolder)`
   - Opens OPC package → extracts all parts → GZip-decompresses `.Compressed` files
2. `CliquidPackage.ReadHeader(inputFile)` → reads `PackageProperties` back into `PackageHeader`
3. Instrument name is parsed from `PackageHeader.Category` (`MASS_SPECTROMETER=...`)
4. `ContentTable` (`PackageHeader.ContentTable` = OPC `Subject`) tells the importer what content to expect
5. `Engine.GetCliquidCompoundSource` reads compounds from the extracted `.sdf`
6. Data is imported/merged into the target LibraryView database

---

## 9. Cliquid (.clq) Source Files

The `CliquidExtension = "clq"` constant and `EngineTypes.Cliquid` / `SourceTypes.Cliquid` indicate that `.clq` files (Cliquid software packages) are also supported as an import source. They use the same OPC container format and are unpacked with the same `CliquidPackage` engine, but the XML content inside follows the `Cliquid.Packaging` schema:

| XML File | Class | Content |
|---|---|---|
| `Compounds.ID.xml` | `CompoundsXml` | List of `CompoundXml` with compound IDs and transitions |
| `SampleColumns.xml` | `SampleColumnsXml` | Column definitions |
| `Test.xml` | `TestXml` | Test/method definitions |
| `ReportStyle.xml` | `ReportStyleXml` | Report template |
| `TestWorkflow.xml` | `WorkflowXml` | Workflow definition |

---

## 10. Physical File Summary

When you open a `.lbp` file as a ZIP archive, you will find:

```
[content_types].xml          ← OPC content types
_rels/.rels                  ← OPC relationships
docProps/core.xml            ← PackageHeader metadata (Dublin Core)
<relative_path>.sdf.Compressed.Part   ← GZip-compressed SqlCE database
<relative_path>/LibViewExportMapper.xml.Part  ← (optional) mapper config
...
```

The `.sdf` itself is a SQL Server Compact Edition 3.5 database (JET-based, proprietary binary format). To read it outside of the packager, extract and decompress the `.sdf.Compressed.Part` file, then open with any SqlCE 3.5-compatible driver.

---

## 11. Key Assembly Reference

| Assembly | Role |
|---|---|
| `CliquidPackage.dll` | OPC container pack/unpack, GZip compression |
| `Clearcore2.LibraryView.Packager.dll` | Export/import orchestration, SqlCE schema creation |
| `Clearcore2.LibraryView.Common.dll` | RC2 encryption (`Crypto`), XY data types |
| `Clearcore2.Library.Core.dll` | Domain entities (Compound, MassSpectrum, Transition, etc.) |
| `System.Data.SqlServerCe.dll` | SqlCE 3.5 database engine |
| `WindowsBase.dll` | `System.IO.Packaging` (OPC/ZIP) |
| `System` | `System.IO.Compression.GZipStream` |
