"""Pure-Python parser for .NET BinaryFormatter streams.

Handles the record types needed for SCIEX qsession RTParts:

*   ``SerializedStreamHeader`` (0x00)
*   ``BinaryLibrary`` (0x0C)
*   ``ClassWithMembersAndTypes`` (0x05)
*   ``SystemClassWithMembersAndTypes`` (0x04)
*   ``ArraySinglePrimitive`` (0x0F)
*   ``BinaryArray`` (0x07)
*   ``MemberPrimitiveTyped`` (0x08)
*   ``BinaryObjectString`` (0x06)
*   ``MemberReference`` (0x09)
*   ``ObjectNull`` (0x0A)
*   ``ObjectNullMultiple256`` (0x0D)
*   ``ObjectNullMultiple`` (0x0E)
*   ``MessageEnd`` (0x0B)

Typed helpers are provided for the specific objects encountered in
MultiQuant RTParts streams (``IntegrationParameters``, ``Hashtable``,
``DateTime``, ``float[]``, ``double[]``, ``int[]``, ``List<int>``).
"""
from __future__ import annotations

import struct
from datetime import datetime, timedelta
from typing import Any
import json
import os
import subprocess
import tempfile
from pathlib import Path


def _find_bf_cli() -> str | None:
    """Locate the ``bf_cli.exe`` C# fallback binary."""
    # Search relative to this module, then in PATH
    module_dir = Path(__file__).resolve().parent
    candidates = [
        module_dir / "bf_cli.exe",
        module_dir.parent / "bf_cli.exe",
        module_dir.parent.parent / "bf_cli.exe",
        Path("bf_cli.exe"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    # Try mono in PATH
    for env_path in os.environ.get("PATH", "").split(os.pathsep):
        p = Path(env_path) / "bf_cli.exe"
        if p.exists():
            return str(p)
    return None


def parse_bf_with_fallback(data: bytes, offset: int = 0) -> tuple[Any, int]:
    """Parse a BinaryFormatter blob, falling back to C# ``bf_cli`` on failure.

    Some .NET classes (e.g. ``XicManagerXic``) use a non-standard wire-format
    ``BinaryTypeEnum`` mapping that the pure-Python parser cannot resolve.
    When the pure-Python parser raises :class:`BinaryFormatterError`, this
    function writes the bytes to a temporary file and invokes ``bf_cli.exe``
    (compiled from ``bf_cli.cs``) via Mono.

    Returns ``(object, consumed_bytes)``.  If both parsers fail, the original
    exception is re-raised.
    """
    try:
        return parse_bf_with_consumed(data, offset)
    except (BinaryFormatterError, UnicodeDecodeError, struct.error):
        pass

    cli = _find_bf_cli()
    if cli is None:
        raise BinaryFormatterError(
            "bf_cli.exe not found; cannot fall back to C# parser"
        ) from None

    with tempfile.NamedTemporaryFile(delete=False, suffix=".bin") as tmp:
        tmp.write(data[offset:])
        tmp_path = tmp.name

    try:
        cmd = ["mono", cli, tmp_path, "0"]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            raise BinaryFormatterError(
                f"bf_cli failed: {result.stderr or result.stdout}"
            )
        output = json.loads(result.stdout)
        if not output.get("success"):
            raise BinaryFormatterError(
                f"bf_cli parse error: {output.get('error', 'unknown')}"
            )
        return output["data"], output["consumed"]
    finally:
        os.unlink(tmp_path)



class BinaryFormatterError(ValueError):
    """Raised when the stream contains an unexpected record or is truncated."""

# Record-type constants
# ──────────────────────────────────────────────────────────
class _Rec:
    Header = 0x00
    RefTypeObject = 0x01
    ClassWithId = 0x02
    SystemClassWithMembers = 0x03
    SystemClassWithMembersAndTypes = 0x04
    ClassWithMembersAndTypes = 0x05
    BinaryObjectString = 0x06
    BinaryArray = 0x07
    MemberPrimitiveTyped = 0x08
    MemberReference = 0x09
    ObjectNull = 0x0A
    MessageEnd = 0x0B
    BinaryLibrary = 0x0C
    ObjectNullMultiple256 = 0x0D
    ObjectNullMultiple = 0x0E
    ArraySinglePrimitive = 0x0F
    ArraySingleObject = 0x10
    ArraySingleString = 0x11


# ──────────────────────────────────────────────────────────
# Primitive-type enum
# ──────────────────────────────────────────────────────────
class _Prim:
    Boolean = 1
    Byte = 2
    Char = 3
    Decimal = 5
    Double = 6
    Int16 = 7
    Int32 = 8
    Int64 = 9
    SByte = 10
    Single = 11
    TimeSpan = 12
    DateTime = 13
    UInt16 = 14
    UInt32 = 15
    UInt64 = 16
    Null = 17
    String = 18


# ──────────────────────────────────────────────────────────
# BinaryType enum
# ──────────────────────────────────────────────────────────
class _BAType:
    """BinaryArrayTypeEnum – distinct from BinaryTypeEnum."""
    Single = 0
    Double = 1
    Decimal = 2
    Boolean = 3
    Int16 = 4
    Int32 = 5
    Int64 = 6
    String = 7
    StringArray = 8
    ObjectArray = 9

class _Ref:
    """Lazy forward-reference placeholder for MS-NRBF MemberReference."""
    __slots__ = ("ref_id",)

    def __init__(self, ref_id: int) -> None:
        self.ref_id = ref_id

    def resolve(self, objects: dict[int, Any]) -> Any:
        return objects.get(self.ref_id)


class _BType:
    Primitive = 0
    String = 1
    Object = 2
    SystemClass = 3
    Class = 4
    ObjectArray = 5
    StringArray = 6
    PrimitiveArray = 7


def _resolve_refs(obj: Any, objects: dict[int, Any]) -> Any:
    """Walk *obj* recursively and replace any ``_Ref`` with its resolved value."""
    if isinstance(obj, _Ref):
        resolved = obj.resolve(objects)
        return _resolve_refs(resolved, objects)
    if isinstance(obj, dict):
        return {k: _resolve_refs(v, objects) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve_refs(v, objects) for v in obj]
    return obj


class _BfReader:
    """Low-level BinaryFormatter byte reader with object reference tracking."""

    def __init__(self, data: bytes, offset: int = 0) -> None:
        self.raw = data
        self.pos = offset
        self.objects: dict[int, Any] = {}
        self.libraries: dict[int, str] = {}
        self.class_info: dict[int, dict[str, Any]] = {}
    def _u8(self) -> int:
        b = self.raw[self.pos]
        self.pos += 1
        return b

    def _i16(self) -> int:
        v = struct.unpack_from("<h", self.raw, self.pos)[0]
        self.pos += 2
        return v

    def _i32(self) -> int:
        v = struct.unpack_from("<i", self.raw, self.pos)[0]
        self.pos += 4
        return v

    def _i64(self) -> int:
        v = struct.unpack_from("<q", self.raw, self.pos)[0]
        self.pos += 8
        return v

    def _f32(self) -> float:
        v = struct.unpack_from("<f", self.raw, self.pos)[0]
        self.pos += 4
        return v

    def _f64(self) -> float:
        v = struct.unpack_from("<d", self.raw, self.pos)[0]
        self.pos += 8
        return v

    def _bool(self) -> bool:
        return self._u8() != 0

    def _7bit_int(self) -> int:
        result = 0
        shift = 0
        while True:
            b = self._u8()
            result |= (b & 0x7F) << shift
            if not (b & 0x80):
                break
            shift += 7
        return result

    def _string(self) -> str:
        length = self._7bit_int()
        s = self.raw[self.pos : self.pos + length].decode("utf-8")
        self.pos += length
        return s

    def _decimal(self) -> float:
        """.NET Decimal → float (lossy but sufficient for metadata)."""
        lo = self._i32()
        mid = self._i32()
        hi = self._i32()
        flags = self._i32()
        sign = -1 if (flags >> 31) else 1
        scale = (flags >> 16) & 0xFF
        value = (hi << 64) | (mid << 32) | lo
        return sign * value / (10 ** scale)

    def _datetime(self) -> datetime:
        """.NET DateTime (Int64 ticks + 2-bit kind) → Python datetime."""
        raw = self._i64()
        ticks = raw & 0x3FFFFFFFFFFFFFFF
        # kind = (raw >> 62) & 0x3  # 0=unspecified, 1=UTC, 2=local
        epoch_offset = 621355968000000000  # ticks between 0001-01-01 and 1970-01-01
        python_ticks = ticks - epoch_offset
        seconds = python_ticks / 10_000_000
        microseconds = (python_ticks % 10_000_000) // 10
        return datetime(1970, 1, 1) + timedelta(seconds=seconds, microseconds=microseconds)

    def _primitive(self, ptype: int) -> Any:
        if ptype == _Prim.Boolean:
            return self._bool()
        if ptype == _Prim.Byte:
            return self._u8()
        if ptype == _Prim.Char:
            return chr(self._i16())
        if ptype == _Prim.Decimal:
            return self._decimal()
        if ptype == _Prim.Double:
            return self._f64()
        if ptype == _Prim.Int16:
            return self._i16()
        if ptype == _Prim.Int32:
            return self._i32()
        if ptype == _Prim.Int64:
            return self._i64()
        if ptype == _Prim.SByte:
            return struct.unpack_from("<b", self.raw, self.pos)[0]
        if ptype == _Prim.Single:
            return self._f32()
        if ptype == _Prim.TimeSpan:
            return self._i64()
        if ptype == _Prim.DateTime:
            return self._datetime()
        if ptype == _Prim.UInt16:
            return struct.unpack_from("<H", self.raw, self.pos)[0]
        if ptype == _Prim.UInt32:
            return struct.unpack_from("<I", self.raw, self.pos)[0]
        if ptype == _Prim.UInt64:
            return struct.unpack_from("<Q", self.raw, self.pos)[0]
        if ptype == _Prim.Null:
            return None
        if ptype == _Prim.String:
            return self._string()
        raise BinaryFormatterError(
            f"Unknown primitive type {ptype} at position {self.pos}"
        )

    # ── record readers ──
    def _read_class_info(self) -> tuple[int, str, int, list[str]]:
        """ClassInfo: obj_id, name, member_count, member_names[]."""
        obj_id = self._i32()
        name = self._string()
        member_count = self._i32()
        member_names = [self._string() for _ in range(member_count)]
        return obj_id, name, member_count, member_names

    def _read_member_type_info(self, member_count: int) -> tuple[list[int], list[Any]]:
        """MemberTypeInfo: binary_types[], type_info[].

        For ``Class`` (bt=4) the additional info is a ``ClassTypeInfo``
        containing a type-name string *and* a 4-byte ``LibraryId``.
        For ``SystemClass`` (bt=3) only the type-name string is present.
        """
        binary_types = [self._u8() for _ in range(member_count)]
        type_infos: list[Any] = []
        for bt in binary_types:
            if bt == _BType.Primitive or bt == _BType.PrimitiveArray:
                type_infos.append(self._u8())
            elif bt in (_BType.Class, _BType.SystemClass):
                type_infos.append(self._string())
                if bt == _BType.Class:
                    self._i32()  # consume LibraryId
            else:
                type_infos.append(None)
        return binary_types, type_infos

    def _read_class_with_members_and_types(self, is_system: bool) -> dict[str, Any]:
        obj_id, class_name, member_count, member_names = self._read_class_info()
        binary_types, type_infos = self._read_member_type_info(member_count)
        library_id = None
        if not is_system:
            library_id = self._i32()
        result: dict[str, Any] = {"__class": class_name}
        for mname, btype, tinfo in zip(member_names, binary_types, type_infos):
            result[mname] = self._read_member_value(btype, tinfo)
        self.objects[obj_id] = result
        self.class_info[obj_id] = {
            "class_name": class_name,
            "member_names": member_names,
            "binary_types": binary_types,
            "type_infos": type_infos,
            "library_id": library_id,
        }
        return result
    def _read_class_with_id(self) -> dict[str, Any]:
        """ClassWithId (0x02) – references a previous ClassInfo by metadata_id."""
        obj_id = self._i32()
        metadata_id = self._i32()
        meta = self.class_info.get(metadata_id)
        if meta is None:
            raise BinaryFormatterError(
                f"ClassWithId references unknown metadata_id {metadata_id}"
            )
        result: dict[str, Any] = {"__class": meta["class_name"]}
        for mname, btype, tinfo in zip(
            meta["member_names"], meta["binary_types"], meta["type_infos"]
        ):
            result[mname] = self._read_member_value(btype, tinfo)
        self.objects[obj_id] = result
        return result

    def _read_array_single_object(self) -> list[Any]:
        """ArraySingleObject (0x10)."""
        obj_id = self._i32()
        length = self._i32()
        result = [self._read_record() for _ in range(length)]
        self.objects[obj_id] = result
        return result
    def _read_system_class_with_members(self) -> dict[str, Any]:
        """SystemClassWithMembers (0x03) – no MemberTypeInfo; read values as records."""
        obj_id, class_name, member_count, member_names = self._read_class_info()
        result: dict[str, Any] = {"__class": class_name}
        for mname in member_names:
            result[mname] = self._read_record()
        self.objects[obj_id] = result
        self.class_info[obj_id] = {
            "class_name": class_name,
            "member_names": member_names,
            "binary_types": [],
            "type_infos": [],
            "library_id": None,
        }
        return result
    def _read_member_value(self, btype: int, tinfo: Any) -> Any:
        # For non-primitive types, ObjectNull (0x01),
        # MemberPrimitiveTyped (0x08), and MemberReference (0x09)
        # markers may appear before the actual value.
        # Primitive types skip marker checks to avoid false positives
        # (e.g. Boolean True = 0x01 collides with ObjectNull).
        if btype not in (_BType.Primitive, _BType.PrimitiveArray):
            rec = self.raw[self.pos]
            if rec == _Rec.ObjectNull:
                self.pos += 1
                return None
            if rec == _Rec.MemberPrimitiveTyped:
                self.pos += 1
                return self._primitive(self._u8())
            if rec == _Rec.MemberReference:
                self.pos += 1
                ref_id = self._i32()
                if ref_id not in self.objects:
                    return _Ref(ref_id)
                return self.objects[ref_id]
        if btype == _BType.Primitive:
            return self._primitive(tinfo)
        if btype == _BType.String:
            return self._read_string_value()
        if btype == _BType.Object:
            return self._read_record()
        if btype == _BType.PrimitiveArray:
            return self._read_record()
        if btype == _BType.ObjectArray:
            return self._read_record()
        if btype == _BType.StringArray:
            return self._read_record()
        if btype in (_BType.Class, _BType.SystemClass):
            return self._read_record()
        raise BinaryFormatterError(
            f"Unsupported BinaryTypeEnum {btype} at position {self.pos}"
        )

    def _read_string_value(self) -> str | None:
        rec = self.raw[self.pos]
        if rec == _Rec.ObjectNull:
            self.pos += 1
            return None
        if rec == _Rec.BinaryObjectString:
            self.pos += 1
            return self._read_binary_object_string()
        if rec == _Rec.MemberReference:
            self.pos += 1
            ref_id = self._i32()
            return self.objects.get(ref_id)
        # Inline record (e.g. ClassWithMembersAndTypes for StringBuilder)
        return self._read_record()

    def _read_binary_object_string(self) -> str:
        obj_id = self._i32()
        val = self._string()
        self.objects[obj_id] = val
        return val

    def _read_array_single_primitive(self) -> list[Any]:
        obj_id = self._i32()
        length = self._i32()
        ptype = self._u8()

        if ptype == _Prim.Single:
            data = struct.unpack_from(f"<{length}f", self.raw, self.pos)
            self.pos += length * 4
        elif ptype == _Prim.Double:
            data = struct.unpack_from(f"<{length}d", self.raw, self.pos)
            self.pos += length * 8
        elif ptype == _Prim.Int32:
            data = struct.unpack_from(f"<{length}i", self.raw, self.pos)
            self.pos += length * 4
        elif ptype == _Prim.Int16:
            data = struct.unpack_from(f"<{length}h", self.raw, self.pos)
            self.pos += length * 2
        elif ptype == _Prim.Int64:
            data = struct.unpack_from(f"<{length}q", self.raw, self.pos)
            self.pos += length * 8
        elif ptype == _Prim.Boolean:
            data = [self.raw[self.pos + i] != 0 for i in range(length)]
            self.pos += length
        elif ptype == _Prim.Byte:
            data = list(self.raw[self.pos : self.pos + length])
            self.pos += length
        else:
            raise BinaryFormatterError(
                f"Unsupported array primitive type {ptype} at {self.pos}"
            )
        data = list(data)
        self.objects[obj_id] = data
        return data

    def _read_binary_array(self) -> list[Any]:
        """BinaryArray (0x07) – MS-NRBF §2.5.5."""
        obj_id = self._i32()
        _array_type = self._u8()  # BinaryArrayTypeEnum (unused by reader)
        rank = self._i32()
        lengths = [self._i32() for _ in range(rank)]
        # LowerBounds are optional; skip them (always zero for 1-D arrays).

        # Element type: BinaryTypeEnum byte followed by additional info.
        element_bt = self._u8()
        element_info: Any = None
        if element_bt in (_BType.Primitive, _BType.PrimitiveArray):
            element_info = self._u8()
        elif element_bt in (_BType.Class, _BType.SystemClass):
            element_info = self._string()
            if element_bt == _BType.Class:
                self._i32()  # LibraryId

        total_length = 1
        for ln in lengths:
            total_length *= ln

        result: list[Any] = []
        for _ in range(total_length):
            if element_bt == _BType.Primitive:
                result.append(self._primitive(element_info))
            elif element_bt in (_BType.Class, _BType.SystemClass, _BType.Object):
                result.append(self._read_record())
            elif element_bt == _BType.String:
                result.append(self._read_string_value())
            else:
                result.append(self._read_record())

        self.objects[obj_id] = result
        return result

    def _read_hashtable(self) -> dict[Any, Any]:
        """Read a System.Collections.Hashtable from BinaryFormatter."""
        # Skip Header (0x00) and BinaryLibrary (0x0C) records
        while self.pos < len(self.raw):
            rec = self.raw[self.pos]
            if rec in (_Rec.SystemClassWithMembersAndTypes, _Rec.ClassWithMembersAndTypes):
                break
            if rec == _Rec.Header:
                self.pos += 17
                continue
            if rec == _Rec.BinaryLibrary:
                self.pos += 1
                lib_id = self._i32()
                lib_name = self._string()
                self.libraries[lib_id] = lib_name
                continue
            raise BinaryFormatterError(
                f"Expected Hashtable class record, got 0x{rec:02x} at {self.pos}"
            )
        if self.pos >= len(self.raw):
            raise BinaryFormatterError("Unexpected end of stream looking for Hashtable")

        is_system = self.raw[self.pos] == _Rec.SystemClassWithMembersAndTypes
        self.pos += 1
        obj = self._read_class_with_members_and_types(is_system)
        # If the object has a direct ``buckets`` member it is the older
        # field-based format; otherwise we try to reconstruct from the
        # ``m_info`` SerializationInfo entries.
        if "buckets" in obj:
            # Old bucket-based format – iterate bucket array and collect
            # non-null key/value pairs.
            result: dict[Any, Any] = {}
            for bucket in obj.get("buckets", []):
                if bucket is not None and isinstance(bucket, dict):
                    key = bucket.get("key")
                    val = bucket.get("val")
                    if key is not None:
                        result[key] = val
            return result

        # SerializationInfo-style (``m_info`` member)
        if "m_info" in obj:
            # The m_info object is itself a ClassWithMembersAndTypes
            # with ``Data`` (array of key/value pairs).
            info = obj.get("m_info", {})
            data = info.get("Data", []) if isinstance(info, dict) else []
            result = {}
            for entry in data:
                if isinstance(entry, dict):
                    key = entry.get("Name")
                    val = entry.get("Value")
                    if key is not None:
                        result[key] = val
            return result

        # Fallback: just return the raw dict if we cannot interpret it.
        return obj

    def _read_record(self) -> Any:
        if self.pos >= len(self.raw):
            raise BinaryFormatterError("Unexpected end of stream")

        rec = self._u8()

        if rec == _Rec.Header:
            struct.unpack_from("<HHI", self.raw, self.pos)  # major, minor, header_id
            self.pos += 8
            struct.unpack_from("<II", self.raw, self.pos)  # top_id, header_ref
            self.pos += 8
            return None  # header has no value

        if rec == _Rec.SystemClassWithMembersAndTypes:
            return self._read_class_with_members_and_types(is_system=True)
        if rec == _Rec.ClassWithMembersAndTypes:
            return self._read_class_with_members_and_types(is_system=False)
        if rec == _Rec.ClassWithId:
            return self._read_class_with_id()
        if rec == _Rec.SystemClassWithMembers:
            return self._read_system_class_with_members()
        if rec == _Rec.BinaryObjectString:
            return self._read_binary_object_string()

        if rec == _Rec.MemberPrimitiveTyped:
            ptype = self._u8()
            return self._primitive(ptype)

        if rec == _Rec.MemberReference:
            ref_id = self._i32()
            return self.objects.get(ref_id)

        if rec == _Rec.ObjectNull or rec == 0x01:
            return None

        if rec == _Rec.MessageEnd:
            return None

        if rec == _Rec.BinaryLibrary:
            lib_id = self._i32()
            lib_name = self._string()
            self.libraries[lib_id] = lib_name
            return None

        if rec == _Rec.ObjectNullMultiple256:
            count = self._u8()
            return [None] * count

        if rec == _Rec.ObjectNullMultiple:
            count = self._i32()
            return [None] * count

        if rec == _Rec.ArraySinglePrimitive:
            return self._read_array_single_primitive()

        if rec == _Rec.BinaryArray:
            return self._read_binary_array()

        if rec == _Rec.ArraySingleString:
            return self._read_array_single_string()
        if rec == _Rec.ArraySingleObject:
            return self._read_array_single_object()

        raise BinaryFormatterError(
            f"Unsupported record type 0x{rec:02x} at position {self.pos - 1}"
        )

    def _read_array_single_string(self) -> list[str | None]:
        obj_id = self._i32()
        length = self._i32()
        result = [self._read_string_value() for _ in range(length)]
        self.objects[obj_id] = result
        return result
    def parse(self) -> Any:
        """Parse the stream and return the top-level deserialized object."""
        result = None
        while self.pos < len(self.raw):
            rec = self.raw[self.pos]
            if rec == _Rec.MessageEnd:
                self.pos += 1
                break
            obj = self._read_record()
            if result is None and obj is not None:
                result = obj
        if result is not None:
            result = _resolve_refs(result, self.objects)
        return result
    @property
    def bytes_consumed(self) -> int:
        return self.pos
# Public API
# ═══════════════════════════════════════════════════════════

def parse_bf(data: bytes, offset: int = 0) -> Any:
    """Parse a BinaryFormatter blob and return the top-level object.

    Returns a ``dict`` for class instances, a ``list`` for arrays, or a
    primitive value for simple types.
    """
    reader = _BfReader(data, offset)
    return reader.parse()


def parse_bf_with_consumed(data: bytes, offset: int = 0) -> tuple[Any, int]:
    """Like ``parse_bf`` but also returns the number of bytes consumed."""
    reader = _BfReader(data, offset)
    obj = reader.parse()
    return obj, reader.bytes_consumed


def parse_integration_parameters(data: bytes, offset: int = 0) -> tuple[dict[str, Any], int]:
    """Parse a BinaryFormatter ``IntegrationParameters`` blob.

    Returns ``(dict, consumed_bytes)`` where *consumed_bytes* is the
    number of bytes consumed from *data* starting at *offset*.
    """
    reader = _BfReader(data, offset)
    obj = reader.parse()
    if isinstance(obj, dict) and "__class" in obj:
        obj.pop("__class", None)
    return (obj if isinstance(obj, dict) else {}, reader.pos - offset)


def parse_hashtable(data: bytes, offset: int = 0) -> dict[Any, Any]:
    """Parse a BinaryFormatter ``System.Collections.Hashtable`` blob."""
    reader = _BfReader(data, offset)
    return reader._read_hashtable()


def parse_datetime(data: bytes, offset: int = 0) -> datetime | None:
    """Parse a BinaryFormatter ``System.DateTime`` blob.

    Expects a ``MemberPrimitiveTyped`` (0x08) record with primitive
    type ``DateTime`` (0x0D).
    """
    if offset >= len(data):
        return None
    reader = _BfReader(data, offset)
    if reader.raw[reader.pos] == _Rec.MemberPrimitiveTyped:
        reader.pos += 1
        ptype = reader._u8()
        if ptype == _Prim.DateTime:
            return reader._datetime()
    # Fallback: try full parse
    reader.pos = offset
    obj = reader.parse()
    if isinstance(obj, datetime):
        return obj
    return None

def parse_float_array(data: bytes, offset: int = 0) -> tuple[list[float] | None, int]:
    """Parse a BinaryFormatter ``float[]`` (``ArraySinglePrimitive``).

    Returns ``(list, consumed_bytes)``.
    """
    if offset >= len(data):
        return None, 0
    reader = _BfReader(data, offset)
    if reader.raw[reader.pos] == _Rec.ArraySinglePrimitive:
        reader.pos += 1
        _ = reader._i32()  # obj_id – tracked internally
        length = reader._i32()
        ptype = reader._u8()
        if ptype == _Prim.Single and length > 0:
            arr = struct.unpack_from(f"<{length}f", reader.raw, reader.pos)
            reader.pos += length * 4
            return list(arr), reader.pos - offset
    reader.pos = offset
    obj = reader.parse()
    if isinstance(obj, list):
        return [float(x) for x in obj], reader.pos - offset
    return None, 0

def parse_double_array(data: bytes, offset: int = 0) -> tuple[list[float] | None, int]:
    """Parse a BinaryFormatter ``double[]`` (``ArraySinglePrimitive``).

    Returns ``(list, consumed_bytes)``.
    """
    if offset >= len(data):
        return None, 0
    reader = _BfReader(data, offset)
    if reader.raw[reader.pos] == _Rec.ArraySinglePrimitive:
        reader.pos += 1
        _ = reader._i32()  # obj_id – tracked internally
        length = reader._i32()
        ptype = reader._u8()
        if ptype == _Prim.Double and length > 0:
            arr = struct.unpack_from(f"<{length}d", reader.raw, reader.pos)
            reader.pos += length * 8
            return list(arr), reader.pos - offset
    reader.pos = offset
    obj = reader.parse()
    if isinstance(obj, list):
        return [float(x) for x in obj], reader.pos - offset
    return None, 0


def parse_int_array(data: bytes, offset: int = 0) -> tuple[list[int] | None, int]:
    """Parse a BinaryFormatter ``int[]`` (``ArraySinglePrimitive``).

    Returns ``(list, consumed_bytes)``.
    """
    if offset >= len(data):
        return None, 0
    reader = _BfReader(data, offset)
    if reader.raw[reader.pos] == _Rec.ArraySinglePrimitive:
        reader.pos += 1
        _ = reader._i32()  # obj_id – tracked internally
        length = reader._i32()
        ptype = reader._u8()
        if ptype == _Prim.Int32 and length > 0:
            arr = struct.unpack_from(f"<{length}i", reader.raw, reader.pos)
            reader.pos += length * 4
            return list(arr), reader.pos - offset
    reader.pos = offset
    obj = reader.parse()
    if isinstance(obj, list):
        return [int(x) for x in obj], reader.pos - offset
    return None, 0



