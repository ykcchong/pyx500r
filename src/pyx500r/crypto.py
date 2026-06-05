"""Pure-Python decryption of SCIEX encrypted SQLite containers.

SCIEX ``.wiff2`` and ``.timeseries.data`` files are SQLite databases encrypted
with the SQLite Encryption Extension (SEE) using its **AES-128-OFB** cipher.

The layout is:

* Page size 4096 bytes, with 12 reserved bytes at the end of every page.
* The 12 reserved bytes (page offset ``4084:4096``) hold a per-page random
  nonce, stored in clear.
* Each page is encrypted with AES-128 in OFB mode over bytes ``[0:4084]``.
* The AES-128 key is the first 16 bytes of the UTF-8 connection password.
* The per-page 16-byte IV is ``struct.pack("<I", page_number) + nonce`` where
  ``page_number`` is the 1-based SQLite page number.
* Page 1 keeps the first 24 bytes of the SQLite header partly in clear: after
  decrypting, bytes ``[0:16]`` are replaced with the ``"SQLite format 3\\x00"``
  magic and bytes ``[16:24]`` (page-size / reserved-size fields) are taken
  verbatim from the ciphertext.
Only the third-party ``cryptography`` package is required for the AES
primitive.
"""

from __future__ import annotations

import struct
from pathlib import Path

from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

__all__ = [
    "WIFF2_PASSWORD",
    "QSESSION_PASSWORD",
    "PAGE_SIZE",
    "QSESSION_PAGE_SIZE",
    "RESERVED_BYTES",
    "decrypt_page",
    "decrypt_database",
]

# Fixed password used for SCIEX wiff2 containers.
WIFF2_PASSWORD = "F90CA3B4-CC7B-4439-A479-2097CB8AE246"

# Fixed password used for SCIEX qsession containers.
QSESSION_PASSWORD = "PQS1 is not Sirius"

PAGE_SIZE = 4096
QSESSION_PAGE_SIZE = 1024
RESERVED_BYTES = 12
_SQLITE_MAGIC = b"SQLite format 3\x00"


def _aes_key(password: str) -> bytes:
    """Return the AES-128 key: the first 16 bytes of the UTF-8 password."""
    key = password.encode("utf-8")[:16]
    if len(key) < 16:
        raise ValueError("password must be at least 16 bytes when UTF-8 encoded")
    return key


def _decrypt_page_core(
    page_number: int,
    page: bytes,
    key: bytes,
    out: bytearray,
    out_offset: int,
    *,
    page_size: int = PAGE_SIZE,
) -> None:
    """Decrypt a single SQLite page into a pre-allocated bytearray buffer."""
    nonce = page[page_size - RESERVED_BYTES:page_size]
    iv = struct.pack("<I", page_number) + nonce
    region = page[:page_size - RESERVED_BYTES]

    keystream = (
        Cipher(algorithms.AES(key), modes.OFB(iv))
        .encryptor()
        .update(b"\x00" * len(region))
    )

    # Fast vectorised XOR when numpy is available (5-10× faster on 4 KB pages)
    try:
        import numpy as np
        region_arr = np.frombuffer(region, dtype=np.uint8)
        key_arr = np.frombuffer(keystream, dtype=np.uint8)
        plain = (region_arr ^ key_arr).tobytes()
    except Exception:
        # Fallback for environments without numpy
        plain = bytes(a ^ b for a, b in zip(region, keystream))

    out[out_offset:out_offset + page_size - RESERVED_BYTES] = plain
    if page_number == 1:
        out[0:16] = _SQLITE_MAGIC
        out[16:24] = page[16:24]


def decrypt_page(
    page_number: int,
    page: bytes,
    password: str = WIFF2_PASSWORD,
    *,
    page_size: int = PAGE_SIZE,
) -> bytes:
    """Decrypt a single SQLite page.

    ``page_number`` is the 1-based SQLite page number.
    ``page_size`` defaults to 4096 (wiff2); use 1024 for qsession files.
    """
    if len(page) != page_size:
        raise ValueError(f"page must be {page_size} bytes, got {len(page)}")

    key = _aes_key(password)
    out = bytearray(page)
    _decrypt_page_core(page_number, page, key, out, 0, page_size=page_size)
    return bytes(out)


def decrypt_database(
    source: str | Path | bytes,
    password: str = WIFF2_PASSWORD,
    *,
    page_size: int = PAGE_SIZE,
) -> bytes:
    """Decrypt a whole encrypted SQLite container into plaintext SQLite bytes.

    ``source`` may be a path or the raw encrypted bytes. The returned bytes are
    a standard SQLite database openable with the :mod:`sqlite3` stdlib module.

    ``page_size`` defaults to 4096 (wiff2); use 1024 for qsession files.
    """
    data = source if isinstance(source, (bytes, bytearray)) else Path(source).read_bytes()
    if len(data) % page_size != 0:
        raise ValueError(
            f"file size {len(data)} is not a multiple of page size {page_size}"
        )

    key = _aes_key(password)
    num_pages = len(data) // page_size
    out = bytearray(len(data))

    for index in range(num_pages):
        start = index * page_size
        page = data[start:start + page_size]
        _decrypt_page_core(index + 1, page, key, out, start, page_size=page_size)

    return bytes(out)
