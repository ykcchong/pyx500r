"""Unit test suite for pyx500r.

Run with::

    python -m unittest tests.test_suite -v

Coverage includes:

    * Crypto: page-level encryption/decryption, password derivation
    * TOF codec: unit tests against hand-crafted streams, edge cases
    * TOF calibration: bin↔mass round-trips, parameter validation
    * API surface: dataclass fields, public exports
"""

from __future__ import annotations

import struct
import unittest

# ---------------------------------------------------------------------------
# Imports from the package under test
# ---------------------------------------------------------------------------
from pyx500r.crypto import (
    PAGE_SIZE,
    RESERVED_BYTES,
    WIFF2_PASSWORD,
    decrypt_database,
    decrypt_page,
)
from pyx500r.tof import TofCalibration, decompress_tof
from pyx500r.models import (
    Chromatogram,
    ExperimentInfo,
    InstrumentInfo,
    SampleInfo,
    SpectrumData,
    SpectrumMetadata,
)

# ===================================================================
# 1. Crypto — encryption/decryption primitives
# ===================================================================


class TestCrypto(unittest.TestCase):
    """Unit tests for the SQLite SEE AES-128-OFB page cipher."""

    # -- password / key --------------------------------------------------

    def test_password_is_expected_format(self) -> None:
        self.assertIsInstance(WIFF2_PASSWORD, str)
        self.assertEqual(len(WIFF2_PASSWORD), 36)  # UUID format
        self.assertEqual(len(WIFF2_PASSWORD.encode("utf-8")), 36)  # all ASCII

    def test_key_is_first_16_bytes(self) -> None:
        key = WIFF2_PASSWORD.encode("utf-8")[:16]
        self.assertEqual(len(key), 16)
        self.assertEqual(key[:15], b"F90CA3B4-CC7B-4")  # first 15 bytes

    # -- constants -------------------------------------------------------

    def test_page_size_constant(self) -> None:
        self.assertEqual(PAGE_SIZE, 4096)

    def test_reserved_bytes_constant(self) -> None:
        self.assertEqual(RESERVED_BYTES, 12)

    # -- decrypt_page ----------------------------------------------------

    def test_decrypt_page_rejects_wrong_size(self) -> None:
        with self.assertRaises(ValueError):
            decrypt_page(1, b"too short")

    def test_decrypt_page_1_restores_sqlite_magic(self) -> None:
        """Page 1 must produce the SQLite header magic after decryption."""
        import os
        from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

        key = WIFF2_PASSWORD.encode("utf-8")[:16]
        nonce = os.urandom(12)
        iv = struct.pack("<I", 1) + nonce
        # Encrypt 4084 zero bytes with OFB
        plain_region = bytes(4084)
        keystream = (
            Cipher(algorithms.AES(key), modes.OFB(iv))
            .encryptor()
            .update(b"\x00" * 4084)
        )
        cipher_region = bytes(a ^ b for a, b in zip(plain_region, keystream))
        # Insert SQLite header bookkeeping bytes at [16:24]
        # Bytes 16:20 = page size (uint16 BE), 20:21 = reserved (uint8) + unused
        header = struct.pack(">H", 4096) + struct.pack(">H", 12) + b"\x00" * 4
        cipher_region = bytearray(cipher_region)
        cipher_region[16:24] = header
        page = bytes(cipher_region) + nonce

        result = decrypt_page(1, page)
        self.assertEqual(result[:16], b"SQLite format 3\x00")
        self.assertEqual(result[16:20], header[:4])

    def test_decrypt_page_round_trip(self) -> None:
        """Encrypt then decrypt should recover original plaintext."""
        import os
        from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

        key = WIFF2_PASSWORD.encode("utf-8")[:16]
        original = os.urandom(4084)
        nonce = os.urandom(12)
        pageno = 42

        iv = struct.pack("<I", pageno) + nonce
        keystream = (
            Cipher(algorithms.AES(key), modes.OFB(iv))
            .encryptor()
            .update(b"\x00" * 4084)
        )
        ciphertext = bytes(a ^ b for a, b in zip(original, keystream))
        page = ciphertext + nonce

        result = decrypt_page(pageno, page)
        self.assertEqual(result[:4084], original)

    def test_decrypt_database_rejects_non_page_aligned_data(self) -> None:
        with self.assertRaises(ValueError):
            decrypt_database(b"not a multiple of 4096")


# ===================================================================
# 2. TOF codec — decompression unit tests
# ===================================================================


class TestTofCodecUnit(unittest.TestCase):
    """Hand-crafted bitstreams to validate every codec path."""

    # -- basic intensity emission ----------------------------------------

    def test_single_byte_intensity(self) -> None:
        bins, ints = decompress_tof(b"\x0a\x07")
        self.assertEqual(bins, [1, 2])
        self.assertEqual(ints, [10, 7])

    def test_zero_intensity_is_skipped(self) -> None:
        bins, ints = decompress_tof(b"\x00\x05")
        self.assertEqual(bins, [2])
        self.assertEqual(ints, [5])

    # -- zero-run / gap tokens -------------------------------------------

    def test_zero_run_advances_bins(self) -> None:
        # 0x83 = high bit set + value 3 → skip 3 bins, then intensity 5
        bins, ints = decompress_tof(b"\x83\x05")
        self.assertEqual(bins, [4])  # start=1, skip 3 → 4
        self.assertEqual(ints, [5])

    def test_consecutive_zero_runs(self) -> None:
        bins, ints = decompress_tof(b"\x82\x81\x05")
        self.assertEqual(bins, [4])  # skip 2, skip 1 → bin 4
        self.assertEqual(ints, [5])

    # -- multi-byte value encodings --------------------------------------

    def test_one_byte_value(self) -> None:
        # token 124 → next byte is value
        bins, ints = decompress_tof(b"\x7c\xc8")
        self.assertEqual(ints, [200])

    def test_two_byte_value(self) -> None:
        # token 125 → next two bytes little-endian
        stream = b"\x7d" + struct.pack("<H", 5000)
        bins, ints = decompress_tof(stream)
        self.assertEqual(ints, [5000])

    def test_four_byte_value(self) -> None:
        # token 126 → next four bytes little-endian
        stream = b"\x7e" + struct.pack("<I", 70000)
        bins, ints = decompress_tof(stream)
        self.assertEqual(ints, [70000])

    def test_mixed_value_sizes(self) -> None:
        # 1-byte val (10), 2-byte val (300), 4-byte val (50000)
        stream = b"\x0a" + b"\x7d" + struct.pack("<H", 300) + b"\x7e" + struct.pack("<I", 50000)
        bins, ints = decompress_tof(stream)
        self.assertEqual(bins, [1, 2, 3])
        self.assertEqual(ints, [10, 300, 50000])

    # -- fixed-bin marker ------------------------------------------------

    def test_fixed_bin_marker_sets_start_bin(self) -> None:
        stream = b"\xff\xff\xff\xff" + struct.pack("<I", 1000) + b"\x09"
        bins, ints = decompress_tof(stream, number_of_time_bins_to_sum=4)
        self.assertEqual(bins, [1000])
        self.assertEqual(ints, [9])

    def test_fixed_bin_marker_uses_step(self) -> None:
        stream = (
            b"\xff\xff\xff\xff"
            + struct.pack("<I", 1000)
            + b"\x09\x09"
        )
        bins, ints = decompress_tof(stream, number_of_time_bins_to_sum=4)
        self.assertEqual(bins, [1000, 1004])

    # -- stop marker & edge cases ----------------------------------------

    def test_stop_marker_terminates(self) -> None:
        bins, ints = decompress_tof(b"\x05\xff\x06")
        self.assertEqual(ints, [5])  # 0x06 after 0xff is ignored
        self.assertEqual(len(bins), 1)

    def test_stop_marker_as_first_byte(self) -> None:
        bins, ints = decompress_tof(b"\xff")
        self.assertEqual(bins, [])
        self.assertEqual(ints, [])

    def test_empty_stream(self) -> None:
        bins, ints = decompress_tof(b"")
        self.assertEqual(bins, [])
        self.assertEqual(ints, [])

    # -- zero-run value sizes --------------------------------------------

    def test_zero_run_one_byte_value(self) -> None:
        # high bit + token 124 → high-bit zero-run, one-byte length
        bins, ints = decompress_tof(b"\xfc\x40\x05")
        # 0xfc: high bit set, token=124 → next byte(0x40=64) is gap length
        # skip 64, then intensity 5
        self.assertEqual(bins, [65])
        self.assertEqual(ints, [5])

    def test_zero_run_two_byte_value(self) -> None:
        # high bit + token 125 → two-byte gap length
        stream = b"\xfd" + struct.pack("<H", 200) + b"\x0a"
        bins, ints = decompress_tof(stream)
        self.assertEqual(bins, [201])
        self.assertEqual(ints, [10])

    # -- min_bin filtering -----------------------------------------------

    def test_min_bin_skips_low_bins(self) -> None:
        bins, ints = decompress_tof(b"\x0a\x0a\x0a", min_bin=2)
        self.assertEqual(bins, [2, 3])
        self.assertEqual(ints, [10, 10])


# ===================================================================
# 3. TOF calibration — bin↔mass conversion
# ===================================================================


class TestTofCalibration(unittest.TestCase):
    """Validate the quadratic calibration formula and round-trip stability."""

    CAL_A = 0.000490142680921281
    CAL_T0 = 6.80890033434671
    TIME_RES = 0.0260416666666667

    def _make_cal(self) -> TofCalibration:
        return TofCalibration(
            cal_a=self.CAL_A, cal_t0=self.CAL_T0, time_resolution=self.TIME_RES
        )

    def test_bin_to_mass_idempotent(self) -> None:
        cal = self._make_cal()
        m1 = cal.bin_to_mass(554908)
        m2 = cal.bin_to_mass(554908)
        self.assertAlmostEqual(m1, m2, places=12)

    def test_mass_to_bin_round_trip(self) -> None:
        cal = self._make_cal()
        mass = cal.bin_to_mass(554908)
        recovered = cal.bin_to_mass(cal.mass_to_bin(mass))
        self.assertAlmostEqual(mass, recovered, places=6)

    def test_bins_to_masses_bulk(self) -> None:
        cal = self._make_cal()
        bins = [554908, 554912, 554916]
        masses = cal.bins_to_masses(bins)
        self.assertEqual(len(masses), 3)
        for b, m in zip(bins, masses):
            self.assertAlmostEqual(cal.bin_to_mass(b), m, places=12)

    def test_mass_to_bin_zero_mass(self) -> None:
        cal = self._make_cal()
        self.assertEqual(cal.mass_to_bin(0.0), 0.0)

    def test_mass_to_bin_negative_mass(self) -> None:
        cal = self._make_cal()
        self.assertEqual(cal.mass_to_bin(-1.0), 0.0)

    # -- validation ------------------------------------------------------

    def test_rejects_zero_cal_a(self) -> None:
        with self.assertRaises(ValueError):
            TofCalibration(cal_a=0.0, cal_t0=1.0, time_resolution=1.0)

    def test_rejects_negative_cal_a(self) -> None:
        with self.assertRaises(ValueError):
            TofCalibration(cal_a=-0.1, cal_t0=1.0, time_resolution=1.0)

    def test_rejects_zero_time_resolution(self) -> None:
        with self.assertRaises(ValueError):
            TofCalibration(cal_a=1.0, cal_t0=1.0, time_resolution=0.0)


# ===================================================================
# 4. Calibrated decompression — TOF + calibration fused
# ===================================================================


class TestCalibratedDecompression(unittest.TestCase):
    """Verify that the fused JIT kernel matches manual calibration."""

    def test_fused_matches_manual(self) -> None:
        stream = b"\xff\xff\xff\xff" + struct.pack("<I", 1000) + b"\x09\x09"
        cal_a = 0.000490142680921281
        cal_t0 = 6.80890033434671
        tdc_res = 0.0260416666666667

        # Un-calibrated path
        bins, ints = decompress_tof(stream, number_of_time_bins_to_sum=4)
        cal = TofCalibration(cal_a=cal_a, cal_t0=cal_t0, time_resolution=tdc_res)
        manual_mz = cal.bins_to_masses(bins)

        # Calibrated (fused) path
        mz, ints2 = decompress_tof(
            stream,
            number_of_time_bins_to_sum=4,
            cal_a=cal_a,
            cal_t0=cal_t0,
            time_resolution=tdc_res,
        )

        self.assertEqual(ints, ints2)
        for a, b in zip(mz, manual_mz):
            self.assertAlmostEqual(a, b, places=10)


# ===================================================================
# 5. API surface — dataclasses & public exports
# ===================================================================


class TestApiSurface(unittest.TestCase):
    """Validate that public dataclasses and exports are correct."""

    def test_sample_info_fields(self) -> None:
        s = SampleInfo(
            index=0, sample_id="abc", name="test",
            source="file.wiff2", start_timestamp="2025-01-01T00:00:00Z",
        )
        self.assertEqual(s.index, 0)
        self.assertEqual(s.sample_id, "abc")
        self.assertTrue(s.__dataclass_fields__)  # is a dataclass

    def test_experiment_info_fields(self) -> None:
        e = ExperimentInfo(
            index=0, experiment_id="1", scan_type="TOFMS",
            ms_level=1, polarity="positive", cycle_count=100,
        )
        self.assertEqual(e.scan_type, "TOFMS")

    def test_instrument_info_fields(self) -> None:
        i = InstrumentInfo(
            sample_index=0, instrument_index=0, device_type=0,
            device_name="MS", model_name="X500", serial_number="SN123",
            is_mass_spectrometer=True,
        )
        self.assertTrue(i.is_mass_spectrometer)

    def test_chromatogram_fields(self) -> None:
        c = Chromatogram(times=[1.0, 2.0], intensities=[100.0, 200.0])
        self.assertEqual(len(c.times), 2)

    def test_spectrum_metadata_fields(self) -> None:
        sm = SpectrumMetadata(
            sample_index=0, experiment_index=0, cycle_index=0,
            scan_time=1.0, scan_type="TOFMS", ms_level=1,
            polarity="positive", point_count=100,
        )
        self.assertEqual(sm.point_count, 100)

    def test_spectrum_data_fields(self) -> None:
        sd = SpectrumData(
            sample_index=0, experiment_index=0, cycle_index=0,
            scan_time=1.0, mz=[50.0, 100.0], intensities=[10, 20],
            centroided=False,
        )
        self.assertEqual(len(sd.mz), 2)
        self.assertFalse(sd.centroided)

    def test_all_models_are_frozen(self) -> None:
        """All model dataclasses should be frozen (immutable)."""
        for cls in (
            SampleInfo, ExperimentInfo, InstrumentInfo,
            Chromatogram, SpectrumMetadata, SpectrumData,
        ):
            with self.subTest(cls=cls.__name__):
                self.assertTrue(
                    cls.__dataclass_params__.frozen,
                    f"{cls.__name__} is not frozen",
                )

    def test_public_exports(self) -> None:
        """All major symbols are exported from pyx500r."""
        import pyx500r
        expected = [
            "WiffReader", "open_wiff2",
            "TofCalibration", "decompress_tof",
            "WIFF2_PASSWORD", "PAGE_SIZE", "RESERVED_BYTES",
            "decrypt_database", "decrypt_page",
            "SampleInfo", "ExperimentInfo", "InstrumentInfo",
            "Chromatogram", "SpectrumData", "SpectrumMetadata",
        ]
        for name in expected:
            self.assertTrue(
                hasattr(pyx500r, name),
                f"pyx500r.{name} is not exported",
            )


# ===================================================================
# main
# ===================================================================

if __name__ == "__main__":
    unittest.main(verbosity=2)
