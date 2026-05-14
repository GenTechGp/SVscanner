#!/usr/bin/env python3
"""
Unit tests for src/retrotransposition_detection.py

Tests cover detect_retrotp() and load_config() in isolation using lightweight
mock VCF records — no real VCF files or pysam VariantFile I/O required.
"""
import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from retrotransposition_detection import detect_retrotp, load_config


# ---------------------------------------------------------------------------
# Minimal mock of a pysam VariantRecord for testing detect_retrotp()
# ---------------------------------------------------------------------------
class MockInfo:
    def __init__(self, data):
        self._data = data

    def __getitem__(self, key):
        if key not in self._data:
            raise KeyError(key)
        return self._data[key]

    def __contains__(self, key):
        return key in self._data


class MockRec:
    def __init__(self, chrom, pos, svtype, svlen, rm_subfamily, rm_coverage, rec_id=None):
        info = {}
        if svtype is not None:
            info["SVTYPE"] = svtype
        if svlen is not None:
            info["SVLEN"] = svlen
        if rm_subfamily is not None:
            info["RM_SUBFAMILY"] = tuple(rm_subfamily)
        if rm_coverage is not None:
            info["RM_SV_COVERAGE"] = tuple(rm_coverage)
        self.info = MockInfo(info)
        self.chrom = chrom
        self.pos = pos
        self.id = rec_id


# ---------------------------------------------------------------------------
# Minimal config fixture
# ---------------------------------------------------------------------------
CONFIG = {
    "AluYa5":  {"class": "SINE",       "confidence": "HIGH", "min_rm_sv_coverage": 0.70, "svlen_min": 50,  "svlen_max": 350},
    "AluSx":   {"class": "SINE",       "confidence": "LOW",  "min_rm_sv_coverage": 0.70, "svlen_min": 50,  "svlen_max": 350},
    "L1Hs":    {"class": "LINE",       "confidence": "HIGH", "min_rm_sv_coverage": 0.50, "svlen_min": 100, "svlen_max": 6500},
    "L1PA3":   {"class": "LINE",       "confidence": "LOW",  "min_rm_sv_coverage": 0.50, "svlen_min": 100, "svlen_max": 6500},
    "SVA_F":   {"class": "Retroposon", "confidence": "HIGH", "min_rm_sv_coverage": 0.60, "svlen_min": 700, "svlen_max": 4500},
}


class TestDetectRetrtp(unittest.TestCase):

    # --- basic positive cases ---

    def test_single_alu_high_confidence(self):
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5"], [0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])
        self.assertEqual(confs, ["HIGH"])

    def test_single_l1hs_high_confidence(self):
        rec = MockRec("chr1", 100, "INS", 3000, ["L1Hs"], [0.75])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["L1Hs"])
        self.assertEqual(confs, ["HIGH"])

    def test_sva_f_high_confidence(self):
        rec = MockRec("chr1", 100, "INS", 2500, ["SVA_F"], [0.65])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["SVA_F"])
        self.assertEqual(confs, ["HIGH"])

    def test_low_confidence_subfamily(self):
        rec = MockRec("chr1", 100, "INS", 280, ["AluSx"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluSx"])
        self.assertEqual(confs, ["LOW"])

    # --- multi-hit cases ---

    def test_multi_hit_same_subfamily(self):
        """Two RM hits of the same active subfamily (RM fragmentation) — both reported, no deduplication."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluYa5"], [0.75, 0.72])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5", "AluYa5"])
        self.assertEqual(confs, ["HIGH", "HIGH"])

    def test_multi_hit_mixed_subfamily_one_active(self):
        """One active hit, one inactive — active position qualifies, inactive gets '.'."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluJb"], [0.75, 0.20])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5", "."])
        self.assertEqual(confs, ["HIGH", "."])

    def test_multi_hit_mixed_confidence(self):
        """Two active subfamilies of different confidence tiers — both reported in order."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluSx"], [0.75, 0.72])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5", "AluSx"])
        self.assertEqual(confs, ["HIGH", "LOW"])

    def test_multi_hit_mixed_class(self):
        """Active hits from two different classes, both meeting thresholds — both reported."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "L1Hs"], [0.75, 0.55])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5", "L1Hs"])
        self.assertEqual(confs, ["HIGH", "HIGH"])

    # --- '.' placeholder and list-length parity ---

    def test_dot_placeholder_coverage_fail(self):
        """Active subfamily but coverage below threshold — gets '.' not omitted."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluSx"], [0.75, 0.50])
        subs, confs = detect_retrotp(rec, CONFIG)
        # AluSx cov 0.50 < 0.70 threshold → '.'
        self.assertEqual(subs, ["AluYa5", "."])
        self.assertEqual(confs, ["HIGH", "."])

    def test_dot_placeholder_svlen_fail(self):
        """Active subfamily but SVLEN outside its range — gets '.' not omitted."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "SVA_F"], [0.75, 0.65])
        subs, confs = detect_retrotp(rec, CONFIG)
        # SVA_F requires svlen in [700, 4500]; 280 is outside → '.'
        self.assertEqual(subs, ["AluYa5", "."])
        self.assertEqual(confs, ["HIGH", "."])

    def test_list_length_parity(self):
        """Output lists are always the same length as RM_SUBFAMILY."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluJb", "AluSx"], [0.75, 0.85, 0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(len(subs), 3)
        self.assertEqual(len(confs), 3)
        self.assertEqual(subs[0], "AluYa5")
        self.assertEqual(subs[1], ".")    # AluJb not in config
        self.assertEqual(subs[2], "AluSx")

    def test_all_inactive_returns_none(self):
        """All RM hits from inactive subfamilies — tag absent (None, None), not a list of '.'."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluJb", "AluJr"], [0.85, 0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)
        self.assertIsNone(confs)

    # --- no-call cases ---

    def test_non_ins_svtype(self):
        rec = MockRec("chr1", 100, "DEL", -280, ["AluYa5"], [0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)
        self.assertIsNone(confs)

    def test_svtype_missing(self):
        rec = MockRec("chr1", 100, None, 280, ["AluYa5"], [0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_svlen_missing(self):
        rec = MockRec("chr1", 100, "INS", None, ["AluYa5"], [0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_no_rm_subfamily(self):
        rec = MockRec("chr1", 100, "INS", 280, None, None)
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_subfamily_not_in_config(self):
        """Subfamily present but not in active list — no call."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluJb"], [0.85])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    # --- threshold boundary tests ---

    def test_coverage_exactly_at_threshold(self):
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5"], [0.70])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])

    def test_coverage_just_below_threshold(self):
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5"], [0.699])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_svlen_at_minimum(self):
        rec = MockRec("chr1", 100, "INS", 50, ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])

    def test_svlen_at_maximum(self):
        rec = MockRec("chr1", 100, "INS", 350, ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])

    def test_svlen_below_minimum(self):
        rec = MockRec("chr1", 100, "INS", 49, ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_svlen_above_maximum(self):
        rec = MockRec("chr1", 100, "INS", 351, ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertIsNone(subs)

    def test_negative_svlen_treated_as_abs(self):
        """Negative SVLEN (as sometimes written for INS) should use abs value."""
        rec = MockRec("chr1", 100, "INS", -280, ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])

    def test_svlen_as_tuple(self):
        """pysam may return SVLEN as a single-element tuple for some VCF variants."""
        rec = MockRec("chr1", 100, "INS", (280,), ["AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        self.assertEqual(subs, ["AluYa5"])

    # --- list padding edge case ---

    def test_coverage_list_shorter_than_subfamily_list(self):
        """Mismatched list lengths are padded with 0.0 — first passes, second gets '.'."""
        rec = MockRec("chr1", 100, "INS", 280, ["AluYa5", "AluYa5"], [0.80])
        subs, confs = detect_retrotp(rec, CONFIG)
        # First hit passes (cov 0.80), second hit gets padded cov 0.0 → '.'
        self.assertEqual(subs, ["AluYa5", "."])
        self.assertEqual(confs, ["HIGH", "."])
        self.assertEqual(len(subs), 2)


class TestLoadConfig(unittest.TestCase):

    def test_load_real_config(self):
        """Integration: load the actual config/retrotp_params.tsv and verify structure."""
        config_path = os.path.join(
            os.path.dirname(__file__), "..", "config", "retrotp_params.tsv"
        )
        if not os.path.exists(config_path):
            self.skipTest("config/retrotp_params.tsv not found")
        cfg = load_config(config_path)
        self.assertGreater(len(cfg), 0)
        for sub, row in cfg.items():
            self.assertIn("confidence", row)
            self.assertIn(row["confidence"], ("HIGH", "LOW"))
            self.assertIn("min_rm_sv_coverage", row)
            self.assertIsInstance(row["min_rm_sv_coverage"], float)
            self.assertIn("svlen_min", row)
            self.assertIn("svlen_max", row)
            self.assertLessEqual(row["svlen_min"], row["svlen_max"])

    def test_load_config_skips_comments(self):
        """Comments and blank lines must not produce config entries."""
        import tempfile
        content = (
            "# comment\n"
            "\n"
            "subfamily\tclass\tconfidence\tmin_rm_sv_coverage\tsvlen_min\tsvlen_max\n"
            "AluYa5\tSINE\tHIGH\t0.70\t50\t350\n"
        )
        with tempfile.NamedTemporaryFile("wt", suffix=".tsv", delete=False) as f:
            f.write(content)
            path = f.name
        try:
            cfg = load_config(path)
            self.assertEqual(list(cfg.keys()), ["AluYa5"])
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
