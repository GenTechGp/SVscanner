"""Tests for POS-outside-contig handling in src/extract_sv.py. Run: python -m unittest test/test_extract_sv_pos_bounds.py"""
import os
import random
import subprocess
import sys
import tempfile
import unittest

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "extract_sv.py")
CHROM_LEN = 1000
INS_SEQ = "GAATG" * 20
HEADER = (
    "##fileformat=VCFv4.2\n"
    f"##contig=<ID=chr1,length={CHROM_LEN}>\n"
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">\n'
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def ins(pos, vid):
    return f"chr1\t{pos}\t{vid}\tG\tG{INS_SEQ}\t.\tPASS\tSVTYPE=INS;SVLEN={len(INS_SEQ)};END={pos}\n"


def dele(pos, vid, svlen=100):
    return f"chr1\t{pos}\t{vid}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-{svlen};END={pos + svlen}\n"


class ExtractSvPosBounds(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.ref = self.path("ref.fa")
        rng = random.Random(1)
        with open(self.ref, "w") as f:
            f.write(">chr1\n" + "".join(rng.choice("ACGT") for _ in range(CHROM_LEN)) + "\n")

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def run_script(self, body):
        vcf = self.path("in.vcf")
        with open(vcf, "w") as f:
            f.write(HEADER + body)
        out, info = self.path("out"), self.path("info.tab")
        r = subprocess.run([sys.executable, SRC, "-v", vcf, "-r", self.ref, "-o", out, "-i", info],
                           capture_output=True, text=True)
        return r, out, info

    def test_out_of_contig_pos_skipped(self):
        # cuteSV writes POS=0/END=0 for an INS at a contig start (cuteSV.INS.9918)
        r, out, info = self.run_script(
            ins(0, "ins_pos0") + dele(0, "del_pos0") + ins(500, "ins_ok") + ins(CHROM_LEN + 1, "ins_n1"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        for vid in ("ins_pos0", "del_pos0", "ins_n1"):
            self.assertIn(f"code: 16) for record (ID:{vid})", r.stdout)
        with open(os.path.join(out, "extract_sv.summary")) as f:
            summary = f.read()
        self.assertIn("3 : POS outside contig (telomeric POS 0 or N+1) (code: 16)", summary)
        self.assertIn("1 : proper SV INS (code: 1)", summary)
        with open(info) as f:
            rows = [l.split("\t") for l in f.read().splitlines()[1:]]
        self.assertEqual([row[7] for row in rows], ["ins_ok"])

    def test_contig_edges_still_processed(self):
        r, _, info = self.run_script(ins(1, "ins_first") + ins(CHROM_LEN, "ins_last"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(info) as f:
            rows = [l.split("\t") for l in f.read().splitlines()[1:]]
        self.assertEqual(sorted(row[7] for row in rows), ["ins_first", "ins_last"])


if __name__ == "__main__":
    unittest.main()
