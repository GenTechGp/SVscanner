"""Tests for src/check_vcf_ids.py. Run: python -m unittest test/test_check_vcf_ids.py"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "check_vcf_ids.py")
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def rec(chrom, pos, vid, svtype="DEL", extra=""):
    return f"{chrom}\t{pos}\t{vid}\tN\t<{svtype}>\t.\tPASS\tSVTYPE={svtype};END={pos + 100}{extra}\n"


class CheckVcfIds(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def write(self, body, name="in.vcf", gz=False):
        p = self.path(name)
        with (gzip.open(p, "wt") if gz else open(p, "wt")) as f:
            f.write(body)
        return p

    def run_script(self, *args):
        r = subprocess.run([sys.executable, SRC, *args], capture_output=True, text=True)
        return r.returncode, r.stdout, r.stderr

    def body_lines(self, p):
        opener = gzip.open if p.endswith(".gz") else open
        with opener(p, "rt") as f:
            return [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]

    def test_clean(self):
        p = self.write(HEADER + rec("chr1", 10, "a") + rec("chr1", 20, "b"))
        self.assertEqual(self.run_script("--vcf", p)[0], 0)

    def test_header_only(self):
        p = self.write(HEADER)
        self.assertEqual(self.run_script("--vcf", p)[0], 0)

    def test_all_missing(self):
        p = self.write(HEADER + rec("chr1", 10, ".") + rec("chr1", 20, "."))
        self.assertEqual(self.run_script("--vcf", p)[0], 1)

    def test_partial_missing_and_duplicates_rewrite(self):
        body = HEADER + rec("chr1", 10, "a") + rec("chr1", 20, ".") + rec("chr1", 30, "a") + rec("chr1", 40, "b", "INS")
        p = self.write(body)
        out, mp = self.path("out.vcf"), self.path("map.tsv")
        code, _, err = self.run_script("--vcf", p, "--mode", "rewrite", "--out", out, "--map", mp)
        self.assertEqual(code, 0)
        self.assertIn("Warning", err)  # duplicate renamed
        ids = [c[2] for c in self.body_lines(out)]
        self.assertEqual(ids, ["a", "SVSCANNER_chr1_20_DEL_2", "SVSCANNER_chr1_30_DEL_3", "b"])
        # everything but the ID column is untouched
        for old, new in zip(self.body_lines(p), self.body_lines(out)):
            self.assertEqual(old[:2] + old[3:], new[:2] + new[3:])
        self.assertEqual(self.run_script("--vcf", out)[0], 0)
        with open(mp) as f:
            self.assertEqual(len(f.readlines()), 3)  # header + 2 changes

    def test_rewrite_clean_is_noop(self):
        p = self.write(HEADER + rec("chr1", 10, "a"))
        out = self.path("out.vcf")
        code, stdout, _ = self.run_script("--vcf", p, "--mode", "rewrite", "--out", out)
        self.assertEqual(code, 0)
        self.assertFalse(os.path.exists(out))
        self.assertIn("nothing rewritten", stdout)

    def test_whitespace_unfixable(self):
        p = self.write(HEADER + rec("chr1", 10, "a b") + rec("chr1", 20, "."))
        out = self.path("out.vcf")
        self.assertEqual(self.run_script("--vcf", p)[0], 2)
        self.assertEqual(self.run_script("--vcf", p, "--mode", "rewrite", "--out", out)[0], 2)
        self.assertFalse(os.path.exists(out))

    def test_mateid_with_bad_ids_unfixable(self):
        body = HEADER + rec("chr1", 10, ".", "BND", ";MATEID=m2") + rec("chr1", 20, ".", "BND", ";MATEID=m1")
        p = self.write(body)
        out = self.path("out.vcf")
        self.assertEqual(self.run_script("--vcf", p, "--mode", "rewrite", "--out", out)[0], 2)
        self.assertFalse(os.path.exists(out))

    def test_mateid_with_good_ids_ok(self):
        body = HEADER + rec("chr1", 10, "m1", "BND", ";MATEID=m2") + rec("chr1", 20, "m2", "BND", ";MATEID=m1")
        self.assertEqual(self.run_script("--vcf", self.write(body))[0], 0)

    def test_gzip_in_and_out(self):
        p = self.write(HEADER + rec("chr1", 10, ".") + rec("chr1", 20, "."), "in.vcf.gz", gz=True)
        out = self.path("out.vcf.gz")
        self.assertEqual(self.run_script("--vcf", p, "--mode", "rewrite", "--out", out)[0], 0)
        self.assertEqual(self.run_script("--vcf", out)[0], 0)
        self.assertEqual(len({c[2] for c in self.body_lines(out)}), 2)

    def test_generated_id_collision(self):
        body = HEADER + rec("chr1", 10, "SVSCANNER_chr1_20_DEL_2") + rec("chr1", 20, ".")
        p = self.write(body)
        out = self.path("out.vcf")
        self.assertEqual(self.run_script("--vcf", p, "--mode", "rewrite", "--out", out)[0], 0)
        ids = [c[2] for c in self.body_lines(out)]
        self.assertEqual(len(set(ids)), 2)

    def test_missing_file(self):
        self.assertEqual(self.run_script("--vcf", self.path("nope.vcf"))[0], 3)


if __name__ == "__main__":
    unittest.main()
