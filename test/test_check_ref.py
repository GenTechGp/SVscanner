"""Tests for src/check_ref.py. Run from test/: python -m unittest test_check_ref"""
import gzip
import os
import random
import subprocess
import sys
import tempfile
import unittest

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "check_ref.py")
CONTIG_LEN = 20000


def random_seq(seed, n=CONTIG_LEN):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class CheckRef(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.chr1, self.chr2 = random_seq(1), random_seq(2)
        self.ref = self.write_fasta({"chr1": self.chr1, "chr2": self.chr2})

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def write_fasta(self, contigs, name="ref.fa"):
        p = self.path(name)
        with open(p, "wt") as f:
            for cid, seq in contigs.items():
                f.write(f">{cid}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")
        return p

    def write_vcf(self, records, contigs=(("chr1", CONTIG_LEN), ("chr2", CONTIG_LEN)), gz=False):
        p = self.path("in.vcf.gz" if gz else "in.vcf")
        header = "##fileformat=VCFv4.2\n"
        header += "".join(f"##contig=<ID={c},length={n}>\n" for c, n in contigs)
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        with (gzip.open(p, "wt") if gz else open(p, "wt")) as f:
            f.write(header + "".join(records))
        return p

    def dels(self, seq, n=30, chrom="chr1", length=60):
        """Sequence-resolved DELs whose REF is taken from seq."""
        out = []
        for i in range(n):
            pos = 100 + i * 500
            ref = seq[pos - 1:pos - 1 + length]
            out.append(f"{chrom}\t{pos}\tDEL{i}\t{ref}\t{ref[0]}\t.\tPASS\tSVTYPE=DEL;SVLEN=-{length - 1}\n")
        return out

    def symbolic(self, seq, n=30, chrom="chr1"):
        """Symbolic <INS> records with only the single anchor base as REF."""
        return [
            f"{chrom}\t{100 + i * 500}\tINS{i}\t{seq[99 + i * 500]}\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=300\n"
            for i in range(n)
        ]

    def run_script(self, vcf, ref=None):
        r = subprocess.run([sys.executable, SRC, "--vcf", vcf, "--ref", ref or self.ref], capture_output=True, text=True)
        return r.returncode, r.stdout + r.stderr

    def test_matching_reference(self):
        rc, out = self.run_script(self.write_vcf(self.dels(self.chr1)))
        self.assertEqual(rc, 0, out)
        self.assertIn("0/30 sampled multi-base REF", out)

    def test_matching_reference_gzipped_vcf(self):
        rc, out = self.run_script(self.write_vcf(self.dels(self.chr1), gz=True))
        self.assertEqual(rc, 0, out)

    def test_wrong_reference_same_contigs(self):
        other = self.write_fasta({"chr1": random_seq(9), "chr2": self.chr2}, name="other.fa")
        rc, out = self.run_script(self.write_vcf(self.dels(self.chr1)), other)
        self.assertEqual(rc, 1, out)
        self.assertIn("does NOT match", out)
        self.assertIn("DEL0 chr1:100", out)

    def test_header_length_mismatch(self):
        vcf = self.write_vcf(self.dels(self.chr1), contigs=(("chr1", CONTIG_LEN + 5),))
        rc, out = self.run_script(vcf)
        self.assertEqual(rc, 1, out)
        self.assertIn(f"chr1 (VCF header {CONTIG_LEN + 5}, reference {CONTIG_LEN})", out)

    def test_no_contig_header_lines(self):
        rc, out = self.run_script(self.write_vcf(self.dels(self.chr1), contigs=()))
        self.assertEqual(rc, 0, out)

    def test_header_contig_absent_from_fasta_but_unused(self):
        vcf = self.write_vcf(self.dels(self.chr1), contigs=(("chr1", CONTIG_LEN), ("chrUn", 5)))
        rc, out = self.run_script(vcf)
        self.assertEqual(rc, 0, out)

    def test_used_contig_missing_with_naming_hint(self):
        recs = [r.replace("chr1\t", "1\t", 1) for r in self.dels(self.chr1)]
        rc, out = self.run_script(self.write_vcf(recs, contigs=()))
        self.assertEqual(rc, 1, out)
        self.assertIn("NOT IN REFERENCE:  1", out)
        self.assertIn("different contig names", out)

    def test_bnd_mate_contig_missing(self):
        recs = self.dels(self.chr1) + ["chr1\t5000\tBND1\tN\tN[chr9:100[\t.\tPASS\tSVTYPE=BND\n"]
        rc, out = self.run_script(self.write_vcf(recs))
        self.assertEqual(rc, 1, out)
        self.assertIn("NOT IN REFERENCE:  chr9", out)

    def test_single_base_refs_used_when_no_long_refs(self):
        vcf = self.write_vcf(self.symbolic(self.chr1))
        rc, out = self.run_script(vcf)
        self.assertEqual(rc, 0, out)
        self.assertIn("mostly single-base", out)
        other = self.write_fasta({"chr1": random_seq(9), "chr2": self.chr2}, name="other.fa")
        rc, out = self.run_script(vcf, other)
        self.assertEqual(rc, 1, out)

    def test_n_refs_are_skipped(self):
        recs = [f"chr1\t{100 + i}\tDEL{i}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-500\n" for i in range(10)]
        rc, out = self.run_script(self.write_vcf(recs))
        self.assertEqual(rc, 0, out)
        self.assertIn("skipped", out)

    def test_lowercase_reference_matches(self):
        soft = self.write_fasta({"chr1": self.chr1.lower(), "chr2": self.chr2}, name="soft.fa")
        rc, out = self.run_script(self.write_vcf(self.dels(self.chr1)), soft)
        self.assertEqual(rc, 0, out)

    def test_few_mismatches_tolerated(self):
        recs = self.dels(self.chr1, n=40)
        recs[0] = recs[0].replace("\tDEL0\t", "\tDEL0\tTTTTTTTTTT", 1)   # REF no longer matches
        rc, out = self.run_script(self.write_vcf(recs))
        self.assertEqual(rc, 0, out)
        self.assertIn("1/40", out)

    def test_missing_input(self):
        rc, out = self.run_script(self.path("nope.vcf"))
        self.assertEqual(rc, 3, out)


if __name__ == "__main__":
    unittest.main()
