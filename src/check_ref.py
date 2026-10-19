#!/usr/bin/env python3
"""
Check that a reference FASTA is the one an SV VCF was called against.

SVscanner fetches DEL, DUP, INV and BND sequences, and the flanks of every SV, from the
reference. Given the wrong reference the run still completes, but those records are
annotated from unrelated sequence. Two checks, in one pass over the VCF:

  contigs  every contig the records use (CHROM, and BND mates in ALT) must be in the
           FASTA, and every ##contig length in the header must match the FASTA length
  REF      a random sample of REF alleles must match the FASTA at their position.
           Multi-base REFs (e.g. sequence-resolved DELs) are used when the VCF has
           enough of them; otherwise single-base REFs are added. 'N' never mismatches.

Records are read as text, never through pysam/htslib; pysam is used only for the FASTA.

Exit codes:
  0  the reference matches the VCF (or the VCF has nothing the REF check can use)
  1  the reference does not match the VCF
  3  the input could not be read
"""

import argparse
import os
import random
import re
import sys
from typing import Dict, List, Optional, Tuple

from check_vcf_ids import MAX_EXAMPLES, examples, iter_records, open_text

EXIT_OK = 0
EXIT_MISMATCH = 1
EXIT_IO = 3

SAMPLE_SIZE = 200
MIN_LONG_REFS = 20          # fewer multi-base REFs than this: add single-base REFs
MAX_MISMATCH_FRACTION = 0.05
SHOW_BASES = 20

CONTIG_ID_RE = re.compile(r"[<,]ID=([^,>]+)")
CONTIG_LENGTH_RE = re.compile(r"[<,]length=(\d+)")
BND_MATE_RE = re.compile(r"[\[\]]([^\[\]:]+):\d+[\[\]]")
REF_RE = re.compile(r"^[ACGTNRYSWKMBDHV]+$")

# (file_line_number, ID, CHROM, POS, REF)
Sample = Tuple[int, str, str, int, str]


def header_contig_lengths(path: str) -> Dict[str, int]:
    lengths = {}
    with open_text(path) as f:
        for line in f:
            if not line.startswith("##"):
                break
            if not line.startswith("##contig=<"):
                continue
            cid, clen = CONTIG_ID_RE.search(line), CONTIG_LENGTH_RE.search(line)
            if cid and clen:
                lengths[cid.group(1)] = int(clen.group(1))
    return lengths


class Reservoir:
    """Uniform random sample of a stream (Algorithm R), seeded so runs are repeatable."""

    def __init__(self, size: int, rng: random.Random):
        self.size, self.rng = size, rng
        self.seen = 0
        self.items: List[Sample] = []

    def add(self, item: Sample) -> None:
        self.seen += 1
        if len(self.items) < self.size:
            self.items.append(item)
        else:
            j = self.rng.randrange(self.seen)
            if j < self.size:
                self.items[j] = item


class Report:
    def __init__(self):
        self.records = 0
        self.used_contigs: Dict[str, int] = {}      # contig -> record count
        self.missing_contigs: List[str] = []        # used, but not in the FASTA
        self.length_mismatches: List[Tuple[str, int, int]] = []   # (contig, header, fasta)
        self.ref_checked = 0
        self.ref_kind = ""
        self.ref_mismatches: List[Tuple[Sample, str]] = []        # (sample, fasta seq)
        self.fasta_contigs: List[str] = []

    @property
    def mismatch_fraction(self) -> float:
        return len(self.ref_mismatches) / self.ref_checked if self.ref_checked else 0.0

    @property
    def status(self) -> int:
        if self.missing_contigs or self.length_mismatches:
            return EXIT_MISMATCH
        if self.mismatch_fraction > MAX_MISMATCH_FRACTION:
            return EXIT_MISMATCH
        return EXIT_OK


def same_bases(vcf_ref: str, fasta_seq: str) -> bool:
    if len(vcf_ref) != len(fasta_seq):
        return False
    return all(a == b or a == "N" or b == "N" for a, b in zip(vcf_ref, fasta_seq))


def check_reference(vcf_path: str, fasta, sample_size: int = SAMPLE_SIZE) -> Report:
    report = Report()
    fasta_lengths = dict(zip(fasta.references, fasta.lengths))
    report.fasta_contigs = list(fasta.references)

    for contig, header_len in header_contig_lengths(vcf_path).items():
        if contig in fasta_lengths and fasta_lengths[contig] != header_len:
            report.length_mismatches.append((contig, header_len, fasta_lengths[contig]))

    rng = random.Random(0)
    long_refs, short_refs = Reservoir(sample_size, rng), Reservoir(sample_size, rng)
    for line_number, _, _, cols in iter_records(vcf_path):
        if len(cols) < 8:
            continue
        report.records += 1
        chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3].upper(), cols[4]
        report.used_contigs[chrom] = report.used_contigs.get(chrom, 0) + 1
        if "[" in alt or "]" in alt:
            for mate in BND_MATE_RE.findall(alt):
                report.used_contigs.setdefault(mate, 0)
        if chrom not in fasta_lengths or not pos.isdigit():
            continue
        if not REF_RE.match(ref) or set(ref) == {"N"}:
            continue
        (long_refs if len(ref) > 1 else short_refs).add((line_number, vid, chrom, int(pos), ref))

    report.missing_contigs = [c for c in report.used_contigs if c not in fasta_lengths]

    if len(long_refs.items) >= MIN_LONG_REFS:
        sample, report.ref_kind = long_refs.items, "multi-base REF alleles"
    else:
        sample = long_refs.items + short_refs.items[: max(0, sample_size - len(long_refs.items))]
        report.ref_kind = "REF alleles (mostly single-base)" if short_refs.items else "multi-base REF alleles"

    for item in sample:
        _, _, chrom, pos, ref = item
        fasta_seq = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
        report.ref_checked += 1
        if not same_bases(ref, fasta_seq):
            report.ref_mismatches.append((item, fasta_seq))
    return report


def naming_hint(report: Report) -> Optional[str]:
    """Suggest a chr-prefix mix-up when the missing contigs exist under the other convention."""
    fasta = set(report.fasta_contigs)
    flipped = [c[3:] if c.startswith("chr") else "chr" + c for c in report.missing_contigs]
    if flipped and all(f in fasta for f in flipped):
        return "the VCF and the reference use different contig names ('chr1' vs '1')"
    return None


def short(seq: str) -> str:
    return seq if len(seq) <= SHOW_BASES else seq[:SHOW_BASES] + f"...({len(seq)} bp)"


def print_report(vcf_path: str, ref_path: str, report: Report) -> None:
    print(f"Reference check: {ref_path}")
    print(f"  VCF:               {vcf_path}")
    print(f"  records:           {report.records}")
    print(f"  contigs used:      {len(report.used_contigs)}")
    if report.missing_contigs:
        print(f"  NOT IN REFERENCE:  {examples(report.missing_contigs, str)}")
        hint = naming_hint(report)
        if hint:
            print(f"  hint:              {hint}")
    if report.length_mismatches:
        print(
            "  LENGTH MISMATCH:   "
            + examples(report.length_mismatches, lambda m: f"{m[0]} (VCF header {m[1]}, reference {m[2]})")
        )
    if report.ref_checked:
        print(
            f"  REF check:         {len(report.ref_mismatches)}/{report.ref_checked} sampled "
            f"{report.ref_kind} differ from the reference ({report.mismatch_fraction:.1%}; "
            f"limit {MAX_MISMATCH_FRACTION:.0%})"
        )
        for (line_number, vid, chrom, pos, ref), fasta_seq in report.ref_mismatches[:MAX_EXAMPLES]:
            print(f"    {vid} {chrom}:{pos} (line {line_number}): VCF {short(ref)}, reference {short(fasta_seq)}")
        if len(report.ref_mismatches) > MAX_EXAMPLES:
            print(f"    ... (+{len(report.ref_mismatches) - MAX_EXAMPLES} more)")
    else:
        print("  REF check:         skipped, no record has a usable REF allele")
    if report.status == EXIT_OK:
        print("  result: the reference matches the VCF")
    else:
        print("  result: the reference does NOT match the VCF")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Check that a reference FASTA matches an SV VCF (see module docstring).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--vcf", required=True, help="Input VCF (plain or gzip/bgzip)")
    p.add_argument("--ref", required=True, help="Reference FASTA (plain or bgzip, indexed or indexable)")
    p.add_argument("--sample", type=int, default=SAMPLE_SIZE, help="REF alleles to check (default: %(default)s)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    for label, path in (("VCF", args.vcf), ("reference", args.ref)):
        if not os.path.isfile(path):
            print(f"Error: {label} not found: {path}", file=sys.stderr)
            return EXIT_IO

    try:
        import pysam
        fasta = pysam.FastaFile(args.ref)
    except (OSError, ValueError) as e:
        print(f"Error: could not open the reference {args.ref}: {e}", file=sys.stderr)
        return EXIT_IO

    try:
        report = check_reference(args.vcf, fasta, args.sample)
    except (OSError, EOFError) as e:
        print(f"Error: could not read {args.vcf}: {e}", file=sys.stderr)
        return EXIT_IO
    finally:
        fasta.close()
    print_report(args.vcf, args.ref, report)
    return report.status


if __name__ == "__main__":
    sys.exit(main())
