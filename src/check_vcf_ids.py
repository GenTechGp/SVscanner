#!/usr/bin/env python3
"""
Check, and optionally repair, the ID column of an SV VCF.

SVscanner joins its annotations back onto the VCF by the ID column, so every record
needs an ID that is present and unique. This script has two modes:

  check    Read the VCF and report missing IDs, duplicated IDs, IDs containing
           whitespace and records carrying MATEID. Never writes the VCF.
  rewrite  Stream the VCF to --out, changing only the ID column: existing unique IDs
           (and the first occurrence of each duplicated ID) are kept, every missing ID
           and every later duplicate gets SVSCANNER_<CHROM>_<POS>_<SVTYPE>_<n>.

Records are handled as text, never through pysam/htslib, so INFO/END and every other
column are passed through verbatim.

Exit codes:
  0  the ID column is usable (check), or the output was written/was not needed (rewrite)
  1  the ID column is not usable but can be repaired with --mode rewrite
  2  the ID column is not usable and must not be repaired: whitespace in an ID, or
     MATEID present - rewriting IDs would leave the MATEID references dangling
  3  the input could not be read or written
"""

import argparse
import gzip
import io
import os
import re
import sys
from typing import Iterator, List, Optional, Set, Tuple

EXIT_OK = 0
EXIT_FIXABLE = 1
EXIT_UNFIXABLE = 2
EXIT_IO = 3

ID_PREFIX = "SVSCANNER"
MAX_EXAMPLES = 5
WHITESPACE_RE = re.compile(r"\s")
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
MATEID_RE = re.compile(r"(?:^|;)MATEID=")


def is_gzipped(path: str) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def open_text(path: str):
    return gzip.open(path, "rt") if is_gzipped(path) else open(path, "rt")


def open_output(path: str):
    """bgzip when the output name asks for compression, so downstream tools can index it."""
    if not path.endswith(".gz"):
        return open(path, "wt")
    try:
        import pysam
        raw = pysam.BGZFile(path, "wb")
    except (ImportError, AttributeError):
        raw = gzip.open(path, "wb")
    return io.TextIOWrapper(raw, encoding="utf-8", newline="")


def iter_records(path: str) -> Iterator[Tuple[int, int, str, List[str]]]:
    """Yield (file_line_number, record_index, raw_line, columns) for each data line.

    record_index is 1-based over data lines only. Malformed lines (fewer than 8
    columns) are yielded with whatever columns they have so callers can count them."""
    record_index = 0
    with open_text(path) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line.strip():
                continue
            record_index += 1
            yield line_number, record_index, line, line.split("\t")


def is_missing(vid: str) -> bool:
    return vid == "." or vid == ""


def svtype_of(info: str) -> str:
    m = SVTYPE_RE.search(info)
    return m.group(1) if m else "SV"


class Report:
    """Everything the check learns about the ID column in one pass."""

    def __init__(self):
        self.records = 0
        self.malformed = 0
        self.missing: List[int] = []        # file line numbers
        self.duplicates: List[Tuple[int, str]] = []   # (line, id) of later occurrences
        self.whitespace: List[Tuple[int, str]] = []
        self.mateid: List[int] = []
        self.ids: Set[str] = set()

    @property
    def bad_ids(self) -> bool:
        return bool(self.missing or self.duplicates)

    @property
    def unfixable_reasons(self) -> List[str]:
        reasons = []
        if self.whitespace:
            reasons.append("an ID contains whitespace, which a VCF does not allow")
        if self.bad_ids and self.mateid:
            reasons.append(
                "missing or duplicated IDs in a VCF that uses MATEID - changing IDs "
                "would leave the MATEID references pointing at the wrong records"
            )
        return reasons

    @property
    def status(self) -> int:
        if self.unfixable_reasons:
            return EXIT_UNFIXABLE
        return EXIT_FIXABLE if self.bad_ids else EXIT_OK


def check_vcf(path: str) -> Report:
    report = Report()
    seen: Set[str] = set()
    for line_number, _, _, cols in iter_records(path):
        report.records += 1
        if len(cols) < 8:
            report.malformed += 1
            continue
        vid, info = cols[2], cols[7]
        if MATEID_RE.search(info):
            report.mateid.append(line_number)
        if is_missing(vid):
            report.missing.append(line_number)
            continue
        if WHITESPACE_RE.search(vid):
            report.whitespace.append((line_number, vid))
        if vid in seen:
            report.duplicates.append((line_number, vid))
        else:
            seen.add(vid)
    report.ids = seen
    return report


def examples(items, fmt) -> str:
    shown = ", ".join(fmt(i) for i in items[:MAX_EXAMPLES])
    more = len(items) - MAX_EXAMPLES
    return shown + (f", ... (+{more} more)" if more > 0 else "")


def print_report(path: str, report: Report) -> None:
    print(f"ID check: {path}")
    print(f"  records:           {report.records}")
    print(f"  missing IDs:       {len(report.missing)}")
    print(f"  duplicated IDs:    {len(report.duplicates)} (later occurrences)")
    print(f"  IDs with spaces:   {len(report.whitespace)}")
    print(f"  records w/ MATEID: {len(report.mateid)}")
    if report.malformed:
        print(f"  malformed lines:   {report.malformed} (fewer than 8 columns; not checked)")
    if report.missing:
        print(f"  missing at lines:  {examples(report.missing, str)}")
    if report.duplicates:
        print(f"  duplicates:        {examples(report.duplicates, lambda d: f'{d[1]} (line {d[0]})')}")
    if report.whitespace:
        print(f"  whitespace IDs:    {examples(report.whitespace, lambda d: f'{d[1]!r} (line {d[0]})')}")
    for reason in report.unfixable_reasons:
        print(f"  CANNOT REPAIR: {reason}", file=sys.stderr)
    if report.status == EXIT_OK:
        print("  result: ID column is usable")
    elif report.status == EXIT_FIXABLE:
        print("  result: ID column is not usable; it can be repaired (--mode rewrite)")
    else:
        print("  result: ID column is not usable and cannot be repaired")


def make_id(chrom: str, pos: str, info: str, record_index: int, taken: Set[str]) -> str:
    """Build SVSCANNER_<CHROM>_<POS>_<SVTYPE>_<n>. <n> is the record index, so IDs are
    unique among themselves; the loop only guards against a pre-existing ID of the
    same shape."""
    base = f"{ID_PREFIX}_{chrom}_{pos}_{svtype_of(info)}_{record_index}"
    candidate, bump = base, 0
    while candidate in taken:
        bump += 1
        candidate = f"{base}_{bump}"
    return candidate


def rewrite_vcf(path: str, out_path: str, map_path: Optional[str], report: Report) -> int:
    taken = set(report.ids)
    seen: Set[str] = set()
    changes: List[Tuple[int, str, str]] = []
    renamed = filled = 0

    out = open_output(out_path)
    try:
        with open_text(path) as f:
            record_index = 0
            for line in f:
                if line.startswith("#") or not line.strip():
                    out.write(line)
                    continue
                record_index += 1
                had_newline = line.endswith("\n")
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 8:
                    vid = cols[2]
                    if is_missing(vid) or vid in seen:
                        new_id = make_id(cols[0], cols[1], cols[7], record_index, taken)
                        taken.add(new_id)
                        if is_missing(vid):
                            filled += 1
                        else:
                            renamed += 1
                        changes.append((record_index, vid, new_id))
                        cols[2] = new_id
                    seen.add(cols[2])
                out.write("\t".join(cols) + ("\n" if had_newline else ""))
    finally:
        out.close()

    if map_path:
        with open(map_path, "wt") as m:
            m.write("RECORD_INDEX\tOLD_ID\tNEW_ID\n")
            for record_index, old, new in changes:
                m.write(f"{record_index}\t{old}\t{new}\n")

    print(f"ID rewrite: filled {filled} missing ID(s), renamed {renamed} duplicated ID(s)")
    if renamed:
        print(
            f"Warning: {renamed} duplicated ID(s) were renamed (first occurrence kept): "
            + examples(
                [c for c in changes if c[1] != "." and c[1] != ""],
                lambda c: f"{c[1]} -> {c[2]}",
            ),
            file=sys.stderr,
        )
    print(f"  output:  {out_path}")
    if map_path:
        print(f"  mapping: {map_path} ({len(changes)} change(s))")
    return len(changes)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Check, or repair, the ID column of an SV VCF (see module docstring).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--vcf", required=True, help="Input VCF (plain or gzip/bgzip)")
    p.add_argument("--mode", choices=["check", "rewrite"], default="check")
    p.add_argument("--out", help="rewrite mode: output VCF; a .gz name is bgzip-compressed")
    p.add_argument("--map", help="rewrite mode: write an old-ID -> new-ID table here")
    args = p.parse_args(argv)
    if args.mode == "rewrite" and not args.out:
        p.error("--mode rewrite requires --out")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    if not os.path.isfile(args.vcf):
        print(f"Error: VCF not found: {args.vcf}", file=sys.stderr)
        return EXIT_IO

    try:
        report = check_vcf(args.vcf)
    except (OSError, EOFError) as e:
        print(f"Error: could not read {args.vcf}: {e}", file=sys.stderr)
        return EXIT_IO
    print_report(args.vcf, report)

    if args.mode == "check":
        return report.status

    if report.status == EXIT_UNFIXABLE:
        return EXIT_UNFIXABLE
    if report.status == EXIT_OK:
        print("ID column is already usable; nothing rewritten")
        return EXIT_OK

    try:
        rewrite_vcf(args.vcf, args.out, args.map, report)
    except OSError as e:
        print(f"Error: could not write {args.out}: {e}", file=sys.stderr)
        return EXIT_IO
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
