#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import sys
from typing import Dict, List, Tuple

try:
    import pysam
except ImportError:
    sys.stderr.write("Error: pysam is required. Install with: pip install pysam\n")
    sys.exit(1)

INFO_HEADER_RE = re.compile(r'^##INFO=<ID=(?P<ID>[^,]+),')


def parse_args():
    ap = argparse.ArgumentParser(
        description="Annotate a VCF with repeat-based INFO tags and optional BND mate annotations."
    )
    ap.add_argument("--vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument(
        "--header",
        required=True,
        help=(
            "File containing INFO header lines to append (e.g., repeat_info.txt). "
            "Ensure it defines BND_MATE_INFO with Number=1,Type=String."
        ),
    )
    ap.add_argument(
        "--tsv",
        required=True,
        help="Primary annotations TSV (columns include ID and tag values)",
    )
    ap.add_argument(
        "--mate",
        required=False,
        default="",
        help=(
            "BND mate annotations TSV (optional, can be empty). Values will be serialized into "
            "BND_MATE_INFO as a single scalar pipe-separated blob: key=v1,v2|key2=u1,u2|..."
        ),
    )
    ap.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output VCF path. Must end with .vcf.",
    )
    ap.add_argument(
        "--write_method",
        choices=["pysam", "pysam_str", "manual_str"],
        default="manual_str",
        help=(
            "Method to write VCF records. "
            "'pysam' uses VariantFile.write(rec) — may drop INFO/END for redundant records; "
            "'pysam_str' modifies rec via pysam then writes str(rec) — END may still be dropped after mutation; "
            "'manual_str' processes records as raw text without pysam parsing — fully preserves INFO/END."
        ),
    )
    return ap.parse_args()


def open_text(path: str):
    if not path or not os.path.exists(path):
        return None
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_info_ids_from_header_file(header_file: str) -> List[str]:
    ids: List[str] = []
    with open(header_file, "rt") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("##INFO="):
                continue
            m = INFO_HEADER_RE.match(line)
            if m:
                ids.append(m.group("ID"))
    return ids


def write_header_only_vcf(input_vcf: str, header_only_path: str) -> None:
    """
    Create a header-only VCF (.vcf) by letting pysam write the header.
    No variant records are written.
    """
    with pysam.VariantFile(input_vcf, "r") as invcf:
        hdr = invcf.header.copy()
    with pysam.VariantFile(header_only_path, "w", header=hdr):
        pass


def strip_overridden_info_ids_from_header_text(
    header_vcf_path: str, ids_to_remove: List[str], cleaned_path: str
) -> None:
    """
    Read header-only VCF as text, remove lines for INFO IDs in ids_to_remove,
    and write to cleaned_path.
    """
    with open(header_vcf_path, "rt") as f:
        lines = f.readlines()

    id_set = set(ids_to_remove)
    kept: List[str] = []
    for line in lines:
        ls = line.strip()
        if ls.startswith("##INFO="):
            m = INFO_HEADER_RE.match(ls)
            if m:
                info_id = m.group("ID")
                if info_id in id_set:
                    sys.stderr.write(
                        f"Warning: INFO ID '{info_id}' overridden by provided header\n"
                    )
                    continue
        kept.append(line)

    with open(cleaned_path, "wt") as out:
        out.writelines(kept)


def create_temp_vcf_with_new_header(
    input_vcf: str, header_file: str, output_path: str
) -> Tuple[str, pysam.VariantHeader]:
    """
    Header-editing pipeline:
    1) Write a header-only VCF from the input header to {output}.temp0.vcf.
    2) Remove overriding INFO IDs to {output}.temp1.vcf.
    3) Load cleaned header via pysam, then add the new INFO lines.
    4) Write header to {output}.temp2.vcf, then append raw input records as text
       (no pysam parsing) so INFO/END is fully preserved.
    Returns (temp2_vcf_path, out_hdr).
    """
    header_only_path = f"{output_path}.temp0.vcf"
    cleaned_hdr_path = f"{output_path}.temp1.vcf"
    temp2_vcf = f"{output_path}.temp2.vcf"

    # Step 1: header-only file from input header
    write_header_only_vcf(input_vcf, header_only_path)

    # Step 2: remove overridden INFO IDs
    override_ids = parse_info_ids_from_header_file(header_file)
    strip_overridden_info_ids_from_header_text(
        header_only_path, override_ids, cleaned_hdr_path
    )

    # Step 3: load cleaned header via pysam and add new INFO lines
    with pysam.VariantFile(cleaned_hdr_path, "r") as cleaned_vcf:
        out_hdr = cleaned_vcf.header.copy()
    with open(header_file, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    out_hdr.add_line(line)
                except ValueError as e:
                    sys.stderr.write(
                        f"Warning: failed to add header line '{line}': {e}\n"
                    )

    # Step 4: write header via pysam, then append raw records as text
    with pysam.VariantFile(temp2_vcf, "w", header=out_hdr):
        pass

    opener = gzip.open if input_vcf.endswith(".gz") else open
    with opener(input_vcf, "rt") as in_fh, open(temp2_vcf, "a") as out_fh:
        for line in in_fh:
            if line.startswith("#"):
                continue
            out_fh.write(line)

    # Cleanup header-only intermediates
    for p in (header_only_path, cleaned_hdr_path):
        try:
            os.remove(p)
        except OSError:
            pass

    return temp2_vcf, out_hdr


def load_tsv_primary(path: str) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """
    Load primary annotations TSV into a dict: id -> {tag: value (string raw from TSV)}.
    Returns (map, tag_order) where tag_order are column names excluding CHROM, POS, ID.
    """
    m: Dict[str, Dict[str, str]] = {}
    tags_order: List[str] = []
    if not path:
        return m, tags_order
    f = open_text(path)
    if f is None:
        return m, tags_order
    header = None
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            tags_order = [c for c in header if c not in ("CHROM", "POS", "ID")]
            continue
        rec = dict(zip(header, cols))
        vid = rec.get("ID")
        if not vid:
            continue
        tagvals = {tag: rec.get(tag, "") for tag in tags_order}
        m[vid] = tagvals
    f.close()
    return m, tags_order


def load_tsv_bnd_mate_blob_scalar(path: str) -> Dict[str, str]:
    """
    Build a single scalar blob for BND_MATE_INFO:
    "key=v1,v2|key2=u1,u2|..."
    """
    m: Dict[str, str] = {}
    if not path:
        return m
    f = open_text(path)
    if f is None:
        return m
    header = None
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            continue
        rec = dict(zip(header, cols))
        vid = rec.get("ID")
        if not vid:
            continue
        groups: List[str] = []
        for k, v in rec.items():
            if k in ("CHROM", "POS", "ID"):
                continue
            if v == "." or v == "":
                continue
            groups.append(f"{k}={v}")
        if groups:
            m[vid] = "|".join(groups)
    f.close()
    return m


def cast_info_value_by_effective_header(
    tag: str, value: str, hdr: pysam.VariantHeader
):
    """
    Cast values using the effective header attached to the current VariantFile.
    Return None for '.' or empty.
    """
    if value == "." or value == "":
        return None

    if tag not in hdr.info:
        return None

    info_def = hdr.info[tag]
    num = info_def.number
    typ = info_def.type  # 'Integer', 'Float', 'String', 'Flag'

    def cast_one(x: str):
        if typ == "Float":
            try:
                return float(x)
            except ValueError:
                return x
        elif typ == "Integer":
            try:
                return int(x)
            except ValueError:
                return x
        else:
            return x

    if num == ".":
        parts = value.split(",")
        return tuple(cast_one(p) for p in parts if p != "")
    else:
        return cast_one(value)


def annotate_pysam_record(
    rec,
    primary_map: Dict[str, Dict[str, str]],
    primary_tags: List[str],
    bnd_scalar_map: Dict[str, str],
    out_hdr: pysam.VariantHeader,
):
    """Annotate a pysam VariantRecord in-place with primary and BND mate tags."""
    vid = rec.id
    if vid and vid in primary_map:
        tagvals = primary_map[vid]
        for tag in primary_tags:
            raw_val = tagvals.get(tag, "")
            val = cast_info_value_by_effective_header(tag, raw_val, out_hdr)
            if val is None:
                continue
            try:
                rec.info[tag] = val
            except KeyError:
                sys.stderr.write(
                    f"Warning: skipping tag '{tag}' not present in header\n"
                )

    if vid and vid in bnd_scalar_map:
        blob = bnd_scalar_map[vid]
        if blob:
            try:
                rec.info["BND_MATE_INFO"] = blob
            except KeyError:
                sys.stderr.write(
                    "Warning: BND_MATE_INFO not present in header; "
                    "ensure Number=1,Type=String in --header\n"
                )


def annotate_vcf(
    vcf_path: str,
    info_header_path: str,
    annotations_tsv: str,
    bnd_mate_tsv: str,
    output_path: str,
    write_method: str = "manual_str",
):
    """
    Two-pass approach:
    - Pass 1: create a temp VCF with updated INFO header (overrides applied via
      text header edit) and raw-text records (END preserved).
    - Pass 2: reopen the temp VCF, assign annotations, and stream records to the
      final output VCF using the chosen write method.
    """
    # Pass 1
    temp_vcf, out_hdr = create_temp_vcf_with_new_header(
        vcf_path, info_header_path, output_path
    )

    # Load annotation maps
    primary_map, primary_tags = load_tsv_primary(annotations_tsv)
    bnd_scalar_map = (
        load_tsv_bnd_mate_blob_scalar(bnd_mate_tsv) if bnd_mate_tsv else {}
    )

    # Pass 2 — each branch writes its own header + records
    if write_method == "pysam":
        with pysam.VariantFile(temp_vcf, "r") as invcf, \
             pysam.VariantFile(output_path, "w", header=out_hdr) as outvcf:
            for rec in invcf:
                annotate_pysam_record(
                    rec, primary_map, primary_tags, bnd_scalar_map, out_hdr
                )
                outvcf.write(rec)

    elif write_method == "pysam_str":
        # Write header via pysam, then append string-serialised records
        with pysam.VariantFile(output_path, "w", header=out_hdr):
            pass
        with pysam.VariantFile(temp_vcf, "r") as invcf, \
             open(output_path, "a") as outvcf:
            for rec in invcf:
                annotate_pysam_record(
                    rec, primary_map, primary_tags, bnd_scalar_map, out_hdr
                )
                # WARNING: pysam_str still goes through htslib's internal
                # serialization. If rec.info is modified (as we do here), END
                # may still be dropped for records where
                # END == POS + len(REF) - 1. Only manual_str fully avoids this.
                # str(rec) from pysam already includes trailing newline.
                outvcf.write(str(rec))

    elif write_method == "manual_str":
        # Write header via pysam, then append text-processed records.
        # Records are never parsed by pysam/htslib, so INFO/END is fully preserved.
        with pysam.VariantFile(output_path, "w", header=out_hdr):
            pass
        with open(temp_vcf, "rt") as invcf, \
             open(output_path, "a") as outvcf:
            for line in invcf:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\n")
                cols = line.split("\t")
                if len(cols) < 8:
                    print("Warning: skipping malformed VCF line (expected at least 8 columns): " + line, file=sys.stderr)
                    continue
                chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
                rest = cols[8:] if len(cols) > 8 else []

                # Determine which tags need stripping (dedup on re-run)
                tags_to_strip = set()
                if vid and vid in primary_map:
                    tags_to_strip.update(primary_tags)
                if vid and vid in bnd_scalar_map:
                    tags_to_strip.add("BND_MATE_INFO")

                if tags_to_strip:
                    info_fields = info.split(";")
                    info_fields = [
                        f
                        for f in info_fields
                        if f.split("=", 1)[0] not in tags_to_strip
                    ]
                    info = ";".join(info_fields)

                # Append primary annotations
                if vid and vid in primary_map:
                    tagvals = primary_map[vid]
                    for tag in primary_tags:
                        raw_val = tagvals.get(tag, "")
                        if raw_val == "" or raw_val == ".":
                            continue
                        if tag not in out_hdr.info:
                            continue
                        if info == "" or info == ".":
                            info = f"{tag}={raw_val}"
                        else:
                            info += f";{tag}={raw_val}"

                # Append BND mate annotation
                if vid and vid in bnd_scalar_map:
                    blob = bnd_scalar_map[vid]
                    if blob:
                        if info == "" or info == ".":
                            info = f"BND_MATE_INFO={blob}"
                        else:
                            info += f";BND_MATE_INFO={blob}"

                record_str = "\t".join(
                    [chrom, pos, vid, ref, alt, qual, filt, info] + rest
                )
                outvcf.write(record_str + "\n")

    # Cleanup temp2
    try:
        os.remove(temp_vcf)
    except OSError:
        pass


def main():
    args = parse_args()
    if not args.output.endswith(".vcf"):
        sys.stderr.write("Error: Output path must end with .vcf\n")
        sys.exit(2)
    print(f"Info: VCF File: {args.vcf}")
    print(f"Info: Header File: {args.header}")
    print(f"Info: Annotations TSV: {args.tsv}")
    print(f"Info: BND Mate TSV: {args.mate}")
    print(f"Info: Output VCF: {args.output}")
    print(f"Info: Write method: {args.write_method}")
    annotate_vcf(
        vcf_path=args.vcf,
        info_header_path=args.header,
        annotations_tsv=args.tsv,
        bnd_mate_tsv=args.mate,
        output_path=args.output,
        write_method=args.write_method,
    )
    sys.stderr.write(f"Written annotated VCF: {args.output}\n")


if __name__ == "__main__":
    main()