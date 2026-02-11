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
        help="Output VCF path. Must end with .vcf or .vcf.gz to control compression.",
    )
    return ap.parse_args()


def ensure_valid_output_suffix(path: str):
    if path.endswith(".vcf") or path.endswith(".vcf.gz"):
        return
    sys.stderr.write("Error: Output path must end with .vcf or .vcf.gz\n")
    sys.exit(2)


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
    # Opening with header writes the header; closing without writing any records yields a header-only VCF
    with pysam.VariantFile(header_only_path, "w", header=hdr):
        pass


def strip_overridden_info_ids_from_header_text(header_vcf_path: str, ids_to_remove: List[str], cleaned_path: str) -> None:
    """
    Read header-only VCF as text, remove lines for INFO IDs in ids_to_remove, and write to cleaned_path.
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
                    sys.stderr.write(f"Warning: INFO ID '{info_id}' overridden by provided header\n")
                    continue  # skip this line
        kept.append(line)

    with open(cleaned_path, "wt") as out:
        out.writelines(kept)


def create_temp_vcf_with_new_header(input_vcf: str, header_file: str, output_path: str) -> str:
    """
    Header-editing pipeline without manual copying:
    1) Write a header-only VCF from the input header to {output}.temp0.vcf.
    2) Remove overriding INFO IDs to {output}.temp1.vcf.
    3) Load cleaned header via pysam, then add the new INFO lines.
    4) Use the resulting header to write {output}.temp2.(vcf|vcf.gz) duplicating input records unchanged.
    Returns the path of temp2 VCF.
    """
    header_only_path = f"{output_path}.temp0.vcf"
    cleaned_hdr_path = f"{output_path}.temp1.vcf"

    # Step 1: header-only file from input header
    write_header_only_vcf(input_vcf, header_only_path)

    # Step 2: remove overridden INFO IDs
    override_ids = parse_info_ids_from_header_file(header_file)
    strip_overridden_info_ids_from_header_text(header_only_path, override_ids, cleaned_hdr_path)

    # Step 3: load cleaned header via pysam and add new INFO lines
    with pysam.VariantFile(cleaned_hdr_path, "r") as cleaned_vcf:
        out_hdr = cleaned_vcf.header.copy()
    with open(header_file, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                out_hdr.add_line(line)

    # Step 4: write temp2 VCF duplicating records with the new header
    temp2_vcf = f"{output_path}.temp2.vcf"
    with pysam.VariantFile(temp2_vcf, "w", header=out_hdr) as outvcf:
        pass
    with pysam.VariantFile(input_vcf, "r") as invcf:
        with open(temp2_vcf, "a") as outvcf:
            for rec in invcf:
                outvcf.write(str(rec)) # to preserve INFO/END filed which is dropped when writing rec directly with pysam; see https://github.com/pysam-developers/pysam/issues/973
                # outvcf.write(rec)

    # Cleanup header-only intermediates
    for p in (header_only_path, cleaned_hdr_path):
        try:
            os.remove(p)
        except OSError:
            pass

    return temp2_vcf


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
        # Group all non-pos/id keys. Values may contain commas; we keep them as-is.
        for k, v in rec.items():
            if k in ("CHROM", "POS", "ID"):
                continue
            if v == "." or v == "":
                continue
            # Use pipe to separate groups to avoid VCF semicolon/comma parsing issues.
            groups.append(f"{k}={v}")
        if groups:
            m[vid] = "|".join(groups)
    f.close()
    return m


def cast_info_value_by_effective_header(tag: str, value: str, hdr: pysam.VariantHeader):
    """
    Cast values using the effective header attached to the current VariantFile.
    Return None for '.' or empty.
    """
    if value == "." or value == "":
        return None

    if tag not in hdr.info:
        # Skip tags not present in header to avoid KeyError
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
        # Number=1 scalar: commas inside are fine for String type
        return cast_one(value)


def annotate_vcf(
    vcf_path: str,
    info_header_path: str,
    annotations_tsv: str,
    bnd_mate_tsv: str,
    output_path: str,
):
    """
    Two-pass approach:
    - Pass 1: create a temp VCF with updated INFO header (overrides applied via text header edit).
    - Pass 2: reopen the temp VCF (effective header attached), assign annotations directly on rec.info,
              and stream records to the final output VCF.
    """
    ensure_valid_output_suffix(output_path)

    # Pass 1
    temp_vcf = create_temp_vcf_with_new_header(vcf_path, info_header_path, output_path)

    # Load maps
    primary_map, primary_tags = load_tsv_primary(annotations_tsv)
    bnd_scalar_map = load_tsv_bnd_mate_blob_scalar(bnd_mate_tsv) if bnd_mate_tsv else {}

    # Pass 2
    with pysam.VariantFile(temp_vcf, "r") as invcf:
        with pysam.VariantFile(output_path, "w", header=invcf.header) as outvcf:
            pass
        with open(output_path, "a") as outvcf:
            for rec in invcf:
                vid = rec.id
                # if vid != "0_Sniffles2.DEL.5624S0":
                #     continue
                if vid and vid in primary_map:
                    tagvals = primary_map[vid]
                    for tag in primary_tags:
                        raw_val = tagvals.get(tag, "")
                        val = cast_info_value_by_effective_header(tag, raw_val, invcf.header)
                        if val is None:
                            continue
                        try:
                            rec.info[tag] = val
                        except KeyError:
                            sys.stderr.write(f"Warning: skipping tag '{tag}' not present in header\n")

                if vid and vid in bnd_scalar_map:
                    blob = bnd_scalar_map[vid]  # single scalar string using pipe delimiter
                    if blob:
                        try:
                            rec.info["BND_MATE_INFO"] = blob
                        except KeyError:
                            sys.stderr.write(
                                "Warning: BND_MATE_INFO not present in header; "
                                "ensure Number=1,Type=String in --header\n"
                            )
                # outvcf.write(rec)
                outvcf.write(str(rec)) # to preserve INFO/END filed which is dropped when writing rec directly with pysam; see https://github.com/pysam-developers/pysam/issues/973

    # Cleanup temp2
    try:
        os.remove(temp_vcf)
    except OSError:
        pass


def main():
    args = parse_args()
    # if ouptut path doesn't end with .vcf or .vcf.gz, exit with error
    if not (args.output.endswith(".vcf")):
        sys.stderr.write("Error: Output path must end with .vcf\n")
        sys.exit(2)
    annotate_vcf(
        vcf_path=args.vcf,
        info_header_path=args.header,
        annotations_tsv=args.tsv,
        bnd_mate_tsv=args.mate,
        output_path=args.output,
    )
    sys.stderr.write(f"Written annotated VCF: {args.output}\n")


if __name__ == "__main__":
    main()