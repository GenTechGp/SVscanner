#!/usr/bin/env python3
import argparse
import sys

try:
    import pysam
except ImportError:
    sys.stderr.write("Error: pysam is required. Install with: pip install pysam\n")
    sys.exit(1)


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Detect retrotransposition candidate insertion SVs from a SVscanner-annotated VCF. "
            "Reads RM_CLASSIFICATION, RM_SUBFAMILY, and RM_SV_COVERAGE tags produced by SVscanner "
            "and flags INS records where a known-active retrotransposon subfamily is present at "
            "sufficient coverage and size. Outputs a TSV compatible with src/annotate_vcf.py."
        )
    )
    ap.add_argument("--vcf", required=True, help="SVscanner-annotated VCF (.vcf or .vcf.gz)")
    ap.add_argument(
        "--config", required=True,
        help="Retrotransposition params TSV (config/retrotp_params.tsv)"
    )
    ap.add_argument("--out", required=True, help="Output annotations TSV path")
    return ap.parse_args()


def load_config(path):
    """
    Load retrotp_params.tsv.
    Returns dict: subfamily_name -> {class, confidence, min_rm_sv_coverage, svlen_min, svlen_max}.
    Lines beginning with '#' and the header row are skipped.
    """
    config = {}
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if cols[0] == "subfamily":
                continue
            if len(cols) < 6:
                sys.stderr.write(f"Warning: skipping malformed config row: {line}\n")
                continue
            subfamily, cls, confidence, min_cov, svlen_min, svlen_max = cols[:6]
            try:
                config[subfamily] = {
                    "class": cls,
                    "confidence": confidence,
                    "min_rm_sv_coverage": float(min_cov),
                    "svlen_min": int(svlen_min),
                    "svlen_max": int(svlen_max),
                }
            except ValueError:
                sys.stderr.write(f"Warning: skipping config row with non-numeric threshold: {line}\n")
    return config


def _get_info_list(rec, tag):
    """Return list of values for a Number=. INFO tag, or empty list if absent/None."""
    try:
        val = rec.info[tag]
    except KeyError:
        return []
    if val is None:
        return []
    return list(val)


def _get_svlen(rec):
    """Return abs(SVLEN) as int, or None if absent or unparseable."""
    try:
        raw = rec.info["SVLEN"]
    except KeyError:
        return None
    if isinstance(raw, (tuple, list)):
        raw = raw[0]
    try:
        return abs(int(raw))
    except (TypeError, ValueError):
        return None


def detect_retrotp(rec, config):
    """
    Evaluate one pysam VariantRecord for retrotransposition candidacy.

    Algorithm:
      1. SVTYPE must be INS.
      2. SVLEN must be present.
      3. For each i in zip(RM_SUBFAMILY, RM_SV_COVERAGE):
           - RM_SUBFAMILY[i] must be in config (known-active subfamily).
           - RM_SV_COVERAGE[i] >= config row's min_rm_sv_coverage.
           - abs(SVLEN) within [svlen_min, svlen_max] for that row.
           Emit (subfamily, confidence) if all pass, else emit ('.', '.').
      4. Return parallel lists of the same length as RM_SUBFAMILY, with '.'
         at positions that did not qualify. Return (None, None) if no position
         qualified (all entries would be '.').

    Output lists are parallel to RM_SUBFAMILY — positionally aligned with all
    other RM_* tags in the SVscanner VCF.
    """
    # Gate: insertions only
    try:
        svtype = rec.info["SVTYPE"]
    except KeyError:
        return None, None
    if svtype != "INS":
        return None, None

    svlen = _get_svlen(rec)
    if svlen is None:
        return None, None

    rm_subfamily = _get_info_list(rec, "RM_SUBFAMILY")
    rm_coverage = _get_info_list(rec, "RM_SV_COVERAGE")

    if not rm_subfamily:
        return None, None

    # Pad to equal length so zip is safe (lists should already be equal-length
    # in well-formed SVscanner output, but guard defensively)
    n = max(len(rm_subfamily), len(rm_coverage))
    rm_subfamily = rm_subfamily + [""] * (n - len(rm_subfamily))
    rm_coverage = list(rm_coverage) + [0.0] * (n - len(rm_coverage))

    out_subfamily = []
    out_confidence = []
    any_hit = False

    for sub, cov in zip(rm_subfamily, rm_coverage):
        row = config.get(sub)
        if row is None:
            out_subfamily.append(".")
            out_confidence.append(".")
            continue
        try:
            cov_f = float(cov)
        except (TypeError, ValueError):
            out_subfamily.append(".")
            out_confidence.append(".")
            continue
        if cov_f < row["min_rm_sv_coverage"] or not (row["svlen_min"] <= svlen <= row["svlen_max"]):
            out_subfamily.append(".")
            out_confidence.append(".")
            continue
        out_subfamily.append(sub)
        out_confidence.append(row["confidence"])
        any_hit = True

    if not any_hit:
        return None, None
    return out_subfamily, out_confidence


def main():
    args = parse_args()
    config = load_config(args.config)
    sys.stderr.write(f"Info: loaded {len(config)} active subfamilies from {args.config}\n")

    n_total = 0
    n_candidates = 0

    with pysam.VariantFile(args.vcf, "r") as vcf, open(args.out, "wt") as out_f:
        out_f.write("ID\tRETROTP_CANDIDATE\tRETROTP_ELEMENT\n")
        for rec in vcf:
            n_total += 1
            subfamilies, confidences = detect_retrotp(rec, config)
            if subfamilies is None:
                continue
            n_candidates += 1
            vid = rec.id if rec.id else f"{rec.chrom}_{rec.pos}"
            out_f.write(
                f"{vid}\t{','.join(confidences)}\t{','.join(subfamilies)}\n"
            )

    sys.stderr.write(f"Info: {n_total} records processed\n")
    sys.stderr.write(f"Info: {n_candidates} retrotransposition candidates written to {args.out}\n")


if __name__ == "__main__":
    main()
