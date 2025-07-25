import sys
import re
from collections import defaultdict

# Mapping from simulated sv_type to preferred VCF SVTYPE(s)
preferred_svtype_map = {
    "tandem repeat contraction": "DEL",
    "deletion": "DEL",
    "approximate tandem repetition": "INS",
    "tandem repeat expansion": "INS",
    "perfect tandem repetition": "INS",
    "tandem duplication": ["DUP", "INS"],
    "inverted tandem duplication": ["INV", "INS"],  # prefer INV
    "inversion": "INV",                             # only INV
    "insertion": "INS"
}

# Load BED regions: REGION_ID → list of variant_ids
bed_regions = defaultdict(list)
with open("matched_regions.bed") as bed_file:
    for line in bed_file:
        fields = line.strip().split("\t")
        if len(fields) < 4:
            continue
        variant_id = fields[0]
        chrom, start, end = fields[1], fields[2], fields[3]
        region_id = f"{chrom}:{start}-{end}"
        bed_regions[region_id].append(variant_id)

# Build region → list of annotated variants
region_annotations = defaultdict(list)
with open("base_ref/simulated_svtypes.tsv") as tsv_file:
    header = tsv_file.readline()
    for line in tsv_file:
        fields = line.strip().split("\t")
        if len(fields) < 13:
            continue

        region_id = f"{fields[9]}:{fields[10]}-{fields[11]}"
        sim_id = fields[0]
        name_field = fields[4]
        sv_type = fields[12]
        svlen = fields[7]
        frac = fields[8]

        # Classify repeat and extract subfamily
        if any(rep in name_field for rep in ["LINE", "SINE", "DNA", "LTR", "ERV", "Retroposon"]):
            classification = "ME"
            match_subfamily = re.search(r'_name=([^_@]+)', name_field)
            subfamily = match_subfamily.group(1) if match_subfamily else "NA"
            match_class = re.search(r'#([^_]+)', name_field)
            repeat_detail = match_class.group(1) if match_class else "NA"
        elif ":" in name_field and name_field.count(":") == 1 and "#" not in name_field:
            classification = "TR"
            repeat_detail = name_field.split(":")[0]
            subfamily = "NA"
        else:
            classification = "Non-repetitive"
            repeat_detail = "NA"
            subfamily = "NA"

        for variant_id in bed_regions.get(region_id, []):
            region_annotations[region_id].append((
                variant_id, region_id, sim_id, sv_type,
                classification, repeat_detail, svlen, frac, subfamily
            ))

# Scoring function with strict VCF SVTYPE match
def score(record):
    variant_id, _, _, sv_type, _, _, _, _, _ = record
    caller_type = variant_id.split(".")[1]  # e.g., DEL, INS, INV

    expected = preferred_svtype_map.get(sv_type)
    valid_call = False

    if isinstance(expected, list):
        if caller_type in expected:
            preference = expected.index(caller_type)
            valid_call = True
    else:
        if caller_type == expected:
            preference = 0
            valid_call = True

    # Assign high penalty if not valid
    if not valid_call:
        preference = 999

    return (preference, variant_id)

# Output one valid variant per region
with open("matched_annotated.tsv", "w") as matched_out:
    for region_id, records in region_annotations.items():
        sorted_records = sorted(records, key=score)
        best = sorted_records[0]
        if score(best)[0] == 999:
            continue  # skip if no valid match
        variant_id, region_id, sim_id, sv_type, classification, repeat_detail, svlen, frac, subfamily = best
        matched_out.write(f"{variant_id}\t{region_id}\t{sim_id}\t{sv_type}\t{svlen}\t{frac}\t{classification}\t{repeat_detail}\t{subfamily}\n")

