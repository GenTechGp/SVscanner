import argparse
import os
import gzip
import time
import pandas as pd

def get_sniffles_record(record_id, input_path):
    """Extract a record from sniffles.vcf.gz and write to outfile."""
    found = False
    record = None
    with gzip.open(input_path, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] == record_id:
                record = line
                found = True
                break

    if not found:
        print(f"[sniffles] Record ID '{record_id}' not found in {input_path}")
    return record

def get_diagram_record(record_id, input_path):
    """Extract the diagram for the given ID from diagram.txt and write to outfile."""
    target_marker = f"({record_id})"
    found = False
    buffer = ""
    blank_count = 0

    with open(input_path, "r") as infile:
        for line in infile:
            if not found:
                if target_marker in line:
                    buffer += line
                    found = True
                continue
            else:
                buffer += line
                if line.strip() == "":
                    blank_count += 1
                    # If we see two blank lines in a row, assume the diagram block ended
                    if blank_count >= 2:
                        break
                else:
                    blank_count = 0

    if not found:
        print(f"[diagram.txt] Diagram for record ID '{record_id}' not found in {input_path}")
    return buffer

def search_sim_record(sim_records, sim_record_hits, vcf_pos, vcf_end, chr, called_svlen, threshold=100):
        """Search for a simulated record in simulated_svtypes.tsv."""
        i = 0
        count = 0
        hit_line_num = 0
        hit_line = None
        for sim_record in sim_records:
            score = 0
            sim_record = sim_record.split('\t')
            sim_pos = int(sim_record[10])
            sim_end = int(sim_record[11])
            diff_pos = abs(sim_pos - vcf_pos)
            diff_end = abs(sim_end - vcf_end)
            diff = min(diff_pos, diff_end)
            sim_chr = sim_record[1]
            sim_svlen = sim_record[7]

            if sim_chr == chr and diff <= threshold and abs(int(sim_svlen)) + threshold > abs(int(called_svlen)):
                score = 1
            elif sim_chr == chr and sim_pos <= vcf_pos <= sim_end and sim_pos <= vcf_end <= sim_end:
                score = 1

            if score > 0:
                count += 1
                hit_line_num = i
                hit_line = sim_record
                sim_record_hits[i] += 1
            i += 1
            
        return hit_line, count, hit_line_num

def search_file(input_path, vcf_pos, vcf_end, chr, threshold=100):
    found = False
    record = None
    lines = ""
    hit_count = 0
    line_num = 0
    hit_line_num = 0
    with open(input_path, "r") as infile:
        # Skip the header line
        next(infile)
        for line in infile:
            record = line.strip().split('\t')
            sim_pos = int(record[10])
            sim_end = int(record[11])
            sim_chr = record[1]
            # diff_pos = abs(sim_pos - vcf_pos)
            diff_end = abs(sim_end - vcf_end)
            if sim_chr == chr and diff_end <= threshold:
                lines += line
                found = True
                hit_count += 1
                hit_line_num = line_num
            line_num += 1

    if not found:
        print(f"[Position {vcf_pos}] No records found within {threshold} bp in {input_path}")
    return lines, hit_count, hit_line_num

def get_tsv_record(record_id, input_path, all=False):
    """Extract the line with the given ID from a .tsv file and write to outfile."""
    found = False
    record = None
    records = []
    with open(input_path, "r") as infile:
        # Skip the header line
        next(infile)
        for line in infile:
            if all:
                records.append(line)
            else:
                if record_id in line:
                    record = line
                    found = True
                    break

    if not all and not found:
        print(f"Record ID '{record_id}' not found in {input_path}")
    return records if all else record

def get_fasta_record(record_id, input_path):
    """Extract a FASTA record by ID from a FASTA file where each sequence is on a single line."""
    found = False
    record = None
    with open(input_path, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()
                if record_id in header:
                    seq = next(infile, "").strip()  # The sequence is on the next line
                    record = f"{header}\n{seq}\n"
                    found = True
                    break

    if not found:
        print(f"[FASTA record ID '{record_id}' not found in {input_path}")
    return record

def fetch_for_id(args, output_path):
    with open(output_path, "w") as outfile:
        input_path = os.path.join(args.input_dir, "svclassifier/simulated_info.tab")
        record = get_tsv_record(args.id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")
        record = record.strip().split('\t')
        relative_id = record[6]
        sniffles_id = record[7]
        vcf_pos = record[3]
        vcf_end = record[4]

        input_path = os.path.join(args.input_dir, "sniffles.vcf.gz")
        record = get_sniffles_record(sniffles_id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")
        record = record.strip().split('\t')

        chr = record[0]
        input_path = os.path.join(args.input_dir, "base_ref/simulated_svtypes.tsv")
        record, hit_count, line_num = search_file(input_path, int(vcf_pos), int(vcf_end), chr)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")
        
        # input_path = os.path.join(args.input_dir, "svclassifier/simulated_annotated.vcf.gz")
        # record = get_sniffles_record(sniffles_id, input_path)
        # outfile.write(f"## {input_path}\n")
        # outfile.write(f"{record}\n")
        
        input_path = os.path.join(args.input_dir, "svclassifier/annotations_out/diagram.txt")
        record = get_diagram_record(sniffles_id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")

        input_path = os.path.join(args.input_dir, "svclassifier/annotations_out/plot_annotate.tsv")
        record = get_tsv_record(sniffles_id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")

        input_path = os.path.join(args.input_dir, "svclassifier/extract_sv_flanks_out/extract_sv.summary")
        record = get_tsv_record(relative_id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")
        
        fasta_file_number = record.split(',')[3].strip()
        input_path = os.path.join(args.input_dir, "svclassifier/extract_sv_flanks_out", f"{fasta_file_number}.fa")
        record = get_fasta_record(relative_id, input_path)
        outfile.write(f"## {input_path}\n")
        outfile.write(f"{record}\n")

def determine_classification(sim_record, svtype_map):
    element_type = sim_record[4]
    if element_type.startswith("DF"):
        # DF000000010.4#LINE/L1_name=L1MC4_3end_@Eutheria_[S:55]
        repeat_classification = element_type.split('_')[0].split('#')[1].split('/')[0]
        return repeat_classification
    else:
        # AGCCCAGG:13
        sv_type = sim_record[12]
        classification = svtype_map[sv_type]
        return classification

def fetch_all(args):
    output_path = os.path.join(args.output_dir, "all_records.txt")
    input_path = os.path.join(args.input_dir, "svclassifier/simulated_info.tab")
    header = ["chrom", "vcf_pos", "relative_ID", "sniffles_ID", "classification", "reciprocal", "sim_record_count", "sim_ID", "sim_ref_pos", "element_details", "simulated_svlen/called_svlen", "simulated_svtype/called_svtype", "sim_classification/called_classification", "sim_reciprocal/called_reciprocal"]
    sv_records = get_tsv_record(args.id, input_path, all=True)
    
    input_path = os.path.join(args.input_dir, "base_ref/simulated_svtypes.tsv")
    sim_records = []
    with open(input_path, "r") as f:
        sim_records = [line.strip() for line in f.readlines()[1:]]
    sim_record_hits = [0] * len(sim_records)
  
    prev_sim_line_num = 0
    no_sim_record_count = 0
    no_sim_record_list = []

    simulated_sv_labels = ["simulated_sv", "not_simulated_sv"]
    called_sv_labels = ["called_sv", "not_called_sv"]
    conf_matrix_sim_called_sv = pd.DataFrame(0, index=simulated_sv_labels, columns=called_sv_labels)
    
    svtype_map = {"deletion": "DEL",
                    "inversion": "INV",
                    "tandem duplication": "TD",
                    "inverted tandem duplication": "ITD",
                    "insertion": "INS",
                    "tandem repeat expansion": "TRE",
                    "tandem repeat contraction": "TRC",
                    "perfect tandem repetition": "PTR",
                    "approximate tandem repetition": "ATR"
                    }

    sim_svtype_labels = list(svtype_map.values())
    called_svtype_labels = ["DEL", "INV", "DUP", "BND", "INS"]
    conf_matrix_svtypes = pd.DataFrame(0, index=sim_svtype_labels, columns=called_svtype_labels)

    sim_classification_labels = ["LTR", "LINE", "SINE", "DNA", "Retroposon", "Random", "TRE", "TRC", "PTR", "ATR"]
    called_classification_labels = ["LTR", "LINE", "SINE", "DNA", "Retroposon", "NON_REPETITIVE", "HOMO", "STR", "TR"]
    conf_matrix_classification = pd.DataFrame(0, index=sim_classification_labels, columns=called_classification_labels)

    sim_reciprocal_labels = []
    for classificatin in ["LTR", "LINE", "SINE", "DNA", "Retroposon"]:
        for percentage in ["1.0", "0.75", "0.5", "0.25"]:
            sim_reciprocal_labels.append(f"{classificatin}_{percentage}")
    called_reciprocal_labels = ["Full", "Partial"]
    conf_matrix_reciprocal = pd.DataFrame(0, index=sim_reciprocal_labels, columns=called_reciprocal_labels)


    with open(output_path, "w") as outfile:
        outfile.write("\t".join(header) + "\n")
        for sv_record in sv_records:
            tsv_record = ""
            sv_record = sv_record.strip().split('\t')
            relative_id = sv_record[6]
            sniffles_id = sv_record[7]
            vcf_pos = sv_record[3]
            vcf_end = sv_record[4]
            chr = sv_record[0]
            called_svtype = relative_id.split('.')[0]
            called_svlen = sv_record[5]

            # print(relative_id)
            # if relative_id != "DUP.630":
            #     continue

            input_path = os.path.join(args.input_dir, "svclassifier/annotations_out/plot_annotate.tsv")
            plot_record = get_tsv_record(sniffles_id, input_path)
            called_classification = plot_record.strip().split('\t')[3]
            called_reciprocal = plot_record.strip().split('\t')[4]

            tsv_record = f"{chr}\t{vcf_pos}\t{relative_id}\t{sniffles_id}\t{called_classification}\t{called_reciprocal}"

            sim_record, count, line_num = search_sim_record(sim_records, sim_record_hits, int(vcf_pos), int(vcf_end), chr, called_svlen)
            if count > 0:
                if count > 1:
                    print(f"Warning: Multiple simulation records found for {relative_id} at position {vcf_pos} with count {count}")
                if line_num < prev_sim_line_num:
                    print("Line number in simulated_svtypes.tsv should be non-decreasing")
                prev_sim_line_num = line_num
            else:
                no_sim_record_count += 1
                no_sim_record_list.append(sv_record)
            tsv_record += f"\t{count}"
            
            if count == 1:
                # print(sim_record)
                sim_id = sim_record[0]
                sim_ref_pos = sim_record[2]
                element_type = sim_record[4]
                sv_type = sim_record[12]
                simulated_svtype = svtype_map[sv_type]
                simulated_svlen = sim_record[7]
                tsv_record += f"\t{sim_id}\t{sim_ref_pos}\t{element_type}"
                
                conf_matrix_svtypes.loc[simulated_svtype, called_svtype] += 1
                
                sim_classification = determine_classification(sim_record, svtype_map)
                conf_matrix_classification.loc[sim_classification, called_classification] += 1

                tsv_record += f"\t{simulated_svlen}/{called_svlen}\t{simulated_svtype}/{called_svtype}\t{sim_classification}/{called_classification}"
                if called_reciprocal != "NA" and sim_classification in ["LTR", "LINE", "SINE", "DNA", "Retroposon"]:
                    percentage = sim_record[8]
                    sim_reciprocal = f"{sim_classification}_{percentage}"
                    conf_matrix_reciprocal.loc[sim_reciprocal, called_reciprocal] += 1
                    tsv_record += f"\t{sim_reciprocal}/{called_reciprocal}"

            tsv_record += "\n"
            outfile.write(tsv_record)
    
    
    no_sv_record_list = []
    for i, hit in enumerate(sim_record_hits):
        if hit == 0:
            no_sv_record_list.append(sim_records[i])
    
    called_sv_sim_hit_freq_table = pd.Series(sim_record_hits).value_counts().sort_index()
    hit_freqs = ""
    for i, hit in enumerate(called_sv_sim_hit_freq_table.index):
        if hit > 0:
            hit_freqs += f"{hit}:{called_sv_sim_hit_freq_table[hit]},"

    conf_matrix_sim_called_sv.loc["simulated_sv", "called_sv"] = f"{len(sim_records)-called_sv_sim_hit_freq_table.get(0, 0)}, {len(sv_records)-no_sim_record_count} ({hit_freqs.strip(',')})"
    conf_matrix_sim_called_sv.loc["simulated_sv", "not_called_sv"] = called_sv_sim_hit_freq_table.get(0, 0)
    conf_matrix_sim_called_sv.loc["not_simulated_sv", "called_sv"] = no_sim_record_count
    conf_matrix_sim_called_sv.loc["not_simulated_sv", "not_called_sv"] = "NA"
    
    print(f"info: SV count in svclassifier/simulated_info.tab: {len(sv_records)}")
    print(f"info: Simulation records count in base_ref/simulated_svtypes.tsv: {len(sim_records)}")
    print(f"info: Count of SVs with no simulation record: {no_sim_record_count}")
    print(f"info: Count of Simulation records with no SV: {called_sv_sim_hit_freq_table.get(0, 0)}")

    with open(os.path.join(args.output_dir, "no_sim_record.txt"), "w") as f:
        INFO_COLS = ["chrom", "querystart", "queryend", "pos", "end", "svlen", "relativeID", "callerID", "ref", "alt"]
        f.write("\t".join(INFO_COLS) + "\n")
        for record in no_sim_record_list:
            f.write("\t".join(record) + "\n")
    with open(os.path.join(args.output_dir, "no_sv_record.txt"), "w") as f:
        SIM_COLS = ["SIM_ID", "CHROM", "POINTER", "INJ_LEN", "NAME", "INJECT_TO_REF", "INJ_SEQ", "SVLEN", "FRAC", "BED_RECORD"]
        f.write("\t".join(SIM_COLS) + "\n")
        for record in no_sv_record_list:
            f.write(f"{record}\n")
    print(f"info: No simulation record written to {os.path.join(args.output_dir, 'no_sim_record.txt')}")
    print(f"info: No SV record written to {os.path.join(args.output_dir, 'no_sv_record.txt')}")\
    
    conf_matrix_sim_called_sv.to_csv(os.path.join(args.output_dir, "conf_matrix_sim_called_sv.txt"), sep='\t')
    conf_matrix_svtypes.to_csv(os.path.join(args.output_dir, "conf_matrix_svtypes.txt"), sep='\t')
    conf_matrix_classification.to_csv(os.path.join(args.output_dir, "conf_matrix_classification.txt"), sep='\t')
    conf_matrix_reciprocal.to_csv(os.path.join(args.output_dir, "conf_matrix_reciprocal.txt"), sep='\t')

def argparser():
    def positive_int(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is an invalid. It must be greater than zero.")
        return ivalue
    parser = argparse.ArgumentParser(
        description="Summarise simulation results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('--input', dest='input_dir', type=str, required=True, help="Input directory")
    required_args.add_argument('--output', dest='output_dir', type=str, required=True, help="Output directory")
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--id', required=False, type=str, help="Record ID to search for")
    return parser

def main():
    start = time.time()
    
    parser = argparser()
    args = parser.parse_args()

    print(f"info: Input Directory: {args.input_dir}")
    print(f"info: Output Directory: {args.output_dir}")
    if args.id:
        print(f"info: Record ID: {args.id}")

    os.makedirs(args.output_dir, exist_ok=True)

    if args.id:
        output_path = os.path.join(args.output_dir, f"{args.id}.txt")
        print(f"info: Output file: {output_path}")
        fetch_for_id(args, output_path)
    else:
        print("No record ID provided. Writing all records to output directory.")
        fetch_all(args)
    
    end = time.time()
    print(f"Run time: {end - start:.3f} seconds")

if __name__ == "__main__":
    main()
