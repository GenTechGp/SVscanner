import subprocess
import pysam
import sys
import argparse
import os
import time
import re

SAMTOOLS="./samtools-1.21/samtools"
BCFTOOLS="./bcftools-1.21/bcftools"
TABIX="./htslib-1.21/tabix"

NO_RETURN_CODES=15

INFO_COLS = ["chrom", "querystart", "queryend", "pos", "end", "svlen", "relativeID", "callerID", "ref", "alt"]

def get_f_len(args, svlen) :
    return min(args.flen, svlen * args.ffac)

def check_program(tool):
    if not (os.path.exists(tool) and os.path.isfile(tool)):
        print("Error: {} does not exist. Run ./scripts/install_tools.sh".format(tool))
        exit(1)

def extract_chromosome_lengths(vcf):
    # Create a dictionary to store chromosome lengths
    chrom_lengths = {}
    
    # Look for the "##contig" lines in the header
    for line in vcf.header.records:
        if line.key == 'contig':  # 'contig' lines contain chromosome length info
            chrom_name = line['ID']  # Chromosome name (e.g., 'chr20')
            chrom_length = int(line['length'])  # Chromosome length (e.g., 62435964)
            chrom_lengths[chrom_name] = chrom_length
    
    return chrom_lengths

def run_subprocess(command, prefix):
    result = None
    try:
        result = subprocess.run(command, check=True, text=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        # if result.stderr:
        #     for line in result.stderr.splitlines():
        #         print(f"[{prefix}  STDERR]: {line}")
    except subprocess.CalledProcessError as e:
        if e.stderr:
            for line in e.stderr.splitlines():
                print(f"[{prefix}  STDERR]: {line}")
        else:
            print(f"Error running {prefix}: {e}")
        sys.exit(1)
    # print(result)
    return result

def get_svlen(record):
    """
    Extracts the SVLEN value from the VCF record.
    - If SVLEN is a tuple, return the first value.
    - If SVLEN is an int, return it directly.
    - If SVLEN is missing, return None.
    """
    svlen = record.info.get("SVLEN", None)
    
    if isinstance(svlen, tuple):  
        return svlen[0]  # Take the first value if it's a tuple
    return svlen  # Return directly if it's an int or None

def process_fasta(input_fasta, output_fasta, new_title):
    """
    Reads a single-record FASTA file, changes the title, calculates sequence length, and writes back.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to save the modified FASTA file.
        new_title (str): New title for the sequence.
    """
    with open(input_fasta, "r") as f:
        lines = f.readlines()

    if not lines or lines[0][0] != ">":
        raise ValueError("Invalid FASTA format: Missing header line.")

    sequence = "".join(line.strip() for line in lines[1:])  # Join sequence lines
    seq_length = len(sequence)

    with open(output_fasta, "w") as f:
        f.write(f">{new_title}\n{sequence}\n")

    # print(f"{new_title}: seq_len = {seq_length}")
    return seq_length, sequence

#bcftools consensus method
def handle_vcf_types_using_bcftools(args, vcf, fasta, record, chrom_lengths, i, f_len):
    svlen = get_svlen(record)
    output_dir = args.out

    # Write a new VCF file with only the ith record
    single_record_vcf = f"{output_dir}/temp_record_{i}.vcf.gz"
    with pysam.VariantFile(single_record_vcf, "w", header=vcf.header) as out_vcf:
        out_vcf.write(record)

    # Index the temporary VCF file using tabix
    tabix_command = [TABIX, "-p", "vcf", single_record_vcf]
    run_subprocess(tabix_command, "tabix")
    
    svtype = record.info.get("SVTYPE", None)

    # Determine chrom:start-end values
    chrom = record.chrom
    start = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end = record.stop+f_len-1 # End position of the REF allele
    chrom_length = chrom_lengths.get(chrom, 0)
    end = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length
    
    # Call the samtools_command to create the reference sequence
    reference_fasta = f"{output_dir}/record_{i}_ref.fa"
    samtools_command = [SAMTOOLS, "faidx", args.ref, f"{chrom}:{start}-{end}", "--output", reference_fasta]
    run_subprocess(samtools_command, "samtools faidx")

    svID = f'{svtype}.{i}'
    # Call the bcftools_command to create the consensus sequence
    temp_fasta = f"{output_dir}/temp_{svID}.fa"
    bcftools_command = [BCFTOOLS, "consensus", "-f", reference_fasta,  single_record_vcf, "--output", temp_fasta]
    run_subprocess(bcftools_command, "bcftoools consensus")

    output_fasta = f"{output_dir}/{svID}.fa"
    seq_len, seq = process_fasta(temp_fasta, output_fasta, svID)

    # print(f"Info: Successfully created consensus sequence for record {i} at {output_fasta}\n")
    
    # Remove temp files
    os.remove(single_record_vcf)
    os.remove(f"{single_record_vcf}.tbi")
    os.remove(reference_fasta)
    os.remove(temp_fasta)

    return [svID, svlen, seq_len, seq, start, end]

#INS
def handle_vcf_types_ins(args, vcf, fasta, record, chrom_lengths, i):
    # svlen = get_svlen(record) #commented because sniffles svlen for ins is incorrect
    seq = record.alts[0]
    svlen = len(seq)
    assert svlen > 0
    f_len = get_f_len(args, svlen)

    output_dir = args.out
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'

    assert record.pos == record.stop

    # Determine chrom:start-end values
    chrom = record.chrom
    start_fl = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end_fl = record.pos-1 # End position of the REF allele

    start_fr = record.stop+1
    chrom_length = chrom_lengths.get(chrom, 0)
    end = record.stop+f_len-1 # End position of the REF allele
    end_fr = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    try:
        seq_fl = fasta.fetch(region=f"{chrom}:{start_fl}-{end_fl}")
        seq_fr = fasta.fetch(region=f"{chrom}:{start_fr}-{end_fr}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise
    # print(len(seq_fl))
    # print(len(seq_fr))
    # print(seq_fl)
    # print(seq)
    # print(seq_fr)

    if svlen != len(seq):
        print(f"Error: record (at {record.chrom}:{record.pos}) svlen ({svlen}) != len(ALT[0]) ({len(seq)})")

    query_seq = seq_fl + seq + seq_fr

    if args.debug:
        _,_,_,seq_bcf,_,_ = handle_vcf_types_using_bcftools(args, vcf, fasta, record, chrom_lengths, i, f_len)
        flag_homozygous_for_ref_allele = record.samples[list(vcf.header.samples)[0]]["GT"] == (0,0)
        if query_seq.lower() != seq_bcf.lower() and not flag_homozygous_for_ref_allele:
            print(f"svlen ({get_svlen(record)}) != len(ALT[0]) ({len(seq)})")
            print(len(seq_fl))
            print(len(seq_fr))
            print(f"Error: record (ID:{record.id}) query_seq != seq_bcf")
            print(">query_seq")
            print(query_seq)
            print(">seq_bcf")
            print(seq_bcf)
            exit(1)

    seq_len = len(query_seq)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{query_seq}\n")
    
    return [svID, svlen, seq_len, start_fl, end_fr]

#DEL
def handle_vcf_types_del(args, vcf, fasta, record, chrom_lengths, i):
    svlen = get_svlen(record) 
    assert svlen < 0
    f_len = get_f_len(args, abs(svlen))

    output_dir = args.out
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'

    # Determine chrom:start-end values
    chrom = record.chrom
    start_f = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end = record.stop+f_len-1 # End position of the REF allele
    chrom_length = chrom_lengths.get(chrom, 0)
    end_f = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    try:
        seq = fasta.fetch(region=f"{chrom}:{start_f}-{end_f}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise    
    query_seq = seq
    seq_len = len(query_seq)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{query_seq}\n")
    
    return [svID, svlen, seq_len, start_f, end_f]

# https://github.com/samtools/bcftools/issues/1778
# DUP
def handle_vcf_types_dup(args, vcf, fasta, record, chrom_lengths, i):
    svlen = get_svlen(record) 
    assert svlen > 0
    f_len = get_f_len(args, svlen)

    output_dir = args.out
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'


    # Determine chrom:start-end values
    chrom = record.chrom
    start_fl = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end_fl = record.pos-1 # End position of the REF allele

    start_fr = record.stop
    chrom_length = chrom_lengths.get(chrom, 0)
    end = record.stop+f_len-1 # End position of the REF allele
    end_fr = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    try:
        seq_fl = fasta.fetch(region=f"{chrom}:{start_fl}-{end_fl}")
        seq_fr = fasta.fetch(region=f"{chrom}:{start_fr}-{end_fr}")
        seq = fasta.fetch(region=f"{chrom}:{record.pos}-{record.stop-1}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise

    # print(len(seq_fl))
    # print(len(seq_fr))
    assert len(seq) == svlen

    query_seq = seq_fl + seq + seq + seq_fr # assuming tandem dups
    seq_len = len(query_seq)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{query_seq}\n")
    
    return [svID, svlen, seq_len, start_fl, end_fr]

# INV
def handle_vcf_types_inv(args, vcf, fasta, record, chrom_lengths, i):
    svlen = get_svlen(record) 
    assert svlen > 0
    f_len = get_f_len(args, svlen)

    output_dir = args.out
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'

    # Determine chrom:start-end values
    chrom = record.chrom
    start_f = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end = record.stop+f_len-1 # End position of the REF allele
    chrom_length = chrom_lengths.get(chrom, 0)
    end_f = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    try:
        seq = fasta.fetch(region=f"{chrom}:{start_f}-{end_f}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise

    query_seq = seq
    seq_len = len(query_seq)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{query_seq}\n")
    
    return [svID, svlen, seq_len, start_f, end_f]

def extract_bnd_target(alt):
    # Match formats like ]chr20:12345], [X:999999[
    match = re.search(r'[\[\]]([\w\.]+:\d+)[\[\]]', alt)
    if match:
        chrom = match.group(1).split(":")[0]  # Extract the chromosome part
        pos = match.group(1).split(":")[1]  # Extract the position part
        return chrom, pos
    else:
        return None, None

#BND
def handle_vcf_types_bnd(args, vcf, fasta, record, chrom_lengths, i):
    f_len = args.flen
    output_dir = args.out
    
    svtype = record.info.get("SVTYPE", None)
    assert svtype in {"BND", "TRA"}, f"Error: SVTYPE {svtype} is not supported for BND processing"
    svtype = "BND"
    svID_0 = f'{svtype}.{i}.0'
    svID_1 = f'{svtype}.{i}.1'

    # BND records have two parts, one for each end of the breakend
    # The first part is the reference sequence for the first end of the breakend
    chrom_0 = record.chrom
    start_f_0 = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end_0 = record.pos+f_len-1 # End position of the REF allele
    chrom_length = chrom_lengths.get(chrom_0, 0)
    end_f_0 = min(end_0, chrom_length)  # Ensure end doesn't exceed chromosome length
    query_seq_0 = ""
    try:
        query_seq_0 = fasta.fetch(region=f"{chrom_0}:{start_f_0}-{end_f_0}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise
    seq_len_0 = len(query_seq_0)

    # The second part is the reference sequence for the second end of the breakend
    chrom_1, pos_1 = extract_bnd_target(record.alts[0])  # Extract the chromosome from the ALT allele
    if chrom_1 is None or pos_1 is None:
        chrom_1 = record.info.get("CHR2") if "CHR2" in record.info else None 
        pos_1 = record.stop
    if chrom_1 is None or pos_1 is None:
        print(f"Error: Unable to extract target chromosome from ALT allele {record.alts[0]} or INFO tags in record {record.id}")
        exit(1)
    pos_1 = int(pos_1)  # Convert position to integer
    start_f_1 = max(pos_1-f_len, 1) # Ensure start is >= 1 (1-based)
    end_1 = pos_1+f_len-1 # End position of the REF allele
    chrom_length = chrom_lengths.get(chrom_1, 0)
    end_f_1 = min(end_1, chrom_length)  # Ensure end doesn't exceed chromosome length
    query_seq_1 = ""
    try:
        query_seq_1 = fasta.fetch(region=f"{chrom_1}:{start_f_1}-{end_f_1}")
    except Exception as e:
        print(f"Error fetching sequence for ({record.id}): {e}")
        raise
    seq_len_1 = len(query_seq_1)

    output_fasta = f"{output_dir}/{svID_0}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID_0}\n{query_seq_0}\n")
    output_fasta = f"{output_dir}/{svID_1}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID_1}\n{query_seq_1}\n")

    return [svID_0, -1, seq_len_0, start_f_0, end_f_0], [svID_1, -1, seq_len_1, start_f_1, end_f_1, chrom_1, pos_1, end_1]

def is_valid_vcf_record(record, args, chrom_lengths, warnings_dict):
    """
    Perform checks to validate a VCF record.

    Return 
    look at RETURN_CODE_DESCRIPTIONS
    """


    # Check if the record is multi-allelic (more than one ALT allele)
    if record.alts is not None and len(record.alts) > 1 and warnings_dict["multi-allelic"] < args.warning_count:
        print(f"Warning: Multi-allelic record (ID:{record.id}) with {len(record.alts)} ALT alleles")
        warnings_dict["multi-allelic"] += 1
        # print(f"Warning: Multi-allelic record at {record.chrom}:{record.pos} (REF: {record.ref}, ALTs: {record.alts})")
        if warnings_dict["multi-allelic"] == args.warning_count:
            print("Warning: Suppressing further multi-allelic warnings")

    # Check if the record has multiple samples (more than one sample with genotype data)
    sample_count = len(record.samples)
    if sample_count > 1 and warnings_dict["multi-sample"] < args.warning_count:
        print(f"Warning: Multi-sample record (ID:{record.id}) with {sample_count} samples")
        warnings_dict["multi-sample"] += 1
        if warnings_dict["multi-sample"] == args.warning_count:
            print("Warning: Suppressing further multi-sample warnings")
        
    if not record.info:
        print(f"Error: record (ID:{record.id}) does not have any INFO keys. Please check the VCF file header")
        return 15

    chrom = record.chrom
    chrom_1 = record.info.get("CHR2") if "CHR2" in record.info else None 
    if chrom_1 is not None:
        chrom = chrom_1
    chrom_length = chrom_lengths.get(chrom, 0)
    if record.stop > chrom_length:
        print(f"Error: record (ID:{record.id}) stop ({record.stop}) exceeds chromosome length ({chrom_length}) for {chrom}")
        return 15

    # Extract SVTYPE from INFO field
    svtype = record.info.get("SVTYPE") if "SVTYPE" in record.info else None
    
    # Jasmine introduces SVLEN=0 for BND
    if svtype in {"BND", "TRA"}:
        return 5
    
    if "SVLEN" in record.info:
        svlen = abs(get_svlen(record))
        if svlen < args.min:
            return 6
        if svlen > args.max:
            return 7
        if svtype not in {"DEL", "INS"} and svlen != record.stop - record.pos:
            return 13
    
    if svtype == "INS" and record.pos != record.stop:
        print(f"Error: record (ID:{record.id}) is an INS with end ({record.end}) not equal to pos ({record.pos})")
        return 15

    # Check if the ALT allele is symbolic (e.g., "<INS>", "<DEL>")
    has_symbolic_alt = any(alt.startswith("<") and alt.endswith(">") for alt in record.alts) if record.alts else False

    # Categorize based on SVTYPE and symbolic alleles
    if svtype == "INS":
        if has_symbolic_alt:
            return 12
        if "SVLEN" in record.info:
            return 1
        else:
            return 8
    elif svtype == "DEL":
        if "SVLEN" in record.info:
            return 2
        else:
            return 9
    elif svtype == "DUP":
        if "SVLEN" in record.info:
            return 3
        else:
            return 10
    elif svtype == "INV":
        if "SVLEN" in record.info:
            return 4
        else:
            return 11
    elif svtype in {"BND", "TRA"}:
        return 5
    else:
        # print(f"Warning: unknown SVTYPE ({svtype}) record at {record.chrom}:{record.pos}")
        return 14
    
    # 3. Further checks can be added here (e.g., check if REF matches the reference, check for missing info, etc.)
    # Example: Check if REF allele matches the reference genome
    # (This can be added if you have the reference FASTA loaded and want to verify the REF allele)

    # Example further check: Make sure QUAL (quality) is above a threshold (e.g., 20)
    # if record.qual < 20:
    #     print(f"Warning: Low quality score at {record.chrom}:{record.pos} (QUAL: {record.qual})")

    # Return True if the record is valid (after all checks), or False otherwise
    # (you could return False if you want to skip invalid records)

def print_vcf_summary(args, vcf_summary):
    RETURN_CODE_DESCRIPTIONS = ["",
                            "proper SV INS", # 1
                            "proper SV DEL", # 2
                            "proper SV DUP", # 3
                            "proper SV INV", # 4
                            "proper SV BND", # 5
                            f"SVLEN >= args.min ({args.min}) is not satisfied", # 6
                            f"SVLEN <= args.max ({args.max}) is not satisfied", # 7
                            "SVLEN missing for INS", # 8
                            "SVLEN missing for DEL", # 9
                            "SVLEN missing for DUP", # 10
                            "SVLEN missing for INV", # 11
                            "INS has symbolic ALT", # 12
                            "malformed record SVLEN != END - POS", # 13
                            "unknown SVTYPE", #14
                            "other"] # 15

    code = 1
    summary_out = f"{args.out}/extract_sv.summary"
    with open(summary_out, "w") as f:
        print("records stats (count : description)")
        f.write(f"records stats (count : description)\n")
        for count in vcf_summary[1:]:
            # print("code {}: record count: {} {}".format(code, count, RETURN_CODE_DESCRIPTIONS[code]))
            if count > 0:
                print(f"{count} : {RETURN_CODE_DESCRIPTIONS[code]} (code: {code})")
            f.write(f"{count} : {RETURN_CODE_DESCRIPTIONS[code]} (code: {code})\n")
            code += 1

def print_record_stats(args, record_stats_arr, balanced_seq_bins):
    summary_out = f"{args.out}/extract_sv.summary"
    with open(summary_out, "a") as f:
        f.write(f".fa, no.of.seqs, total.seq.len\n")
        for i, b in enumerate(balanced_seq_bins):
            f.write(f"{i}.fa, {len(b)}, {sum(x[2] for x in b)}\n")

        f.write(f"('SV_ID', svlen, query_seq_len, i.fa)\n")
        for record in record_stats_arr:
            f.write(f"{record[0]},{record[1]},{record[2]},{record[5]}\n")

def balance_bins(items, n_bins):
    """
    Distribute items into n bins such that the total sum of value2 in each bin is as balanced as possible.

    :param items: List of tuples (string_id, value1, value2)
    :param n_bins: Number of bins
    :return: List of bins, where each bin is a list of tuples
    """
    # Ensure the number of items is greater than or equal to the number of bins
    if len(items) < n_bins:
        # raise ValueError("The number of items must be greater than or equal to the number of bins.")
        n_bins = len(items)
    
    # Sort items by value2 in descending order for greedy balancing
    if n_bins > 1:
        items = sorted(items, key=lambda x: x[2], reverse=True)

    # Initialize bins
    bins = [[] for _ in range(n_bins)]
    bin_sums = [0] * n_bins  # Track sum of value2 in each bin

    # Distribute items using a greedy approach
    for item in items:
        # Find the bin with the minimum total sum
        min_bin_index = bin_sums.index(min(bin_sums))
        bins[min_bin_index].append(item)
        bin_sums[min_bin_index] += item[2]  # Add value2 to bin sum
        item.append(min_bin_index)  # Append bin index to item for tracking

    return bins

def concat_fasta(args, balanced_seq_bins):
    output_dir = args.out
    for i, b in enumerate(balanced_seq_bins):
        output_fasta = f"{output_dir}/{i}.fa"
        with open(output_fasta, "w") as output: 
            for x in b:
                input_fasta = f"{output_dir}/{x[0]}.fa"
                with open(input_fasta, "r") as input:
                    lines = input.readlines()
                    svID = lines[0]
                    query_seq = lines[1]
                    output.write(f"{svID}{query_seq}")
                os.remove(input_fasta)

def write_to_id_file(args, vcf, record, record_stats):
    with open(args.info, "a") as f:
        # info="${chr}\t${startFlank}\t${endFlank}\t${pos}\t${end}\t${len}\t${id}\t${callerID}\t${ref}\t${alt}"
        info="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(record.chrom, record_stats[3], record_stats[3]+record_stats[2], record.pos, record.stop, record_stats[1], record_stats[0], record.id, record.ref, record.alts[0])
        # print(info)
        f.write(f"{info}\n")

def process_vcf_records(args):
    # Load the VCF file
    vcf = pysam.VariantFile(args.vcf)
    print(f"\nWarning: Please make sure the input ref ({args.ref}) is same as the reference in the vcf header ({args.vcf})\n")

    fasta = pysam.FastaFile(args.ref)

    # Extract chromosome lengths
    chrom_lengths = extract_chromosome_lengths(vcf)

    vcf_summary = [0]*(NO_RETURN_CODES+1)

    record_stats_arr = []

    warnings_dict = {"multi-allelic": 0, "multi-sample": 0}
    
    # Iterate over each variant in the VCF file
    for i, record in enumerate(vcf):
        # if record.id != "Sniffles2.INS.BS0":
        #     continue

        # Perform error checks (you can add your custom checks here)
        ret = is_valid_vcf_record(record, args, chrom_lengths, warnings_dict)
        vcf_summary[ret] += 1
        if ret == 1:
            record_stats = handle_vcf_types_ins(args, vcf, fasta, record, chrom_lengths, i)
        elif ret == 2:
            record_stats = handle_vcf_types_del(args, vcf, fasta, record, chrom_lengths, i)
        elif ret == 3:
            record_stats = handle_vcf_types_dup(args, vcf, fasta, record, chrom_lengths, i)
        elif ret == 4:
            record_stats = handle_vcf_types_inv(args, vcf, fasta, record, chrom_lengths, i)
        elif ret == 5:
            record_stats, second = handle_vcf_types_bnd(args, vcf, fasta, record, chrom_lengths, i)
            with open(args.info, "a") as f:
                # info="${chr}\t${startFlank}\t${endFlank}\t${pos}\t${end}\t${len}\t${id}\t${callerID}\t${ref}\t${alt}"
                info="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(second[5], second[3], second[3]+second[2], second[6], second[7], second[1], second[0], record.id, record.ref, record.alts[0])
                f.write(f"{info}\n")
            record_stats_arr.append(second[:5])
        else:
            print(f"Skipped. vcf check failed (code: {ret}) for record (ID:{record.id})")
            continue
        write_to_id_file(args, vcf, record, record_stats)
        # if i == 100:
        #     break
        record_stats_arr.append(record_stats)
    
    balanced_seq_bins = balance_bins(record_stats_arr, args.n)
    concat_fasta(args, balanced_seq_bins)
    
    print_vcf_summary(args, vcf_summary)
    print_record_stats(args, record_stats_arr, balanced_seq_bins)

def argparser():

    def positive_int(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is an invalid. It must be greater than zero.")
        return ivalue


    parser = argparse.ArgumentParser(
        description="Process VCF file and create alternative sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('-v', '--vcf', required=True, type=str, help="Path to the input VCF file (compressed or uncompressed)")
    required_args.add_argument('-r', '--ref', required=True, type=str, help="Path to the reference FASTA file")
    required_args.add_argument('-o', '--out', required=True, type=str, help="Path to the output directory")
    required_args.add_argument('-i', '--info', required=True, type=str, help="Path to write the info file")
    
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--min', required=False, type=positive_int, default=50, help="The minimum length of the SV")
    optional_args.add_argument('--max', required=False, type=positive_int, default=50000, help="The maximum length of the SV")
    optional_args.add_argument('--flen', required=False, type=positive_int, default=2000, help="The maximum detectable period size supported by TRF to determine the length of flanking sequences")
    optional_args.add_argument('--ffac', required=False, type=positive_int, default=10, help="Multiplication factor for SVLEN to determine the length of flanking sequences")
    optional_args.add_argument('-n', required=False, type=positive_int, default=1, help="Number of equal-sized output fasta files")
    optional_args.add_argument('--debug', required=False, action='store_true', help="Debug mode")
    optional_args.add_argument('--warning_count', required=False, type=positive_int, default=10, help="The maximum warning count for each type of warning")
    optional_args.add_argument('-h', '--help', action='help', help="Show this help message and exit")

    return parser

if __name__ == "__main__":
    start = time.time()
    
    parser = argparser()
    args = parser.parse_args()
    
    print(f"Info: VCF File: {args.vcf}")
    print(f"Info: Reference FASTA: {args.ref}")
    print(f"info: Output Directory: {args.out}")
    print(f"info: SV info file: {args.info}")
    print(f"info: Min SV Length: {args.min}")
    print(f"Info: Max SV Length: {args.max}")
    print(f"Info: Number of output fasta files: {args.n}")
    print(f"Info: flen: {args.flen}")
    print(f"Info: ffac: {args.ffac}")

    if args.debug:
        print(f"Info: Debug mode: {args.debug}")
        check_program(SAMTOOLS)
        check_program(BCFTOOLS)
        check_program(TABIX)

    if not os.path.exists(args.out):
        os.mkdir(args.out)
    else:
        print("Error: {} output dir already exists.".format(args.out))
        exit(1)
    
    with open(args.info, "w") as f:
        f.write("\t".join(INFO_COLS) + "\n")

    process_vcf_records(args)

    end = time.time()
    print(f"Run time: {end - start:.3f} seconds")
