import subprocess
import pysam
import sys
import argparse
import os
from Bio.Seq import Seq

SAMTOOLS="samtools-1.21/samtools"
BCFTOOLS="bcftools-1.21/bcftools"
TABIX="htslib-1.21/tabix"

NO_RETURN_CODES=12

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
    return seq_length

def handle_vcf_types_0(args, vcf, record, chrom_lengths, i, f_len):
    output_dir = args.out
    # Write a new VCF file with only the ith record
    single_record_vcf = f"{output_dir}/temp_record_{i}.vcf.gz"
    with pysam.VariantFile(single_record_vcf, "w", header=vcf.header) as out_vcf:
        out_vcf.write(record)

    # Index the temporary VCF file using tabix
    tabix_command = [TABIX, "-p", "vcf", single_record_vcf]
    run_subprocess(tabix_command, "tabix")
    
    svtype = record.info.get("SVTYPE", None)
    svlen = get_svlen(record) #todo:check multi-allelic scenario
    if abs(svlen) < f_len:
        f_len += f_len

    # Determine chrom:start-end values
    chrom = record.chrom
    start = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end = record.stop + f_len # End position of the REF allele
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
    seq_len = process_fasta(temp_fasta, output_fasta, svID)

    # print(f"Info: Successfully created consensus sequence for record {i} at {output_fasta}\n")
    
    # Remove temp files
    os.remove(single_record_vcf)
    os.remove(f"{single_record_vcf}.tbi")
    os.remove(reference_fasta)
    os.remove(temp_fasta)

    return svID, svlen, seq_len

# https://github.com/samtools/bcftools/issues/1778
# DUP
def handle_vcf_types_1(args, vcf, record, chrom_lengths, i, f_len):
    output_dir = args.out
    fasta = pysam.FastaFile(args.ref)
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'

    svlen = get_svlen(record) #todo:check multi-allelic scenario
    if abs(svlen) < f_len:
        f_len += f_len

    # Determine chrom:start-end values
    chrom = record.chrom
    start_fl = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end_fl = record.pos # End position of the REF allele

    start_fr = record.stop
    chrom_length = chrom_lengths.get(chrom, 0)
    end = record.stop + f_len # End position of the REF allele
    end_fr = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    seq_fl = fasta.fetch(region=f"{chrom}:{start_fl}-{end_fl}")
    seq_fr = fasta.fetch(region=f"{chrom}:{start_fr}-{end_fr}")

    seq = fasta.fetch(region=f"{chrom}:{record.pos}-{record.stop}")

    seq_consensus = seq_fl + seq + seq + seq_fr # assuming tandem dups
    seq_len = len(seq_consensus)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{seq_consensus}\n")
    
    return svID, svlen, seq_len

# INV
def handle_vcf_types_2(args, vcf, record, chrom_lengths, i, f_len):
    output_dir = args.out
    fasta = pysam.FastaFile(args.ref)
    
    svtype = record.info.get("SVTYPE", None)
    svID = f'{svtype}.{i}'

    svlen = get_svlen(record) #todo:check multi-allelic scenario
    if abs(svlen) < f_len:
        f_len += f_len

    # Determine chrom:start-end values
    chrom = record.chrom
    start_fl = max(record.pos-f_len, 1) # Ensure start is >= 1 (1-based)
    end_fl = record.pos # End position of the REF allele

    start_fr = record.stop
    chrom_length = chrom_lengths.get(chrom, 0)
    end = record.stop + f_len # End position of the REF allele
    end_fr = min(end, chrom_length)  # Ensure end doesn't exceed chromosome length

    seq_fl = fasta.fetch(region=f"{chrom}:{start_fl}-{end_fl}")
    seq_fr = fasta.fetch(region=f"{chrom}:{start_fr}-{end_fr}")

    seq = fasta.fetch(region=f"{chrom}:{record.pos}-{record.stop}")

    seq_rc = Seq(seq).reverse_complement()

    seq_consensus = seq_fl + seq_rc + seq_fr # assuming tandem dups
    seq_len = len(seq_consensus)

    output_fasta = f"{output_dir}/{svID}.fa"
    with open(output_fasta, "w") as f:
        f.write(f">{svID}\n{seq_consensus}\n")
    
    return svID, svlen, seq_len

def handle_vcf_types_3(args, vcf, record, chrom_lengths, i, f_len):
    output_dir = args.out
    return 1

def is_valid_vcf_record(record, args):
    """
    Perform checks to validate a VCF record. Checks include:
    - Multi-allelic records
    - Multi-sample records

    Return 
    look at RETURN_CODE_DESCRIPTIONS
    """

    # 1. Check if the record is multi-allelic (more than one ALT allele)
    if record.alts is not None and len(record.alts) > 1:
        print(f"Warning: Multi-allelic record at {record.chrom}:{record.pos}")
        # print(f"Warning: Multi-allelic record at {record.chrom}:{record.pos} (REF: {record.ref}, ALTs: {record.alts})")

    # 2. Check if the record has multiple samples (more than one sample with genotype data)
    sample_count = len(record.samples)
    if sample_count > 1:
        print(f"Warning: Multi-sample record at {record.chrom}:{record.pos} with {sample_count} samples")
    
    # Extract SVTYPE from INFO field
    svtype = record.info.get("SVTYPE", None)

    # 3. Check if the ALT allele is symbolic (e.g., "<INS>", "<DEL>")
    has_symbolic_alt = any(alt.startswith("<") and alt.endswith(">") for alt in record.alts) if record.alts else False

    if "SVLEN" in record.info:
        svlen = abs(get_svlen(record))
        if svlen < args.min or svlen > args.max:
            return 1
    
    if svtype == "INS" and record.pos != record.stop:
        print(f"Error: record at {record.chrom}:{record.pos} is an INS with end ({record.end}) not equal to pos")
        return -1

    # Categorize based on SVTYPE and symbolic alleles
    if svtype in {"DEL"}:
        if "SVLEN" in record.info:
            return 2
        else:
            return 3
    elif svtype == "INS":
        if has_symbolic_alt:
            return 4
        if "SVLEN" in record.info:
            return 5
        else:
            return 6
    elif svtype in {"DUP"}:
        if "SVLEN" in record.info:
            return 7
        else:
            return 8
    elif svtype in {"INV"}:
        if "SVLEN" in record.info:
            return 9
        else:
            return 10
    elif svtype == "BND":
        return 11
    else:
        print(f"Warning: unknown SVTYPE ({svtype}) record at {record.chrom}:{record.pos}")
        return 12
    
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
                            f"args.min ({args.min}) <= SVLEN <= args.max ({args.max}) not satisfied",
                            "proper SV DEL",
                            "SVLEN missing for DEL",
                            "INS has symbolic ALT",
                            "proper SV INS",
                            "SVLEN missing for INS",
                            "proper SV DUP",
                            "SVLEN missing for DUP",
                            "proper SV INV",
                            "SVLEN missing for INV",
                            "proper SV BND",
                            "unknown SVTYPE"]

    code = 1
    with open(args.summary, "w") as f:
        print("records stats (count : description)")
        f.write(f"records stats (count : description)\n")
        for count in vcf_summary[1:]:
            # print("code {}: record count: {} {}".format(code, count, RETURN_CODE_DESCRIPTIONS[code]))
            print(f"{count} : {RETURN_CODE_DESCRIPTIONS[code]}")
            f.write(f"{count} : {RETURN_CODE_DESCRIPTIONS[code]}\n")
            code += 1

def print_record_stats(args, record_stats_arr, balanced_seq_bins):
    with open(args.summary, "a") as f:
        f.write(f"('SV_ID', svlen, consensus_seq_len)\n")
        for record in record_stats_arr:
            f.write(f"{record}\n")
        for i, b in enumerate(balanced_seq_bins):
            f.write(f"Bin {i+1}: {b}, Total sum: {sum(x[2] for x in b)}")


def balance_bins(items, n_bins):
    """
    Distribute items into n bins such that the total sum of value2 in each bin is as balanced as possible.

    :param items: List of tuples (string_id, value1, value2)
    :param n_bins: Number of bins
    :return: List of bins, where each bin is a list of tuples
    """
    # Sort items by value2 in descending order for greedy balancing
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
                    seq_consensus = lines[1]
                    output.write(f">{svID}{seq_consensus}")
                os.remove(input_fasta)

def process_vcf_records(args):
    # Load the VCF file
    vcf = pysam.VariantFile(args.vcf)
    print(f"\nWarning: Please make sure the input ref ({args.ref}) is same as the reference in the vcf header ({args.vcf})\n")

    # Extract chromosome lengths
    f_len = args.flank
    chrom_lengths = extract_chromosome_lengths(vcf)

    vcf_summary = [0]*(NO_RETURN_CODES+1)

    record_stats_arr = []
    
    # Iterate over each variant in the VCF file
    for i, record in enumerate(vcf):
        # Perform error checks (you can add your custom checks here)
        ret = is_valid_vcf_record(record, args)
        vcf_summary[ret] += 1
        if ret == 2 or ret == 5:
            record_stats = handle_vcf_types_0(args, vcf, record, chrom_lengths, i, f_len)
        elif ret == 7:
            record_stats = handle_vcf_types_1(args, vcf, record, chrom_lengths, i, f_len)
        elif ret == 9:
            record_stats = handle_vcf_types_2(args, vcf, record, chrom_lengths, i, f_len)
        elif ret == 11:
            record_stats = handle_vcf_types_3(args, vcf, record, chrom_lengths, i, f_len)
        else:
            print(f"Skipping invalid record at {record.chrom}:{record.pos}")
            continue
        # if i == 100:
        #     break
        record_stats_arr.append(record_stats)
    
    balanced_seq_bins = balance_bins(record_stats_arr, 10)
    concat_fasta(args, balanced_seq_bins)
    
    print_vcf_summary(args, vcf_summary)
    print_record_stats(args, record_stats_arr, balanced_seq_bins)


    

def check_program(tool):
    command = [tool, "--version"]
    result = run_subprocess(command, tool)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process VCF file and create alternative sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('-v', '--vcf', required=True, type=str, help="Path to the input VCF file (compressed or uncompressed)")
    required_args.add_argument('-r', '--ref', required=True, type=str, help="Path to the reference FASTA file")
    required_args.add_argument('-o', '--out', required=True, type=str, help="Path to the output directory")
    
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--min', required=False, type=int, default=50, help="The minimum length of the SV")
    optional_args.add_argument('--max', required=False, type=int, default=50000, help="The maximum length of the SV")
    optional_args.add_argument('--flank', required=False, type=int, default=2000, help="The maximum detectable period size supported by TRF to determine the flanking sequences")
    optional_args.add_argument('--summary', required=False, type=str, default="extract_sv.summary", help="Path to write the summary file")
    optional_args.add_argument('-h', '--help', action='help', help="Show this help message and exit")

    args = parser.parse_args()

    print(f"Info: VCF File: {args.vcf}")
    print(f"Info: Reference FASTA: {args.ref}")
    print(f"info: Output Directory: {args.out}")
    print(f"info: Min SV Length: {args.min}")
    print(f"Info: Max SV Length: {args.max}")

    if not os.path.exists(args.out):
        os.mkdir(args.out)
    else:
        print("Error: {} output dir already exists.".format(args.out))
        exit(1)

    check_program(SAMTOOLS)
    check_program(BCFTOOLS)
    check_program(TABIX)

    process_vcf_records(args)
