import subprocess
import pysam
import sys
import argparse
import os
from Bio.Seq import Seq
import time
import random

def generate_random_dna(length):
    """Generate a random DNA sequence of given length with a fixed seed."""
    return ''.join(random.choices('ACGT', k=length))

def read_Mob_file(args):
    data = []
    with open(args.mob, "r") as file:
        name = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name:  # Store previous sequence
                    seq = "".join(sequence)
                    seq_len = len(seq)
                    if seq_len > args.min and seq_len < args.max:
                        data.append((name, seq))
                name = line[1:].strip().replace(" ", "_")  
                # name = line.split()[0][1:]  # Extract name (without '>')
                sequence = []
            else:
                sequence.append(line)
        seq = "".join(sequence)
        seq_len = len(seq)
        if seq_len > args.min and seq_len < args.max:
            data.append((name, seq))
    
    return data, len(data)

def read_Rep_file(args):
    data = []
    with open(args.rep) as file:
        for line in file:
            parts = line.strip().split("\t")
            count, sequence = parts[3].split("x")
            seq_len = int(count) * len(sequence)
            if seq_len > args.min and seq_len < args.max:
                data.append((int(count), sequence))
    
    return data, len(data)


# DEL - deletion (delete [start,end])
# chr1	10000000	11000000	deletion	None	0
# INV - inversion (invert [start,end])
# chr1	20000000	21000000	inversion	None	0
# TD - tandem duplication (duplicate [start,end]) following e.g. the resulting ref will have 2 copies in total
# chr1	30000000	31000000	tandem duplication	2	0
# ITD - inverted tandem duplication (duplicate [start,end] and invert duplicated segment) following e.g. the resulting ref will have 2 copies in total where one is inverted
# chr1	40000000	41000000	inverted tandem duplication	2	0
# INS - insertion (insert sequence immediately after end)
# chr1	50999999	51000000	insertion	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	0
# TRE - tandem repeat expansion (expand a pre-existent microsatellite, with start and end coordinates (chr1	50481	50513	16xGT) following e.g. the resulting ref will have extra 200 TAs. start position and the existing sequence do not matter.
# chr1	503138	503172	tandem repeat expansion	TA:200	0
# TRC - tandem repeat contraction (contract a pre-existent microsatellite, as above) (chr1	1773676	1773730	27xTG) following e.g. the resulting ref will remove 20 TGs.
# chr1	1773676	1773730	tandem repeat contraction	TG:20	0
# PTR - perfect tandem repetition (insert a new perfect tandem repetition immediately after end). similar to insertion
# chr1	60999999	61000000	perfect tandem repetition	AT:30	0
# ATR - approximate tandem repetition (insert a new approximate tandem repetition immediately after end). similar to PTR but with some error bases injected here and there
# chr1	70999999	71000000	approximate tandem repetition	AT:30:3	0


def decide_sv_type(mob_seqs, rep_seqs, sv_types):
    col_5 = random.randint(0, 10)
    col_6 = random.randint(0, 10)
    inject_to_ref = True
    sv_len = 0
    pointer = 0
    bed_tuple = ()
    seq_type = ""
    sv_type = ""

    if args.debug:
        seq_type = "repeat"
        sv_type = "approximate tandem repetition"

    if seq_type == "":
        # get keys of sv_types dict
        # randomly choose one of the keys
        seq_type = random.choice(list(sv_types.keys()))
    
    if seq_type == "mobile":
        
        if sv_type == "":
            sv_type = random.choice(sv_types["mobile"])
        inj_mob_seq = random.choice(mob_seqs)
        inj_seq = inj_mob_seq[1]
        inj_len = len(inj_seq)
        name = f"{inj_mob_seq[0]}"

        end = pointer + inj_len

        if sv_type == "deletion":
            bed_tuple = (pointer+1, end, sv_type, "None", 0)
            sv_len = inj_len
        elif sv_type == "inversion":
            bed_tuple = (pointer+1, end, sv_type, "None", 0)
            sv_len = inj_len
        elif sv_type == "tandem duplication":
            bed_tuple = (pointer+1, end, sv_type, str(col_5), 0)
            sv_len = inj_len*col_5-1
        elif sv_type == "inverted tandem duplication":
            bed_tuple = (pointer+1, end, sv_type, str(col_5), 0)
            sv_len = inj_len*col_5-1
        elif sv_type == "insertion":
            inject_to_ref = False
            bed_tuple = (pointer-1, pointer, sv_type, inj_seq, 0)
            sv_len = inj_len
        else:
            print(f"Error: incorrect seq_type ({seq_type}) and sv_type ({sv_type}) combination")
            exit()
        
        return (inject_to_ref, inj_seq, name, bed_tuple, sv_len)

    elif seq_type == "repeat":
        
        if sv_type == "":
            sv_type = random.choice(sv_types["repeat"])
        inj_rep_seq = random.choice(rep_seqs)
        rep_count = inj_rep_seq[0]
        motif = inj_rep_seq[1]
        inj_seq = motif*rep_count
        inj_len = len(inj_seq)
        name = f"{motif}:{rep_count}"

        end = pointer + inj_len

        if sv_type == "tandem repeat expansion":
            max_rep_count = min(inj_len*2, args.max)/len(motif)
            rep_count = random.randint(rep_count+1, int(max_rep_count))
            TRE = f"{motif}:{rep_count}"
            bed_tuple = (pointer, end, sv_type, TRE, 0)
            sv_len = len(motif)*rep_count
        elif sv_type == "tandem repeat contraction":
            rep_count = random.randint(1, rep_count-1)
            TRC = f"{motif}:{rep_count}"        
            bed_tuple = (pointer, end, sv_type, TRC, 0)
            sv_len = len(motif)*rep_count
        elif sv_type == "perfect tandem repetition":
            inject_to_ref = False
            PTR = f"{motif}:{rep_count}"
            bed_tuple = (pointer-1, pointer, sv_type, PTR, 0)
            sv_len = inj_len
        elif sv_type == "approximate tandem repetition": #inj_seq is incorrect in this case and only at VISOR HACk stage we know the inj_seq
            inject_to_ref = False
            error = random.randint(1,10)
            ATR = f"{motif}:{rep_count}:{error}"
            bed_tuple = (pointer-1, pointer, sv_type, ATR, 0)
            sv_len = len(motif)*rep_count+error
        else:
            print(f"Error: incorrect seq_type ({seq_type}) and sv_type ({sv_type}) combination")
            exit()
        return (inject_to_ref, inj_seq, name, bed_tuple, sv_len)

def get_f_len(args, svlen):
    return min(args.flen, svlen * args.ffac)

def create_ref(args, sv_types):
    mob_seqs, mob_count = read_Mob_file(args)
    rep_seqs, rep_count = read_Rep_file(args)
    print(f"Info: Filtered Mobile Elements count: {mob_count}")
    print(f"Info: Filtered Repeats count: {rep_count}")
    
    sv_count = args.n

    chrom = "ref0"
    ref = ""
    pointer = 0

    sv_list = [] #(inject_to_ref, inj_seq, name, bed_record, sv_len)
    flen_array = []
    for sv_i in range(0, sv_count):
        sv_i = decide_sv_type(mob_seqs, rep_seqs, sv_types)
        sv_list.append(sv_i)
        sv_len = sv_i[4]
        f_len = args.flen # get_f_len(args, sv_len)
        flen_array.append(f_len)

    f_len = 0
    ref_tsv_out = f"{args.out}/simulated_svtypes.tsv"
    hack_bed_out = f"{args.out}/visor_hack.bed"
    sv_i = 0
    sv_dict = {}
    #initialize sv_dict with all sv types
    for sv_type in sv_types:
        for sv_subtype in sv_types[sv_type]:
            sv_dict[sv_subtype] = 0
    # create a dict to store the frequency of each SV type
    with open(ref_tsv_out, "w") as f1, open(hack_bed_out, "w") as f2:
        f1.write("CHROM\tPOINTER\tINJ_LEN\tNAME\tINJECT_TO_REF\tINJ_SEQ\tBED_RECORD\n")
        
        while sv_i < sv_count:
        # for i in range(0, sv_count):
            f_len = flen_array[sv_i]
            ran_len = random.randint(f_len, 1.5*f_len)
            ran_seq = generate_random_dna(ran_len)
            ref = ref + ran_seq
            pointer += ran_len
            
            inject_to_ref, inj_seq, name, bt, sv_len = sv_list[sv_i][0], sv_list[sv_i][1], sv_list[sv_i][2], sv_list[sv_i][3], sv_list[sv_i][4]
            inj_len = len(inj_seq)
            # (pointer+1, end, sv_type, "None", 0)
            bed_record = f"{chrom}\t{bt[0]+pointer}\t{bt[1]+pointer}\t{bt[2]}\t{bt[3]}\t{bt[4]}"
            
            f1.write(f"{chrom}\t{pointer}\t{inj_len}\t{name}\t{inject_to_ref}\t{inj_seq}\t{bed_record}\n") 
            f2.write(f"{bed_record}\n")
            if inject_to_ref:
                ref = ref + inj_seq
                pointer += inj_len
            if pointer > args.len:
                print(f"Warning: reference length reached {args.len} bases after simulating {sv_i} SVs. Finishing simulation...")
                break
            sv_type_key = sv_list[sv_i][3][2]
            sv_dict[sv_type_key] += 1
            sv_i += 1

    ran_len = random.randint(f_len, 1.5*f_len)
    ran_seq = generate_random_dna(ran_len)
    ref = ref + ran_seq
    pointer += ran_len

    ref_out = f"{args.out}/base_ref.fa"
    with open(ref_out, "w") as f:
        f.write(f">ref0\n")
        f.write(f"{ref}\n") 
    
    print(f"Info: Number of SVs simulated: {sv_i}")
    print(f"Info: SV types simulated: {sv_dict}")
    ref_length = len(ref)
    print(f"Info: Reference length: {ref_length} bases")
    print(f"Info: Reference file created: {ref_out}")
    print(f"Info: Visor hack file created: {hack_bed_out}")

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
    required_args.add_argument('-m', '--mob', required=True, type=str, help="Path to the Mobile Elements file")
    required_args.add_argument('-r', '--rep', required=True, type=str, help="Path to the Repeats file")
    required_args.add_argument('-o', '--out', required=True, type=str, help="Path to the output directory")
    
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--seed', required=False, type=positive_int, default=42, help="The seed for random generator")
    optional_args.add_argument('--len', required=False, type=positive_int, default=100000000, help="The length of the simulated base reference (before editing with SVs)")
    optional_args.add_argument('-n', required=False, type=positive_int, default=100, help="Number of SVs to simulate")
    optional_args.add_argument('--min', required=False, type=positive_int, default=50, help="Minimum length of SV")
    optional_args.add_argument('--max', required=False, type=positive_int, default=50000, help="Maximum length of SV")
    optional_args.add_argument('--flen', required=False, type=positive_int, default=2000, help="The maximum detectable period size supported by TRF to determine the length of flanking sequences")
    optional_args.add_argument('--ffac', required=False, type=positive_int, default=10, help="Multiplication factor for SVLEN to determine the length of flanking sequences")
    optional_args.add_argument('--svtypes', required=False, type=str, default="", help="File that contains the SV types to simulate. If not provided, all types will be used. Format: check docs/SV_simulation.md")
    optional_args.add_argument('--debug', required=False, action='store_true', help="Debug mode")
    optional_args.add_argument('-h', '--help', action='help', help="Show this help message and exit")

    return parser

if __name__ == "__main__":
    start = time.time()
    parser = argparser()
    args = parser.parse_args()
    random.seed(args.seed)

    
    print(f"Info: Mobile Elements File: {args.mob}")
    print(f"Info: Repeats File: {args.rep}")
    print(f"info: seed: {args.seed}")
    print(f"info: base reference length: {args.len}")
    print(f"info: number of SVs to simulate: {args.n}")
    print(f"info: base reference length (not gauranteed, fewer SVs result in a shorter length): {args.len}")
    print(f"info: Output Directory: {args.out}")
    print(f"Info: Max SV Length (not gauranteed, e.g. duplications, expansions can exceed this limit): {args.max}")
    print(f"Info: Min SV Length: {args.min}")

    if args.debug:
        print(f"Info: Debug mode: {args.debug}")

    if not os.path.exists(args.out):
        os.mkdir(args.out)
    else:
        print("Error: {} output dir already exists.".format(args.out))
        exit(1)

    sv_types = {}
    if args.svtypes:
        # read the svtypes file and filter the SV types
        print(f"Info: SV types file: {args.svtypes}")
        with open(args.svtypes, "r") as f:
            # read line split by : first. put first element as key and second as value
            for line in f:
                parts = line.strip().split(":")
                if len(parts) > 1:
                    result = [x for x in parts[1].split(",") if x]
                    # result is not empty
                    if result:
                        sv_types[parts[0]] = result
    else:
        sv_types = {"mobile": ["deletion", "inversion", "tandem duplication", "inverted tandem duplication", "insertion"],
                "repeat": ["tandem repeat expansion", "tandem repeat contraction", "perfect tandem repetition", "approximate tandem repetition"]}

    #check if sv_types is empty
    if sv_types == {}:
        print("Error: No SV types provided")
        exit(1)

    #print the SV types
    for sv_type in sv_types:
        print(f"Info: SV types for {sv_type}: {sv_types[sv_type]}")
        
    create_ref(args, sv_types)

    end = time.time()
    print(f"Run time: {end - start:.3f} seconds")
