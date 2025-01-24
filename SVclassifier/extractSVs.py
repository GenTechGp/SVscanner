import sys
import argparse
import pysam

def getFlankingCoordinates(pos, end, flanking, chrom_length):
    """
    Calculates the start and end flanking coordinates ensuring it is not less than and 
    does not exceed the chromosome size
    """
    startFlanking = max(0, pos - flanking - 1)
    endFlanking = min(chrom_length - 1, end + flanking - 1)

    return startFlanking, endFlanking

def getSequences(fasta, chrom, pos, end, startFlanking, endFlanking):
    """
    Fetches the sequences 
    """
    # Start position to be 0-based and half-open (end exclusive)
    sequence = fasta.fetch(chrom, pos, end-1) if pos-1 >= 0 else ''
    flankingSequence = fasta.fetch(chrom, startFlanking, endFlanking)

    return sequence, flankingSequence

def get_flanking_regions(vcf_file, fasta_file, output_file, min_size, max_size, flanking_factor):
    """
    Extract flanking regions and sequences for structural variants (SVs) from a VCF file.

    The function extracts flanking regions (scaled by the 'flanking factor') for each SV based on its type and writes to a file 
    - Insertions: SV sequence is obtained from the ALT tag in the VCF and concatenated between the flanking regions to create the query
    - Other: SV sequence is fetched from the reference based on the positions

    Input: 
        - vcf_file (str): Path to the input VCF file containing structural variants.
        - fasta_file (str): Path to the FASTA file for the reference genome.
        - output_file (str): Path to the output file where extracted data will be written.
        - min_size (int): Minimum size of SVs to include in the output.
        - max_size (int): Maximum size of SVs to include in the output.
        - flanking_factor (float): Scaling factor to determine the size of flanking regions 
                                 based on the SV length.
    """
    vcf = pysam.VariantFile(vcf_file)
    fasta = pysam.FastaFile(fasta_file)
    total_sv = 0
    filtered_min = 0
    filtered_max = 0
    symbolic_ins = 0

    # Count for print of filtering 
    sv_count = {'INS' : 0, 'DEL' : 0, 'INV' : 0, 'DUP' : 0, 'BND' : 0}
    min_count = {'INS' : 0, 'DEL' : 0, 'INV' : 0, 'DUP' : 0, 'BND' : 0}
    max_count = {'INS' : 0, 'DEL' : 0, 'INV' : 0, 'DUP' : 0, 'BND' : 0}
    filtered_count = {'INS' : 0, 'DEL' : 0, 'INV' : 0, 'DUP' : 0, 'BND' : 0}

    with open(output_file, 'w') as out:
        for i, record in enumerate(vcf): 
            chrom = record.chrom
            pos = record.pos # pysam converts 1-base VCF to 0-base
            ref = record.ref
            alt = record.alts[0]
            svtype = record.info.get('SVTYPE')
            end = record.stop
            callerID = record.id

            # Generate a unique ID
            svID = f'{svtype}.{i}'
            sv_count[svtype] += 1
            # Get the length of the chromosome sequence
            chrom_length = fasta.lengths[fasta.references.index(chrom)]
              
            # INSERTION 
            if svtype == 'INS':  
                if alt == '<INS>':
                    symbolic_ins += 1
                    continue

                length = len(alt)
                startFlanking = max(0, pos - flanking_factor * length)
                endFlanking = min(chrom_length - 1, pos + flanking_factor * length)

                sequence = alt[len(ref):]
                flanking_region = fasta.fetch(chrom, startFlanking, endFlanking)

                # Find the position within the flanking sequence to insert the new sequence
                if startFlanking < pos < endFlanking:
                    # Join the flanking sequences and inserted sequence (between them) to make the ALT allele as a FASTA sequence
                    left_flank_end = pos - 1
                    right_flank_start = pos - 1
                    
                    if left_flank_end < startFlanking:
                        left_flank = ''
                    else:
                        left_flank = fasta.fetch(chrom, startFlanking, left_flank_end)
                    
                    if right_flank_start > endFlanking:
                        right_flank = ''
                    else:
                        right_flank = fasta.fetch(chrom, right_flank_start + 1, endFlanking)
                    
                    # Concatenate the sequence between flanking regions 
                    flankingSequence = left_flank + sequence + right_flank
                else:
                    flankingSequence = flanking_region
            
            # DELETION
            elif svtype == 'DEL' and end > pos:  
                if alt == '<DEL>' or ref == '<DEL>':
                    length = abs(end - pos) 

                length = len(ref)
                startFlanking, endFlanking = getFlankingCoordinates(pos, pos+length, flanking_factor*length, chrom_length)
                sequence, flankingSequence = getSequences(fasta, chrom, pos, pos+length, startFlanking, endFlanking)

            # DUPLICATION            # INVERSION  
            elif svtype == 'DUP' or svtype == "INV":
                length = end - pos + 1 
                startFlanking = max(0, pos - flanking_factor * length)
                endFlanking = min(chrom_length - 1, end + flanking_factor * length - 1)
                startFlanking, endFlanking = getFlankingCoordinates(pos, end, flanking_factor*length, chrom_length)
                sequence, flankingSequence = getSequences(fasta, chrom, pos, pos+length, startFlanking, endFlanking)

            else:
                continue

            # Filter out any SV smaller than the min and greater than max -> output only within 
            if length >= min_size and length <= max_size:
                out.write(f'{chrom}\t{startFlanking}\t{endFlanking}\t{pos}\t{end}\t{length}\t{svID}\t{flankingSequence}\t{sequence}\t{svtype}\t{callerID}\n')
                filtered_count[svtype] += 1
                total_sv += 1
            elif length < min_size:
                min_count[svtype] += 1
                filtered_min += 1
            elif length > max_size:
                max_count[svtype] +=1
                filtered_max += 1


    #------------------------ Command Line Output -----------------------------#           
    print(f'Total SV:           \t{total_sv:>6}')
    print(f'Filtered (<{min_size}):   \t{filtered_min:>6}')
    print(f'Filtered (>{max_size}):   \t{filtered_max:>6}')

    # Header with SV types
    print("-" * 75)
    print(f"{'SV Type Breakdown':<20} {'INS':>10} {'DEL':>10} {'INV':>10} {'DUP':>10} {'BND':>10}")
    print("-" * 75)

    # For each filtering stage, print the corresponding values for DEL, INS, INV, DUP, and BND
    for stage, count_dict in [('Original', sv_count), ('Min Filtered', min_count), ('Max Filtered', max_count), ('After Filtering', filtered_count)]:
        print(f"{stage:<20} {count_dict['INS']:>10} {count_dict['DEL']:>10} {count_dict['INV']:>10} {count_dict['DUP']:>10} {count_dict['BND']:>10}")
    print("-" * 35)
    # Symbolic Insertions
    print(f"{'Symbolic Insertions':<10} {symbolic_ins:>10}")
    print("-" * 35)


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="Extract and format sequences including variants from a VCF file.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-vcf",
        help="Path to the input VCF file (compressed or uncompressed)", required=True)
    parser.add_argument("-fa",
        help="Path to the reference FASTA file.", required=True)
    parser.add_argument("-out",
        help="Path to the output file", required=True)
    parser.add_argument("-min",
        help="The minimum length of the SV", default=50, type=int)
    parser.add_argument("-max",
        help="The maximum length of the SV", default=50000, type=int)
    parser.add_argument("-flank",
        help="Flanking factor surrounding the SV ", default=10, type=int)

    args = parser.parse_args()

    get_flanking_regions(args.vcf, args.fa, args.out, args.min, args.max, args.flank)
