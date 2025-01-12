import pysam
import csv
from statistics import mean
from tqdm import tqdm

def read_bed_file(bed_file):
    duplications = []
    duplicationDepth = {}
    with open(bed_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) >= 3:  # Ensure there are enough columns
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                dupID = row[3] if len(row) > 3 else f"dup_{len(duplications)}"
                duplications.append((chrom, start, end, dupID))
                duplicationDepth.update({dupID: {'chr': chrom, 'start': start, 'end': end, 'dupDepth': 0, 'leftDepth': 0, 'rightDepth': 0}})
    return duplications, duplicationDepth

# def calculate_region_depth(mybam, chrom, start, end):
#     depths = []

#     for pileupcolumn in mybam.pileup(chrom, start, end, truncate=True):
#         dup_del = 0
#         # print(f'Pileup Column {start}-{end}')
#         pileupcolumn.set_min_base_quality(0)
#         # Adjust depth calculation based on whether the base is deleted or not
#         for read in pileupcolumn.pileups:
#             if read.is_del or read.is_refskip: 
#                 dup_del += 1
#         depth = pileupcolumn.nsegments - dup_del
#         # print(f'{pileupcolumn.pos}: {pileupcolumn.nsegments} - {dup_del} = {depth}')
#         depths.append(depth)

#     return int(mean(depths)) if depths else None

def calculate_region_depth(mybam, chrom, start, end):
    depths = []
    
    # Iterate over the region from start to end with a progress bar
    depth_array = mybam.count_coverage(chrom, start, end, quality_threshold = 0)

    depths = [sum(position) for position in zip(*depth_array)]

    return int(mean(depths)) if depths else None

def calculate_read_depth(duplicationDepth, bam_file, ref_file, duplications):
    mybam = pysam.AlignmentFile(bam_file, "rb")
    reference = pysam.FastaFile(ref_file)

    for chrom, start, end, dupID in tqdm(duplications, desc="Processing duplications"):
        flanking = 2 * (end - start)
        chrom_length = reference.get_reference_length(chrom)

        leftFlankingStart = max(0, start - flanking - 1)  # Set to 0 if less than 0
        rightFlankingEnd = min(end + flanking, chrom_length)

        # Calculate depths for the duplication region
        dup_depth = calculate_region_depth(mybam, chrom, start-1, end-1)
        left_depth = calculate_region_depth(mybam, chrom, leftFlankingStart, start)
        right_depth = calculate_region_depth(mybam, chrom, end-1, rightFlankingEnd)

        if dup_depth:
            duplicationDepth[dupID]['dupDepth'] = dup_depth
        if left_depth:
            duplicationDepth[dupID]['leftDepth'] = left_depth
        if right_depth:
            duplicationDepth[dupID]['rightDepth'] = right_depth

    mybam.close()
    return duplicationDepth


if __name__ == "__main__":
    import argparse
    import sys
    import csv

    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="duplication depth",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-bed",
        help="BED file containing duplications")
    
    parser.add_argument("-bam",
        help="BAM file containing aligned reads")
    
    parser.add_argument("-o",
        help="Output file for calculated depths")
    
    parser.add_argument("-ref",
                    help="ref file ")

    args = parser.parse_args()
    bam_file = args.bam
    ref_file = args.ref

    duplications, duplicationDepth = read_bed_file(args.bed)

    # Update calculate_read_depth call to exclude args.ref
    duplicationDepth = calculate_read_depth(duplicationDepth, bam_file, ref_file, duplications)

    # Write results to the output file
    with open(args.o, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        # Write rows without header
        for duplication in duplicationDepth:
            chr_value = duplicationDepth[duplication]['chr']
            start_value = duplicationDepth[duplication]['start']
            end_value = duplicationDepth[duplication]['end']
            meanLeft = duplicationDepth[duplication]['leftDepth'] if 'leftDepth' in duplicationDepth[duplication] else ''
            meanDup = duplicationDepth[duplication]['dupDepth'] if 'dupDepth' in duplicationDepth[duplication] else ''
            meanRight = duplicationDepth[duplication]['rightDepth'] if 'rightDepth' in duplicationDepth[duplication] else ''
            writer.writerow([chr_value, start_value, end_value, duplication, meanDup, meanLeft, meanRight])
