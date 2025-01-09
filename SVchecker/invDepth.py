import pysam
import pysamstats
import sys
import argparse
import csv
from tqdm import tqdm

def read_bed_file(bed_file):
    inversions = []
    inversionDepth = {}
    with open(bed_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) >= 3:  # Ensure there are enough columns
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                invID = row[3]
                inversions.append((chrom, start, end, invID))
                inversionDepth.update({invID: {'chr' : chrom, 'start': start, 'end' : end, 'startDepth' : 0, 'endDepth' : 0}})
    return inversions, inversionDepth

def calculate_read_depth(inversionDepth, bam_file, inversions, fa_file):
    mybam = pysam.AlignmentFile(bam_file, "rb")
    
    for chrom, start, end, invID in tqdm(inversions, desc="inversions"):
        for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=start-1, end=start, truncate=True, fafile=fa_file):
            inversionDepth[invID]['startDepth'] = stat['reads_all']-stat['deletions']
        for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=end-1, end=end, truncate=True, fafile=fa_file):
            inversionDepth[invID]['endDepth'] = stat['reads_all']-stat['deletions']
    mybam.close()
    return inversionDepth


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="inversion depth",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-bed",
        help="bed file")
    
    parser.add_argument("-bam",
        help="bam file ")

    parser.add_argument("-ref",
        help="ref file ")

    parser.add_argument("-o",
        help="output file")

    args = parser.parse_args()
    bam_file = args.bam

    inversions, inversionDepth = read_bed_file(args.bed)

    inversionDepth = calculate_read_depth(inversionDepth, bam_file, inversions, args.ref)

    # Open the CSV file in write mode with tab delimiter and no header
    with open(args.o, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        # Write rows without header
        for inversion in inversionDepth:
            chr_value = inversionDepth[inversion]['chr']
            start_value = inversionDepth[inversion]['start']
            end_value = inversionDepth[inversion]['end']
            start_depth = inversionDepth[inversion]['startDepth'] if 'startDepth' in inversionDepth[inversion] else ''
            end_depth = inversionDepth[inversion]['endDepth'] if 'endDepth' in inversionDepth[inversion] else ''
            writer.writerow([chr_value, start_value, end_value, inversion, start_depth, end_depth])



# def calculate_read_depth(inversionDepth, bam_file, inversions, fa_file):
#     mybam = pysam.AlignmentFile(bam_file, "rb")
#     # depth_results = []
    
#     for chrom, start, end, invID in inversions:
#         for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=start-1, end=start, truncate=True, fafile=fa_file):
#             # print(stat)

#             # print(f"{invID}\t{stat['reads_all']}\t{stat['deletions']}\t={stat['reads_all']-stat['deletions']-stat['mismatches']}")
#             # print(stat)
#             inversionDepth[invID]['startDepth'] = stat['reads_all']-stat['deletions']
#         for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=end-1, end=end, truncate=True, fafile=fa_file):
#             # print(f"{invID}\t{stat['reads_all']}\t{stat['deletions']}\t={stat['reads_all']-stat['deletions']-stat['mismatches']}")
#             # print(stat)
#             inversionDepth[invID]['endDepth'] = stat['reads_all']-stat['deletions']
#     mybam.close()
#     return inversionDepth