import argparse
import csv
import sys
from tqdm import tqdm
import pysam

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
                invID = row[3] if len(row) > 3 else f"{chrom}:{start}-{end}"  # Handle missing ID
                inversions.append((chrom, start, end, invID))
                inversionDepth[invID] = {
                    'chr': chrom, 'start': start, 'end': end,
                    'startDepth': 0, 'endDepth': 0
                }
    return inversions, inversionDepth

def calculate_region_depth(mybam, chrom, start, end):
    inv_del = 0
    depth = 0
    for pileupcolumn in mybam.pileup(chrom, start, end, truncate=True, min_base_quality=0):
        # Adjust depth calculation based on whether the base is deleted or not
        for read in pileupcolumn.pileups:
            if read.is_del or read.is_refskip: 
                inv_del += 1
        depth = pileupcolumn.nsegments - inv_del
    return depth

def calculate_read_depth(inversionDepth, bam_file, inversions):
    mybam = pysam.AlignmentFile(bam_file, "rb")
    
    for chrom, start, end, invID in tqdm(inversions, desc="Calculating read depth for inversions"):
        # Version 1 
        inversionDepth[invID]['startDepth'] = calculate_region_depth(mybam, chrom, start-1, start)
        inversionDepth[invID]['endDepth'] = calculate_region_depth(mybam, chrom, end-1, end)

        # Version 2 
        # print(f"{inversionDepth[invID]['startDepth']} {inversionDepth[invID]['endDepth']}")
        # start_depth_array = mybam.count_coverage(chrom, start-1, start, quality_threshold = 0)
        # end_depth_array = mybam.count_coverage(chrom, end-1, end, quality_threshold = 0)
        # start_depth = sum(arr[0] for arr in start_depth_array)
        # end_depth = sum(arr[0] for arr in end_depth_array)
        # inversionDepth[invID]['startDepth'] = start_depth
        # inversionDepth[invID]['endDepth'] = end_depth

    mybam.close()
    return inversionDepth


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="Inversion depth calculation",
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-bed",
                        help="BED file containing inversions", required=True)
    parser.add_argument("-bam",
                        help="BAM file containing aligned reads", required=True)
    parser.add_argument("-o",
                        help="Output file to store calculated depths", required=True)


    args = parser.parse_args()
    bam_file = args.bam

    inversions, inversionDepth = read_bed_file(args.bed)

    inversionDepth = calculate_read_depth(inversionDepth, bam_file, inversions)

    # Write results to the output file
    with open(args.o, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        # Write rows without a header
        for inversion in inversionDepth:
            chr_value = inversionDepth[inversion]['chr']
            start_value = inversionDepth[inversion]['start']
            end_value = inversionDepth[inversion]['end']
            start_depth = inversionDepth[inversion]['startDepth']
            end_depth = inversionDepth[inversion]['endDepth']
            writer.writerow([chr_value, start_value, end_value, inversion, start_depth, end_depth])
