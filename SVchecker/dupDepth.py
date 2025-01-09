import pysam
import pysamstats
import sys
import argparse
import csv
from statistics import mean 
from tqdm import tqdm

# Calculate coverage within duplication (start-end)
# Calculate left flanking region (2*size of duplication)
# Calculate right flanking region (2*size of duplication)

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
                dupID = row[3]
                duplications.append((chrom, start, end, dupID))
                duplicationDepth.update({dupID: {'chr' : chrom, 'start': start, 'end' : end, 'dupDepth' : 0, 'leftDepth' : 0, 'rightDepth' : 0}})
    return duplications, duplicationDepth

def calculate_read_depth(duplicationDepth, bam_file, duplications, fa_file):
    mybam = pysam.AlignmentFile(bam_file, "rb")
    
    for chrom, start, end, dupID in tqdm(duplications, desc="duplications"):
        flanking = 2*(end - start) 
        dupDepths = []
        leftFlankingDepths = []
        rightFlankingDepths = []

        leftFlankingStart = max(0, start - flanking - 1) # Set to 0 if less than 0
        rightFlankingEnd = end + flanking
        for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=start-1, end=end-1, truncate=True, fafile=fa_file):
            dupDepths.append(stat['reads_all']-stat['deletions'])
        for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=leftFlankingStart, end=start, truncate=True, fafile=fa_file):
            leftFlankingDepths.append(stat['reads_all']-stat['deletions'])
        for stat in pysamstats.stat_variation(mybam, chrom=chrom, start=end-1, end=rightFlankingEnd, truncate=True, fafile=fa_file):
            rightFlankingDepths.append(stat['reads_all']-stat['deletions'])
        
        if leftFlankingDepths != []:
            leftDepth = int(mean(leftFlankingDepths))
            duplicationDepth[dupID]['leftDepth'] = leftDepth
        if dupDepths != []:
            dupDepth = int(mean(dupDepths))
            duplicationDepth[dupID]['dupDepth'] = dupDepth
        if rightFlankingDepths != []:
            rightDepth = int(mean(rightFlankingDepths))
            duplicationDepth[dupID]['rightDepth'] = rightDepth

    mybam.close()
    return duplicationDepth


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="duplication depth",
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

    duplications, duplicationDepth = read_bed_file(args.bed)
    # print(duplications)

    duplicationDepth = calculate_read_depth(duplicationDepth, bam_file, duplications, args.ref)

    # Open the CSV file in write mode with tab delimiter and no header
    with open(args.o, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        # Write rows without header
        for duplication in duplicationDepth:
            chr_value = duplicationDepth[duplication]['chr']
            start_value = duplicationDepth[duplication]['start']
            end_value = duplicationDepth[duplication]['end']
            meanLeft = duplicationDepth[duplication]['leftDepth'] if 'leftDepth' in duplicationDepth[duplication] else ''
            meanDup =  duplicationDepth[duplication]['dupDepth'] if 'dupDepth' in duplicationDepth[duplication] else ''
            meanRight = duplicationDepth[duplication]['rightDepth'] if 'rightDepth' in duplicationDepth[duplication] else ''
            writer.writerow([chr_value, start_value, end_value, duplication, meanDup, meanLeft, meanRight])
