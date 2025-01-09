import pysam
import argparse
import re

def calculateEndPos(read):
    """
    Calculate the end position of a read based on the CIGAR string.
    """
    if read.cigartuples is None:
        # Return None if CIGAR string is not available
        return None
    
    end_pos = read.reference_start  # Start with the reference start position
    for cigar_op, length in read.cigartuples:
        # CIGAR operations that consume reference positions
        if cigar_op in {0, 2, 3, 7, 8}:  # M, D, N, =, X
            end_pos += length
    return end_pos

def simplifyCIGAR(cigar, end, start):
    """
    Create CIGAR string with left and right clipping (and match)
    """
    left_match = re.match(r'^(\d+)([SH])', cigar)
    right_match = re.search(r'(\d+)([SH])$', cigar)
    
    left = f"{left_match.group(1)}{left_match.group(2)}" if left_match else ""
    right = f"{right_match.group(1)}{right_match.group(2)}" if right_match else ""
    
    cigar = f"{left}{end-start}M{right}"

    return cigar

def process_read(read):
    if read.is_unmapped:
        return 

    # Get the start position
    start = read.reference_start  # 0-based position for BED format
    # Get the end position based on the CIGAR string
    end = calculateEndPos(read)

    if not end:
        return 
    
    strand = "+" if not read.is_reverse else "-"
    flag = read.flag
    mismatch = read.get_tag("NM") if read.has_tag("NM") else None
    cigar = simplifyCIGAR(read.cigarstring, end, start)
    # Check for SA CIGAR strings 
    if read.has_tag("SA"): 
        alignments = read.get_tag("SA") 
        # Create SA tag for current alignment
        tag = f'{read.reference_name},{int(start)+1},{strand},{cigar},{read.mapping_quality},{mismatch}'
        # Get CIGAR string of the primary alignment (first in SA list for Sniffles)
        if read.is_supplementary: 
            primary = re.match(r'^([^;]+)', alignments).group(1)
        # Read is the primary alignment 
        else:
            primary = tag
        
        # Add the current alignment to list of all supplementary alignments
        alignments = f'{tag};{alignments}'
    else: 
        # No supplementary alignments -> primary and SA are the same
        alignments = f'{read.reference_name},{int(start)+1},{strand},{cigar},{read.mapping_quality},{mismatch}'
        primary = alignments
    
    if read.mapping_quality > 0:
        # Print read information in BED format
        print(f"{read.reference_name}\t{start}\t{end}\t{read.query_name}\t{read.mapping_quality}\t{strand}\t{flag}\t{mismatch}\t{primary}\t{alignments}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="inversion", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="BAM file")
    args = parser.parse_args()

    if not args.i:
        parser.error("Please specify the input BAM file with -i")

    # Open the BAM file
    bam_file = args.i
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # Iterate through each read in the BAM file
    for read in samfile:
        process_read(read)
    # Close the BAM file
    samfile.close()

