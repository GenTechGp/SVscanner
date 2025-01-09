import sys
import pandas as pd
import argparse
import re

def getReadSuffix(row):
    """
    Determine which breakpoint of the read intersects the SV
    """
    svStart, svEnd, readStart, readEnd = row['svStart'], row['svEnd'], row['readStart'], row['readEnd']
    if svStart <= readStart <= svEnd:
        return 'S', readStart
    elif svStart <= readEnd <= svEnd:
        return 'E', readEnd
    else:
        return None, None


def readOverlaps(overlaps):
    """
    Parses the overlaps file determining whether the start/end of the alignment intersects with the SV
    """
    # Readname(7) ReadDirectionCode(e.g. LR) strand(9) position(start-5) quality(8)
    columns = ["svChr", "svStart", "svEnd", "svID", "readStart", "readEnd", "readID", "quality", "strand", "flag", "NM", "primary_tag", "SA_tag"]
    df = pd.read_csv(overlaps, sep='\s+', header=None, names=columns)

    ## Generate intersect type
    # Get the suffix - overlap information
    df[['readSuffix', 'readPos']] = df.apply(getReadSuffix, axis=1, result_type="expand")

    # Filter out read alignments where start/end are not within buffer region
    df = df[df['readSuffix'].notna()]
    sv_suffix = df['svID'].str[-1]
    # Remove the suffix from the svID and create column
    df['svID'] = df['svID'].str.replace(r'[.][LR]$', '', regex=True)

    # Combine the suffix with the read suffix to create the intersect column
    df['intersect'] = sv_suffix + df['readSuffix']

    # Calculate the distance between read and breakpoint
    # Original breakpoint (midpoint between the buffer value) minus the read start/end position 
    df['distance'] = abs(((df['svStart'] + df['svEnd']) // 2) - df['readPos']).astype(int)

    # Select the relevant columns
    df = df[['svID', 'svStart', 'svEnd', 'readID', 'readStart', 'readEnd', 'intersect', 'quality', 'strand', 'distance', 'svChr', 'NM', 'primary_tag', 'SA_tag']] 
    
    return df

def getClipping(cigar, strand):
    """
    Converts the left and right clipping of the CIGAR string to list, standardising the orientation 
    """
    left_match = re.match(r'^(\d+)([SH])', cigar)
    right_match = re.search(r'(\d+)([SH])$', cigar)
    
    # Extract values or default to 0 if not found
    left = int(left_match.group(1)) if left_match else 0
    right = int(right_match.group(1)) if right_match else 0

    if strand == '+':
        clipping = [left, right]
    elif strand == '-':
        clipping = [right,left]

    return clipping

def parseAlignment(alignment):
    alignmentInfo = alignment.split(',')
    return {
        "chr": alignmentInfo[0],
        "start": int(alignmentInfo[1]) - 1,  # Adjust start position by 1
        "strand": alignmentInfo[2],
        "cigar": alignmentInfo[3],
        "mapq": alignmentInfo[4],
        "mismatch": alignmentInfo[5] # Number of mismatches as part of unique identifier
    }

def getOrder(group):
    """
    Sorts the supplementary alignment tags based on CIGAR string, returning each alignment with its associated position 
    input 
        group (pandas) 
            svID 11111 22222 readID 11100 11122 + ... SA1;SA2
            svID 11111 22222 readID 22211 22233 + ... SA1;SA2
    output
        order (dict)    {SA2 : 1, SA1 : 2}
    """
    # Convert the supplementary alignment entries of the group (same SV and read) from pandas object to list
    # [SA1;SA2, SA1;SA2, SA1;SA2]
    SA_tags = group['SA_tag'].tolist()

    # Combine (if double tags come from overlap)
    SA_tags = list(set(SA_tags))
    order = None

    # Convert to list of supplementary alignments [SA1, SA2]
    suppAlignments = list(set(tag for tag in SA_tags[0].split(';') if tag)) # ensure no empty 

    alignmentsOrder = []
    
    # Convert supplementary alignments to read order
    if len(suppAlignments) > 1:
        for alignment in suppAlignments:
            info = parseAlignment(alignment)
            clipping = getClipping(info['cigar'], info['strand'])
            alignmentID = f'{info["chr"]},{info["start"]},{info["strand"]},{info["mapq"]},{info["mismatch"]}'
            alignmentsOrder.append((alignmentID, clipping))

        # Sort alignments based on the smallest left clip (first) and the smallest right clip (last)
        alignmentsOrder.sort(key=lambda x: (x[1][0], x[1][1]))
        order = { alignmentID: i+1 for i, (alignmentID, _) in enumerate(alignmentsOrder) }
    
    # Single Alignment - no read order change
    else:
        info = parseAlignment(suppAlignments[0])
        alignmentID = f'{info["chr"]},{info["start"]},{info["strand"]},{info["mapq"]},{info["mismatch"]}'
        order = {alignmentID : 1}

    return order

def outputOrdered(outputFile, result_df):
    """
    Outputs reads with the read order information 
        File 1: (args.o) outputs the SV details, read details with associated read order as integer
        File 2: (args.o.SA_tag) outputs the read followed by an ordered string of the supplementary tags
    """
    # Print the DataFrame with only the specified columns, and without the index
    columns = ['svID', 'svStart', 'svEnd', 'readID', 'readStart', 'readEnd', 'intersect', 'quality', 'strand', 'distance', 'svChr', 'primary_tag', 'readOrder']    
    result_df.to_csv(outputFile, sep='\t', columns=columns, index=False)

    # Print the Supplementary alignments in a separate file
    suppAlignmentsOrderColumns = ['readID', 'primary_tag', 'SA']
    unique_df = result_df.drop_duplicates(subset='readID')
    unique_df.to_csv(f'{outputFile}.SA_tag', sep='\t', columns=suppAlignmentsOrderColumns, index=False)


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="Order Supplementary Alignments ",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",
        help="path to input of intersected reads")
        # chr svStart svEnd svID.L/R readStart readEnd readID qual strand flag mismatch primary_tag SA_tags 
    parser.add_argument("-o",
        help="path to output file of ordered intersected reads [stdout]")
        # svID svStart svEnd readID readStart readEnd intersect[LS,RS,LE,RE] qual strand distance dupChr primary_tag readOrder

    args = parser.parse_args()

    # 1. Read all the overlaps into a dataframe
    df = readOverlaps(args.i)

    order_col = []

    # 2. Group alignments based on the SV and readID
    readGroups = df.groupby(['svID', 'readID'])

    # 3. Get read order for each read group 
    for (svID, read), group in readGroups:
        order = getOrder(group)
        orderedTags = ";".join(order.keys())

        # Add the order to the corresponding DataFrame rows
        if order:
            group['readOrder'] = group.apply(
                # Create key - based on coordinates, quality and mismatch to fetch order assocaited with SA tag
                lambda row: order.get(
                    (f"{row['svChr']},{row['readStart']},{row['strand']},{row['quality']},{row['NM']}")
                ), 
            axis=1)
        else:
            group['readOrder'] = 1

        order_col.append(group)

        # Order the original supplementary alignments tag string
        group['SA'] = orderedTags
    
    # 4. Add the ordered supplementary alignment numbers to the original dataframe
    result_df = pd.concat(order_col, ignore_index=True)

    # 5. Output the updated read order
    outputOrdered(args.o, result_df)

