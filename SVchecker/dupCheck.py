import sys
import pandas as pd
import argparse
import csv


################################################################################
#                                 Parse Inputs                                 #
################################################################################
def readIntersects(intersects):
    """
    Reads in the ordered intersects file into a pandas dataframe
    """
    columns = []

    columns = ['duplicationID', 'dupStart', 'dupEnd', 'readID', 'readStart', 'readEnd', 'intersect', 'quality', 'strand', 'distance', 'dupChr','primary_tag', 'readOrder']
    df = pd.read_csv(intersects, sep='\t', names=columns, dtype=str, header=0)

    return df

def readDepth(depthFile):
    depth = {}
    with open(depthFile, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            key = row[3]
            meanSV = int(row[4])
            left = int(row[5])
            right = int(row[6])

            if left != 0:
                leftRatio = round(meanSV / left, 2)
            else:
                leftRatio = float('0')  # or some other value, e.g., 0 or a placeholder
            if right != 0:
                rightRatio = round(meanSV / right, 2)
            else:
                rightRatio = float('0')   
            depth[key] = {
                'mean': int(row[4]),
                'left': left,
                'right': right,
                'left_ratio': leftRatio, 
                'right_ratio': rightRatio,
                'chr': row[0],
                'start': row[1],
                'end': row[2],
            }
    return depth

################################################################################
#                           ASSESS SUPPORTING READS                            #
################################################################################

#------------------------------ Quality Checks---------------------------------#
def filter_reads(group, base_split, max_split, min_len, mapq):
    '''
    Filters out any reads that exceed the max number of splits or reads with single low quality alignment less
    input
    '''
    isSupport = True
    flag = ''

    # determine alignment lengths 
    readStart = group['readStart'].astype(int)
    readEnd = group['readEnd'].astype(int)
    alignment_lengths = readEnd - readStart
    
    # calculate splits
    max_read = int(group['readEnd'].max())
    min_read = int(group['readStart'].min())
    total_length =  max_read - min_read
    additional_splits = int((total_length / 1000) * max_split)
    allowed_splits = base_split + additional_splits
    num_splits = group['readOrder'].nunique() - 1

    # determine MAPQ
    group['quality'] = group['quality'].astype(int)  # Ensure quality is integer
    filtered_quality = group[group['quality'] >= mapq]
    total_reads = group['readOrder'].nunique()    
    num_reads = filtered_quality['readOrder'].nunique()    

    # All alignments below min --> reject
    if (alignment_lengths < min_len).all():
        isSupport = False
        flag = 'ALIGNMENT_LENGTH'
    # Alignments greater than max split --> reject
    elif num_splits > allowed_splits:
        isSupport = False
        flag = 'NUM_SPLIT'
    # Alignments have low MAPQ --> reject
    elif num_reads < 2 and total_reads > 1:
        isSupport = False
        flag = 'MAPQ'
    
    return isSupport, flag


#------------------------------- Creating Pairs -------------------------------#
def getAlignments(group):
    # Extract readOrder, intersect, and strand as tuples
    alignments = [(row['readOrder'], row['intersect'], row['strand']) for index, row in group.iterrows()]
    # Sort the list of tuples by readOrder
    alignments = sorted(alignments, key=lambda x: x[0])
    # Separate alignments into LS and RE lists
    LS = [(readOrder, strand) for readOrder, intersect, strand in alignments if intersect == 'LS']
    RE = [(readOrder, strand) for readOrder, intersect, strand in alignments if intersect == 'RE']

    return LS, RE


def getPair(set1, set2):
    # Get index of set 1
    # Switched strand --> list is empty 
    if not set1 or set1 == []:
        return None, set1, set2

    readOrder1 = set1[0][0]
    strand1 = set1[0][1]

    j = 0
    pair = None
    while j < len(set2) and (set2[j][0] <= readOrder1 or set2[j][1] != strand1):
        j += 1
    if j < len(set2):
        pair = ([readOrder1, set2[j][0]], strand1)
        set1.remove((readOrder1, strand1))
        set2.remove((set2[j][0], set2[j][1]))
    return pair, set1, set2

def getPairings(L, R):
    """
    Gets the complete set of pairings. Switches starting intersect when strand direction changes

    """
    pairs = []
    set1 = L.copy()
    set2 = R.copy()
    currentStrand = L[0][1]
    for readOrder, strand in L:
        # Change in strand direction -> switch pair direction
        if currentStrand != strand:
            pair, set1, set2 = getPair(set2, set1)
        
        # Same strand
        else:
            pair, set1, set2 = getPair(set1, set2)
        
        if pair: 
            pairs.append(pair)

        currentStrand = strand
    unused = set1 + set2
    return pairs, unused

def isConsecutive(readOrders):
    """
    Determiens if read order of alignments is consecutive or has breaks
 
    Input:
        readOrders (list): 
    Returns:
        (string):
    """
    # Sort read orders to handle unordered inputs
    readOrders = sorted(readOrders)

    # Check if each read order is same or consecutive e.g. (1, 2) or (1, 1, 2, 2)
    for i in range(len(readOrders) - 1):
        currReadOrder = int(readOrders[i])
        nextReadOrder = int(readOrders[i + 1])
        if currReadOrder == nextReadOrder or currReadOrder + 1 == nextReadOrder:
            continue
        else:
            return 'INTERSPERSED'
    return "TANDEM"

#--------------------------------- Main Cases ---------------------------------#
def isSupporting(group):
    """
    Determines if the alignments and breakpoints of the read support a duplication 
    1. Finds pairings
    2. Determines tandem/interspersed and duplication/repeat

    Input:
        readOrders (list): 
    Returns:
        (string):
    """
    support = False
    flag = ''

    # Convert groups into list of tuples (readOrder, strand)
        # L = [(1, '+'), (4, '-'), (5, '+')]
        # R = [(2, '+'), (3, '-'), (6, '+')]
    LS, RE = getAlignments(group)
    if not LS or not RE:
        # Other pairing e.g. LR, RS 
        flag = 'MISSING_INTERSECT'
    else: 
        sorted_group = group.sort_values(by='readOrder')
        # print(sorted_group[['readID', 'readStart', 'readEnd', 'quality', 'intersect', 'readOrder']].to_string(index=False))
        strand = group['strand'].iloc[0]
        if strand == '+':
            # Start with RE
            pairs, unused = getPairings(RE, LS)
        elif strand =='-':
            # Start with LS
            pairs, unused = getPairings(LS, RE)

        readOrders = []
        strands = set()
        for (reads, strand) in pairs:
            # List of read order e.g. [1, 2, 3, 3]
            readOrders.extend(reads)
            # Strand direction  ['+'], ['-'] or ['+', '-']
            strands.add(strand)

        # 1. Assess Read Order (Tandem/Interspersed)
        readOrderType = isConsecutive(readOrders)
        flag = readOrderType
        
        # 2. Count pairs
        numPairs = len(pairs)
        if numPairs == 0:
            flag = 'MISSING_PAIRS'
        elif numPairs == 1:
            flag = f'{flag}_DUPLICATION'
            support = True

        # 3. Check Duplex
        elif numPairs == 2 and '+' in strands and '-' in strands:
            flag = f'{flag}_DUPLICATION_DUPLEX'
            support = True
        else:
            support = True
            flag = f'{flag}_REPEAT_{len(pairs)}'

            firstRead = pairs[0][0][0]
            lastRead = pairs[-1][0][-1]
            unusedFiltered = [read for read in unused if firstRead < read[0] < lastRead]
            if unusedFiltered:
                maximal = len(pairs) + len(unusedFiltered)
                flag = f'{flag}.{maximal}'

        # print(f'{strand} {flag}')
    # print()
    return support, flag

################################################################################
#                               FLAG INFORMATION                               #
################################################################################
def countDuplicationFlags(details):
    reasonsRejected = details['FFLAG']
    passedFeature = details['SFLAG']

    dupFlagCount = {
        'DUPLEX': 0, 
        'TANDEM_DUPLICATION' :0,
        'TANDEM_REPEAT' :0,
        'INTERSPERSED_DUPLICATION':0, 
        'INTERSPERSED_REPEAT':0,
        'MISSING_PAIRS':0, 
        'MISSING_INTERSECT': 0,
        'NUM_SPLIT' : 0,
        'ALIGNMENT_LENGTH': 0, 
        'MAPQ' : 0,
    }

    # List of passed and rejected flags
    passedFlags = ['TANDEM_DUPLICATION', 'TANDEM_REPEAT', 'INTERSPERSED_DUPLICATION', 'INTERSPERSED_REPEAT', 'DUPLEX']
    rejectedFlags = ['MISSING_INTERSECT', 'MISSING_PAIRS', 'NUM_SPLIT', 'ALIGNMENT_LENGTH', 'MAPQ']

    # Count occurrences of passed flags
    for item in passedFeature:
        for reason in passedFlags:
            if reason in item:
                dupFlagCount[reason] += 1

    # Count occurrences in rejectedFlags
    for item in reasonsRejected:
        for reason in rejectedFlags:
            if reason in item:
                dupFlagCount[reason] += 1

    return dupFlagCount

################################################################################
#                                 OUTPUT FILES                                 #
################################################################################
def outputSummary(duplications, depth, summaryFile):
    for dupID, details in duplications.items():
        dupFlag = countDuplicationFlags(details)
        
        support = details['SUPPORT']
        rejected = details['NOTSUPPORT']
        
        # Coordinates
        chr= depth[dupID]['chr']
        start = depth[dupID]['start']
        end = depth[dupID]['end']

        # Depth Information 
        mean = depth[dupID]['mean']
        left = depth[dupID]['left']
        right = depth[dupID]['right']
        depthRatioLeft = depth[dupID]['left_ratio']
        depthRatioRight = depth[dupID]['right_ratio']

        rnames = ','.join(details['RNAMES']) if details['RNAMES'] else 'NA'

        rows.append([dupID, support, rnames, 
                    chr, start, end, mean, left, right, depthRatioLeft, depthRatioRight, 
                    dupFlag['TANDEM_DUPLICATION'], dupFlag['TANDEM_REPEAT'], dupFlag['INTERSPERSED_DUPLICATION'], dupFlag['INTERSPERSED_REPEAT'], 
                    rejected, dupFlag['MISSING_INTERSECT'], dupFlag['MISSING_PAIRS']])

    # Create and write to TSV
    df = pd.DataFrame(rows, columns=['ID', 'TOTAL_PASSED', 'RNAMES', 
                                     'chr', 'start', 'end', 'depth', 'leftDepth', 'rightDepth', 'leftDepthRatio', 'rightDepthRatio', 
                                     'TANDEM_DUPLICATION', 'TANDEM_REPEAT', 'INTERSPERSED_DUPLICATION', 'INTERSPERSED_REPEAT', 
                                     'TOTAL_REJECTED', 'MISSING_INTERSECT', 'MISSING_PAIRS'])
    if summaryFile:
        df.to_csv(summaryFile, sep='\t', header=False, index=False)
    else:
        print(df)

def outputReadDetails(duplications, readsFile):
    if readsFile:
        with open(readsFile, "w") as f:
            f.write("dupID\tread\tfilter\n")
            for dupID, details in duplications.items():
                readsPassed = details['RNAMES']
                readsFailed = details['RFAILED']
                reasonsSupport = details['SFLAG']
                reasonsRejected = details['FFLAG']
                for i, read in enumerate(readsPassed):
                    f.write(f"{dupID}\t{read}\t{reasonsSupport[i]}\n")
                for i, read in enumerate(readsFailed):
                    f.write(f"{dupID}\t{read}\t{reasonsRejected[i]}\n")
    
if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="duplication",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",
        help="path to input of intersected reads")
    parser.add_argument("-o",
        help="path to output file name [stdout]")
    parser.add_argument("-r",
        help="path to output file name for read details [stdout]")
    parser.add_argument("-d",
        help="depth bed file")
    
    # Caller default filtering parameters
    parser.add_argument("-mapq",
        help="Alignments with mapping quality lower than this value will be ignored", default=25, type=float)
        # CuteSV: Ignores reads that only report alignments with not longer than bp. - 20
    parser.add_argument("-min-alignment-length",
        help="Reads with alignments shorter than this length (in bp) will be ignored", default=1000, type=float)
        # CuteSV: Ignores reads that only report alignments with not longer than bp - 500
    parser.add_argument("-max-splits-kb",
        help="Additional number of splits per kilobase read sequence allowed before reads are ignored", default=0.1, type=float)
    parser.add_argument("-max-splits-base",
        help="Base number of splits allowed before reads ignored", default=3, type=float)

    args = parser.parse_args()
    
    # Read depth into dictionary 
    depth = readDepth(args.d)
    
    # Read all the overlaps into a dataframe
    df = readIntersects(args.i)
    # df['quality'] = df['quality'].astype(int)
    # df = df[df['quality'] >= args.mapq]


    keys = df['duplicationID']
    duplications = {key: {'RNAMES' : [], 'SFLAG' : [], 'RFAILED' : [], 'FFLAG' : [], 'SUPPORT' : 0, 'NOTSUPPORT' : 0} for key in keys.unique()}
    
    numDuplications = 0
    numReads = 0
    numPassed = 0
    numRejected = 0

    # Check duplications and their supporting reads 
    readGroups = df.groupby(['duplicationID', 'readID'])
    for (dupID, read), group in readGroups:
        # print(f'{dupID} {read}')
        numReads += 1
        isSupport, flag = filter_reads(group, args.max_splits_base, args.max_splits_kb, args.min_alignment_length, args.mapq)

        if isSupport:
            isSupport, flag = isSupporting(group)
        if isSupport:
            duplications[dupID]['RNAMES'].append(read)
            duplications[dupID]['SUPPORT'] += 1
            duplications[dupID]['SFLAG'].append(flag)
            numPassed += 1
        else: 
            duplications[dupID]['RFAILED'].append(read)
            duplications[dupID]['FFLAG'].append(flag)
            duplications[dupID]['NOTSUPPORT'] += 1
            numRejected += 1
        
    rows = []

    # Output files  
    numDuplications =  len(duplications)
    outputSummary(duplications, depth, args.o)
    outputReadDetails(duplications, args.r) if args.r else None
    print(f"Duplications Checked:      {numDuplications:>10}")
    print(f"Read Groups Checked:       {numReads:>10}")
    print(f"    Passed:                {numPassed:>10}")
    print(f"    Rejected:              {numRejected:>10}")


