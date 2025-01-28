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

    columns = ['invID', 'invStart', 'invEnd', 'readID', 'readStart', 'readEnd', 'intersect', 'quality', 'strand', 'distance', 'invChr','primary_tag', 'readOrder']
    df = pd.read_csv(intersects, sep='\t', names=columns, dtype=str, header=0)

    return df

def readDepth(depthFile):
    depth = {}
    with open(depthFile, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            key = row[3]
            left_bp = int(row[4])
            right_bp = int(row[5])

            depth[key] = {
                'left_bp': left_bp,
                'right_bp': right_bp,
                'chr': row[0],
                'start': row[1],
                'end': row[2],
            }
    return depth

################################################################################
#                         Check Supporting Breakpoints                         #
################################################################################

#------------------------------ Quality Checks---------------------------------#
def filter_reads(group, base_split, max_split, min_len):

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
        
    # All alignments below min --> reject
    if (alignment_lengths < min_len).all():
        isSupport = False
        flag = 'ALIGNMENT_LENGTH'
    # Alignments greater than max split --> reject
    elif num_splits > allowed_splits:
        isSupport = False
        flag = 'NUM_SPLIT'

    return isSupport, flag
    


def checkPrimary(read):
    """
    1. Alignments are in the same region as the primary chromosome 
    2. Primary chromosome has MAPQ > 25
    """
    primaryPass = True
    # Check if 'invChr' and 'primaryChr' are the same for each row
    # print(read)
    primary_tag = read['primary_tag']
    primary_tag = read['primary_tag'].iloc[0].split(',')
    primaryChr = primary_tag[0]
    quality = primary_tag[-2]

    readChr = read['invChr'].iloc[0]
    # print(f'Read {readChr}, Primary {primaryChr}')

    if readChr != primaryChr: 
        primaryPass = False
    elif int(quality) < 25:
        primaryPass = False
    
    # print(isMatchPrimary)
    return primaryPass

def getAlignmentInfo(read):
    numAlignment = read['numAlignment'].iloc[0]
    chromosomes = read['alignedChr'].iloc[0]

    numAlignedChr = len(chromosomes.split(','))

    return numAlignment, numAlignedChr


#------------------------------- Pairs & Opposite------------------------------#
# Function to check inversion support for each group
def getStrands(group):
    """
        input:
            group       pandas.core.frame.DataFrame
                e.g.
                          invID                                readID intersect  quality strand  distance
                Sniffles2.INV.119BCS1  55b85f4e-1caa-5365-83b5-96cdb0795fb1        LE       60      -         1
                Sniffles2.INV.119BCS1  55b85f4e-1caa-5365-83b5-96cdb0795fb1        LE       60      +         5
                Sniffles2.INV.119BCS1  55b85f4e-1caa-5365-83b5-96cdb0795fb1        LS       60      +         1
                Sniffles2.INV.119BCS1  55b85f4e-1caa-5365-83b5-96cdb0795fb1        RE       60      +         1
                Sniffles2.INV.119BCS1  55b85f4e-1caa-5365-83b5-96cdb0795fb1        RS       60      -         1       
        output: 
            intersects  dict  e.g. {'LE': ['-', '+'], 'LS': ['+'], 'RS': ['-'], 'RE': ['+']}
    """
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)  # Display all columns
    pd.set_option('display.width', 1000)  # Set the width to a large enough value to prevent wrapping
    pd.set_option('display.max_colwidth', None)  # Set max column width to prevent truncation

    pd.set_option('display.max_columns', None)
    # print(group.to_string(index=False))
    intersects = {}
    intersectCombo = ['LE', 'LS', 'RS', 'RE']
    multiIntersects = []
    isMulti = None

    for combo in intersectCombo:
        if any(group['intersect'] == combo):
            filtered_df = group[group['intersect'] == combo]
            # Filter for check of the strand for the highest quality alignment first            
            filtered_df = filtered_df.sort_values(by='quality', ascending=False)

            strands = filtered_df['strand'].tolist()
            starts = filtered_df['readStart'].tolist()
            ends = filtered_df['readEnd'].tolist()
            invStarts = filtered_df['invStart'].tolist()
            invEnds = filtered_df['invEnd'].tolist()

            # strand = sorted_df.iloc[0]['strand']

            intersects[combo] = {'strand': strands, 'invStart' : invStarts, 'invEnd' : invEnds, 'readStart': starts, 'readEnd': ends}
    
            if len(strands) > 1:
                multiIntersects.append(combo)

    if multiIntersects != []:
        isMulti = "MULTI-" + ",".join(multiIntersects) 


    return intersects, isMulti


def hasOppositeStrand(intersect1, intersect2):
    check = {'+' : '-', '-' : '+'}
    opposite = False
    for strand in intersect1:
        if check[strand] in intersect2:
            opposite = True
            break
    return opposite

def isRightSupporting(LE, RE):
    """
    Buffer         |---|                |---|
    LE        -------|
    RE                               -----| 		
    """        
    support = False 
    flag = None
    for i in range(len(LE['readStart'])):
        if LE['readStart'][i] < LE['invStart'][i]:
            support = hasOppositeStrand(LE['strand'], RE['strand'])
            if not support:
                flag = 'ORIENTATION'
                break
        else: 
            flag = 'ERROR'
    
    return support, flag


def isLeftSupporting(RS, LS):
    """
    Buffer         |---|                |---|
    LS                |------- 					
    RS                                     |------- 	
    """
    support = False 
    flag = None
    for i in range(len(RS['readEnd'])):
        if RS['readEnd'][i] > RS['invEnd'][i]:
            support = hasOppositeStrand(RS['strand'], LS['strand'])
            if not support:
                flag = 'ORIENTATION'
                break
        else: 
            flag = 'ERROR'

    return support, flag

def isLSREsame(LS, RE):
    same = True
    if LS['readStart'] != RE['readStart'] or LS['readEnd'] != RE['readEnd']:
        same = False
    return same


#---------------------------- Additional Annotations --------------------------#

def isREExtension(RE, LE):
    """
    Buffer         |---|                |---|
    RE      ------------------------------| 					
    LE         ------|                            	
    """
    isExtension = False
    # Check RE (end of middle alignment) extends left beyond inversion start (RE_start < SV_left)

    for i in range(len(RE['readStart'])):
        if RE['readStart'][i] < LE['invStart'][0]:
            isExtension = True
            break
    return isExtension

def isLSExtension(LS, RS):
    """
    Buffer         |---|                |---|
    LS               |---------------------------------- 					
    RS                                    |----                                   	
    """
    isExtension = False
    # Check if LS (start of middle alignment) extends right beyond inversion end (LS_end > SV_right)
    for i in range(len(LS['readEnd'])):
        if LS['readEnd'][i] > RS['invEnd'][0]:
            isExtension = True
            break
    return isExtension

#--------------------------------- Main Cases ---------------------------------#
def isSupporting(group):
    support = False
    supportLeft = False
    supportRight = False
    flag = None
    flagLeft = None
    flagRight = None

    # ## FILTER OUT READS
    if not checkPrimary(group):
        return False, 'PRIMARY_CHR'
    
    # Get the strand direction of the overlapping reads 
    intersects, isMulti = getStrands(group)

    if len(intersects) <= 1:
        return False, 'MISSING_BOTH'
    
    ## CASE 1: BOTH BP - Read end overlaps left breakpoint and right breakpoint
    # if 'LE' in intersects and 'RE' in intersects and 'LS' in intersects and 'RS' in intersects: 
    if set(intersects.keys()) == {'LE', 'RE', 'LS', 'RS'}:
        supportRight, flagRight = isRightSupporting(intersects['LE'], intersects['RE'])
        supportLeft, flagLeft = isLeftSupporting(intersects['RS'], intersects['LS'])

        # Check whether any reads extend beyond BP
        if isREExtension(intersects['RE'], intersects['LE']):
            flagRight = 'EXTENSION_RE'
        if isLSExtension(intersects['LS'], intersects['RS']):
            flagLeft = 'EXTENSION_LS'

        if supportLeft and supportRight:
            support = True
        # Create Flag
        if flagLeft is None and flagRight is None:
            flag = 'BOTH'
        elif flagLeft is not None and flagRight is None:
            flag = flagLeft
        elif flagLeft is None and flagRight is not None:
            flag = flagRight
        else:
            flag = flagLeft + ',' + flagRight

    
    ## CASE 2: RIGHT / END OF READS
    elif set(intersects.keys()) == {'LE', 'RE'} or set(intersects.keys()) == {'LE', 'RE', 'LS'}:
        support, flagRight = isRightSupporting(intersects['LE'], intersects['RE'])
    
        if 'LS' in intersects:
            if isLSREsame(intersects['LS'], intersects['RE']): 
                flagRight = 'EXTENSION_RE2LS'
            else: 
                support = False
                flagRight = 'MISSING_RIGHT_BP'
        else:
            if isREExtension(intersects['RE'], intersects['LE']):
                flagRight = 'EXTENSION_RE'

        # Add additional Info
        if support and flagRight is None:
            flag = 'RIGHT' # No additional info 
        else:
            flag = flagRight
               
    ## CASE 3: LEFT BP / START OF READS
    elif set(intersects.keys()) == {'LS', 'RS'} or set(intersects.keys()) == {'LS', 'RS', 'RE'}:
        support, flagLeft = isLeftSupporting(intersects['RS'], intersects['LS'])
        
        if 'RE' in intersects:
            if isLSREsame(intersects['LS'], intersects['RE']): 
                flagLeft = 'EXTENSION_LS2RE'
            else: 
                support = False
                flagLeft = 'MISSING_LEFT_BP'
        else:
            if isLSExtension(intersects['LS'], intersects['RS']):
                flagLeft = 'EXTENSION_LS'

        # Add additional Info
        if support and flagLeft is None:
            flag = 'LEFT' # No additional info 
        else:
            flag = flagLeft
    

    ## OTHER
    elif set(intersects.keys()) == {'RE', 'RS'} or set(intersects.keys()) == {'LE', 'RE', 'RS'}:
        """        
        Buffer         |---|                |---|
        RE                             -------| 					
        RS                                    |-------
        """ 	
        flag = 'MISSING_LEFT_BP'
    elif set(intersects.keys()) == {'LE', 'LS'} or set(intersects.keys()) == {'LE', 'LS', 'RS'}:
        """       
        Buffer         |---|                |---|
        LE        -------|
        LS               |-----
        """
        flag = 'MISSING_RIGHT_BP'
    elif set(intersects.keys()) == {'LS', 'RE'} or set(intersects.keys()) == {'LE', 'RS'}:
        flag = 'MISSING_BOTH'

    if flag == None:
        print(group)
        print(support)
        print("NO FLAG")

    if isMulti is not None:
        flag = f"{flag},{isMulti}"
    return support, flag

################################################################################
#                               FLAG INFORMATION                               #
################################################################################
def addReadReason(readFlagsCount, flag):
    if flag:
        if flag.startswith("MISSING"):
            readFlagsCount['MISSING'] += 1
        elif flag.startswith("ORIENTATION"):
            readFlagsCount['ORIENTATION'] += 1
        elif flag.startswith("NUM_SPLIT"):
            readFlagsCount['NUM_SPLIT'] += 1
        elif flag.startswith("ALIGNMENT_LENGTH"):
            readFlagsCount['ALIGNMENT_LENGTH'] += 1
        else:
            readFlagsCount[flag] += 1
    return readFlagsCount

def countInversionFlags(details):
    reasonsRejected = details['FFLAG']
    passedFeature = details['SFLAG']
    
    inversionFlagCount = {
        'MISSING_LEFT_BP': 0,
        'MISSING_RIGHT_BP': 0,
        'MISSING_BOTH': 0,
        'EXTENSION': 0,
        'MULTI': 0,
        'ORIENTATION': 0,
        'NUM_SPLIT' : 0,
        'ALIGNMENT_LENGTH' : 0, 
    }

    # List of passed and rejected flags
    passedFlags = ['EXTENSION', 'MULTI']
    rejectedFlags = ['MISSING_LEFT_BP', 'MISSING_RIGHT_BP', 'MISSING_BOTH', 'ORIENTATION', 'NUM_SPLIT', 'ALIGNMENT_LENGTH']

    # Count occurrences of passed flags
    for item in passedFeature:
        for reason in passedFlags:
            if reason in item:
                inversionFlagCount[reason] += 1

    # Count occurrences in rejectedFlags
    for item in reasonsRejected:
        for reason in rejectedFlags:
            if reason in item:
                inversionFlagCount[reason] += 1

    return inversionFlagCount

################################################################################
#                                 OUTPUT FILES                                 #
################################################################################

def outputSummary(inversions, depth, summaryFile):
    for invID, details in inversions.items():
        invFlag = countInversionFlags(details)

        support = details['SUPPORT']
        rejected = details['NOTSUPPORT']
        
        # Coordinates
        chr= depth[invID]['chr']
        start = depth[invID]['start']
        end = depth[invID]['end']

        # Depth Information 
        left = depth[invID]['left_bp']
        right = depth[invID]['right_bp']

        rnames = ','.join(details['RNAMES']) if details['RNAMES'] else 'NA'

        rows.append([invID, support, rnames, 
                    chr, start, end, left, right,
                    invFlag['EXTENSION'], invFlag['MULTI'], 
                    rejected, invFlag['MISSING_LEFT_BP'], invFlag['MISSING_RIGHT_BP'], invFlag['MISSING_BOTH'], invFlag['ORIENTATION']])

        # -->ID	chr	start	end	startDepth	endDepth	CALSER_SUPPOER TOTAL_PASSED EXTENTION MULTI TOTAL_REJECTED LEFT RIGHT BOTH ORIENTATION FILTER

    # Create and write to TSV
    df = pd.DataFrame(rows, columns=['ID', 'TOTAL_PASSED', 'RNAMES', 
                                     'chr', 'start', 'end', 'leftDepth', 'rightDepth',
                                     'EXTENSION', 'MULTI', 
                                     'TOTAL_REJECTED', 'MISSING_LEFT_BP', 'MISSING_RIGHT_BP', 'MISSING_BOTH', 'ORIENTATION'])
    if summaryFile:
        df.to_csv(summaryFile, sep='\t', header=False, index=False)
    else:
        print(df)

def outputReadDetails(inversions, readsFile):
    if readsFile:
        with open(readsFile, "w") as f:
            f.write("invID\tread\tfilter\n")
            for invID, details in inversions.items():
                readsPassed = details['RNAMES']
                readsFailed = details['RFAILED']
                reasonsSupport = details['SFLAG']
                reasonsRejected = details['FFLAG']
                for i, read in enumerate(readsPassed):
                    f.write(f"{invID}\t{read}\t{reasonsSupport[i]}\n")
                for i, read in enumerate(readsFailed):
                    f.write(f"{invID}\t{read}\t{reasonsRejected[i]}\n")


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="inversion",
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
        help="Alignments with mapping quality lower than this value will be ignored", default=25)
        # CuteSV: Minimum mapping quality value of alignment to be taken into account.	 - 20
    parser.add_argument("-min-alignment-length",
        help="Reads with alignments shorter than this length (in bp) will be ignored", default=1000)
        # CuteSV: Ignores reads that only report alignments with not longer than bp - 500
    parser.add_argument("-max-splits-kb",
        help="Additional number of splits per kilobase read sequence allowed before reads are ignored", default=0.1)
        # CuteSV: Maximum number of split segments a read may be aligned before it is ignored.
    parser.add_argument("-max-splits-base",
        help="Base number of splits allowed before reads ignored", default=3)
        # CuteSV: Maximum number of split segments a read may be aligned before it is ignored. 7
    

    args = parser.parse_args()

    # Read depth into dictionary 
    depth = readDepth(args.d)

    # Read all the overlaps into a dataframe
    df = readIntersects(args.i)
    numOverlaps = len(df)

    # # Remove reads with quality below 25
    df['quality'] = df['quality'].astype(int)
    df = df[df['quality'] >= args.mapq]

    pd.set_option('display.max_rows', None)


    numOverlapsFiltered = len(df)

    numInversions = 0
    numReads = 0
    numPassed = 0
    numRejected = 0
    readFlagsCount = {'ORIENTATION': 0, 
            'MISSING' : 0,
            # 'EXTENSION' : 0, 
            'PRIMARY_CHR' :0, 
            'NUM_SPLIT' : 0, 
            'ALIGNMENT_LENGTH' : 0}
    
    keys = df['invID']
    inversions = {key: {'RNAMES' : [], 'SFLAG' : [], 'RFAILED' : [], 'FFLAG' : [], 'SUPPORT' : 0, 'NOTSUPPORT' : 0} for key in keys.unique()}


    # Check each read group (primary & supplementary) against inversion it overlaps 
    readGroups = df.groupby(['invID', 'readID'])
    
    for (invID, read), group in readGroups:
        numReads += 1
        # If the read supports an inversion
        isSupport, flag = filter_reads(group, args.max_splits_base, args.max_splits_kb, args.min_alignment_length)
        if isSupport:
            isSupport, flag = isSupporting(group)
        if isSupport:
            flag = f"PASSED_{flag}"
            # Add the read name 
            inversions[invID]['RNAMES'].append(read)
            inversions[invID]['SUPPORT'] += 1
            inversions[invID]['SFLAG'].append(flag)
            numPassed += 1
        else: 
            readFlagsCount = addReadReason(readFlagsCount, flag)
            # readFlagsCount[reasonFailed] += 1
            inversions[invID]['RFAILED'].append(read)
            inversions[invID]['FFLAG'].append(flag)
            inversions[invID]['NOTSUPPORT'] += 1
            numRejected += 1
        
    rows = []
    outputSummary(inversions, depth, args.o)
    

    numInversions = len(inversions)
    # Printing with aligned formatting
    print(f"Inversion Checked:         {numInversions:>10}")
    print(f"Reads Overlapping          {numOverlaps:>10}")
    print(f"Reads Overlapping (MAPQ>25){numOverlapsFiltered:>10}")
    print(f"Read Groups Checked:       {numReads:>10}")
    print(f"    Passed:                {numPassed:>10}")
    print(f"    Rejected:              {numRejected:>10}")
    for reason, count in readFlagsCount.items():
        print(f"        {reason:<22}{count:>7}")




    #                 IL                  IR
    #                  |-------------------|
    # Buffer         |---|                |---|
    # Read LE   -------|
    # Read LS          |-------
    # Read RE                          -----|
    # Read RS                               |--------

    # LE    +   -
    # LS    -   +
    # RE    -   +
    # RS    +   -
