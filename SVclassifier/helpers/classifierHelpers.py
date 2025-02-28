import csv

def read_sv_info(sv_file):
    sv_info = {}
    with open(sv_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for sv in reader:
            # if sv[6] != "INS.36":
            #     continue
            # Extract values from the columns
            chrom = sv[0]
            query_start = int(sv[1])
            query_end = int(sv[2])
            pos = int(sv[3])
            end = int(sv[4])
            length = int(sv[5])
            sv_id = sv[6]
            callerID = sv[7]
                    
            rel_start = pos - query_start
            rel_end = rel_start + length 
            
            # Add the values to the dictionary
            sv_info[sv_id] = {                
                'ID': sv_id, 
                'chrom': chrom,
                'position': pos,
                # 'end' : end, 
                'length': length, 
                'rel_start' : rel_start, 
                'rel_end' : rel_end,
                'start' : 0, 
                'end' : query_end - query_start,
                'callerID' : callerID,  
                'header' : f'{sv_id}\t({callerID})\t{chrom}:{query_start}-{query_end}\t{chrom}:{pos}-{end}\t{length}\t',
            }

    return sv_info

def calculate_intersection_length(repeat_start, repeat_end, rel_start, rel_end):
    overlap_start = max(rel_start, repeat_start)
    overlap_end = min(rel_end, repeat_end)
    
    if overlap_start < overlap_end:
        intersection = overlap_end - overlap_start
    else:
        intersection = 0
    
    return intersection

def get_RM_classification(family, isVisualise):
    if isVisualise:
        classifications = ['SINE', 'LINE', 'LTR', 'DNA', 'Simple_repeat', 'RNA', 'Retroposon', 'Low_complexity', 'Other']
    else: 
        classifications = ['SINE', 'LINE', 'LTR', 'DNA', 'Retroposon']
    
    for classification in classifications:
        if classification in family:
            return classification

    return None

def get_TRF_classification(period_size):
    # Classify as HOMO, STR, or TR based on period size 
    if period_size == 1: 
        key = 'HOMO'
    elif period_size >= 2 and period_size <= 12:
        key = 'STR'
    elif period_size > 12:
        key = 'TR'
    return key 



def parse_TRF_lines(lines, sv, min_intersect):
    """
    Determines which tandem repeats entries intersect with the SV from calculating the intersection 
    """
    tandem_repeats = []

    for trf_line in lines:
        trf_data = trf_line.split()

        repeat_start = int(trf_data[0])
        repeat_end = int(trf_data[1])
        period_size = int(trf_data[2])
        copy_number = float(trf_data[3])
        consensus_repeat = trf_data[13]
        
        intersection = calculate_intersection_length(repeat_start, repeat_end, sv['rel_start'], sv['rel_end'])
        intersect_fraction = round(intersection/ sv['length'], 6)

        repeat_info = {
            'repeat_start': repeat_start,
            'repeat_end': repeat_end,
            'period_size': period_size,
            'motif': consensus_repeat,

            'copy_number': copy_number,
            'key': get_TRF_classification(period_size),
        }

        if intersect_fraction and float(intersect_fraction) > min_intersect:
            repeat_info.update({'intersection' : intersect_fraction})
            tandem_repeats.append(repeat_info)
    
    return tandem_repeats


def read_trf(sv_info, trf_file, min_intersect):
    """
    Processes Tandem Repeat Finder output (.dat), adding relevant entries to sv_info
    
    Input:
    - sv_info (dict)            Dictionary containing the strucutral variants
    - trf_file (str):           Path to TRF output file
    - min_intersect (float):    Minimum intersection fraction between repeat and SV

    Output: 
    - sv_info (dict):           Updated sv_info with TRF annotations
    """

    trf_lines = None

    # @INS.1
    # 1040 1111 23 3.4 23 67 25 68 9 18 50 22 1.76 GAGGGCGTCTGGTCGTCCTGAGG GAGGGCGTCTGGTCGTCCTGAGGGAGGGCCGGTGTTGGTGAGGGCATCTGGTCGTCCTGAGGGAGGGGGTCT GGTGAGAGACGCTGCCGCAGAGCCGCCCGAGAGGGAGGGTCAGTGTTGGT TCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATTTTTTATTTT
    # 1020 1106 39 2.2 39 93 0 147 11 14 51 21 1.74 GAGGGAGGGCCAGTGTTGGTGAGGGCATCTGGTCGTCCT GAGGGAGGGTCAGTGTTGGTGAGGGCGTCTGGTCGTCCTGAGGGAGGGCCGGTGTTGGTGAGGGCATCTGGTCGTCCTGAGGGAGGG CACAGAGGGAAACAAGGGGAGGTGAGAGACGCTGCCGCAGAGCCGCCCGA GGTCTTCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATTTTTT
    # 1153 1200 5 9.8 5 95 4 89 18 0 0 81 0.70 TTTTA TTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTATTTT AGGGGGTCTTCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATT GAAACGGAGTCTCACTCTTGCCTAGGCTGGAGTGCAATGGCGCAATCTCG
    
    with open(trf_file, 'r') as file:
        trf_text = file.read()

    # Split each SV into list
    entries = trf_text.strip().split('@')[1:]   
   
    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        trf_lines = lines[1:]

        sv_id = header.lstrip('@')

        # Get the relevant start, end ... for the sv
        if sv_id in sv_info:
            sv = sv_info[sv_id]
            # Parse and add the TRFs that overlap 
            trf_result = parse_TRF_lines(trf_lines, sv, min_intersect)

            sv_info[sv_id]['TRF'] = trf_result
    
    return sv_info


def positionInRepeatFraction(strand, repeat_begin, repeat_end, repeat_left, intersection):
    """
    Calculates the sv in relation to the repeat
    - Proportion (entire query)
    - fraction (sv only)
    
    """
    # $12           #13         $14 
    # repeat_begin  repeat_end  repeat_left

    proportion = None
    fraction = None 
    if strand == '+':
        # ($13-$12)/(($13-$12)+$12+$14)
        proportion = (repeat_end - repeat_begin) / ((repeat_end - repeat_begin) + repeat_left)
        fraction = intersection / ((repeat_end - repeat_begin) + repeat_left)
    elif strand == 'C':    
        # ($13-$14)/(($13-$14)+$14+$12)
        proportion = (repeat_end - repeat_left) / ((repeat_begin + repeat_end) - repeat_left)
        fraction = intersection / ((repeat_begin + repeat_end) - repeat_left)

    return proportion, fraction 

def read_rm(sv_info, rm_file, min_intersect, isVisualise):
    """
    Processes RepeatMasker output (.out), adding relevant entries to sv_info
    
    Input:
    - sv_info (dict)            Dictionary containing the strucutral variants
    - rm_file (str):            Path to RepeatMasker output file
    - min_intersect (float):    Minimum intersection fraction between repeat and SV
    - isVisualise (bool):       True (visualisation)/ False (annotation)

    Output: 
    - sv_info (dict):           Updated sv_info with repeatMasker annotations
    """
    with open(rm_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        for row in reader:
            sv_id = row[4]
            te_start = int(row[5])
            te_end = int(row[6])
            strand = row[8]
            repeat = row[9]
            family = row[10]
            te_id = row[14]

            classification = get_RM_classification(family, isVisualise)

            # Add element if ['SINE', 'LINE', 'LTR', 'DNA', 'Retroposon']
            if classification and sv_id in sv_info:
                
                # Remove the parentheses around the repeat
                repeat_begin = int(row[11].strip("()"))
                repeat_end = int(row[12].strip("()"))
                repeat_left = int(row[13].strip("()"))
               
                sv = sv_info[sv_id]
                
                # Add if the transposable element overlaps with the SV
                intersection = calculate_intersection_length(te_start, te_end, sv['rel_start'], sv['rel_end'])
                intersect_fraction = round(intersection / sv['length'], 6)
                element_proportion, element_fraction = positionInRepeatFraction(strand, repeat_begin, repeat_end, repeat_left, intersection)
                
            
                # Add if the intersection is > min_intersect (e.g. 5%)
                if intersect_fraction and intersect_fraction > min_intersect:
                    element = {
                        'te_start': te_start,
                        'te_end': te_end,
                        'family': family,
                        'repeat': repeat,
                        'element_coverage' : element_fraction,
                        'element_proportion' : element_proportion,
                        'class': classification, 
                        'intersection' : intersect_fraction, 
                        'te_id' : te_id, 
                    }

                    # First element for SV - add the 'RM'
                    if 'RM' not in sv_info[sv_id]:
                        sv_info[sv_id]['RM'] = {}
                    # First element of class for SV
                    if classification not in sv_info[sv_id]['RM']:
                        sv_info[sv_id]['RM'][classification] = []
                    
                    sv_info[sv_id]['RM'][classification].append(element)

    return sv_info