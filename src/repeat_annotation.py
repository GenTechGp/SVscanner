import pysam
import sys
import argparse
import os
from Bio.Seq import Seq
import time
import csv
import pandas as pd
import math
import copy

BND_LEN_THRESHOLD=0.25

ANNOTATE_NEW_TAGS_HEADER="\
##INFO=<ID=RM_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of repeat class covering the SV [SINE,LINE,LTR,DNA,Retroposon or NON-REPETITIVE]\">\n\
##INFO=<ID=RM_ELEMENTS_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the transposable element covered by the SV\">\n\
##INFO=<ID=RM_ELEMENT_PROPORTION,Number=1,Type=String,Description=\"Proportion of the query sequence (includes flanking region) found in the transposable element\">\n\
##INFO=<ID=RM_RECIPROCAL,Number=1,Type=String,Description=\"Type of transposition [Full/Partial/Minimal]\">\n\
##INFO=<ID=RM_SV_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the SV by the transposable element\">\n\
##INFO=<ID=RM_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of the SV covered by transposable elements\">\n\
##INFO=<ID=TRF_CLASSIFICATION,Number=1,Type=String,Description=\"Classification(s) of tandem repeat class covering the SV [HOMO,STR,TR or NON-REPETITIVE]\">\n\
##INFO=<ID=TRF_SV_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the SV by the tandem repeat(s)\">\n\
##INFO=<ID=TRF_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of the SV covered by tandem repeats\">\n\
##INFO=<ID=TRF_PERIOD_SIZE,Number=1,Type=String,Description=\"Period size of the repeat(s)\">\n\
##INFO=<ID=TRF_COPY_NUMBER,Number=1,Type=String,Description=\"Copy number of the repeat(s)\">\n\
##INFO=<ID=CONSENSUS_REPEAT,Number=1,Type=String,Description=\"Motif of repeat(s) found by Tandem Repeat Finder\">\n\
##INFO=<ID=FINAL_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of SV as repetitive element based on TRF and RepeatMasker results\">\n\
##INFO=<ID=DISEASE_GENE,Number=1,Type=String,Description=\"STR disease associated with gene\">\n\
##INFO=<ID=STRCHIVE_MOTIF,Number=1,Type=String,Description=\"Is consensus repeat a version (rotation/complement) of pathogenic motif(s) annotated by STRchive \">\n\
##INFO=<ID=PATHOGENIC_MIN,Number=1,Type=String,Description=\"Minimum pathogenic number\">\n\
"

ANNOTATE_COLS=["CHROM","POS","ID","REF","ALT","RM_CLASSIFICATION","RM_ELEMENTS_COVERAGE","RM_ELEMENT_PROPORTION",\
               "RM_SV_COVERAGE","RM_TOTAL_SV_COVERAGE","RM_RECIPROCAL","TRF_CLASSIFICATION","TRF_PERIOD_SIZE",\
                "TRF_COPY_NUMBER","CONSENSUS_REPEAT","TRF_SV_COVERAGE","TRF_TOTAL_SV_COVERAGE","FINAL_CLASSIFICATION",\
                    "DISEASE_GENE","STRCHIVE_MOTIF","PATHOGENIC_MIN"]

PLOT_COLS=["SVTYPE","SVLEN","CLASSIFICATION","RECIPROCAL"]

def read_sv_info(sv_file):
    print(f"Info: Reading {args.info} file...")
    sv_info = {}
    with open(sv_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for sv in reader:
            # if sv[6] != "BND.1":
            #     continue
            # Extract values from the columns
            chrom = sv[0]
            query_start = int(sv[1])
            query_end = int(sv[2])
            pos = int(sv[3])
            end = int(sv[4])
            sv_len = int(sv[5])
            sv_id = sv[6]
            callerID = sv[7]
            ref_seq = sv[8]
            alt_seq = sv[9]

            if "BND" in sv_id:
                assert sv_len == -1
                rel_start = pos - query_start - (pos - query_start)*BND_LEN_THRESHOLD
                rel_end = pos - query_start + (query_end - pos)*BND_LEN_THRESHOLD
                sv_len = rel_end - rel_start
            else:
                if "DEL" in sv_id:
                    assert sv_len < 0
                    sv_len = abs(sv_len)
                rel_start = pos - query_start
                rel_end = rel_start + sv_len 
            
            # Add the values to the dictionary
            sv_info[sv_id] = {                
                'ID': sv_id, 
                'chrom': chrom,
                'position': pos,
                'sv_len': sv_len, 
                'rel_start' : rel_start, 
                'rel_end' : rel_end,
                'start' : 0, 
                'end' : query_end - query_start,
                'callerID' : callerID,  
                'header' : f'{sv_id}\t({callerID})\t{chrom}:{query_start}-{query_end}\t{chrom}:{pos}-{end}\t{sv_len}\t',
                'REF': ref_seq,
                'ALT': alt_seq,
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

def read_rm(args, sv_info, isVisualise):
    """
    Processes RepeatMasker output (.out), adding relevant entries to sv_info
    
    Input:
    - sv_info (dict)            Dictionary containing the strucutral variants
    - isVisualise (bool):       True (visualisation)/ False (annotation)

    Output: 
    - sv_info (dict):           Updated sv_info with repeatMasker annotations
    """
    def get_RM_classification(family, isVisualise):
        if isVisualise:
            classifications = ['SINE', 'LINE', 'LTR', 'DNA', 'Simple_repeat', 'RNA', 'Retroposon', 'Low_complexity', 'Other']
        else: 
            classifications = ['SINE', 'LINE', 'LTR', 'DNA', 'Retroposon']
        
        for classification in classifications:
            if classification in family:
                return classification

        return None
    
    def positionInRepeatFraction(strand, repeat_begin, repeat_end, repeat_left, intersection):
        """
        Calculates the sv in relation to the repeat
        - Proportion (entire query)
        - fraction (sv only)
        
        """

        repeat_length = repeat_end + repeat_left
        element_proportion = (repeat_end - repeat_begin) / repeat_length
        element_coverage = intersection / repeat_length
        return element_proportion, element_coverage

        # $12           #13         $14 
        # repeat_begin  repeat_end  repeat_left
        # if strand == '+':
            # proportion = (repeat_end - repeat_begin) / ((repeat_end - repeat_begin) + repeat_left)
            # element_coverage = intersection / ((repeat_end - repeat_begin) + repeat_left)
        # elif strand == 'C':
            # proportion = (repeat_end - repeat_left) / ((repeat_begin + repeat_end) - repeat_left)
            # element_coverage = intersection / ((repeat_begin + repeat_end) - repeat_left)

    print(f"Info: Reading {args.rm} file with isVisualise:{isVisualise}...")
    with open(args.rm, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        row_count = 0
        for row in reader:
            row_count += 1
            sv_id = row[4]
            te_start = int(row[5])
            te_end = int(row[6])
            strand = row[8]
            repeat = row[9]
            family = row[10]
            te_id = row[14]

            # print("{}\t{}".format(repeat, family))
            assert te_start < te_end, f"TE start {te_start} should be less than TE end {te_end} for sv_id {sv_id}."

            classification = get_RM_classification(family, isVisualise)

            # Add element if ['SINE', 'LINE', 'LTR', 'DNA', 'Retroposon']
            if classification and sv_id in sv_info:
                
                if strand == '+':
                    # Remove the parentheses around the repeat
                    repeat_begin = int(row[11].strip("()"))
                    repeat_end = int(row[12].strip("()"))
                    repeat_left = int(row[13].strip("()"))
                elif strand == 'C':
                    # Remove the parentheses around the repeat
                    repeat_begin = int(row[13].strip("()"))
                    repeat_end = int(row[12].strip("()"))
                    repeat_left = int(row[11].strip("()"))

                sv = sv_info[sv_id]
                # assert repeat_begin < repeat_end, f"row {row_count}: sv_id {sv_id} Repeat begin {repeat_begin} should be less than repeat end {repeat_end} for strand {strand}."
                if not repeat_begin < repeat_end:
                    intersection = calculate_intersection_length(te_start, te_end, sv['rel_start'], sv['rel_end'])
                    sv_coverage = round(intersection / sv['sv_len'], 6)
                    print(f"{intersection},{sv_coverage} row {row_count}: sv_id {sv_id},{sv_info[sv_id]['callerID']} Repeat begin {repeat_begin} should be less than repeat end {repeat_end} for strand {strand}.")
                    continue

                # print("query's matching len:{}, repeat's matching len:{}".format(te_end-te_start,repeat_end-repeat_begin))
                
                # Add if the transposable element overlaps with the SV
                intersection = calculate_intersection_length(te_start, te_end, sv['rel_start'], sv['rel_end'])
                sv_coverage = round(intersection / sv['sv_len'], 6)
                element_proportion, element_coverage = positionInRepeatFraction(strand, repeat_begin, repeat_end, repeat_left, intersection)
            
                # Add if the intersection is > min_intersect (e.g. 5%)
                if sv_coverage and sv_coverage > args.min_sv_coverage:
                    element = {
                        'te_start': te_start,
                        'te_end': te_end,
                        'family': family,
                        'repeat': repeat,
                        'element_coverage' : element_coverage,
                        'element_proportion' : element_proportion,
                        'classification': classification, 
                        'sv_coverage' : sv_coverage, 
                        'te_id' : te_id, 
                    }

                    # First element for SV - add the 'RM'
                    if 'RM' not in sv_info[sv_id]:
                        sv_info[sv_id]['RM'] = {}
                    # First element of class for SV
                    if classification not in sv_info[sv_id]['RM']:
                        sv_info[sv_id]['RM'][classification] = []
                    
                    sv_info[sv_id]['RM'][classification].append(element)

                    # print(sv_info[sv_id]['RM'])

    return sv_info

def read_trf(args, sv_info):
    """
    Processes Tandem Repeat Finder output (.dat), adding relevant entries to sv_info
    
    Input:
    - sv_info (dict)            Dictionary containing the strucutral variants

    Output: 
    - sv_info (dict):           Updated sv_info with TRF annotations
    """
    def get_TRF_classification(period_size):
        # Classify as HOMO, STR, or TR based on period size 
        if period_size == 1: 
            key = 'HOMO'
        elif period_size >= 2 and period_size <= 12:
            key = 'STR'
        elif period_size > 12:
            key = 'TR'
        return key 
    
    def parse_TRF_lines(args, lines, sv):
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
            sv_coverage = round(intersection/ sv['sv_len'], 6)

            repeat_info = {
                'repeat_start': repeat_start,
                'repeat_end': repeat_end,
                'period_size': period_size,
                'motif': consensus_repeat,

                'copy_number': copy_number,
                'key': get_TRF_classification(period_size),
            }

            if sv_coverage and float(sv_coverage) > args.min_sv_coverage:
                repeat_info.update({'sv_coverage' : sv_coverage})
                tandem_repeats.append(repeat_info)
        
        return tandem_repeats

    trf_lines = None

    # @INS.1
    # 1040 1111 23 3.4 23 67 25 68 9 18 50 22 1.76 GAGGGCGTCTGGTCGTCCTGAGG GAGGGCGTCTGGTCGTCCTGAGGGAGGGCCGGTGTTGGTGAGGGCATCTGGTCGTCCTGAGGGAGGGGGTCT GGTGAGAGACGCTGCCGCAGAGCCGCCCGAGAGGGAGGGTCAGTGTTGGT TCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATTTTTTATTTT
    # 1020 1106 39 2.2 39 93 0 147 11 14 51 21 1.74 GAGGGAGGGCCAGTGTTGGTGAGGGCATCTGGTCGTCCT GAGGGAGGGTCAGTGTTGGTGAGGGCGTCTGGTCGTCCTGAGGGAGGGCCGGTGTTGGTGAGGGCATCTGGTCGTCCTGAGGGAGGG CACAGAGGGAAACAAGGGGAGGTGAGAGACGCTGCCGCAGAGCCGCCCGA GGTCTTCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATTTTTT
    # 1153 1200 5 9.8 5 95 4 89 18 0 0 81 0.70 TTTTA TTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTATTTT AGGGGGTCTTCTTCACATTCTCACCTCATTTCTTTTCACTCAGCAGGATT GAAACGGAGTCTCACTCTTGCCTAGGCTGGAGTGCAATGGCGCAATCTCG
    print(f"Info: Reading {args.trf} file...")
    with open(args.trf, 'r') as file:
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
            trf_result = parse_TRF_lines(args, trf_lines, sv)

            sv_info[sv_id]['TRF'] = trf_result
    
    return sv_info

def load_strchive(str_file):
    """
    Reads tab delimited STRchive bed file, and organises data into dict with the 
    key as chromosome names and values as lists of the STRs. 
    It determines rotations and reverse complement of the pathogenic motifs. 
    Input: 
        - str_file (str): The path to the STRchive file.
    Output: 
        - data_dict (dict): Dictionary containing STRchive entries
    """
    def rotate_reverse(repeat):
        """
        Takes a sequence and finds the rotations and reverse complements of that sequence 
        Input: 
            repeat:             string                  e.g. CAT
        Output: 
            motifs
                repeat_motifs:      list of strings         e.g. repeat ['CAT', 'ATC', 'TCA'] reverse ['ATG', 'GAT', 'TGA']
        """
        motifs = []

        if repeat != None:
            n = len(repeat)
            for i in range(n):
                rotation = repeat[i:]+repeat[:i]
                seq = Seq(rotation)
                reverse_rotation = str(seq.reverse_complement())
                motifs.append(rotation)
                motifs.append(reverse_rotation)
        
        return motifs
    
    data_dict = {}
    required_headers = {'#chrom', 'start', 'stop', 'id', 
                        'pathogenic_motif_reference_orientation', 'pathogenic_min'}

    # Open the file
    with open(str_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')  # Assuming tab-delimited file

        # CHECK: Correct file from STRchive bed file 
        headers = set(reader.fieldnames)
        if not headers.issuperset(required_headers):
            missing_headers = required_headers - headers
            print(f"Could not annotate STRchive as file is incorrectly formatted and/or missing required headers: {', '.join(missing_headers)}.", file=sys.stderr)
            return None

        for row in reader:
            pathogenic_motifs = []
            # Extract relevant columns, including start and stop coordinates
            motifs = row['pathogenic_motif_reference_orientation'].split(',')
            for repeat in motifs:
                pathogenic_motifs.extend(rotate_reverse(repeat))
            
            chrom = row['#chrom']
            entry = {
                'start': int(row['start']),
                'stop': int(row['stop']),
                'id': row['id'],
                'pathogenic_motif_reference_orientation': pathogenic_motifs,
                'pathogenic_min': row['pathogenic_min']  
            }

            # Store the data in the dictionary under the appropriate chromosome
            if chrom not in data_dict:
                data_dict[chrom] = []
            data_dict[chrom].append(entry)
    
    return data_dict

def calculate_divisor(bucket_percentage):
    if bucket_percentage <= 0 or bucket_percentage >= 1:
        # Set to default 
        bucket_percentage = 0.05
    divisor = int(1/bucket_percentage)
    return divisor

def filter_rm(sv_info):
    """
    Determines non-overlapping transposable elements within each element type e.g. SINE
    Input:
        - sv_info (dict): A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
    Output
        - rm_elements (dict): A dictionary containing the SV ID updated with the list of non-overlapping transposable elements 
    
    """

    rm_elements = {}
    for sv_id in sv_info:
        rm_elements[sv_id] = sv_info[sv_id]
        
        if 'RM' in sv_info[sv_id]:
            all_non_overlapping = []  # Flattened list for non-overlapping elements
            class_reciprocal_map = {}
            for classification in sv_info[sv_id]['RM']:
                transposition = "Minimal"

                # 1. Sum the coverage for the element class
                element_list = sv_info[sv_id]['RM'][classification]
                # pprint(element_list)

                # 2. Sort elements in decreasing order of intersect length
                sorted_elements = sorted(element_list, key=lambda x: x.get('sv_coverage', 0), reverse=True)

                # 3. Initialize an empty list to hold non-overlapping TEs
                non_overlapping = []

                # Function to check overlap
                def is_overlapping(te1, te2):
                    assert te1['te_start'] < te1['te_end'], "Invalid interval for te1"
                    assert te2['te_start'] < te2['te_end'], "Invalid interval for te2"
                    return not (te1['te_end'] <= te2['te_start'] or te2['te_end'] <= te1['te_start'])

                # 4. Iterate through sorted TEs and select non-overlapping ones
                for te in sorted_elements:
                    if all(not is_overlapping(te, existing_te) for existing_te in non_overlapping):
                        non_overlapping.append(te)
                
                # Determine if transposition is Full/Partial/Minimal
                total_sv_coverage = 0
                element_coverage_complete = True
                for te in non_overlapping:
                    total_sv_coverage += te['sv_coverage']
                    if te['element_coverage'] < 0.75:
                        element_coverage_complete = False
                
                if total_sv_coverage >= 0.75 and element_coverage_complete:
                    transposition = "Full"
                elif total_sv_coverage >= 0.75 and not element_coverage_complete:
                    transposition = "Partial"

                # Append non-overlapping elements from this element group to the flat list
                all_non_overlapping.extend(non_overlapping)
                class_reciprocal_map[classification] = transposition
            
            # 5. Order by start position 
            all_non_overlapping = sorted(all_non_overlapping, key=lambda x: x['te_start'])

            # Assign non-overlapping list back to the new dictionary
            rm_elements[sv_id].update({'RM': all_non_overlapping})
            rm_elements[sv_id].update({'transposition' : class_reciprocal_map})

    return rm_elements

def filter_trf(sv_info, interval_divisor):
    """
    Determines non-overlapping tandem repeats (based on intersection and period size)
    Input:
        - sv_info (dict): A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
        - internal_divisor (int): Divisor used to group intersection fractions into interval e.g. 20 interval_divisor -> 5% intervals
    Output
        - trf_elements (dict): A dictionary containing the SV ID updated with the list of non-overlapping trfs in 'TRF' and total_fraction
    """
    def sum_fractions(intersections):
        total_fraction = 0
        for entry in intersections: 
            total_fraction += entry['sv_coverage'] 
        return round(total_fraction,2)
    
    for sv_id in sv_info:
        if 'TRF' in sv_info[sv_id]:
            trf_list = sv_info[sv_id]['TRF']
            
            # Custom key function that groups intersections by 0.05 bins (5%)
            def custom_sort_key(trf):
                # Prioritises largest intersection length in intervals
                # Priotisies lowest period size within the intervals
                intersection_bin = math.floor(trf['sv_coverage'] * interval_divisor) / interval_divisor
                return (-intersection_bin, trf['period_size'])
            
            # 2: Sort TRFs by the custom key
            sorted_trfs = sorted(trf_list, key=custom_sort_key)
            
            # 3: Initialize an empty list to hold non-overlapping TRFs
            non_overlapping = []

            # Function to check overlap
            def is_overlapping(trf1, trf2):
                return not (trf1['repeat_end'] <= trf2['repeat_start'] or trf2['repeat_end'] <= trf1['repeat_start'])

            # 4: Iterate through sorted TRFs and select non-overlapping ones
            for trf in sorted_trfs:
                if all(not is_overlapping(trf, existing_trf) for existing_trf in non_overlapping):
                    non_overlapping.append(trf)

            # 5. Sort by start position 
            non_overlapping = sorted(non_overlapping, key=lambda x: x['repeat_start'])

            sv_info[sv_id]['TRF'] = non_overlapping
    
            total_fraction = sum_fractions(non_overlapping)
            sv_info[sv_id]['total_fraction'] = str(total_fraction)

    return sv_info

def create_rm_tsv_record(sv_info, sv_id, tsv_out):
    """
    Create and write a TSV record for RepeatMasker (RM) annotations of a given structural variant (SV).
    - sv_info (dict):        A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
    - sv_id (str):           The ID of the structural variant for which to create the TSV record.
    - tsv_out (file object): A file object to which the TSV record will be written.
    """
    sv_data = sv_info[sv_id]
    chrom = sv_data['chrom']
    pos = sv_data['position']
    sv_len = sv_data['sv_len']
    sv_type = sv_id.split('.')[0]
    if 'RM' in sv_data:
        rm_entries = sv_data['RM']
        for entry in rm_entries:
            sv_coverage = round(entry['sv_coverage'], 2)
            element_coverage = round(entry['element_coverage'], 2)
            element_proportion = round(entry['element_proportion'],2)
            classification = entry['classification']
            repeat_type = entry['repeat']

            tsv_out.write(f'{sv_id}\t{chrom}\t{pos}\t{sv_len}\t{sv_type}\t{sv_coverage}\t{element_coverage}\t{element_proportion}\t{classification}\t{repeat_type}\n')

def create_trf_tsv_record(sv_info, sv_id, tsv_out):
    """
    Create and write a TSV record for TandemRepeatFinder annotations of a given structural variant (SV).
    - sv_info (dict):        A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
    - sv_id (str):           The ID of the structural variant for which to create the TSV record.
    - tsv_out (file object): A file object to which the TSV record will be written.
    """
    sv_data = sv_info[sv_id]
    chrom = sv_data['chrom']
    pos = sv_data['position']
    sv_len = sv_data['sv_len']
    sv_type = sv_id.split('.')[0]
    if 'TRF' in sv_data:
        trf_entries = sv_data['TRF']
        for entry in trf_entries:
            copy_number = entry['copy_number']
            period_size = entry['period_size']
            sv_coverage = round(entry['sv_coverage'], 2)
            classification = entry['key']
            consensus_repeat = entry['motif']
            tsv_out.write(f'{sv_id}\t{chrom}\t{pos}\t{sv_len}\t{sv_type}\t{sv_coverage}\t{copy_number}\t{period_size}\t{classification}\t{consensus_repeat}\n')

def print_results(title, repeat_count, total, output_total):
    """    
    Prints a formatted summary of the number of repeats (HOMO, STR, TR, SINE, LINE, LTR, DNA, Retroposon, Non-repetitive)
    """
    print(f"\n{title}:")
    print("-" * 35)
    for key, value in repeat_count.items():
        percentage = (value / total) * 100
        print(f"{key:<15}: {value:>5}  {percentage:>6.2f}%")
    if output_total:
        print("-" * 20)
        print(f"{'Total':<15}: {total:>5}")

def get_annot_info(args, sv_info, sv_id, strchive):
    """
    Create and write a new VCF record, including annotations from RepeatMasker and TRF data.
    It determines the final classification based on the coverage of the SV by
    repeats (TRF or mobile element) and writes the new VCF record to the provided VCF file.
    Input
    - vcf_file (VCF file object):   The file object where the new VCF record will be written.
    - sv_vcf (VCF file object):     The VCF file containing the original SV data to retrieve relevant records.
    - sv_info (dict):               A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
    - sv_id (str):                  The ID of the SV for which to create a VCF record.

    -> Writes entry to the VCF

    Output
    - tuple:                        A tuple containing the final classification and transposition type of the SV for summary info
    """
    def prepare_rm_data(rm_data):
        """Extract repeat data fields from RepeatMasker records."""
        sv_coverage = []
        element_coverage = []
        classifications = []
        element_proportion = []
        
        for element in rm_data:
            sv_coverage.append(str(round(element['sv_coverage'], 2)))
            element_coverage.append(str(round(element['element_coverage'], 2)))
            classifications.append(element['classification'])
            element_proportion.append(str(round(element['element_proportion'], 2)))

        return {
            'RM_CLASSIFICATION': ','.join(classifications) if classifications else 'NA',
            'RM_ELEMENTS_COVERAGE': ','.join(element_coverage) if element_coverage else 'NA',
            'RM_ELEMENT_PROPORTION': ','.join(element_proportion) if element_proportion else 'NA',
            'RM_SV_COVERAGE': ','.join(sv_coverage) if sv_coverage else 'NA',
        }
        
    def prepare_trf_data(trf_data):
        """Extract repeat data fields from TRF records."""
        sv_coverage = []
        consensus_repeats = []
        classifications = []
        period_sizes = []
        copy_number = []

        for repeat in trf_data:
            sv_coverage.append(str(round(repeat['sv_coverage'], 2)))
            consensus_repeats.append(repeat['motif'])
            classifications.append(repeat['key'])
            period_sizes.append(str(repeat['period_size']))
            copy_number.append(str(repeat['copy_number']))

        return {
            'TRF_CLASSIFICATION': ','.join(classifications) if classifications else 'NA',
            'TRF_PERIOD_SIZE': ','.join(period_sizes) if period_sizes else 'NA',
            'TRF_COPY_NUMBER': ','.join(copy_number) if copy_number else 'NA',
            'CONSENSUS_REPEAT': ','.join(consensus_repeats) if consensus_repeats else 'NA',
            'TRF_SV_COVERAGE': ','.join(sv_coverage) if sv_coverage else 'NA',
        }

    def calculate_total_coverage(sv_start, sv_end, rm_data):
        """
        Determines the total coverage of the SV by tranposable elements
        The function computes the proportion of the SV length that is covered by non-overlapping 
        RepeatMasker entries. Overlapping TEs are handled to ensure no double-counting 
        of coverage.
        Input:
            - sv_start (int):   The start position of the SV.
            - sv_end (int):     The end position of the SV.
            - rm_data (list of dict) list of filtered repeatMasker entries, where each dict contains
                    - 'te_start' (int): Start position of the repeat
                    - 'te_end' (int):   End position of the repeat
        Output
            - str: A string representation of the total coverage as a proportion (rounded to 2 decimal places).
        """
        sv_length = sv_end - sv_start
        total_length = 0
    
        # Initialize the end of the last repeat to track overlaps
        current_end = 0

        for element in rm_data:
            te_start = element['te_start']
            te_end = element['te_end']

            if te_start < sv_start:
                te_start = sv_start
            if te_end > sv_end:
                te_end = sv_end
            
            # If the te starts after the current_end, add the full length
            if te_start > current_end:
                total_length += te_end - te_start
            # If the te overlaps with the previous one, only add the non-overlapping part
            else:
                total_length += te_end - current_end
            
            # Update current_end to the maximum of the current te's end and previous end
            current_end = max(current_end, te_end)

        return str(round(total_length / sv_length, 2))

    def highest_fraction(args, repeats, fractions, reciprocal_map):
        """
        Initialize a dictionary to store the sum of fractions for each repeat type
        Calculates the total sum of fractions for each unique repeat type (e.g. SINE, LINE / STR, TR)
        and identifies the repeat type with the highest accumulated fraction.
        Inputs:
            - repeats (list of str):        List of repeat types (e.g., ["SINE", "LINE", "LTR"]).
            - fractions (list of float):    List of fractions corresponding to the repeat
        Outputs:
            - max_repeat (str):             Repeat type with hightest total 
            - fraction (float):             Corresponding highest total_fraction value
        """
        fraction_sum = {}
        # Sum the fractions for each repeat type
        for repeat, fraction in zip(repeats, fractions):
            if repeat in fraction_sum:
                fraction_sum[repeat] += float(fraction)
            else:
                fraction_sum[repeat] = float(fraction)
        
        # Find the repeat type with the highest total fraction
        max_repeat, max_fraction = max(fraction_sum.items(), key=lambda x: x[1])
        if max_fraction < args.min_class_sv_coverage:
            # print(f"WARNING: {max_repeat} has a total fraction of {max_fraction} which is below the threshold of {args.min_class_sv_coverage}.")
            return 'NON_REPETITIVE'
        else:
            if max_repeat in reciprocal_map:
                transposition = reciprocal_map[max_repeat]
                if transposition in ['Full', 'Partial']:
                    return max_repeat
                else:
                    return 'NON_REPETITIVE'
            else:
                return max_repeat
    
    def get_final_classification(args, trf_class, rm_class, total_trf_coverage, total_rm_coverage, trf_data, rm_data, reciprocal_map):
        """
        Function that determines the final classification of the SV --> Tandem Repeat or Transposable element
        When multiple repeat types are present it compares the coverage and selects the repeat with the highest coverage fraction
        Input:
        - trf_class (str): The classification of the SV based on tandem repeat data. Multiple classes are comma-separated.
        - rm_class (str): The classification of the SV based on transposable element data. Multiple classes are comma-separated.
        - total_trf_coverage (float): The total coverage fraction of tandem repeats for the SV.
        - total_rm_coverage (float): The total coverage fraction of transposable elements for the SV.
        - trf_data (dict): Additional TRF-related data, including TRF coverage fractions for each type.
                                Example: {'TRF_SV_COVERAGE': '0.5,0.3,0.2'}
        - rm_data (dict): Additional RM-related data, including RM coverage fractions for each type.
                            Example: {'RM_SV_COVERAGE': '0.4,0.6'}
        """
        trf_classes = trf_class.split(',')
        rm_classes = rm_class.split(',')
        trf_coverages = trf_data['TRF_SV_COVERAGE'].split(',')
        rm_coverages = rm_data['RM_SV_COVERAGE'].split(',')

        if (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
            # TRF has higher coverage
            if total_trf_coverage > total_rm_coverage:
                classification = highest_fraction(args, trf_classes, trf_coverages, reciprocal_map)
            # Mobile element from RepeatMasker has higher coverage
            else:
                classification = highest_fraction(args, rm_classes, rm_coverages, reciprocal_map)
        # TRF (no repeatMasker)
        elif (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class in ['NA', 'NON_REPETITIVE']):
            classification = highest_fraction(args, trf_classes, trf_coverages, reciprocal_map)
        # RM (no TRF intersect)
        elif (trf_class in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
            classification = highest_fraction(args, rm_classes, rm_coverages, reciprocal_map)
        # None
        else: 
            classification = 'NON_REPETITIVE'
        
        return classification

    def pick_reciprocal(classification, reciprocal_map):
        """
        Function that determines the transposition type based on the classification
        Input:
        - classification (str): The final classification of the SV based on repeat data.
        - reciprocal_map (dict): A dictionary mapping classifications to transposition types.
        """
        if classification in reciprocal_map:
            return reciprocal_map[classification]
        else:
            return 'NA'

    sv_data = sv_info[sv_id]
    chrom = sv_data['chrom']
    pos = sv_data['position']
    end = pos + sv_data['sv_len']
    callerID = sv_data['callerID']

    ####### 1. Get annotations for RepeatMasker and TRF ######

    rm_data = prepare_rm_data(sv_data.get('RM', []))
    trf_data = prepare_trf_data(sv_data.get('TRF', []))
    
    ### Check total coverage of the SV by RepeatMasker and TRF > 50%
    # --> Replaces any <50% as non-repetitive

    # RepeatMasker
    total_rm_coverage = calculate_total_coverage(sv_data['rel_start'], sv_data['rel_end'], sv_data['RM']) if 'RM' in sv_data else 'NA'
    if total_rm_coverage != 'NA' and float(total_rm_coverage) < args.min_total_sv_coverage:
        for key in rm_data:
            rm_data[key] = 'NA'
        rm_data['RM_CLASSIFICATION'] = 'NON_REPETITIVE'
    rm_data['RM_TOTAL_SV_COVERAGE'] = total_rm_coverage

    # TRF
    total_trf_coverage = sv_data.get('total_fraction', 'NA')
    if total_trf_coverage != 'NA' and float(total_trf_coverage) < args.min_total_sv_coverage:
        for key in trf_data:
            trf_data[key] = 'NA'
        trf_data['TRF_CLASSIFICATION'] = 'NON_REPETITIVE'
    trf_data['TRF_TOTAL_SV_COVERAGE'] = total_trf_coverage

    ####### 2. Determine Final Classification ######
    # Both TRF and RM
    trf_class = trf_data['TRF_CLASSIFICATION']
    rm_class = rm_data['RM_CLASSIFICATION']
    reciprocal_map = sv_data['transposition'] if 'transposition' in sv_data else 'NA'
    final_classification = get_final_classification(args, trf_class, rm_class, total_trf_coverage, total_rm_coverage, trf_data, rm_data, reciprocal_map)
    transposition = pick_reciprocal(final_classification, reciprocal_map)

    # Unpacks annotations
    info = {
        'CHROM' : chrom,
        'POS' : pos,
        'ID' : callerID,
        'RM_RECIPROCAL' : transposition,
        'FINAL_CLASSIFICATION' : final_classification,
        'DISEASE_GENE' : '',
        'STRCHIVE_MOTIF' : '',
        'PATHOGENIC_MIN' : '',
        **rm_data,  # Unpack repeat data from RepeatMasker annotations
        **trf_data, # Unpack repeat data from TRF annotations
        'REF' : sv_data['REF'],
        'ALT' : sv_data['ALT'],
    }

    ####### Check if it intersects any known STR sites #####
    if strchive and chrom in strchive:
        for entry in strchive[chrom]:
            if end >= entry['start'] and pos <= entry['stop']:
                info.update({'DISEASE_GENE' : entry['id']})
                pathogenic_motifs = entry['pathogenic_motif_reference_orientation']
                consensus_repeat = trf_data['CONSENSUS_REPEAT'].split(',')
                strchive_motif_list = []
                for repeat in consensus_repeat:
                    if repeat in pathogenic_motifs:
                        strchive_motif_list.append('YES')
                    else: 
                        strchive_motif_list.append('NO')
                strchive_motif = ','.join(strchive_motif_list)
                info.update({'STRCHIVE_MOTIF' : strchive_motif })
                info.update({'PATHOGENIC_MIN' : entry['pathogenic_min']})
                break

    return info
    ####### 3. Obtain info from original VCF, and add annotations ######
    # Get record from SV VCF 
    # records = sv_vcf.fetch(chrom, pos-1, pos+1)
    # for record in records:
    #     # Check SnifflesID matches
    #     if callerID == record.id:
    #         record.id=sv_id
    #         # Update info fields with new annotations
    #         for annotation in info:
    #             record.info[annotation] = info[annotation]
    #         # Write the updated record to the VCF file
    #         # vcf_file.write(record)
    #         break

    # return info['FINAL_CLASSIFICATION'], info['RM_TRANSPOSITION']

def output_annotations(args, strchive, sv_info):
    """
    Process structural variants and write to VCF.
    This function reads structural variant information from a provided SV VCF file, processes each variant, 
    and writes the annotated results to new VCF and TSV files. It classifies repeats and structural variants 
    based on repeat annotations and covers the following processes:
    1. Creation of a VCF record with repeat classification.
    2. Writing repeat classification data to TSV files for RepeatMasker (RM) and TRF data.
    3. Generates summary statistics of repeat types and SV classifications.
    Input
    - sv_info (dict): A dictionary containing the structural variant information, with SV IDs as keys.
    - sv_vcf_file (str): Path to the input SV VCF file containing structural variants.
    - vcf_output (str): Path to the output VCF file where the processed results will be written.
    - tsv_rm_out (str): Path to the output TSV file for RepeatMasker annotations.
    - tsv_trf_out (str): Path to the output TSV file for TRF annotations.
    - min_repetitive (float): Minimum threshold for repeat coverage to consider a repeat as repetitive (e.g., 0.5 for 50%)
    - strchive (dict): STR disease/gene with pathogenic motif and number 
    """
    # Initialise count for breakdown of repeat types
    sv_count = 0
    repeat_count = {'HOMO' : 0, 'STR' : 0, 'TR' : 0, 'Full' : 0, 'Partial' : 0, 'NON_REPETITIVE' : 0}
    count_by_sv = {
        'INS': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : [0,0], 'LINE' : [0,0], 'LTR' : [0,0], 'DNA' : [0,0], 'Retroposon' : [0,0], 'NON_REPETITIVE' : 0},
        'DEL': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : [0,0], 'LINE' : [0,0], 'LTR' : [0,0], 'DNA' : [0,0], 'Retroposon' : [0,0], 'NON_REPETITIVE' : 0},
        'DUP': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : [0,0], 'LINE' : [0,0], 'LTR' : [0,0], 'DNA' : [0,0], 'Retroposon' : [0,0], 'NON_REPETITIVE' : 0},
        'INV': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : [0,0], 'LINE' : [0,0], 'LTR' : [0,0], 'DNA' : [0,0], 'Retroposon' : [0,0], 'NON_REPETITIVE' : 0},
        'BND': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : [0,0], 'LINE' : [0,0], 'LTR' : [0,0], 'DNA' : [0,0], 'Retroposon' : [0,0], 'NON_REPETITIVE' : 0},
    }

    annotate_tsv_out = f"{args.out}/vcf_annotate.tsv"
    plot_annotate_out = f"{args.out}/plot_annotate.tsv"
    with open(annotate_tsv_out, "w") as f, open(plot_annotate_out, "w") as f2:
        # print("\t".join(ANNOTATE_COLS))
        f.write("#")
        f.write("\t".join(ANNOTATE_COLS))
        f.write("\n")

        f2.write("\t".join(PLOT_COLS))
        f2.write("\n")

        # Process each structural variant and write to VCF
        transposition_map = {"Full": 0, "Partial":1}
        for sv_id in sv_info:
            sv_type = sv_id.split('.')[0]
            sv_len = sv_info[sv_id]['sv_len']
            sv_count += 1
            info = get_annot_info(args, sv_info, sv_id, strchive)
            # print("\t".join(str(info[key]) for key in ANNOTATE_COLS))
            f.write("\t".join(str(info[key]) for key in ANNOTATE_COLS))
            f.write("\n")

            classification, transposition = info['FINAL_CLASSIFICATION'], info['RM_RECIPROCAL']
            if classification in repeat_count:
                repeat_count[classification] += 1
                count_by_sv[sv_type][classification] += 1
            # A repeat element 
            else: 
                repeat_count[transposition] += 1
                count_by_sv[sv_type][classification][transposition_map[transposition]] += 1
            
            f2.write(f"{sv_type}\t{sv_len}\t{classification}\t{transposition}\n")

    
    header_out = f"{args.out}/vcf_header.txt"
    with open(header_out, "w") as f:
            # print(ANNOTATE_NEW_TAGS_HEADER)
            f.write(f"{ANNOTATE_NEW_TAGS_HEADER}")

    tsv_rm_out = f"{args.out}/rm_annotate.tsv"
    tsv_trf_out = f"{args.out}/trf_annotate.tsv"
    # Write RepeatMasker and TRF records to separate TSV files 
    with open(tsv_rm_out, 'w') as tsv_rm_out, open(tsv_trf_out, 'w') as tsv_trf_out:
        # Write headers for both files
        tsv_rm_out.write('ID\tCHR\tPOS\tSVLEN\tSV_TYPE\tSV_COVERAGE\tELEMENT_COVERAGE\tELEMENT_PROPORTION\tFAMILY\tREPEAT\n')
        tsv_trf_out.write('ID\tCHR\tPOS\tSVLEN\tSV_TYPE\tSV_COVERAGE\tCOPY_NUMBER\tPERIOD_SIZE\tCLASSIFICATION\tCONSENSUS_REPEAT\n')

        # Iterate once through sv_info and write to both files
        for sv_id in sv_info:
            create_rm_tsv_record(sv_info, sv_id, tsv_rm_out)  # Write to rm file
            create_trf_tsv_record(sv_info, sv_id, tsv_trf_out)  # Write to trf file

    # Print results summary
    print_results("Classification summary", repeat_count, sv_count, True)
    # Print results summary by SV type 
    sv_type_df = pd.DataFrame(count_by_sv).T
        #     HOMO STR TR ...
        # INS
        # DEL
    print(sv_type_df.to_string())

def create_SV(args, sv_info):
    """
    Creates and writes structural variant (SV) and repeat diagrams to an output file.

    The function processes information from structural variants, RepeatMasker (RM), 
    and Tandem Repeat Finder (TRF) data, generating a graphical representation of the 
    SV and its associated repeats, both in the flanking region and within the SV region itself.

    Input
    - sv_info (dict):       A dictionary containing the structural variants
    - diagram_length (int): The length of the diagram to generate for the SV and repeats.
    - output_file (str):    The path to the output file where the diagrams will be written.
    """
    def create_trf(repeat, diagram_length, scale, start):
        """
        Creates the diagram for a TRF entry
        """
        repeat_start = repeat['repeat_start']
        repeat_end = repeat['repeat_end']
        
        repeat_start_scale = int((repeat_start - start) * scale)
        repeat_end_scale = int((repeat_end - start) * scale)

        repeat_start_scaled = max(0, repeat_start_scale)
        repeat_end_scaled = min(diagram_length - 1, repeat_end_scale)

        trf_diagram = [' '] * diagram_length
        period_size = repeat['period_size']

        for i in range(repeat_start_scaled, repeat_end_scaled + 1):
            trf_diagram[i] = '.'
        # Set boundary if within the scale
        if 0 <= repeat_start_scale < diagram_length:
            trf_diagram[repeat_start_scaled] = '['  
        if 0 <= repeat_end_scale < diagram_length:
            trf_diagram[repeat_end_scaled] = ']' 
        trf_diagram = ''.join(trf_diagram)
        trf_diagram = format_trf_info(f' {period_size}', trf_diagram, diagram_length)

        return trf_diagram

    def create_rm(rm_diagram, element, diagram_length, scale, start):
        """
        Creates the diagram for a RepeatMasker entry
        """
        element_start = element['te_start']
        element_end = element['te_end']
        
        element_start_scale = int((element_start - start) * scale)
        element_end_scale = int((element_end - start) * scale)

        element_start_scaled = max(0, element_start_scale)
        element_end_scaled = min(diagram_length - 1, element_end_scale)

        for i in range(element_start_scaled, element_end_scaled + 1):
            rm_diagram[i] = '.'
        # Set boundary if within the scale
        if 0 <= element_start_scale < diagram_length:
            if rm_diagram[element_start_scaled] != '.':
                rm_diagram[element_start_scaled] = '|'  
            else:
                rm_diagram[element_start_scaled] = '['
        if 0 <= element_end_scale < diagram_length:
            if rm_diagram[element_end_scaled] != '.':
                rm_diagram[element_end_scaled] = '|'
            else:
                rm_diagram[element_end_scaled] = ']'
        
        return rm_diagram
    
    def get_repeat_info(sv_data, key):
        if key in sv_data:
            return sv_data[key]
        return []
    
    def create_diagram(diagram_length, start, end):
        """
        Creates the diagram for SV and SV with flanking 
        """
        if start < 0 or start >= diagram_length or end < 0 or end >= diagram_length:
            return None
        
        diagram = ['-'] * diagram_length
        diagram[start] = '['
        diagram[end] = ']'
        for i in range(start+1, end):
            diagram[i] = '#'
        diagram = ''.join(diagram)

        return diagram

    def format_trf_info(info, diagram, diagram_length, field_width=20):
        info_str = f" {info}".ljust(field_width)  # Left-align the info within a fixed-width field
        diagram_str = f"{diagram}".ljust(diagram_length + 1)  # Ensure diagram is aligned by adding a marker (*)
        return f"{info_str}{diagram_str}"

    def add_rm_diagrams(rm, rm_output_flanking, rm_output, diagram_length, flanking_scale, sv_scale, start, sv_start):
        for classification in rm:
            elements = rm[classification]
            # Initialize diagrams for both flanking and SV regions
            rm_diagram_flanking = [' '] * diagram_length
            rm_diagram_sv = [' '] * diagram_length
            
            intersect_percentages = []
            repeat_names = set()

            for element in elements:
                # Process flanking region
                rm_diagram_flanking = create_rm(rm_diagram_flanking, element, diagram_length, flanking_scale, start)
                intersection = round(element['sv_coverage'] * 100, 2)
                intersect_percentages.append(f'{intersection}%')

                # Process SV region
                rm_diagram_sv = create_rm(rm_diagram_sv, element, diagram_length, sv_scale, sv_start)
                repeat_names.add(element['repeat'])

            # Convert diagrams to string format
            formatted_info = f' {classification}'
            rm_diagram_flanking = format_trf_info(formatted_info, ''.join(rm_diagram_flanking), diagram_length)
            rm_diagram_sv = format_trf_info(formatted_info, ''.join(rm_diagram_sv), diagram_length)

            # Append results to respective outputs
            rm_output_flanking.append(f'{rm_diagram_flanking}\t{",".join(intersect_percentages)}\n')
            rm_output.append(f'{rm_diagram_sv}\t{",".join(repeat_names)}\n')

        return rm_output_flanking, rm_output

    def add_trf_diagrams(trf, trf_output_flanking, trf_output, diagram_length, flanking_scale, sv_scale, start, sv_start, motif_length=80):
        for repeat in trf:
            # Process flanking region
            trf_diagram_flanking = create_trf(repeat, diagram_length, flanking_scale, start)
            intersect = round(repeat['sv_coverage'] * 100, 2)
            trf_output_flanking.append(f'{trf_diagram_flanking}\t{intersect}%\n')

            # Process SV region
            motif = repeat['motif']
            if len(motif) > motif_length:
                motif = f'*{motif_length}plus'

            trf_diagram_sv = create_trf(repeat, diagram_length, sv_scale, sv_start)
            trf_output.append(f'{trf_diagram_sv}\t{motif}\n')

        return trf_output_flanking, trf_output

    def write_flanking(f, flanking_diagram, rm_output_flanking, trf_output_flanking):
        """
        Writes the SV and flanking diagram with RepeatMasker
        and TRF entries 
        """
        if rm_output_flanking != [] or trf_output_flanking != []:
            f.write(f"{flanking_diagram}\n")
            if rm_output_flanking:
                f.write(' -RM-\n')
                f.writelines(rm_output_flanking)
            if trf_output_flanking:
                f.write(' -TRF-\n')
                f.writelines(trf_output_flanking)
            f.write('\n')
    
    def write_SV(f, sv_diagram, rm_output, trf_output):
        """
        Writes the SV diagram with RepeatMasker (RM) 
        and TRF entries 
        """
        if rm_output != [] or trf_output != []:

            f.write(f"{sv_diagram}\n")
            if rm_output:
                f.write(' -RM-\n')
                f.writelines(rm_output)
            if trf_output:
                f.write(' -TRF-\n')
                f.writelines(trf_output)
            f.write('\n')

    output_file = f"{args.out}/diagram.txt"
    diagram_length = args.len
    with open(output_file, 'w') as f:  # Open the file for writing
        for sv_id in sv_info:
            sv_data = sv_info[sv_id]
            
            id_str = sv_data['header']
            rm = get_repeat_info(sv_data, 'RM')
            trf = get_repeat_info(sv_data, 'TRF')

            # Check if there is RepeatMasker and TRF
            if rm == [] and trf == []:
                continue

            rm_output_flanking = []
            trf_output_flanking = []
            rm_output = []
            trf_output = []

            start, end = sv_data['start'], sv_data['end']
            sv_start, sv_end, sv_len = sv_data['rel_start'], sv_data['rel_end'], sv_data['sv_len']

            # Calcualte scaling factors for diagrams
            flanking_scale = diagram_length / end
            sv_scale = diagram_length / sv_len
            
            sv_start_scaled = int((sv_start - start) * flanking_scale)
            sv_end_scaled = int((sv_end - start) * flanking_scale)
            
            # Generate the SV w/ flanking diagram 
            # Don't generate the flanking diagram if sv ends after the flanking
            if sv_end > end and sv_end_scaled > 100:
                id_str = id_str.replace('\t', ' ')
                print(f"Unable to generate sv + flanking diagram for {id_str}")
                flanking_diagram = None
            else:
                flanking_diagram = create_diagram(diagram_length, sv_start_scaled, sv_end_scaled)
                flanking_diagram = format_trf_info('SV & flanking', flanking_diagram, diagram_length)
            
            # Generate the SV Diagram 
            sv_diagram = create_diagram(diagram_length, 0, diagram_length-1)
            sv_diagram = format_trf_info('SV Diagram', sv_diagram, diagram_length)

            # Create the diagrams for the RepeatMasker and TRF entries
            rm_output_flanking, rm_output = add_rm_diagrams(rm, rm_output_flanking, rm_output, diagram_length, flanking_scale, sv_scale, start, sv_start)
            trf_output_flanking, trf_output = add_trf_diagrams(trf, trf_output_flanking, trf_output, diagram_length, flanking_scale, sv_scale, start, sv_start, motif_length=80)

            # Output the diagrams
            if rm_output_flanking != [] or trf_output_flanking != [] or rm_output != [] or trf_output != []:
                f.write(f"{id_str}\n")
                if flanking_diagram:
                    write_flanking(f, flanking_diagram, rm_output_flanking, trf_output_flanking)
                
                write_SV(f, sv_diagram, rm_output, trf_output)
                f.write(f"\n")

def argparser():
    def positive_int(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is an invalid. It must be greater than zero.")
        return ivalue
    
    def positive_float(value):
        ivalue = float(value)
        if ivalue <= 0 or ivalue > 1:
            raise argparse.ArgumentTypeError(f"{value} is an invalid. It must be a value between 0 and 1")
        return ivalue

    parser = argparse.ArgumentParser(
        description="Process VCF file and create alternative sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('-v', '--vcf', required=True, type=str, help="Path to the input VCF file (compressed or uncompressed)")
    required_args.add_argument('--rm', required=True, type=str, help="Path to RepeatMasker .out file")
    required_args.add_argument('--trf', required=True, type=str, help="Path to TRF .dat file")
    required_args.add_argument('-i', '--info', required=True, type=str, help="Path to SV info file")
    required_args.add_argument('--str', required=True, type=str, help="Path to the strchive bed file")
    required_args.add_argument('-o', '--out', required=True, type=str, help="Path to the output directory")
    
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--min_sv_coverage', required=False, type=positive_float, default=0.05, help="The minimum intersection between a repeat element and SV (aka sv_coverage) e.g. 0.05 (5%) (0 < min_sv_coverage < 1)")
    optional_args.add_argument('--min_class_sv_coverage', required=False, type=positive_float, default=0.25, help="The minimum class sv coverage by repeat elements to be considered repetitive")
    optional_args.add_argument('--min_total_sv_coverage', required=False, type=positive_float, default=0.75, help="The minimum total sv coverage by repeat elements to be considered repetitive")
    optional_args.add_argument('--div', required=False, type=positive_float, default=0.05, help="The chosen intervals to prioritise period size over intersection (0 < divisor < 1)")
    optional_args.add_argument('-l', '--len', required=False, type=positive_int, default=100, help="Diagram length")
    optional_args.add_argument('--debug', required=False, action='store_true', help="Debug mode")
    optional_args.add_argument('-h', '--help', action='help', help="Show this help message and exit")

    return parser

if __name__ == "__main__":
    start = time.time()
    
    parser = argparser()
    args = parser.parse_args()
    
    print(f"Info: VCF File: {args.vcf}")
    print(f"Info: RepeatMasker: {args.rm}")
    print(f"info: TandemRepeatFinder: {args.trf}")
    print(f"info: SV info file: {args.info}")
    print(f"Info: strachive bed: {args.str}")
    print(f"Info: min_sv_coverage: {args.min_sv_coverage}")
    print(f"Info: min_total_sv_coverage: {args.min_total_sv_coverage}")
    print(f"Info: min_class_sv_coverage: {args.min_class_sv_coverage}")
    print(f"Info: divisor: {args.div}")
    print(f"info: Output Directory: {args.out}")
    print(f"info: Output vcf header file: {args.out}/vcf_header.txt")
    print(f"info: Output vcf tsv file: {args.out}/vcf_annotate.tsv")
    print(f"info: Output plot tsv file: {args.out}/plot_annotate.tsv")
    print(f"info: Output RM annoatations file: {args.out}/rm_annotate.tsv")
    print(f"info: Output TRF annotations file: {args.out}/trf_annotate.tsv")
    print(f"info: SV diagram file: {args.out}/diagram.txt")

    if args.debug:
        print(f"Info: Debug mode: {args.debug}")

    if not os.path.exists(args.out):
        os.mkdir(args.out)
    else:
        print("Error: {} output dir already exists.".format(args.out))
        exit(1)

    # Read in data
    sv_info = read_sv_info(args.info)
    sv_info = read_trf(args, sv_info)
    sv_info_vis = copy.deepcopy(sv_info)
    sv_info_vis = read_rm(args, sv_info_vis, True)
    sv_info = read_rm(args, sv_info, False)

    strchive = None
    if args.str:
        strchive = load_strchive(args.str)

    # Filter overlaps
    divisor = calculate_divisor(args.div)
    filtered_rm_info = filter_rm(sv_info)
    filtered_sv_info = filter_trf(filtered_rm_info, divisor)

    # Write to files 
    output_annotations(args, strchive, filtered_sv_info)

    # Diagrams
    create_SV(args, sv_info_vis)

    end = time.time()
    print(f"Run time: {end - start:.3f} seconds")
