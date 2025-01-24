
import sys
import argparse
import pandas as pd
import pysam
import math
import os
import csv
from helpers.classifierHelpers import read_sv_info, read_trf, read_rm
from Bio.Seq import Seq

################################################################################
#                         ERROR CHECKS AND PREPROCESSING                       #
################################################################################
def check_files(rm_file, trf_file):
    # Check file is not empty
    if not os.path.getsize(rm_file) > 0:
        sys.stderr.write(f"Error (args.rm): {rm_file} is empty \n")
        sys.exit(1)
    if not os.path.getsize(trf_file) > 0:
        sys.stderr.write(f"Error (args.trf): {trf_file} is empty \n")
        sys.exit(1)


def calculate_divisor(bucket_percentage):
    if bucket_percentage <= 0 or bucket_percentage >= 1:
        # Set to default 
        bucket_percentage = 0.05
    divisor = int(1/bucket_percentage)

    return divisor

#------------------------------- STRchive -------------------------------------#
def rotate_reverse(repeat):
    """
    Takes a sequence and finds the rotations and reverse complements of that sequence 
    Input: 
        repeat:             string                  e.g. CAT
    Output: 
        motifs
            repeat_motifs:      list of strings         e.g. ['CAT', 'ATC', 'TCA']
            reverse_motifs:     list of strings         e.g. ['ATG', 'GAT', 'TGA']
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



################################################################################
#                                  FILTERING                                   #
################################################################################

#-----------------------------    RepeatMasker    -----------------------------#

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
        
        transposition = "FRAGMENT"
        if 'RM' in sv_info[sv_id]:
            all_non_overlapping = []  # Flattened list for non-overlapping elements
            for element in sv_info[sv_id]['RM']:

                # 1. Sum the coverage for the element class
                element_list = sv_info[sv_id]['RM'][element]
                # pprint(element_list)

                # 2. Sort elements in decreasing order of intersect length
                sorted_elements = sorted(element_list, key=lambda x: x.get('intersection', 0), reverse=True)

                # 3. Initialize an empty list to hold non-overlapping TEs
                non_overlapping = []

                # Function to check overlap
                def is_overlapping(te1, te2):
                    return not (te1['te_end'] <= te2['te_start'] or te2['te_end'] <= te1['te_start'])

                # 4. Iterate through sorted TEs and select non-overlapping ones
                transposition_fraction = 0
                element_coverage_complete = True
                for te in sorted_elements:
                    if all(not is_overlapping(te, existing_te) for existing_te in non_overlapping):
                        non_overlapping.append(te)
                        # Determine if transposition is complete or not complete
                        transposition_fraction += te['intersection']
                        if te['element_coverage'] < 0.75:
                            element_coverage_complete = False

                if transposition_fraction > 0.75 and element_coverage_complete:
                    transposition = "COMPLETE"

                # Append non-overlapping elements from this element group to the flat list
                all_non_overlapping.extend(non_overlapping)
            
            # 5. Order by start position 
            all_non_overlapping = sorted(all_non_overlapping, key=lambda x: x['te_start'])

            # Assign non-overlapping list back to the new dictionary
            rm_elements[sv_id].update({'RM': all_non_overlapping})
            rm_elements[sv_id].update({'transposition' : transposition})

    return rm_elements

#---------------------------  Tandem Repeat Finder ----------------------------#

def sum_fractions(intersections):
    total_fraction = 0
    for entry in intersections: 
        total_fraction += entry['intersection'] 
    return round(total_fraction,2)

def filter_tr(sv_info, interval_divisor):
    """
    Determines non-overlapping tandem repeats (based on intersection and period size)
    Input:
        - sv_info (dict): A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
        - internal_divisor (int): Divisor used to group intersection fractions into interval e.g. 20 interval_divisor -> 5% intervals
    Output
        - trf_elements (dict): A dictionary containing the SV ID updated with the list of non-overlapping trfs in 'TRF' and total_fraction
    """
    for sv_id in sv_info:
        if 'TRF' in sv_info[sv_id]:
            trf_list = sv_info[sv_id]['TRF']
            
            # Custom key function that groups intersections by 0.05 bins (5%)
            def custom_sort_key(trf):
                # Prioritises largest intersection length in intervals
                # Priotisies lowest period size within the intervals
                intersection_bin = math.floor(trf['intersection'] * interval_divisor) / interval_divisor
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


################################################################################
#                               OUTPUT ANNOTATIONS                             #
################################################################################
#------------------------- Coverage Calculations ------------------------------#
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


def highest_fraction(repeats, fractions):
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
    max_repeat = max(fraction_sum, key=fraction_sum.get)

    return max_repeat


#-------------------------------- Extract Info ---------------------------------
def prepare_rm_data(rm_data):
    """Extract repeat data fields from TRF records."""
    sv_coverage = []
    element_coverage = []
    classifications = []
    parts = []
    
    for repeat in rm_data:
        sv_coverage.append(str(round(repeat['intersection'], 2)))
        element_coverage.append(str(round(repeat['element_coverage'], 2)))
        classifications.append(repeat['class'])
        parts.append(str(round(repeat['element_proportion'], 2)))

    return {
        'RM_CLASSIFICATION': ','.join(classifications) if classifications else 'NA',
        'RM_ELEMENTS_COVERAGE': ','.join(element_coverage) if element_coverage else 'NA',
        'RM_ELEMENT_PROPORTION': ','.join(parts) if parts else 'NA',
        'RM_SV_COVERAGE': ','.join(sv_coverage) if sv_coverage else 'NA',
    }

def prepare_trf_data(trf_data):
    """Extract repeat data fields from RepeatMasker records."""
    sv_coverage = []
    consensus_repeats = []
    classifications = []
    period_sizes = []
    copy_number = []

    for repeat in trf_data:
        sv_coverage.append(str(round(repeat['intersection'], 2)))
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

################################# Output Files #################################

#------------------------------------ VCF -------------------------------------#
def create_vcf_header(sv_vcf):
    """
    Function that returns the VCF header with relevant annotations
    """
    header = sv_vcf.header
    
    header.add_line("##source=Classifier_v1.0")
    header.add_line(f"##command={' '.join(sys.argv)}")
    header.add_line("##INFO=<ID=CALLER_ID,Number=1,Type=String,Description=\"Caller ID for the SV\">")
    ## Header Lines for RepeatMasker Annotations
    header.add_line("##INFO=<ID=RM_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of repeat class covering the SV [SINE,LINE,LTR,DNA,Retroposon or NON-REPETITIVE]\">")
    header.add_line("##INFO=<ID=RM_ELEMENTS_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the transposable element covered by the SV\">")
    header.add_line("##INFO=<ID=RM_ELEMENT_PROPORTION,Number=1,Type=String,Description=\"Proportion of the query sequence (includes flanking region) found in the transposable element\">")
    header.add_line("##INFO=<ID=RM_TRANSPOSITION,Number=1,Type=String,Description=\"Type of transposition [COMPLETE/FRAGMENT]\">")
    header.add_line("##INFO=<ID=RM_SV_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the SV by the transposable element\">")
    header.add_line("##INFO=<ID=RM_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of the SV covered by transposable elements\">")

    ## Header Lines for Tandem Repeat Finder Annotations
    header.add_line("##INFO=<ID=TRF_CLASSIFICATION,Number=1,Type=String,Description=\"Classification(s) of tandem repeat class covering the SV \[HOMO,STR,TR or NON-REPETITIVE\]\">")
    header.add_line("##INFO=<ID=TRF_SV_COVERAGE,Number=1,Type=String,Description=\"Coverage(s) of the SV by the tandem repeat(s)\">")
    header.add_line("##INFO=<ID=TRF_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of the SV covered by tandem repeats\">")
    header.add_line("##INFO=<ID=TRF_PERIOD_SIZE,Number=1,Type=String,Description=\"Period size of the repeat(s)\">")
    header.add_line("##INFO=<ID=TRF_COPY_NUMBER,Number=1,Type=String,Description=\"Copy number of the repeat(s)\">")
    header.add_line("##INFO=<ID=CONSENSUS_REPEAT,Number=1,Type=String,Description=\"Motif of repeat(s) found by Tandem Repeat Finder\">")
    
    header.add_line("##INFO=<ID=FINAL_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of SV as repetitive element based on TRF and RepeatMasker results\">")
    
    ## Header Lines for STRchive Annotations
    header.add_line("##INFO=<ID=DISEASE_GENE,Number=1,Type=String,Description=\"STR disease associated with gene\">")
    header.add_line("##INFO=<ID=STRCHIVE_MOTIF,Number=1,Type=String,Description=\"Is consensus repeat a version (rotation/complement) of pathogenic motif(s) annotated by STRchive \">")
    header.add_line("##INFO=<ID=PATHOGENIC_MIN,Number=1,Type=String,Description=\"Minimum pathogenic number\">")

    return header


def get_final_classification(trf_class, rm_class, total_trf_coverage, total_rm_coverage, trf_data, rm_data):
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
    if (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
        # TRF has higher coverage
        if total_trf_coverage > total_rm_coverage:
            trf_classes = trf_class.split(',')
            # Determine main classification if SV is made up of different repeat types
            if len(trf_classes) > 1 and not all(x == trf_classes[0] for x in trf_classes):
                coverages = trf_data['TRF_SV_COVERAGE'].split(',')
                classification = highest_fraction(trf_classes, coverages)
            else:
                classification = trf_classes[0]
        # Mobile element from RepeatMasker has higher coverage
        else:
            rm_classes = rm_class.split(',')
            if len(rm_classes) > 1 and not all(x == rm_classes[0] for x in rm_classes):
                coverages = rm_data['RM_SV_COVERAGE'].split(',')
                classification = highest_fraction(rm_classes, coverages)
            else: 
                classification = rm_classes[0]
    # TRF (no repeatMasker)
    elif (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class in ['NA', 'NON_REPETITIVE']):
        trf_classes = trf_class.split(',')
        # Determine main classification if SV is made up of different repeat types
        if len(trf_classes) > 1 and not all(x == trf_classes[0] for x in trf_classes):
            coverages = trf_data['TRF_SV_COVERAGE'].split(',')
            classification = highest_fraction(trf_classes, coverages)
        else:
            classification = trf_classes[0]
    # RM (no TRF intersect)
    elif (trf_class in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
        rm_classes = rm_class.split(',')
        if len(rm_classes) > 1 and not all(x == rm_classes[0] for x in rm_classes):
            coverages = rm_data['RM_SV_COVERAGE'].split(',')
            classification = highest_fraction(rm_classes, coverages)
        else: 
            classification = rm_classes[0]
    # None
    else: 
        classification = 'NON_REPETITIVE'
    
    return classification

def create_vcf_record(vcf_file, sv_vcf, sv_info, sv_id, min_repetitive, strchive):
    """
    Create and write a new VCF record, including annotations from RepeatMasker and TRF data.
    It determines the final classification based on the coverage of the SV by
    repeats (TRF or mobile element) and writes the new VCF record to the provided VCF file.
    Input
    - vcf_file (VCF file object):   The file object where the new VCF record will be written.
    - sv_vcf (VCF file object):     The VCF file containing the original SV data to retrieve relevant records.
    - sv_info (dict):               A dictionary containing the SVs with TRF and RM data, keyed by SV ID.
    - sv_id (str):                  The ID of the SV for which to create a VCF record.
    - min_repetitive (float):       The minimum coverage (as a decimal) for a repeat to be considered repetitive (e.g., 0.5 for 50%).

    -> Writes entry to the VCF

    Output
    - tuple:                        A tuple containing the final classification and transposition type of the SV for summary info
    """
    
    sv_data = sv_info[sv_id]
    chrom = sv_data['chrom']
    pos = sv_data['position']
    end = pos + sv_data['length']
    callerID = sv_data['callerID']

    ####### 1. Get annotations for RepeatMasker and TRF ######
    transposition = sv_data['transposition'] if 'transposition' in sv_data else 'NA'

    rm_data = prepare_rm_data(sv_data.get('RM', []))
    trf_data = prepare_trf_data(sv_data.get('TRF', []))
    
    ### Check total coverage of the SV by RepeatMasker and TRF > 50%
    # --> Replaces any <50% as non-repetitive

    # RepeatMasker
    total_rm_coverage = calculate_total_coverage(sv_data['rel_start'], sv_data['rel_end'], sv_data['RM']) if 'RM' in sv_data else 'NA'
    rm_data.update({'RM_TOTAL_SV_COVERAGE': total_rm_coverage})    
    if total_rm_coverage != 'NA' and float(total_rm_coverage) < min_repetitive:
        rm_data.update({key: 'NA' for key in rm_data if key != 'RM_TOTAL_SV_COVERAGE'})
        rm_data['RM_CLASSIFICATION'] = 'NON_REPETITIVE'
   
    # TRF
    total_trf_coverage = sv_data.get('total_fraction', 'NA')
    trf_data.update({'TRF_TOTAL_SV_COVERAGE': total_trf_coverage})
    if total_trf_coverage != 'NA' and float(total_trf_coverage) < min_repetitive:
        trf_data.update({key: 'NA' for key in trf_data if key != 'TRF_TOTAL_SV_COVERAGE'})
        trf_data['TRF_CLASSIFICATION'] = 'NON_REPETITIVE'

    # Unpacks annotations
    info = {
        'CALLER_ID' : callerID,
        **rm_data,  # Unpack repeat data from RepeatMasker annotations
        'RM_TRANSPOSITION' : transposition,
        **trf_data, # Unpack repeat data from TRF annotations
    }

    ####### 2. Determine Final Classification ######
    # Both TRF and RM
    trf_class = trf_data['TRF_CLASSIFICATION']
    rm_class = rm_data['RM_CLASSIFICATION']

    classification = get_final_classification(trf_class, rm_class, total_trf_coverage, total_rm_coverage, trf_data, rm_data)
    
    info.update({'FINAL_CLASSIFICATION' : classification})
    
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

    ####### 3. Obtain info from original VCF, and add annotations ######
    # Get record from SV VCF 
    records = sv_vcf.fetch(chrom, pos-1, pos+1)
    for record in records:
        # Check SnifflesID matches
        if callerID == record.id:
            record.id=sv_id
            # Update info fields with new annotations
            for annotation in info:
                record.info[annotation] = info[annotation]
            # Write the updated record to the VCF file
            vcf_file.write(record)
            break

    return info['FINAL_CLASSIFICATION'], info['RM_TRANSPOSITION'] 


#------------------------------------ TSV -------------------------------------#
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
    length = sv_data['length']
    sv_type = sv_id.split('.')[0]
    if 'RM' in sv_data:
        rm_entries = sv_data['RM']
        for entry in rm_entries:
            sv_coverage = round(entry['intersection'], 2)
            element_coverage = round(entry['element_coverage'], 2)
            element_proportion = round(entry['element_proportion'],2)
            classification = entry['class']
            repeat_type = entry['repeat']

            tsv_out.write(f'{sv_id}\t{chrom}\t{pos}\t{length}\t{sv_type}\t{sv_coverage}\t{element_coverage}\t{element_proportion}\t{classification}\t{repeat_type}\n')

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
    length = sv_data['length']
    sv_type = sv_id.split('.')[0]
    if 'TRF' in sv_data:
        trf_entries = sv_data['TRF']
        for entry in trf_entries:
            copy_number = entry['copy_number']
            period_size = entry['period_size']
            sv_coverage = round(entry['intersection'], 2)
            classification = entry['key']
            consensus_repeat = entry['motif']
            tsv_out.write(f'{sv_id}\t{chrom}\t{pos}\t{length}\t{sv_type}\t{sv_coverage}\t{copy_number}\t{period_size}\t{classification}\t{consensus_repeat}\n')

         
#----------------------------------- STDOUT -----------------------------------#
def print_results(title, repeat_count, total, output_total):
    """    
    Prints a formatted summary of the number of repeats (HOMO, STR, TR, SINE, LINE, LTR, DNA, Retroposon, Non-repetitive)
    """
    print(f"\n{title}:")
    print("-" * 35)
    for key, value in repeat_count.items():
        percentage = (value / total) * 100
        print(f"{key:<15}: {value:>5}  ({percentage:>6.2f}%)")
    if output_total:
        print("-" * 20)
        print(f"{'Total':<15}: {total:>5}")


#-------------------------------- MAIN OUTPUT ---------------------------------#
def output_annotations(sv_info, sv_vcf_file, vcf_output, tsv_rm_output, tsv_trf_output, min_repetitive, strchive):
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
    - tsv_rm_output (str): Path to the output TSV file for RepeatMasker annotations.
    - tsv_trf_output (str): Path to the output TSV file for TRF annotations.
    - min_repetitive (float): Minimum threshold for repeat coverage to consider a repeat as repetitive (e.g., 0.5 for 50%)
    - strchive (dict): STR disease/gene with pathogenic motif and number 
    """
    sv_vcf = pysam.VariantFile(sv_vcf_file)
    header = create_vcf_header(sv_vcf)
    vcf_out =  pysam.VariantFile(vcf_output, 'w', header=header)

    # Initialise count for breakdown of repeat types
    sv_count = 0
    repeat_count = {'HOMO' : 0, 'STR' : 0, 'TR' : 0, 'COMPLETE' : 0, 'FRAGMENT' : 0, 'NON_REPETITIVE' : 0}
    count_by_sv = {
        'INS': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : 0, 'LINE' : 0, 'LTR' : 0, 'DNA' : 0, 'Retroposon' : 0, 'NON_REPETITIVE' : 0},
        'DEL': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : 0, 'LINE' : 0, 'LTR' : 0, 'DNA' : 0, 'Retroposon' : 0, 'NON_REPETITIVE' : 0},
        'INV': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : 0, 'LINE' : 0, 'LTR' : 0, 'DNA' : 0, 'Retroposon' : 0, 'NON_REPETITIVE' : 0},
        'DUP': {'HOMO': 0, 'STR': 0, 'TR': 0, 'SINE' : 0, 'LINE' : 0, 'LTR' : 0, 'DNA' : 0, 'Retroposon' : 0, 'NON_REPETITIVE' : 0},
    }

    # Process each structural variant and write to VCF 
    for sv_id in sv_info:
        sv_type = sv_id.split('.')[0]
        sv_count += 1
        classification, transposition = create_vcf_record(vcf_out, sv_vcf, sv_info, sv_id, min_repetitive, strchive)
        if classification in repeat_count:
            repeat_count[classification] += 1
            count_by_sv[sv_type][classification] += 1
        # A repeat element 
        else: 
            repeat_count[transposition] += 1
            count_by_sv[sv_type][classification] += 1

    # Write RepeatMasker and TRF records to separate TSV files 
    with open(tsv_rm_output, 'w') as tsv_rm_out, open(tsv_trf_output, 'w') as tsv_trf_out:
        # Write headers for both files
        tsv_rm_out.write('ID\tCHR\tPOS\tSVLEN\tSV_TYPE\tSV_COVERAGE\tELEMENT_COVERAGE\tELEMENT_PROPORTION\tFAMILY\tREPEAT\n')
        tsv_trf_out.write('ID\tCHR\tPOS\tSVLEN\tSV_TYPE\tSV_COVERAGE\tCOPY_NUMBER\tPERIOD_SIZE\tCLASSIFICATION\tCONSENSUS_REPEAT\n')

        # Iterate once through sv_info and write to both files
        for sv_id in sv_info:
            create_rm_tsv_record(sv_info, sv_id, tsv_rm_out)  # Write to rm file
            create_trf_tsv_record(sv_info, sv_id, tsv_trf_out)  # Write to trf file

    # Print results summary
    print_results("Repeat Types", repeat_count, sv_count, True)
    # Print results summary by SV type 
    sv_type_df = pd.DataFrame(count_by_sv).T
        #     HOMO STR TR ...
        # INS
        # DEL
    print(sv_type_df)


################################################################################
if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

    parser = MyParser(description="parse RepeatMasker output",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-rm",
        help="Path to processed RepeatMasker .out file")
    parser.add_argument("-trf",
        help="trf file")
    parser.add_argument("-id",
        help="sv id file")
    parser.add_argument("-vcf",
        help="output vcf file")
    parser.add_argument("-trf_tsv",
        help="output tsv file")
    parser.add_argument("-rm_tsv",
        help="output tsv file")
    parser.add_argument("-sv_vcf",
        help="path to Sniffles VCF")
    parser.add_argument("-strchive",
        help="path to strchive bed file")
    parser.add_argument("-min", "--min_intersect",  type=float,
        help="Minimum intersection between repeat and SV e.g. 0.05 (5%) (0 < min_intersect < 1)", default=0.05)
    parser.add_argument("-mr", "--min_repetitive",  type=float,
        help="Minimum coverage of SV by repeats to be considered repetitive", default=0.5)
    parser.add_argument("-div", "--divisor",  type=float,
        help="Chosen intervals to prioritise period size over intersection (0 < divisor < 1)", default=0.05)
    

    args = parser.parse_args()
    sv_info = read_sv_info(args.id)
    # Check files
    check_files(args.rm, args.trf)

    # Read in data
    sv_info = read_rm(sv_info, args.rm, args.min_intersect, False)
    sv_info = read_trf(sv_info, args.trf, args.min_intersect)
    
    if args.strchive:
        strchive = load_strchive(args.strchive)
    else:
        strchive = None

    # Filter overlaps
    divisor = calculate_divisor(args.divisor)
    filtered_rm_info = filter_rm(sv_info)
    filtered_sv_info = filter_tr(filtered_rm_info, divisor)

    # Write to files 
    output_annotations(filtered_sv_info, args.sv_vcf, args.vcf, args.rm_tsv, args.trf_tsv, args.min_repetitive, strchive)
