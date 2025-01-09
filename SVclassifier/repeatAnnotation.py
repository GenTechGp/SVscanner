
import sys
import argparse
import pandas as pd
import pysam
import math
import os
from helpers.classifierHelpers import read_sv_info, read_trf, read_rm

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


################################################################################
#                                  FILTERING                                   #
################################################################################

#-----------------------------    RepeatMasker    -----------------------------#

def filter_rm(sv_info):
    """
    Determines non-overlapping transposable elements within each element type e.g. SINE
    
    """
    rm_elements = {}
    for sv_id in sv_info:
        rm_elements[sv_id] = sv_info[sv_id]
        
        transposition = "FRAGMENT"
        if 'RM' in sv_info[sv_id]:
            all_non_overlapping = []  # Flattened list for non-overlapping elements
            for element in sv_info[sv_id]['RM']:

                # Sum the coverage for the element class
                element_list = sv_info[sv_id]['RM'][element]
                # pprint(element_list)

                # Sort elements in decreasing order of intersect length
                sorted_elements = sorted(element_list, key=lambda x: x.get('intersection', 0), reverse=True)

                # Initialize an empty list to hold non-overlapping TEs
                non_overlapping = []

                # Function to check overlap
                def is_overlapping(te1, te2):
                    return not (te1['te_end'] <= te2['te_start'] or te2['te_end'] <= te1['te_start'])

                # Iterate through sorted TEs and select non-overlapping ones
                transposition_fraction = 0
                element_coverage_complete = True
                for te in sorted_elements:
                    if all(not is_overlapping(te, existing_te) for existing_te in non_overlapping):
                        non_overlapping.append(te)
                        # Determine if transposition is complete or not complete
                        transposition_fraction += te['intersection']
                        if te['element_coverage'] < 0.75:
                            element_coverage_complete = False

                # print(f'{element} : {transposition_fraction}')
                if transposition_fraction > 0.75 and element_coverage_complete:
                    transposition = "COMPLETE"

                # Append non-overlapping elements from this element group to the flat list
                all_non_overlapping.extend(non_overlapping)
            
            # Order by start position 
            all_non_overlapping = sorted(all_non_overlapping, key=lambda x: x['te_start'])
            # all_non_overlapping.update({'transposition' : transposition})

            # Assign the flattened non-overlapping list back to the new dictionary
            rm_elements[sv_id].update({'RM': all_non_overlapping})
            rm_elements[sv_id].update({'transposition' : transposition})

    return rm_elements

#---------------------------  Tandem Repeat Finder ----------------------------#


def sumFractions(intersections):
    total_fraction = 0
    for entry in intersections: 
        total_fraction += entry['intersection'] 
    return round(total_fraction,2)

def filter_tr(sv_info, interval_divisor):
    for sv_id in sv_info:
        if 'TRF' in sv_info[sv_id]:
            trf_list = sv_info[sv_id]['TRF']
            
            # Custom key function that groups intersections by 0.05 bins (5%)
            def custom_sort_key(trf):
                # Prioritises largest intersection length in 5% intervals
                # Priotisies lowest period size within the 5% intervals
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

            sv_info[sv_id]['TRF'] = non_overlapping
    
            total_fraction = sumFractions(non_overlapping)
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
        - sv_start (int): The start position of the SV.
        - sv_end (int): The end position of the SV.
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

    return max_repeat, fraction_sum[max_repeat]

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
    header.add_line('##command="annotateVCF.py [-h] [-rm RM] [-id ID] [-vcf VCF] [-tsv TSV] [-sv_vcf SVVCF]"')
    header.add_line("##INFO=<ID=CALLER_ID,Number=1,Type=String,Description=\"Caller ID for the SV\">")
    ## Header Lines for RepeatMasker Annotations
    header.add_line("##INFO=<ID=RM_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of repeat class covering the SV [SINE,LINE,LTR,DNA,Retroposon or NON-REPETITIVE]\">")
    header.add_line("##INFO=<ID=RM_ELEMENTS_COVERAGE,Number=1,Type=String,Description=\"Fraction of the mobile element covered by the SV\">")
    header.add_line("##INFO=<ID=RM_ELEMENT_PROPORTION,Number=1,Type=String,Description=\"Proportion of the query sequence (includes flanking region) found in the mobile element\">")
    header.add_line("##INFO=<ID=RM_TRANSPOSITION,Number=1,Type=String,Description=\"Type of transposition [COMPLETE/FRAGMENT]\">")
    header.add_line("##INFO=<ID=RM_SV_COVERAGE,Number=1,Type=String,Description=\"Fraction of the SV covered by the mobile element\">")
    header.add_line("##INFO=<ID=RM_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of SV covered by mobile elements\">")

    ## Header Lines for Tandem Repeat Finder Annotations
    header.add_line("##INFO=<ID=TRF_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of tandem repeat class covering the SV [HOMO,STR,TR or NON-REPEPTITIVE]\">")
    header.add_line("##INFO=<ID=TRF_TOTAL_SV_COVERAGE,Number=1,Type=String,Description=\"Total coverage of SV\">")
    header.add_line("##INFO=<ID=TRF_SV_COVERAGE,Number=1,Type=String,Description=\"Coverage of SV\">")
    header.add_line("##INFO=<ID=TRF_PERIOD_SIZE,Number=1,Type=String,Description=\"Period size of the repeat\">")
    header.add_line("##INFO=<ID=TRF_COPY_NUMBER,Number=1,Type=String,Description=\"Total copy number of the repeat\">")
    header.add_line("##INFO=<ID=CONSENSUS_REPEAT,Number=1,Type=String,Description=\"Motif of repeat found by Tandem Repeat Finder\">")
    header.add_line("##INFO=<ID=FINAL_CLASSIFICATION,Number=1,Type=String,Description=\"Classification of SV as repetitive element based on TRF and RepeatMasker results\">")

    return header


def create_vcf_record(vcf_file, sv_vcf, sv_info, sv_id, min_repetitive):
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
    callerID = sv_data['callerID']

    ####### 1. Get annotations for RepeatMasker and TRF ######
    transposition = sv_data['transposition'] if 'transposition' in sv_data else 'NA'

    rm_data = prepare_rm_data(sv_data.get('RM', []))
    repeat_data = prepare_trf_data(sv_data.get('TRF', []))
    
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
    repeat_data.update({'TRF_TOTAL_SV_COVERAGE': total_trf_coverage})
    if total_trf_coverage != 'NA' and float(total_trf_coverage) < min_repetitive:
        repeat_data.update({key: 'NA' for key in repeat_data if key != 'TRF_TOTAL_SV_COVERAGE'})
        repeat_data['TRF_CLASSIFICATION'] = 'NON_REPETITIVE'

    # Unpacks annotations
    info = {
        'CALLER_ID' : callerID,
        **rm_data,  # Unpack repeat data from RepeatMasker annotations
        'RM_TRANSPOSITION' : transposition,
        **repeat_data, # Unpack repeat data from TRF annotations
    }

    ####### 2. Determine Final Classification ######
    # Both TRF and RM
    trf_class = repeat_data['TRF_CLASSIFICATION']
    rm_class = rm_data['RM_CLASSIFICATION']

    if (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
        # TRF has higher coverage
        if total_trf_coverage > total_rm_coverage:
            trf_classes = trf_class.split(',')
            # Determine main classification if SV is made up of different repeat types
            if len(trf_classes) > 1 and not all(x == trf_classes[0] for x in trf_classes):
                coverages = repeat_data['TRF_SV_COVERAGE'].split(',')
                classification, fraction = highest_fraction(trf_classes, coverages)
            else:
                classification = trf_classes[0]
        # Mobile element from RepeatMasker has higher coverage
        else:
            rm_classes = rm_class.split(',')
            if len(rm_classes) > 1 and not all(x == rm_classes[0] for x in rm_classes):
                coverages = rm_data['RM_SV_COVERAGE'].split(',')
                classification, fraction = highest_fraction(rm_classes, coverages)
            else: 
                classification = rm_classes[0]
    # TRF (no repeatMasker)
    elif (trf_class not in ['NA', 'NON_REPETITIVE']) and (rm_class in ['NA', 'NON_REPETITIVE']):
        trf_classes = trf_class.split(',')
        # Determine main classification if SV is made up of different repeat types
        if len(trf_classes) > 1 and not all(x == trf_classes[0] for x in trf_classes):
            coverages = repeat_data['TRF_SV_COVERAGE'].split(',')
            classification, fraction = highest_fraction(trf_classes, coverages)
        else:
            classification = trf_classes[0]
    # RM (no TRF intersect)
    elif (trf_class in ['NA', 'NON_REPETITIVE']) and (rm_class not in ['NA', 'NON_REPETITIVE']):
        rm_classes = rm_class.split(',')
        if len(rm_classes) > 1 and not all(x == rm_classes[0] for x in rm_classes):
            coverages = rm_data['RM_SV_COVERAGE'].split(',')
            classification, fraction = highest_fraction(rm_classes, coverages)
        else: 
            classification = rm_classes[0]
    # None
    else: 
        classification = 'NON_REPETITIVE'
    
    info.update({'FINAL_CLASSIFICATION' : classification})
    
    ####### 3. Obtain info from original VCF, and add annotations ######
    # Get record from SV VCF 
    records = sv_vcf.fetch(chrom, pos-1, pos+1)
    for record in records:
        # Check SnifflesID matches
        if callerID == record.id:
            allele_tuple = (record.ref, record.alts[0])

            new_record = vcf_file.new_record(
                contig=record.chrom,
                start=record.pos,
                stop=record.stop,
                id=sv_id,
                alleles=allele_tuple,
                qual=record.qual,
                filter=record.filter,
                info=record.info,
            )
            # Update info fields with new annotations
            for annotation in info:
                new_record.info[annotation] = info[annotation]
            # Add sample data (e.g., genotypes) from original record
            for sample in record.samples:
                for sample_format in record.format:
                    # Access the value for the current sample and format
                    new_record.samples[sample][sample_format] = record.samples[sample][sample_format]
            # Write the updated record to the VCF file
            vcf_file.write(new_record)

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
def output_annotations(sv_info, sv_vcf_file, vcf_output, tsv_rm_output, tsv_trf_output, min_repetitive):
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
    - min_repetitive (float): Minimum threshold for repeat coverage to consider a repeat as repetitive (e.g., 0.5 for 50%).
    """
    sv_vcf = pysam.VariantFile(sv_vcf_file)
    header = create_vcf_header(sv_vcf)
    vcf_out =  pysam.VariantFile(vcf_output, 'w', header=header)
    
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
        classification, transposition = create_vcf_record(vcf_out, sv_vcf, sv_info, sv_id, min_repetitive)
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

    # create_plot(count_by_sv, vcf_output)

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

    # Filter overlaps
    divisor = calculate_divisor(args.divisor)
    filtered_rm_info = filter_rm(sv_info)
    filtered_sv_info = filter_tr(filtered_rm_info, divisor)

    # Write to files 
    output_annotations(filtered_sv_info, args.sv_vcf, args.vcf, args.rm_tsv, args.trf_tsv, args.min_repetitive)


# ################################################################################
# #                              PLOT CREATION                                   #
# ################################################################################
# def create_plot(data, output_path):
#     print(f"\nStructural Variant Types:")
#     print("-" * 35)
#     directory = os.path.dirname(output_path)

#     # Convert the dictionary to a DataFrame
#     directory = os.path.dirname(output_path)

#     # Convert the dictionary to a DataFrame
#     df = pd.DataFrame(data)

#     # Reset the index to get 'repeat_type' as a column
#     df = df.reset_index().rename(columns={'index': 'repeat_type'})

#     # Melt the dataframe to long format for seaborn plotting
#     df_melt = df.melt(id_vars='repeat_type', var_name='SV_type', value_name='count')

#     # Define shades for each SV type (DEL, DUP, INS, INV) - these will be colors
#     sv_palette = {
#         'INS': '#9BBB59',    # Green
#         'DEL': '#4BACC6',    # Blue
#         'INV': '#8064A2',    # Purple
#         'DUP': '#DC127C',    # Pink
#     }

#     # Filter data for INS and DEL plot
#     df_melt_ins_del = df_melt[df_melt['SV_type'].isin(['INS', 'DEL'])]

#     # Create the first plot for INS and DEL
#     plt.figure(figsize=(12, 6))
#     g1 = sns.barplot(x='repeat_type', y='count', hue='SV_type', data=df_melt_ins_del, palette=sv_palette)
#     plt.title('Counts of Repeats by SV Type (INS and DEL)')
#     plt.xlabel('Repeat Type')
#     plt.ylabel('Count')
#     plt.legend(title='SV Type', loc='upper right', bbox_to_anchor=(1.1, 1))
#     plt.tight_layout()
#     plt.savefig(f'{directory}/count_ins_del.png')  # Save the INS and DEL plot
#     plt.clf()  # Clear the figure for the next plot

#     # Filter data for INV and DUP plot
#     df_melt_inv_dup = df_melt[df_melt['SV_type'].isin(['INV', 'DUP'])]

#     # Create the second plot for INV and DUP
#     plt.figure(figsize=(12, 6))
#     g2 = sns.barplot(x='repeat_type', y='count', hue='SV_type', data=df_melt_inv_dup, palette=sv_palette)
#     plt.title('Counts of Repeats by SV Type (INV and DUP)')
#     plt.xlabel('Repeat Type')
#     plt.ylabel('Count')
#     plt.legend(title='SV Type', loc='upper right', bbox_to_anchor=(1.1, 1))
#     plt.tight_layout()
#     plt.savefig(f'{directory}/count_inv_dup.png')  # Save the INV and DUP plot
#     plt.clf()  # Clear the figure for the next plot
    