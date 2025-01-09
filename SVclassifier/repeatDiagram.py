
import sys
import argparse
from helpers.classifierHelpers import read_sv_info, read_trf, read_rm

def format_trf_info(info, diagram, diagram_length, field_width=20):
    info_str = f" {info}".ljust(field_width)  # Left-align the info within a fixed-width field
    diagram_str = f"{diagram}".ljust(diagram_length + 1)  # Ensure diagram is aligned by adding a marker (*)
    return f"{info_str}{diagram_str}"

################################################################################
#                                CREATE DIAGRAM                                #
################################################################################
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

################################################################################
#                                OUTPUT DIAGRAM                                #
################################################################################
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

################################################################################
#                            MAIN DIAGRAM FUNCTION                             #
################################################################################
            
def get_repeat_info(sv_data, key):
    if key in sv_data:
        return sv_data[key]
    return []

def create_SV(sv_info, diagram_length, output_file):
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
    motif_length = 80
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

            sv_start = sv_data['rel_start']
            sv_end = sv_data['rel_end']
            sv_length = sv_data['length']

            flanking_scale = diagram_length / end
            sv_start_scaled = int((sv_start - start) * flanking_scale)
            sv_end_scaled = int((sv_end - start) * flanking_scale)


            ###################### SV & FLANKING DIAGRAM #######################
            # Don't generate the flanking diagram if sv ends after the flanking
            if sv_end > end and sv_end_scaled > 100:
                id_str = id_str.replace('\t', ' ')
                print(f"Unable to generate sv + flanking diagram for {id_str}")
                # print(f'\tSV_Start {start}   SV_End {end}')
                # print(f'\tRelative Start {sv_start}   Relative {sv_end}')
                # print(f'Relative Start {sv_start}   Relative {sv_end}')
                flanking_diagram = None
            else:
                flanking_diagram = create_diagram(diagram_length, sv_start_scaled, sv_end_scaled)
                flanking_diagram = format_trf_info('SV & flanking', flanking_diagram, diagram_length)
            
            
            ## REPEAT MASKER 
            for classification in rm:
                elements = rm[classification]
                rm_diagram = [' '] * diagram_length
                info = classification
                intersect_percentages = []
                for element in elements:
                    rm_diagram = create_rm(rm_diagram, element, diagram_length, flanking_scale, start)
                    intersection = round(element['intersection']*100, 2)
                    intersect_percentages.append(f'{intersection}%')

                rm_diagram = ''.join(rm_diagram)
                rm_diagram = format_trf_info(f' {info}', rm_diagram, diagram_length)
                intersects = ','.join(intersect_percentages)
                    
                rm_output_flanking.append(f'{rm_diagram}\t{intersects}\n')
            
            ## TANDEM REPEAT FINDER  
            # Write out repeats as fraction of flanking region with SV
            for repeat in trf:                
                trf_diagram = create_trf(repeat, diagram_length, flanking_scale, start)
                intersect = round(repeat['intersection']*100,2)
                trf_output_flanking.append(f'{trf_diagram}\t{intersect}%\n')

            ########## SV DIAGRAM: Write out repeats as fraction of SV #########
            sv_scale = diagram_length / sv_length
            sv_diagram = create_diagram(diagram_length, 0, diagram_length-1)
            sv_diagram = format_trf_info('SV Diagram', sv_diagram, diagram_length)

            ## REPEAT MASKER
            for classification in rm:

                elements = rm[classification]
                repeat_names = []

                rm_diagram = [' '] * diagram_length
                info = classification
                for element in elements:
                    repeat_name = element['repeat']

                    rm_diagram = create_rm(rm_diagram, element, diagram_length, sv_scale, sv_start)
                    if repeat_name not in repeat_names:
                        repeat_names.append(repeat_name)
               
                rm_diagram = ''.join(rm_diagram)
                rm_diagram = format_trf_info(f' {info}', rm_diagram, diagram_length)
                names = ','.join(repeat_names)
                rm_output.append(f'{rm_diagram}\t{names}\n')

            # TANDEM REPEAT FINDER
            for repeat in trf:
                motif = repeat['motif']
                if len(motif) > motif_length:
                    motif = f'*{motif_length}plus'

                trf_diagram = create_trf(repeat, diagram_length, sv_scale, sv_start)
                trf_output.append(f'{trf_diagram}\t{motif}\n')

            # Output the diagrams
            if rm_output_flanking != [] or trf_output_flanking != [] or rm_output != [] or trf_output != []:
                f.write(f"{id_str}\n")
                if flanking_diagram:
                    write_flanking(f, flanking_diagram, rm_output_flanking, trf_output_flanking)
                
                write_SV(f, sv_diagram, rm_output, trf_output)
                f.write(f"\n")


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

    parser = MyParser(description="parse RepeatMasker output",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-rm",
        help="RepeatMasker out file")
    parser.add_argument("-trf",
        help="RepeatMasker out file")
    parser.add_argument("-sv",
        help="SV ID file")
    parser.add_argument("-length",
        help="Diagram length", type=int, default=100)
    parser.add_argument("-out",
        help="output_file")
    parser.add_argument("-min", "--min_intersect",  type=float,
        help="Minimum intersection between repeat and SV e.g. 0.05 (5%) (0 < min_intersect < 1)", default=0.05)

    args = parser.parse_args()

    sv_info = read_sv_info(args.sv)

    if (args.trf):
        sv_info = read_trf(sv_info, args.trf, args.min_intersect)
    
    if (args.rm):
        sv_info = read_rm(sv_info, args.rm, args.min_intersect, True)

    create_SV(sv_info, args.length, args.out)