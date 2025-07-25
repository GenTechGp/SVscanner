import sys
import os
import argparse
import gzip


def parse_info_field(info):
    """Parses the INFO field of a VCF entry into a dictionary."""
    return dict(item.split('=') for item in info.split(';') if '=' in item)

def open_vcf(filepath):
    """Opens a VCF file, supporting both plain text and gzipped formats."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def process_vcf(vcf_file, output_dir):
    output_path = os.path.join(output_dir, 'end_mismatch.tsv')
    with open_vcf(vcf_file) as vcf, open(output_path, 'w') as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\tEND\tSVLEN\n")
        sv_types = {'DEL': 0, 'INS': 0, 'DUP': 0, 'INV': 0, 'BND': 0}
        symbolic_types = {'DEL': 0, 'INS': 0, 'DUP': 0, 'INV': 0, 'BND': 0}
        for line in vcf:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            info_dict = parse_info_field(info)
            sv_type = info_dict.get('SVTYPE')
            sv_types[sv_type] += 1

            if alt.startswith('<') and alt.endswith('>'):
                symbolic_types[sv_type] += 1
                continue  # Skip symbolic ALT entries

            info_end_str = info_dict.get('END')
            svlen_str = info_dict.get('SVLEN')
            pos = int(pos)
            ref_len = len(ref)
            alt_len = len(alt)
            calculated_end = pos + ref_len - 1

            if sv_type not in ['DEL', 'INS']:
                calculated_svlen = calculated_end - pos
            else:
                calculated_svlen = alt_len - ref_len


            if info_end_str is None or int(info_end_str) != calculated_end:
                # Report the mismatch or missing END
                assert int(info_end_str) == pos
                assert ref_len == alt_len
                out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{calculated_end}\t{calculated_svlen}\n")
                # out.write(f"{chrom}\t{pos}\t{vid}\t{sv_type}\t{calculated_end}\t{calculated_svlen}\n")
            else:
                if svlen_str is not None:
                    assert int(svlen_str) == calculated_svlen, (
                        f"Assertion failed: SVLEN {svlen_str} does not match calculated {calculated_svlen} "
                        f"for {chrom}:{pos} ID={fields}")
        # print the summary of SV types
        print("SVTYPE counts:")
        for sv_type, count in sv_types.items():
            print(f"{sv_type}: {count}")
        print(f"Total entries processed: {sum(sv_types.values())}")

        # print the summary of symbolic types
        print("Symbolic SVTYPE counts:")
        for sv_type, count in symbolic_types.items():
            print(f"{sv_type}: {count}")
        print(f"Total symbolic entries processed: {sum(symbolic_types.values())}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check VCF INFO.END field against POS+len(REF)-1')
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be .gz)')
    parser.add_argument('--out', required=True, help='Output directory for the TSV file')

    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    process_vcf(args.vcf, args.out)
