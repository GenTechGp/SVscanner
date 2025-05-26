import os
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Bin records by length with capacity limit and write to single file.")
    parser.add_argument("--input", help="Path to input file")
    parser.add_argument("--output", help="Path to output file")
    parser.add_argument("--bin_interval", type=int, default=10000, help="Length range for each bin (default: 1000)")
    parser.add_argument("--bin_capacity", type=int, default=500, help="Max number of records per bin (default: 1000)")
    return parser.parse_args()

def main():
    args = parse_args()

    bin_counts = defaultdict(int)
    accepted_records = []

    with open(args.input, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # skip malformed lines
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue  # skip lines with non-integer coordinates

            length = end - start
            motif = parts[3]
            motif_count = length // len(motif)
            if length < 50 or motif_count < 2:
                continue  # skip invalid lengths
            
            
            bin_index = (length - 1) // args.bin_interval
            if bin_counts[bin_index] < args.bin_capacity:
                record = f"{parts[0]}\t{start}\t{end}\t{motif_count}x{motif}\n"
                accepted_records.append(record)
                bin_counts[bin_index] += 1

    # Write to a single output file
    with open(args.output, 'w') as out_file:
        out_file.writelines(accepted_records)

    # Print bin summary
    print("Bin summary:")
    for bin_index in sorted(bin_counts):
        low = bin_index * args.bin_interval + 1
        high = (bin_index + 1) * args.bin_interval
        print(f"Bin {bin_index} (length {low}–{high}): {bin_counts[bin_index]} records")

if __name__ == "__main__":
    main()
