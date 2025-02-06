import argparse
import random

# Function to generate a random DNA sequence of a given length
def generate_random_dna(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choices(bases, k=length))

# Function to write sequences to a .fa file
def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">Sequence_{i}\n")  # Write the header
            f.write(f"{seq}\n")          # Write the sequence

# Main function
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate random DNA sequences and write them to a .fa file.")
    parser.add_argument('-l', '--length', type=int, required=True, help="Length of each DNA sequence.")
    parser.add_argument('-n', '--number', type=int, required=True, help="Number of sequences to generate.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file path.")
    parser.add_argument('-s', '--seed', type=int, default=None, help="Seed for random number generation (optional).")
    args = parser.parse_args()

    # Set the seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
        print(f"Using seed: {args.seed}")

    # Generate random DNA sequences
    sequences = [generate_random_dna(args.length) for _ in range(args.number)]

    # Write sequences to the output .fa file
    write_sequences_to_fasta(sequences, args.output)
    print(f"{args.number} sequences of length {args.length} written to {args.output}")

# Run the script
if __name__ == "__main__":
    main()