import argparse
import random

# Function to generate a random DNA sequence of a given length
def generate_random_dna(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choices(bases, k=length))

# Function to generate a sequence with tandem repeats at a random position
def generate_sequence_with_tandem_repeats(total_length, min_repeat_length=1, max_repeat_length=2000, min_repeats=2):
    # Randomly decide the length of the tandem repeat unit
    max_repeat_length = min(total_length, max_repeat_length)
    repeat_unit_length = random.randint(min_repeat_length, max_repeat_length)
    # Randomly generate the tandem repeat unit
    repeat_unit = generate_random_dna(repeat_unit_length)
    # Randomly decide the number of repeats
    max_repeats = int(total_length/repeat_unit_length) #todo: this assumes the total length (including the flanking regions have repeats. verify how often this is correct by looking at real-world SVs)
    num_repeats = random.randint(min_repeats, max_repeats)
    # Generate the tandem repeat region
    tandem_repeat_region = repeat_unit * num_repeats
    # Randomly decide the starting position of the tandem repeat region
    start_pos = random.randint(0, total_length - len(tandem_repeat_region))
    # Generate the rest of the sequence with random DNA
    before_repeat = generate_random_dna(start_pos)
    after_repeat = generate_random_dna(total_length - start_pos - len(tandem_repeat_region))
    # Combine the sequences
    full_sequence = before_repeat + tandem_repeat_region + after_repeat
    return full_sequence, repeat_unit, num_repeats, start_pos

# Function to write sequences to a .fa file
def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for i, (seq, repeat_unit, num_repeats, start_pos) in enumerate(sequences, start=1):
            f.write(f">Sequence_{i} TandemRepeat: {repeat_unit} x{num_repeats} StartPos: {start_pos}\n")  # Write the header
            f.write(f"{seq}\n")  # Write the sequence

# Main function
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate DNA sequences with tandem repeats at random positions and write them to a .fa file.")
    parser.add_argument('-l', '--length', type=int, required=True, help="Length of each DNA sequence.")
    parser.add_argument('-n', '--number', type=int, required=True, help="Number of sequences to generate.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file path.")
    parser.add_argument('-s', '--seed', type=int, default=None, help="Seed for random number generation (optional).")
    args = parser.parse_args()

    # Set the seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
        print(f"Using seed: {args.seed}")

    # Generate sequences with tandem repeats at random positions
    sequences = []
    for _ in range(args.number):
        seq, repeat_unit, num_repeats, start_pos = generate_sequence_with_tandem_repeats(args.length)
        sequences.append((seq, repeat_unit, num_repeats, start_pos))

    # Write sequences to the output .fa file
    write_sequences_to_fasta(sequences, args.output)
    print(f"{args.number} sequences of length {args.length} written to {args.output}")

# Run the script
if __name__ == "__main__":
    main()