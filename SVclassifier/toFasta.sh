#!/bin/bash
#$ -S /bin/bash
#$ -terse
#$ -cwd
#$ -l mem_requested=8G 
#$ -l tmp_requested=60G 
#$ -pe smp 16
#$ -N duplicationChecker
#$ -V

input_file=$1
fasta_file=$2
id_file=$3

# Initialize the fasta file
> $id_file
> $fasta_file
  
while IFS=$'\t' read -r chr startFlank endFlank pos end len id sequence sv_sequence type callerID; do
    # Create the info 
    info="${chr}\t${startFlank}\t${endFlank}\t${pos}\t${end}\t${len}\t${id}\t${callerID}"
    
    # Write to the fasta file in FASTA format
    echo ">${id}" >> "$fasta_file"
    echo "${sequence}" >> "$fasta_file"
    echo -e "${info}" >> "$id_file"
done < "$input_file"