#!/usr/bin/env python3

import sys
import argparse

def process_trf_file(filename):
    """
    Processes TRF (Tandem Repeats Finder) output file.

    Args:
    filename (str): The path to the TRF output file (file ends in .dat).

    Returns:
    list: A list of unique sequences trandem repeats extracted from the TRF.dat file, with each sequence being at least 10 bp long.
    """
    sequences = []
    with open(filename, 'r') as f:
        for line in f:
            # Check if line starts with a number.
            if line[0].isdigit():
                # Split the line by space and extract the sequence.
                sequence = line.split()[-2]
                if len(sequence) >= 10: # Check if sequence is longer than 10 bp.
                    sequences.append(sequence)

    # Remove duplicates.
    sequences = list(dict.fromkeys(sequences))
    
    return sequences

def main():
    """
    Main function to execute the script.
    Processes a TRF file and outputs the sequences in FASTA format.
    """
    parser = argparse.ArgumentParser(description="Extract sequences from TRF output and build a FASTA formated tandem repeat library. The script extracts sequences from the TRF output, ensuring they are unique and at least 10 bp long, and then outputs them in FASTA format with headers based on the specified species.")
    
    # Required argument for the input TRF file.
    parser.add_argument('-in', '--input', required=True, help="Input TRF file (file ends in .dat)")
    
    # Optional argument for specifying the species initials.
    parser.add_argument('-species', '--species', default='Zm', help="Species initials (e.g., 'Zm' for Zea mays). Defaults to 'Zm' if not specified.")

    args = parser.parse_args()

    # Process the TRF file to get sequences.
    sequences = process_trf_file(args.input)
    
    # Output the sequences in FASTA format.
    for idx, seq in enumerate(sequences, start=1):
        print(f">{args.species}_trf_{idx}#Satellite/Satellite")
        print(seq)

if __name__ == "__main__":
    main()
