import argparse
import re

# Parse the GFF file and extract gene features.
def parse_gff(gff_file):
    with open(gff_file) as gff:
        for line in gff:
            # Skip comment lines and empty lines.
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            # Extract chrom ID, gene start position and and gene end position.
            if parts[2] == 'gene':
                yield parts[0], int(parts[3]), int(parts[4])

# Modify a FASTA file based on gene features, minimum length, and hardmasking.
def modify_fasta(fasta_file, genes, min_length, hardmask):
    with open(fasta_file) as fasta:
        current_seq = bytearray()  # Mutable bytearray for the current sequence.
        current_seq_name = ''
        for line in fasta:
            # If line is a header, process the previous sequence.
            if line.startswith('>'):
                if current_seq:
                    yield process_sequence(current_seq, current_seq_name, genes, min_length, hardmask)
                current_seq_name = line[1:].strip()
                current_seq.clear()
            else:
                # Extend the current sequence with the new line.
                current_seq.extend(line.strip().encode())
        # Process the last sequence in the file.
        if current_seq:
            yield process_sequence(current_seq, current_seq_name, genes, min_length, hardmask)

# Process sequence by unmasking genes, applying minlength thresholded for retention of masked bases, and hardmasking.
def process_sequence(current_seq, current_seq_name, genes, min_length, hardmask):
    # Unmask gene regions if gene list is provided.
    if genes:
        for gene in genes:
            if gene[0] == current_seq_name:
                start, end = gene[1] - 1, gene[2]
                current_seq[start:end] = current_seq[start:end].upper()

    # Unmask short streteches of masked sequences based on minlength.
    if min_length > 1:
        regex = re.compile(r'(?<![a-z])[a-z]{1,' + str(min_length-1) + r'}(?![a-z])')
        current_seq = regex.sub(lambda x: x.group().upper(), current_seq.decode()).encode()
    elif min_length == 1:
        current_seq = current_seq.upper()

    # Hardmask the sequence by replacing remaining masked sequences with 'N'.
    if hardmask:
        current_seq = re.sub(rb'[a-z]', b'N', current_seq)

    # Return the modified sequence in FASTA format.
    return f'>{current_seq_name}\n{current_seq.decode()}\n'

# Main function to execute the script logic.
def main(fasta_file, gff_file, min_length, hardmask):
    # Parse the GFF file if provided.
    genes = list(parse_gff(gff_file)) if gff_file else []
    # Modify the FASTA file based on the parsed genes and other parameters.
    modified_fasta = modify_fasta(fasta_file, genes, min_length, hardmask)
    # Print each modified sequence.
    for seq in modified_fasta:
        print(seq, end='')

if __name__ == '__main__':
    # Create an argument parser with usage discription.
    parser = argparse.ArgumentParser(
        description='This tool processes a FASTA file by unmasking repeat-masked sequences that overlap with gene annotations in gff, '
                    'Unmasks sequences shorter than the specified minimum length to help open genespace, '
                    'and optionally, hardmasking the remaining softmasked sequences. '
                    'Use this to prepare genome fasta for gene annotation with MAKER/BRAKER2/etc.'
    )
    # Define command-line arguments with detailed help messages.
    parser.add_argument('-genome', type=str, required=True, help='The path to the FASTA file of the genome.')
    parser.add_argument('-gff', type=str, required=False, help='The path to the GFF file containing gene features. This argument is optional.')
    parser.add_argument('-minlength', type=int, default=500, help='The minimum length of continuous softmasked sequence to be retainined in the output. Sequences shorter than minlength will unmasked. Default:500')
    parser.add_argument('-hardmask', action='store_true', help='Enable hardmasking to convert all remaining softmasked sequence to "N" after processing gene features and minlength sequences.')
    # Parse the arguments.
    args = parser.parse_args()

    # Run the main function with the parsed arguments.
    main(args.genome, args.gff, args.minlength, args.hardmask)
