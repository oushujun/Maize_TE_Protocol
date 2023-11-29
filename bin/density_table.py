import argparse
import subprocess
import os
import sys

def check_bed_files():
    """
    Checks for the presence of .bed files in the current directory.
    If found, exits the script with a warning.
    """
    if any(file.endswith('.bed') for file in os.listdir('.')):
        sys.exit("Error: bed files detected in the working directory. "
                 "Please remove them or navigate to a new directory.")

def reformat_gff(input_gff, output_file):
    """
    Reformats the GFF file to include only a subset of columns.
    """
    with open(input_gff, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                outfile.write(f"{parts[0]}\t{parts[3]}\t{parts[4]}\t{parts[2]}\n")

def create_bed_files(input_file):
    """
    Creates .bed files for each type found in the GFF file.
    """
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            with open(parts[3] + ".bed", 'a') as outfile:
                outfile.write(line)

def create_genome_fai(genome_file):
    """
    Creates a genome index file using samtools.
    """
    subprocess.run(["samtools", "faidx", genome_file], stderr=subprocess.DEVNULL)

def make_windows(fai_file, output_file):
    """
    Generates window sizes across the genome using bedtools.
    """
    modified_fai = "genome.fa.fai2"
    with open(fai_file, 'r') as infile, open(modified_fai, 'w') as outfile:
        for line in infile:
            parts = line.split('\t')
            outfile.write(f"{parts[0]}\t{parts[1]}\n")

    with open(output_file, 'w') as outfile:
        subprocess.run(["bedtools", "makewindows", "-g", modified_fai, "-w", "1000000", "-s", "500000"], stdout=outfile, stderr=subprocess.DEVNULL)

def calculate_density(window_file):
    """
    Calculates the density for each TE type using bedtools coverage.
    """
    for bed_file in os.listdir('.'):
        if bed_file.endswith('.bed') and bed_file != window_file:
            type_name = bed_file.replace('.bed', '')
            density_file = type_name + '.density.bed'
            with open(density_file, 'w') as outfile:
                subprocess.run(["bedtools", "coverage", "-a", window_file, "-b", bed_file], stdout=outfile, stderr=subprocess.DEVNULL)

def merge_and_format():
    """
    Merges and formats the output.
    """
    for density_file in os.listdir('.'):
        if density_file.endswith('.density.bed'):
            type_name = density_file.replace('.density.bed', '')
            with open(density_file, 'r') as infile:
                for line in infile:
                    parts = line.strip().split('\t')
                    output_line = f"{parts[0]}\t{parts[1]}\t{parts[6]}\t{type_name}"
                    print(output_line)

def main():
    """
    Main function to orchestrate the density extraction workflow.
    """
    parser_description = (
        "Extracts density information for each TE superfamily from the EDTA GFF files.\n"
        "Requires python3, bedtools, and samtools. All can be installed with conda.\n"
        "Be sure to run this in a directory where there are not already bed files present."
    )
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument('-genome', type=str, required=True, help='Path to the genome file (genome.fa)')
    parser.add_argument('-gff', type=str, required=True, help='Path to the EDTA GFF file (genome.mod.EDTA.TEanno.gff3)')

    args = parser.parse_args()

    try:
        check_bed_files()
        reformatted_gff = "reformat.gff"
        reformat_gff(args.gff, reformatted_gff)
        create_bed_files(reformatted_gff)
        create_genome_fai(args.genome)
        window_file = "windows.bed.sizes"
        make_windows(args.genome + ".fai", window_file)
        calculate_density(window_file)
        merge_and_format()
    except Exception as e:
        print(f"Error encountered: {e}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
