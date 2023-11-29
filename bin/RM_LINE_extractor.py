import sys
import argparse

def process_fasta(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    counter_l1 = 1
    counter_rte = 1
    new_records = []

    # Boolean to remember if the current sequence should be included.
    keep_seq = False
    current_seq = ""

    for line in lines:
        line = line.strip()
        
        # If it's a header line.
        if line.startswith('>'):
            # If we were keeping the previous sequence, add it to the new records.
            if keep_seq:
                new_records.append(current_seq)
                current_seq = ""
            
            header_clean = line.split(' ')[0]  # Ignore anything after the first space.
            if "#LINE/L1" in header_clean: # Add the necessary RepeatMasker TE suffix.
                new_records.append(f">Zm_LINE_{counter_l1}#LINE/L1")
                counter_l1 += 1
                keep_seq = True
            elif "#LINE/RTE-BovB" in header_clean:
                new_records.append(f">Zm_LINE_{counter_rte}#LINE/RTE")
                counter_rte += 1
                keep_seq = True
            else:
                keep_seq = False
        else:
            if keep_seq:
                current_seq += line

    # Handling the last sequence if it should be kept.
    if keep_seq:
        new_records.append(current_seq)

    return "\n".join(new_records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts LINEs from the RepeatModeler/RepeatMasker outfile and assigns a RepeatMasker-recognizable TE ontology to the header.')
    parser.add_argument('-in', '--input', type=str, required=True, help='Path to input FASTA file (ie, consensi.fa.classified). Its in the RM_* dir.')

    args = parser.parse_args()

    transformed_fasta = process_fasta(args.input)
    sys.stdout.write(transformed_fasta)
