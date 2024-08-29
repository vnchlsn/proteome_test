import os
import csv
import numpy as np
from localcider.sequenceParameters import SequenceParameters

# Function to calculate hydrophobicity
def calculate_hydrophobicity(sequence):
    hydropathy_values = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    return sum(hydropathy_values.get(aa, 0) for aa in sequence) / len(sequence)

# Function to read a FASTA file
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        name = None
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if name:
                    sequences[name] = sequence.replace('*', '')
                name = line[1:].strip()
                sequence = ''
            else:
                sequence += line.strip()
        if name:
            sequences[name] = sequence.replace('*', '')
    return sequences

# Main function
def main(input_fasta, output_csv):
    sequences = read_fasta(input_fasta)
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Protein Name', 'Kappa', 'Length', 'Hydrophobicity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for name, sequence in sequences.items():
            try:
                seq_obj = SequenceParameters(sequence=sequence)
                kappa = seq_obj.get_kappa()
                length = len(sequence)
                hydrophobicity = calculate_hydrophobicity(sequence)
                writer.writerow({
                    'Protein Name': name,
                    'Kappa': kappa,
                    'Length': length,
                    'Hydrophobicity': hydrophobicity
                })
            except Exception as e:
                print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    input_fasta = 'protein - Copy.fasta'  # Change this to your input file
    output_csv = 'Marc_prefiltered_length_out.csv'    # Change this to your desired output file
    main(input_fasta, output_csv)
    