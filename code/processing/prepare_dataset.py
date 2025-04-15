##########################################################################################
# This script crops sequences to fit ESM-2 input limits, creates subsets of sequences from 
# complete FASTA file, and generates scrambled variants.
##########################################################################################

import os
import random
import pandas as pd
from pathlib import Path

#---------------------------------------------------------------------------------------#
# Functions
#---------------------------------------------------------------------------------------#

def crop_sequence(seq, seq_length):
    """
    Crop sequences to specified length by removing positions in the tail.

    Args:
    - seq: Sequence.
    - seq_length: Maximum length of the sequence.

    Returns:
    - Cropped sequence.
    """

    seq_cropped = seq[0:seq_length]
    return seq_cropped

def crop_fasta(input_path, output_path, length_dict, esm_max = 1024):
    """
    Crop sequences in FASTA file to specified length by removing positions in the tail.
    The sequences are cropped from specified length if it exists in the specified dictionary. 
    Otherwise, it will be cropped with respect to the upper limit for ESM-2 (default: 1024).

    Args:  
    - input_path: Input FASTA for sequences.
    - output_path: Resulting FASTA file path.
    - length_dict: Dict with the maximum length for each sequence.
    - esm_max: Upper limit for ESM-2.

    Returns:
    - FASTA file with cropped sequence.
    """

    sequences = []
    with open(input_path, "r") as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    # Crop sequence if the accession is in length_dict
                    accession = header.split("|")[1].strip() 
                    if accession in length_dict:  
                        seq_cropped = crop_sequence("".join(seq), length_dict[accession])
                        sequences.append((header, seq_cropped))
                     # Crop the sequence according to maximal esm input length
                    else:
                        sequences.append((header, "".join(seq)[:esm_max]))
                header, seq = line, []
            else:
                seq.append(line)
        # Process last sequence
        if header:
            if header[1:] in length_dict:
                seq_cropped = crop_sequence("".join(seq), length_dict[accession])
                sequences.append((header, seq_cropped))
            else:
                sequences.append((header, "".join(seq)[:esm_max])) 

    # Write the cropped sequences to the output file
    with open(output_path, "w") as out_f:
        for header, seq in sequences:
            out_f.write(f"{header}\n{seq}\n")

def read_fasta(file_path):
    """
    Read in a FASTA file.

    Args: 
    - infile_path: FASTA file.
    """

    sequences = []
    with open(file_path, "r") as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
        if header:
            sequences.append((header, "".join(seq))) 
    return sequences

def extract_random_sample(input_path, output_path, num_sequences=100):
    """
    Extract random sequences to generate subset.

    Args: 
    - Infile_path: Source FASTA file.
    - Output_path: Subset FASTA file.
    - num_sequences: Number of sequence to sample.

    Returns:
    - A subset FASTA file.
    """

    sequences = read_fasta(input_path)
    random_sequences = random.sample(sequences, min(num_sequences, len(sequences)))

    with open(output_path, "w") as out_f:
        for header, seq in random_sequences:
            out_f.write(f"{header}\n{seq}\n")

def scramble_sequence(sequence):
    """
    Scramble a sequence.

    Args: 
    - sequence: sequence of interest.

    Returns:
    - Scrambled variant of the input sequence. 
    """

    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

def scramble_fasta(infile_path, outfile_dir, outfile_name):
    """
    Scramble all sequences in FASTA file.

     Args: 
    - infile_path: Input FASTA for sequences.
    - outfile_dir: Name of directory for scrambled FASTA file.
    - outfile_name: Name for scrambled FASTA file.

    Returns:
    - Scrambled FASTA file.
    """

    seq = ""
    read_acc = False
    infile = Path(infile_path)

    if not infile.is_file():
        print(f"The input file is invalid. Invalid file: {infile}")
        return
    
    os.makedirs(outfile_dir, exist_ok=True)
    outfile = os.path.join(outfile_dir, outfile_name)

    # Open the input file for reading and the output file for writing
    with open(infile, "r") as infile, open(outfile, "w") as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                # If there's an ongoing sequence, scramble it and write to the output file
                if read_acc:
                    scrambled_seq = scramble_sequence(seq)
                    outfile.write(f">{acc}\n{scrambled_seq}\n")
                    seq = ""    # Reset the sequence after writing it
                # Read the accession number (header)
                acc = line[1:]  # Remove the '>' from the header line
                read_acc = True
            else:
                # Append the sequence data
                seq += line
        
        # Scramble and write the last sequence after the loop
        if seq:
            scrambled_seq = scramble_sequence(seq)
            outfile.write(f">{acc}\n{scrambled_seq}\n")


def extract_accessions_from_fasta(input_path, output_path):
    """
    Extracts accession numbers from a FASTA file and writes them to a text file.
    
    Args:
        input_path: Path to the input FASTA file.
        output_path: Path to the output text file.

    Returns:
    - Text file with accession numbers.
    """
    
    accessions = []
    with input_path.open() as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().split("|")
                if len(parts) >= 2:
                    accessions.append(parts[1])

    with output_path.open("w") as out:
        for acc in accessions:
            out.write(acc + "\n")

#---------------------------------------------------------------------------------------#
# Generate datasets
#---------------------------------------------------------------------------------------#

# Read adjusted sequence lengths
csv_df = pd.read_csv("../../data/complete/additional/seq_lengths_adj.csv")

# Create a dictionary with accession as the key and adjusted length as the value
length_dict = dict(zip(csv_df['accession'], csv_df['length_adj']))

# Specify the file name of the fasta and resulting cropped file
bio_fasta = Path("../../data/complete/sequences/biological_seq/UP000005640_9606_sp.fasta") 
bio_cropped_fasta = Path("../../data/complete/sequences/biological_seq/UP000005640_9606_sp_cropped.fasta")

# Crop sequences
crop_fasta(bio_fasta, bio_cropped_fasta, length_dict, esm_max = 1024)

# Specify the filename of subsets
subset_01000_cropped_fasta = Path("../../data/subset_01000/sequences/biological_seq/seq_01000.fasta")
subset_00100_cropped_fasta = Path("../../data/subset_00100/sequences/biological_seq/seq_00100.fasta")

# Sample random sequences
extract_random_sample(bio_cropped_fasta, subset_01000_cropped_fasta, num_sequences=1000)
extract_random_sample(subset_01000_cropped_fasta, subset_00100_cropped_fasta, num_sequences=100)

# Extract the accesions for the sequences across subsets.
extract_accessions_from_fasta(subset_01000_cropped_fasta, Path("../../data/subset_01000/additional/accessions.txt"))
extract_accessions_from_fasta(subset_00100_cropped_fasta, Path("../../data/subset_00100/additional/accessions.txt"))

# Specify directories for scrambled sequences (- only for subsets)
subset_01000_scramble_seq_dir = Path("../../data/subset_01000/sequences/scrambled_seq/")
subset_00100_scramble_seq_dir = Path("../../data/subset_00100/sequences/scrambled_seq/")

# Scramble sequences
for i in range(1, 6):
    scramble_fasta(subset_01000_cropped_fasta, subset_01000_scramble_seq_dir, f"scramble_{i:02d}.fasta")
    scramble_fasta(subset_00100_cropped_fasta, subset_00100_scramble_seq_dir, f"scramble_{i:02d}.fasta")