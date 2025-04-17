##########################################################################################
# This script contains functions used to prepare sequences before running ESM-2.
##########################################################################################

import os
import random
import pandas as pd
from pathlib import Path

#---------------------------------------------------------------------------------------#
# Functions
#---------------------------------------------------------------------------------------#

def remove_trembl(input_path, output_path):
    """
    Removes sequences from TrEMBL in Uniprot fasta file. 

    Args:  
    - input_path: Input FASTA for sequences.
    - output_path: Resulting FASTA file path.

    Returns:
    - FASTA file without TrEMBL sequence.
    """
    
    infile  = Path(input_path)
    outfile = Path(output_path)
    
    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")

    sp_count = 0         # SWISS-PROT count
    tr_count = 0         # TrEMBL count
    write_fasta = False  # Flag for writing SWISS-PROT sequences

    print("Removing TrEMBL sequences...")
    
    with open(infile, "r") as infile, open(outfile, "w") as outfile:
        for line in infile:
            if line.startswith(">"):                            
                parts = line.split("|")                         
                if len(parts) > 2:
                    accession = parts[1]                        
                    if line.startswith(">sp|"):                 
                        sp_count += 1  
                        outfile.write(line)            
                        write_fasta = True                      
                    elif line.startswith(">tr|"):               
                        tr_count += 1 
                        write_fasta = False                     
                    else:
                        write_fasta = False                     
            else:
                if write_fasta:                                
                    outfile.write(line)

    infile.close()
    outfile.close()

    print("Number of Swiss-Prot sequences (retained):", sp_count)
    print("Number of TrEMBL sequences (removed):", tr_count)

def extract_accessions_from_fasta(input_path, output_path):
    """
    Extracts accession numbers from a FASTA file and writes them to a text file.
    
    Args:
        input_path: Path to the input FASTA file.
        output_path: Path to the output text file.

    Returns:
    - Text file with accession numbers.
    """

    infile = Path(input_path)
    outfile = Path(output_path)

    if not infile.is_file():
        print(f"The input file is invalid. Invalid file: {infile}")
        return
    
    accessions = []

    print("Extracting accessions...")

    with open(infile, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                parts = line.strip().split("|")
                if len(parts) >= 2:
                    accessions.append(parts[1])

    with open(outfile, "w") as outfile:
        for acc in accessions:
            outfile.write(acc + "\n")

    infile.close()
    outfile.close()

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

    infile  = Path(input_path)
    outfile = Path(output_path)
    
    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")

    sequences = []
    rm_count = 0         # SWISS-PROT count

    print("Cropping sequences...")

    with open(infile, "r") as infile:
        header, seq = None, []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    # Crop sequence if the accession is in length_dict
                    accession = header.split("|")[1].strip() 
                    if accession in length_dict:
                        if length_dict[accession] > 0:  
                            seq_cropped = crop_sequence("".join(seq), length_dict[accession])
                            sequences.append((header, seq_cropped))
                        else: 
                            rm_count += 1
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
    
    infile.close()

    # Write the cropped sequences to the output file
    with open(outfile, "w") as outfile:
        for header, seq in sequences:
            outfile.write(f"{header}\n{seq}\n")

    outfile.close()

    print("Number sequences removed due to specified length of 0:", rm_count)

def read_fasta(input_path):
    """
    Read in a FASTA file.

    Args: 
    - infile_path: FASTA file.
    """
    
    infile  = Path(input_path)
    
    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")
    
    sequences = []

    with open(infile, "r") as infile:
        header, seq = None, []
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
        if header:
            sequences.append((header, "".join(seq))) 

    infile.close()

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

    print("Generating subset...")

    sequences = read_fasta(input_path)
    random_sequences = random.sample(sequences, min(num_sequences, len(sequences)))

    outfile = Path(output_path)

    with open(outfile, "w") as outfile:
        for header, seq in random_sequences:
            outfile.write(f"{header}\n{seq}\n")

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
    
    infile = Path(infile_path)

    if not infile.is_file():
        print(f"The input file is invalid. Invalid file: {infile}")
        return
    
    seq = ""
    read_acc = False
    
    print("Scrambling sequences...")
    
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
    
    infile.close()
    outfile.close()