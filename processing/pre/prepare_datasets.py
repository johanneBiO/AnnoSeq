##########################################################################################
# This script crops sequences to fit ESM-2 input limits, creates subsets of sequences from 
# complete FASTA file, and generates scrambled variants.
##########################################################################################

import os
from utils.seq_processors import *

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

# Initialize the random number generator
random.seed(100)

# Read adjusted sequence lengths
csv_df = pd.read_csv("../data/complete/additional/seq_lengths_adj.csv")

# Create a dictionary with accession as the key and adjusted length as the value
length_dict = dict(zip(csv_df['accession'], csv_df['length_adj']))

# Specify the file name of the fasta and resulting cropped file
bio_fasta = "../data/complete/sequences/seq_complete_sp.fasta"
bio_cropped_fasta = "../data/complete/sequences/seq_complete_sp_cropped.fasta"

# Crop sequences
crop_fasta(bio_fasta, bio_cropped_fasta, length_dict, esm_max = 1024)

# Extract the accesion for cropped sequences
extract_accessions_from_fasta(bio_cropped_fasta, "../data/complete/additional/accessions_sp_cropped.txt")

# Specify the filename of subsets
os.makedirs("../data/subset_01000/sequences/biological_seq", exist_ok=True)
os.makedirs("../data/subset_00100/sequences/biological_seq", exist_ok=True)
subset_01000_cropped_fasta = "../data/subset_01000/sequences/biological_seq/seq_01000.fasta"
subset_00100_cropped_fasta = "../data/subset_00100/sequences/biological_seq/seq_00100.fasta"

# Sample random sequences
extract_random_sample(bio_cropped_fasta, subset_01000_cropped_fasta, num_sequences=1000)
extract_random_sample(subset_01000_cropped_fasta, subset_00100_cropped_fasta, num_sequences=100)

# Extract the accesions for the sequences across subsets.
extract_accessions_from_fasta(subset_01000_cropped_fasta, "../data/subset_01000/additional/accessions.txt")
extract_accessions_from_fasta(subset_00100_cropped_fasta, "../data/subset_00100/additional/accessions.txt")

# Specify directories for scrambled sequences (- only for subsets)
subset_01000_scramble_seq_dir = "../data/subset_01000/sequences/scrambled_seq/"
subset_00100_scramble_seq_dir = "../data/subset_00100/sequences/scrambled_seq/"
os.makedirs(subset_01000_scramble_seq_dir, exist_ok=True)
os.makedirs(subset_00100_scramble_seq_dir, exist_ok=True)

# Scramble sequences
for i in range(1, 6):
    scramble_fasta(subset_01000_cropped_fasta, subset_01000_scramble_seq_dir, f"scramble_{i:02d}.fasta")
    scramble_fasta(subset_00100_cropped_fasta, subset_00100_scramble_seq_dir, f"scramble_{i:02d}.fasta")