##########################################################################################
# Removes sequences from TrEMBL in the raw fasta file. 
##########################################################################################

import os
from utils.seq_processors import remove_trembl, extract_accessions_from_fasta

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

# Define file names
fasta_file = "../data/_raw/UP000005640_9606.fasta"
res_dir = "../data/complete"
fasta_sp = os.path.join(res_dir, "sequences/seq_complete_sp.fasta")
txt_acc = os.path.join(res_dir, "additional/accessions_sp.txt")

# Remove TrEMBL entries
remove_trembl(fasta_file, fasta_sp)

# Get the accessions
extract_accessions_from_fasta(fasta_sp, txt_acc)