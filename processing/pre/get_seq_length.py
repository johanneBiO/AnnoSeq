##########################################################################################
# This script outputs csv files descriping the sequence length for cropped sequences.
##########################################################################################

import os
import sys
sys.path.append(os.path.abspath(".."))

from utils.seq_processors import *

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

# Specify the file name of the fasta and the cropped file
fasta = "../../data/complete/sequences/seq_complete_sp.fasta"
fasta_cropped= "../../data/complete/sequences/seq_complete_sp_cropped.fasta"

# Specify the output file names
out_csv = "../../data/complete/additional/seq_length_sp_all.csv"
out_cropped_csv = "../../data/complete/additional/seq_length_sp_cropped.csv"

# Compute the sequence length
fasta_to_length_csv(fasta, out_csv)
fasta_to_length_csv(fasta_cropped, out_cropped_csv)