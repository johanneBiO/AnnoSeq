##########################################################################################
# Johanne B. Overgaard
# Latest update: 13/3 2025

# remove_trembl: 
# Removes sequences from TrEMBL in Uniprot fasta file. 
##########################################################################################

import sys
import re

#---------------------------------------------------------------------------------------#

# Input file name
fsa_name = input("Enter the name of the UniProt FASTA file: ")

# Try to open input and output files
try:
    infile = open(fsa_name, "r")
    sp_fasta_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_sp.fasta", "w")
    sp_acc_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_sp_acc.txt", "w")
    tr_acc_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_tr_acc.txt", "w")
except IOError as error:
    print("Can't open file, reason:", str(error))
    sys.exit(1)

#---------------------------------------------------------------------------------------#

sp_count = 0         # Swiss-Prot count
tr_count = 0         # TrEMBL count
write_fasta = False  # Flag for writing Swiss-Prot sequences

# Process the FASTA file
for line in infile:
    if line.startswith(">"):                            # Header line
        parts = line.split("|")                         # Split by '|'
        if len(parts) > 2:
            accession = parts[1]                        # Extract accession number
            if line.startswith(">sp|"):                 # Swiss-Prot
                sp_count += 1
                sp_acc_outfile.write(accession + "\n")  # Write to accession file
                sp_fasta_outfile.write(line)            # Write header to Swiss-Prot FASTA
                write_fasta = True                      # Enable writing sequence
            elif line.startswith(">tr|"):               # TrEMBL
                tr_count += 1
                tr_acc_outfile.write(accession + "\n")  # Write to accession file
                write_fasta = False                     # Stop writing sequence
            else:
                write_fasta = False                     # Unrecognized format (unlikely)
    else:
        if write_fasta:                                 # Write Swiss-Prot sequence lines
            sp_fasta_outfile.write(line)

# Print summary
print("Number of Swiss-Prot sequences:", sp_count)
print("Number of TrEMBL accessions:", tr_count)

# Close files
infile.close()
sp_fasta_outfile.close()
sp_acc_outfile.close()
tr_acc_outfile.close()