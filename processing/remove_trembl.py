##########################################################################################
# Removes sequences from TrEMBL in Uniprot fasta file. 
##########################################################################################

import re
import sys

#---------------------------------------------------------------------------------------#
# Input file name
#---------------------------------------------------------------------------------------#

fsa_name = input("Enter the name of the UniProt FASTA file: ")

try:
    infile = open(fsa_name, "r")
    sp_fasta_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_sp.fasta", "w")
    sp_acc_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_sp_acc.txt", "w")
    tr_acc_outfile = open(f"{fsa_name.rsplit('.', 1)[0]}_tr_acc.txt", "w")
except IOError as error:
    print("Can't open file, reason:", str(error))
    sys.exit(1)

#---------------------------------------------------------------------------------------#
# Process FASTA file
#---------------------------------------------------------------------------------------#

sp_count = 0         # Swiss-Prot count
tr_count = 0         # TrEMBL count
write_fasta = False  # Flag for writing Swiss-Prot sequences

for line in infile:
    if line.startswith(">"):                            
        parts = line.split("|")                         
        if len(parts) > 2:
            accession = parts[1]                        
            if line.startswith(">sp|"):                 
                sp_count += 1
                sp_acc_outfile.write(accession + "\n")  
                sp_fasta_outfile.write(line)            
                write_fasta = True                      
            elif line.startswith(">tr|"):               
                tr_count += 1
                tr_acc_outfile.write(accession + "\n")  
                write_fasta = False                     
            else:
                write_fasta = False                     
    else:
        if write_fasta:                                
            sp_fasta_outfile.write(line)

print("Number of Swiss-Prot sequences:", sp_count)
print("Number of TrEMBL accessions:", tr_count)

infile.close()
sp_fasta_outfile.close()
sp_acc_outfile.close()
tr_acc_outfile.close()