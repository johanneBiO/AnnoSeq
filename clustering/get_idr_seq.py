##########################################################################################
# This script cut out IDR regions for specific clusters and computes the amino acid 
# composition for each cluster. 
##########################################################################################

import os
import pandas as pd
from collections import defaultdict

def load_fasta_as_dict(fasta_path):
    
    fasta_dict = {}
    
    with open(fasta_path, "r") as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    fasta_dict[header] = "".join(seq_lines)
                header = line.split("|")[1].strip() 
                seq_lines = []
            else:
                seq_lines.append(line)
        # Last record
        if header:
            fasta_dict[header] = "".join(seq_lines)

    return fasta_dict

def compute_aa_composition(sequences):

    aa_counts = defaultdict(int)
    total = 0
    
    for seq in sequences:
        for aa in seq:
            aa_counts[aa] += 1
            total += 1

    # Convert to frequencies
    return {aa: count / total for aa, count in aa_counts.items()}

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

idr_data = pd.read_csv("results/idr_data.csv")

fasta_dict = load_fasta_as_dict("../data/subset_01000/sequences/biological_seq/seq_01000.fasta")

cluster_sequences = defaultdict(list)

for _, row in idr_data.iterrows():

    accession = row["accession"]
    start = int(row["start_position"])
    end = int(row["end_position"])
    cluster = row["cluster_number"]
    
    seq = fasta_dict.get(accession)
    
    if seq:
        sub_seq = seq[start:end]
        cluster_sequences[cluster].append(sub_seq)

cluster_compositions = {cluster: compute_aa_composition(seqs) 
                        for cluster, seqs in cluster_sequences.items()}

cluster_composition_df = pd.DataFrame.from_dict(cluster_compositions, orient='index').fillna(0)

cluster_composition_df.to_csv("results/idr_cluster_aa_composition.csv", index=True)