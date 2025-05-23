##########################################################################################
# This script runs a modified version of ESM-2 to extract attention scores. Both raw and 
# normalized attention scores are aggregated using row means.
# Optionally, the model can be run on pre-generated scrambled sequences, which should be 
# stored in a shared folder.
##########################################################################################

import os
import sys
import time
import torch
import argparse
from torch.utils.data import DataLoader
from utils.esm_processors import SequenceDataset, load_esm_model, collate_fn, get_esm2_output, read_acc_seqs_from_fasta, aggregate_batches

#---------------------------------------------------------------------------------------#
# PREPARATION
#---------------------------------------------------------------------------------------#

# Time
start_time = time.time()

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-seq", "--seq_input", help = "Input sequence fasta file.")
parser.add_argument("-res", "--res_dir", help = "Result folder.")
parser.add_argument("-s", "--scr_run", action = "store_true", help = "Run ESM for scrambled sequences.")
parser.add_argument("-scr", "--scr_dir", help = "Scrambled sequence folder.")
parser.add_argument("-n", "--norm",  action = "store_true", help = "Extract normalized attention scores.")
parser.add_argument("-r", "--row",  action = "store_true", help = "Row-summarize the attention scores.")
parser.add_argument("-q", "--qua",  action = "store_true", help = "Quantiles of the attention scores.")
parser.add_argument("-bs", "--batch_size", type=int, default=20, help = "Batch size for processing.")
args = parser.parse_args()

bio_fasta_file = args.seq_input
res_dir = args.res_dir
scr = args.scr_run
scr_dir = args.scr_dir
norm_attn = args.norm
row_attn = args.row
qua_attn = args.qua
batch_size = args.batch_size

# Create output directories
if scr:
    bio_res_dir = os.path.join(res_dir, "biological_seq")
    scr_res_dir = os.path.join(res_dir, "scrambled_seq")
else:
    bio_res_dir = os.path.join(res_dir)

# Use customized ESM-2 model
sys.path.insert(0, '../esm/')
model, alphabet, batch_converter = load_esm_model()

# Layer indices
layer_indices = list(range(33))

#---------------------------------------------------------------------------------------#
# BIOLOGICAL SEQUENCES
#---------------------------------------------------------------------------------------#

# Prepare data
sequences = read_acc_seqs_from_fasta(bio_fasta_file)  
dataset = SequenceDataset(sequences, batch_converter)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False, collate_fn=lambda b: collate_fn(b, batch_converter))

# Run ESM-2
get_esm2_output(dataloader, model, alphabet, bio_res_dir, layer_indices, norm_attn, row_attn, qua_attn)

# Aggregate results
batch_res_dir = os.path.join(bio_res_dir, "batch_results")
aggregate_batches(batch_res_dir, bio_res_dir, layer_indices, norm_attn)

#---------------------------------------------------------------------------------------#
# SCRAMBLED SEQUENCES
#---------------------------------------------------------------------------------------#

if scr:

    file_count = sum(os.path.isfile(os.path.join(scr_dir, f)) for f in os.listdir(scr_dir))
    file_names = [f for f in os.listdir(scr_dir) if os.path.isfile(os.path.join(scr_dir, f))]

    print(f"Number of fasta files in '{scr_dir}': {file_count}")

    for i in range(file_count):
        
        # Specify the scrambled sequences 
        scr_fasta_file = os.path.join(scr_dir, file_names[i])
        res_dir_i = os.path.join(scr_res_dir, f"scramble_{i+1:02d}")
        
        # Prepare data
        sequences = read_acc_seqs_from_fasta(scr_fasta_file)  
        dataset = SequenceDataset(sequences, batch_converter)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False, collate_fn=lambda b: collate_fn(b, batch_converter)) 

        # Run ESM-2
        get_esm2_output(dataloader, model, alphabet, res_dir_i, layer_indices, norm_attn, row_attn, qua_attn)

        # Aggregate results
        batch_res_dir = os.path.join(res_dir_i, "batch_results")
        aggregate_batches(batch_res_dir, res_dir_i, layer_indices, norm_attn)

# Time 
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script took {elapsed_time/60:.2f} minutes to run.")