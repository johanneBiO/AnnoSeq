##########################################################################################
# This script contains functions used by main.py to extract ESM-2 embeddings and attention 
# scores. 
# Note: attention scores are aggregated using row means.
##########################################################################################

import os
import gc
import esm
import ujson
import torch
import random
import shutil
import numpy as np
from pathlib import Path
from torch.utils.data import Dataset

#---------------------------------------------------------------------------------------#
# Functions
#---------------------------------------------------------------------------------------#

def read_acc_seqs_from_fasta(infile_path):
    """
    Read in a FASTA file and prepare it to be used for the ESM-2 model.

    Args: 
    - infile_path: FASTA file.
    
    Returns: 
    - List of tuples. Contains accessions and sequences, e.g. [(acc, aTHNtem..)..()].
    """
    
    infile = Path(infile_path)
    
    if not infile.is_file():
        print(f"The input file was invalid. Invalid file was {infile}")

    accs = list()
    sequences = list()
    seq = ""
    read_acc = False

    infile = open(infile, "r")
    readfile = infile.readlines()
    infile.close()

    for line in readfile:
        line = line.strip()
        if line.startswith(">"):
            acc = line.split(">")[1]
            if read_acc:
                accs.append(acc)
                sequences.append(seq)
                seq = ""
            else:
                accs.append(acc)
        else:
            seq += line
            read_acc = True

    # Get last sequence
    sequences.append(seq)
    accs_and_sequences = tuple(zip(accs, sequences))
    return accs_and_sequences

class SequenceDataset(Dataset):
    """
    Custom dataset to handle sequences.
    """

    def __init__(self, sequences, batch_converter):
        self.sequences = sequences
        self.batch_converter = batch_converter

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]

def load_esm_model():
    """
    Load the ESM-2 model and return the model and alphabet.
    """

    # Look for GPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # Load ESM-2
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    model.to(device)
    batch_converter = alphabet.get_batch_converter()
    return model, alphabet, batch_converter

def collate_fn(batch, batch_converter):
    """
    Collate function to convert sequences into batch tokens.
    """

    batch_labels, batch_strs, batch_tokens = batch_converter(batch)
    return batch_labels, batch_strs, batch_tokens

def get_esm2_output(dataloader, model, alphabet, res_dir, layer_indices, norm = False, row = True, qua = False):
    """
    Compute ESM-2 embeddings and attention scores for sequences using batching.

    Args:
    - dataloader: PyTorch DataLoader.
    - model: ESM-2 model.
    - alphabet: ESM-2 alphabet.
    - output_dir: Directory to save JSON results.
    - norm: Save the normalized attention scores.
    
    Returns:
    - Embeddings and attention scores as JSON files for each bacth.
    """
    
    # Check for GPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Make result folder for batch results
    batch_res_dir = os.path.join(res_dir, "batch_results")
    os.makedirs(batch_res_dir, exist_ok=True)

    # Define the method for summarizing results (row-wise or column-wise)
    if row:
        sum_dim = 0
    else:
        sum_dim = 1

    # Set quantiles
    quantiles = 0.9 #torch.tensor([0, 0.25, 0.5, 0.75, 1]) 
        
    # Run ESM-2
    for batch_idx, (batch_labels, batch_strs, batch_tokens) in enumerate(dataloader):
        
        batch_tokens = batch_tokens.to(device)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])
        
        # Get embeddings for the final layer
        token_representations = results["representations"][33]

        # Store sequence embeddings after removing padding
        batch_embeddings = []

        for i, tokens_len in enumerate(batch_lens):
            batch_token = batch_tokens[i]
            if batch_token[0] == alphabet.cls_idx and (batch_token[-1] == alphabet.eos_idx or batch_token[-1] == alphabet.padding_idx):
                batch_embeddings.append(token_representations[i, 1 : tokens_len - 1].cpu().numpy().tolist())
            else:
                batch_embeddings.append(token_representations[i, :, :].cpu().numpy().tolist())

        # Save the embeddings of the current batch to a JSON file
        batch_embedding_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_embedding.json")
        with open(batch_embedding_file, "w") as outfile:
            ujson.dump(batch_embeddings, outfile)

        # Extract and store attention scores (raw)
        batch_attention_raw = {j: [] for j in layer_indices}

        for j in layer_indices:
            token_attention_raw = results["attn_raw"][j]
            for i, tokens_len in enumerate(batch_lens):
                batch_token = batch_tokens[i]
                if batch_token[0] == alphabet.cls_idx and (batch_token[-1] == alphabet.eos_idx or batch_token[-1] == alphabet.padding_idx):
                    batch_attention = token_attention_raw[i, 1:tokens_len - 1, 1:tokens_len - 1].cpu()
                else: 
                    batch_attention = token_attention_raw[i, :, :].cpu()
                
                # Summarize the attention signal
                if qua:
                    attention_features = torch.quantile(batch_attention, quantiles, dim=sum_dim).numpy().tolist()
                else:
                    attention_features = batch_attention.mean(dim=sum_dim).numpy().tolist()
                
                batch_attention_raw[j].append(attention_features)
                
        # Save attention scores for the current batch
        batch_attention_raw_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_attn_raw.json")
        with open(batch_attention_raw_file, "w") as outfile:
            ujson.dump(batch_attention_raw, outfile)
        
        # Extract and store normalized attention scores if specified
        if norm:
        
            batch_attention_norm = {j: [] for j in layer_indices}

            for j in layer_indices:
                token_attention_norm = results["attn_norm"][j]
                for i, tokens_len in enumerate(batch_lens):
                    batch_token = batch_tokens[i]
                    if batch_token[0] == alphabet.cls_idx and (batch_token[-1] == alphabet.eos_idx or batch_token[-1] == alphabet.padding_idx):
                        batch_attention = token_attention_norm[i, 1:tokens_len - 1, 1:tokens_len - 1].cpu()
                    else: 
                        batch_attention = token_attention_norm[i, :, :].cpu()
                    
                    # Summarize the attention signal
                    if qua:
                        attention_features = torch.quantile(batch_attention, quantiles, dim=sum_dim).numpy().tolist()
                    else:
                        attention_features = batch_attention.mean(dim=sum_dim).numpy().tolist()

                    batch_attention_norm[j].append(attention_features)
        
            batch_attention_norm_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_attn_norm.json")
            with open(batch_attention_norm_file, "w") as outfile:
                ujson.dump(batch_attention_norm, outfile)
            
            # Clean up
            del token_attention_norm, batch_attention_norm
            
        # Clean up and release memory
        del batch_tokens, results, token_representations, batch_embeddings, token_attention_raw, batch_attention, batch_attention_raw, attention_features
        gc.collect()
        torch.cuda.empty_cache()

def aggregate_batches(batch_res_dir, res_dir, layer_indices, norm = False):
    """
    Aggregate all batch embeddings and attention scores into separate files per layer
    in a memory-efficient way by writing intermediate results incrementally.

    Args:
    - batch_res_dir: Directory containing the batch-wise embeddings and attention JSON files.
    - res_dir: Directory where the aggregated results will be saved.
    - layer_indices: List of indices for layers, needed for attention aggregation.
    - norm: Include normalized attention scores.

    Returns:
    - Aggregated embeddings and attention scores as JSON files.
    """

    # Number of batches
    num_batches = len([f for f in os.listdir(batch_res_dir) if f.endswith("_embedding.json")])

    # Create output directories
    embedding_dir = os.path.join(res_dir, "embedding")
    attention_raw_dir = os.path.join(res_dir, "attention", "raw")
    os.makedirs(embedding_dir, exist_ok=True)
    os.makedirs(attention_raw_dir, exist_ok=True)

    if norm:
        attention_norm_dir = os.path.join(res_dir, "attention", "norm")
        os.makedirs(attention_norm_dir, exist_ok=True)

    # Open a file to write embeddings incrementally
    with open(os.path.join(embedding_dir, "embedding.json"), "w") as embedding_file:
        
        # Open separate files for each layer's attention scores
        attn_raw_files = {j: open(os.path.join(attention_raw_dir, f"layer_{j+1:02d}_attn_raw.json"), "w") for j in layer_indices}
        
        attn_norm_files = {}
        if norm:
            attn_norm_files = {j: open(os.path.join(attention_norm_dir, f"layer_{j+1:02d}_attn_norm.json"), "w") for j in layer_indices}

        for batch_idx in range(num_batches):
            batch_embedding_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_embedding.json")
            batch_attention_raw_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_attn_raw.json")
            batch_attention_norm_file = os.path.join(batch_res_dir, f"batch_{batch_idx+1}_attn_norm.json")

            # Process embeddings
            if os.path.isfile(batch_embedding_file):
                with open(batch_embedding_file, "r") as infile:
                    batch_embeddings = ujson.load(infile)
                ujson.dump(batch_embeddings, embedding_file)
                # Separate batches with a newline
                embedding_file.write("\n")

            # Process raw attention scores
            if os.path.isfile(batch_attention_raw_file):
                with open(batch_attention_raw_file, "r") as infile:
                    batch_attention_raw = ujson.load(infile)
                for j in layer_indices:
                    j_key = str(j)
                    ujson.dump(batch_attention_raw[j_key], attn_raw_files[j])
                    # Separate batches with a newline
                    attn_raw_files[j].write("\n")

            # Process norm attention scores (only if norm = True and file exists)
            if norm and os.path.isfile(batch_attention_norm_file):
                with open(batch_attention_norm_file, "r") as infile:
                    batch_attention_norm = ujson.load(infile)
                for j in layer_indices:
                    j_key = str(j)
                    ujson.dump(batch_attention_norm[j_key], attn_norm_files[j])
                    attn_norm_files[j].write("\n")

        # Close all attention score files
        for f in attn_raw_files.values():
            f.close()
        for f in attn_norm_files.values():
            f.close()

    # Clean up old batch files after aggregation
    if os.path.exists(batch_res_dir):
        shutil.rmtree(batch_res_dir) 