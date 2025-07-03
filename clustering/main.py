##########################################################################################
# This script 
##########################################################################################

import gc
import os
import time
import h5py
import ujson
import pickle
import pyreadr
import argparse
import numpy as np
import pandas as pd
from joblib import load
from pathlib import Path
from scipy.ndimage import gaussian_filter1d 
from utils.get_region_embeddings import load_xgb_models, predict_and_smooth, extract_regions, process_data_with_models, summarize_region

#---------------------------------------------------------------------------------------#
# PREPARATION
#---------------------------------------------------------------------------------------#

# Time
start_time = time.time()

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-att", "--att_data", help = "Processed attention data.")
parser.add_argument("-mod", "--mod_dir", help = "Model folder.")
parser.add_argument("-emb", "--emb_file", help = "ESM-2 embedding file.")
parser.add_argument("-acc", "--acc_file", help = "Accession file.")
parser.add_argument("-res", "--res_dir", help = "Result folder.")
parser.add_argument("-a", "--anno_pre",  action = "store_true", help = "Use annotation file.")
parser.add_argument("-a_rds", "--anno_rds", help = "Name of RDS annotation data file.")
args = parser.parse_args()

attn_file = args.att_data
model_dir = args.mod_dir
emb_file = args.emb_file
acc_file = args.acc_file
res_dir = args.res_dir
predefined_anno = args.anno_pre
anno_file = args.anno_rds

# Create output directories
res_dir = os.path.join(res_dir)

#---------------------------------------------------------------------------------------#
# PREDICT REGIONS
#---------------------------------------------------------------------------------------#

# Read accession list
with open(acc_file, 'r') as f:
    acc = [line.strip() for line in f.readlines()]

if predefined_anno:
    anno = pyreadr.read_r(anno_file)
    anno = list(anno.values())[0]
    region_df = anno[anno['accession'].isin(acc)].reset_index(drop=True)
else:
    df = pd.read_pickle(attn_file)
    df = df.drop(columns=["category"])
    region_df = process_data_with_models(df, model_dir)

#---------------------------------------------------------------------------------------#
# EXTRACT AND SUMMARIZE EMBEDDINGS
#---------------------------------------------------------------------------------------#

# Create dict: accession -> list of (start, end) regions with annotation index
regions_by_accession = {}
for i, row in region_df.iterrows():
    acc_id = row['accession']
    if acc_id not in regions_by_accession:
        regions_by_accession[acc_id] = []
    regions_by_accession[acc_id].append((i, int(row['start_position']), int(row['end_position'])))

# Container for summary vectors and metadata rows
summary_vectors = []
meta_rows = []
index = 0

print("Processing batches...")

with open(emb_file, 'r') as f:
    for line_number, line in enumerate(f, 1):
        batch = ujson.loads(line)
        
        print(f"Processing batch {line_number}, size {len(batch)}")

        for seq in batch:
            seq_id = acc[index]
            if seq_id in regions_by_accession:
                emb = np.array(seq)
                for (anno_idx, start, end) in regions_by_accession[seq_id]:
                     summary_vec = summarize_region(emb, start, end)
                     summary_vectors.append(summary_vec)
                     meta_row = region_df.loc[anno_idx, ['accession', 'feature_type', 'start_position', 'end_position']].to_dict()
                     meta_rows.append(meta_row)
            
            index = index + 1 

        del batch, emb
        gc.collect()

# After processing all batches, convert summary_vectors and metadata to arrays/DataFrames
summary_matrix = np.vstack(summary_vectors)
meta_df = pd.DataFrame(meta_rows).reset_index(drop=True)

with h5py.File(os.path.join(res_dir, "esm_regions.h5"), "w") as f:
    # Save embedding matrix
    f.create_dataset("summary_matrix", data=summary_matrix)

    # Save metadata as datasets in a group
    meta_group = f.create_group("metadata")
    for col in meta_df.columns:
        meta_group.create_dataset(col, data=meta_df[col].astype(str).values)

# Time 
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script took {elapsed_time/60:.2f} minutes to run.")