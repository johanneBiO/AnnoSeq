##########################################################################################
# This script process the attention scores after running ESM-2 to prepare the data to be 
# used as input to a classification model.
##########################################################################################

import os
import csv
import time
import ujson
import pyreadr
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from collections.abc import Iterable
from scipy.ndimage import gaussian_filter1d

#---------------------------------------------------------------------------------------#
# PREPARATION
#---------------------------------------------------------------------------------------#

# Time
start_time = time.time()

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-at", "--attn", help = "Folder of attention scores across layers.")
parser.add_argument("-ac", "--acc", help = "File for accesions.")
parser.add_argument("-an", "--anno", help = "RDS object with annotation.")
parser.add_argument("-pa", "--part", help = "CSV describing the data set partition.")
parser.add_argument("-r", "--res", help = "Resulting processed file.")
parser.add_argument("-q", "--qua",  action = "store_true", help = "Quantiles of the attention scores.")
parser.add_argument("-csv", "--csv_file",  action = "store_true", help = "Save CSV data table.")
parser.add_argument("-r_csv", "--res_csv", help = "Name of CSV data file.")
args = parser.parse_args()

attn_path = args.attn
acc_path = args.acc
anno_path = args.anno
part_path = args.part
res_path = args.res
qua_attn = args.qua
csv_file = args.csv_file
res_csv = args.res_csv

#---------------------------------------------------------------------------------------#
# ACCESSIONS
#---------------------------------------------------------------------------------------#

# Load sequence accessions
with open(acc_path, 'r') as f:
    acc = [line.strip() for line in f.readlines()]

#---------------------------------------------------------------------------------------#
# ATTENTION
#---------------------------------------------------------------------------------------#

# Combine the attention scores from all layers
combined_data = []

for layer_idx in range(1, 34):
    file_path = os.path.join(attn_path, f"layer_{str(layer_idx).zfill(2)}_attn_raw.json")

    # Each line is a batch
    with open(file_path, 'r') as f:
        data = [ujson.loads(line) for line in f]  

    # Flatten batches. Shape: [layer][sequence][5 features][position]
    layer_data = [sequence for batch in data for sequence in batch]
    combined_data.append(layer_data)

# Initialize an empty list to store the rows of the table
table_data = []

# Loop through each layer
if qua_attn:
    for layer_idx, layer in enumerate(combined_data):
        # Loop through each sequence within the layer
        for seq_idx, sequence in enumerate(layer):
            # Loop through each feature in the sequence
            for feature_idx, feature in enumerate(sequence):
                for pos_idx, value in enumerate(feature):
                        row = [acc[seq_idx], layer_idx + 1, feature_idx + 1, pos_idx + 1, value]
                        table_data.append(row)
else:
    for layer_idx, layer in enumerate(combined_data):
    # Loop through each sequence within the layer
        for seq_idx, sequence in enumerate(layer):
            for pos_idx, value in enumerate(sequence):
                row = [acc[seq_idx], layer_idx + 1, 1, pos_idx + 1, value]
                table_data.append(row)
        
# Make the data into a panda dataframe
df = pd.DataFrame(table_data, columns=["accession", "layer", "feature", "position", "value"])

# Create 'layer_feature' column
df["layer_feature"] = df["layer"].map(lambda x: f"layer_{x:02}_feature_") + df["feature"].map(lambda x: f"{x:02}")

# Pivot to wide format
df = df.pivot(index=["accession", "position"], columns="layer_feature", values="value").reset_index()

# Smoothing
def smooth_group(group):
    group = group.sort_values('position')
    for col in layer_columns:
        group[col] = gaussian_filter1d(group[col].values, sigma=1)
    return group

layer_columns = [col for col in df.columns if col.startswith("layer")]

df = df.groupby('accession', group_keys=False).apply(smooth_group)

#---------------------------------------------------------------------------------------#
# ANNOTATIONS
#---------------------------------------------------------------------------------------#

# Read the annotations
anno = pyreadr.read_r(anno_path)

# Assuming only one object in the RDS file
anno_df = list(anno.values())[0]

# Merge the annotations with the attention scores
merged_df = pd.merge(
    anno_df,
    df,
    on=["accession", "position"],
    how="right",  
    sort=False
)

# Fill NAs with unannotated
merged_df[['category', 'feature_type']] = merged_df[['category', 'feature_type']].fillna("unannotated")

# Modify category and feature type column
merged_df['category'] = merged_df['category'].str.replace(' ', '_')
merged_df['feature_type'] = merged_df['feature_type'].str.replace(' ', '_')

#---------------------------------------------------------------------------------------#
# PARTITION
#---------------------------------------------------------------------------------------#

# Get the partition info
part_df = pd.read_csv(part_path)

# Merge the partition info with the table
merged_df = pd.merge(
    part_df,
    merged_df,
    on="accession",
    how="right", 
    sort=False
)

#---------------------------------------------------------------------------------------#
# EXPORT
#---------------------------------------------------------------------------------------#

# Save the final results as pickle files
merged_df.to_pickle(res_path)

# Save as CSV as well if specified
if csv_file:
    merged_df.to_csv(res_csv, index = False)

# Time 
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script took {elapsed_time/60:.2f} minutes to run.")