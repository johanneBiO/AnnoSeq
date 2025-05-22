##########################################################################################
# This script partition the dataset including the subsets into training and test sets.
##########################################################################################

import os
import numpy as np
import pandas as pd

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

# Function for loading accessions 
def load_accessions(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())
    
# Function to save partition file
def save_partition_file(path, train_set, test_set):

    file_path = os.path.join(path, "partition.csv")

    data = (
        [(acc, "train") for acc in sorted(train_set)] +
        [(acc, "test") for acc in sorted(test_set)]
    )

    df = pd.DataFrame(data, columns=["accession", "partition"])
    df.to_csv(file_path, index=False)

complete_acc = load_accessions("../data/complete/additional/accessions_sp_cropped.txt")
subset_01000_acc = load_accessions("../data/subset_01000/additional/accessions.txt")
subset_00100_acc = load_accessions("../data/subset_00100/additional/accessions.txt")

# Split the complete set
# Remove any accession that's in a subset - these should not be included in the final test set
valid_complete_acc = complete_acc - subset_01000_acc - subset_00100_acc
valid_complete_acc = list(valid_complete_acc)

# Split (80/20)
n_complete = len(valid_complete_acc)
test_complete = set(np.random.choice(valid_complete_acc, size=int(0.2 * n_complete), replace=False))
train_complete = complete_acc - test_complete

# Split subset 1000
subset_01000_acc = list(subset_01000_acc)
n_01000 = len(subset_01000_acc)
test_01000 = set(np.random.choice(subset_01000_acc, size=int(0.2 * n_01000), replace=False))
train_01000 = set(subset_01000_acc) - test_01000

# Split subset 100
subset_00100_acc = list(subset_00100_acc)
n_00100 = len(subset_00100_acc)
test_00100 = set(np.random.choice(subset_00100_acc, size=int(0.2 * n_00100), replace=False))
train_00100 = set(subset_00100_acc) - test_00100

# Save partition CSVs
save_partition_file("../data/complete/additional/", train_complete, test_complete)
save_partition_file("../data/subset_01000/additional", train_01000, test_01000)
save_partition_file("../data/subset_00100/additional/", train_00100, test_00100)