##########################################################################################
# This script trains XGBoost models for each annotation type present in the data.
# The final model is applied on a test set to evaluate performance. 
##########################################################################################

import os
import time
import joblib
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter
from xgboost import XGBClassifier
from joblib import Parallel, delayed
from scipy.ndimage import gaussian_filter1d
from imblearn.over_sampling import SMOTENC
from imblearn.under_sampling import RandomUnderSampler
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV, GroupKFold

#---------------------------------------------------------------------------------------#
# SET UP
#---------------------------------------------------------------------------------------#

# Time
start_time = time.time()

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file_input", help = "Input CSV file.")
parser.add_argument("-res", "--res_dir", help = "Result folder.")
parser.add_argument("-s", "--small",  action = "store_true", help = "Run on subset 100.")
parser.add_argument("-c", "--comp",  action = "store_true", help = "Complete grid search.")
parser.add_argument("-n", "--n_size", type=int, default=10, help = "Number of parameter combinations to try in random search.")
args = parser.parse_args()

file_name = args.file_input
res_dir = args.res_dir
small = args.small
comp = args.comp
n_size = args.n_size

# Check if file exists
if not Path(file_name).is_file():
    print(f"The input file was invalid.")

#---------------------------------------------------------------------------------------#
# MODEL TRAINING
#---------------------------------------------------------------------------------------#

# Hyperparameters
param_grid = {
    'n_estimators': [200, 400, 600],
    'max_depth': [4, 6, 8],
    'learning_rate': [0.01, 0.1, 0.2],
    'colsample_bytree': [0.2, 0.5, 0.8],
    "subsample": [0.4, 0.6, 0.8],
    'lambda': [0.1, 1, 5],
    'alpha': [0.1, 1, 5]
}

# Function for handling class imbalance
def handle_class_imbalance(X_train, y_train, target_class, control_class):
    class_counts = Counter(y_train["feature_type"])
    majority_count = class_counts[control_class]
    minority_count = class_counts[target_class]
    
    # If the ratio is greater than or equal to 10:1, apply undersampling
    if majority_count / minority_count >= 10:
        undersampler = RandomUnderSampler(sampling_strategy={control_class: minority_count * 10}, random_state=100)
        oversampler = SMOTENC(sampling_strategy=0.2, random_state=100, categorical_features=[0])
        X_resampled, y_resampled = undersampler.fit_resample(X_train, y_train)
        X_resampled, y_resampled = oversampler.fit_resample(X_resampled, y_resampled)
        return X_resampled, y_resampled, True
    else:
        return X_train, y_train, False

def train_model_for_class(target_class, X_train, y_train, X_test, y_test, param_grid):
    # Skip the control class
    if target_class == control_class:
        return None
    
    # Subset train and test to only keep rows of the target class or the control class
    train_mask = (y_train["feature_type"] == target_class) | (y_train["feature_type"] == control_class)
    test_mask = (y_test["feature_type"] == target_class) | (y_test["feature_type"] == control_class)
    
    X_train_subset = X_train[train_mask]
    y_train_subset = y_train[train_mask]
    y_train_subset = y_train_subset.drop(columns=["accession", "position"])
    
    X_test_subset = X_test[test_mask]
    y_test_subset = y_test[test_mask]
    
    # Apply class imbalance handling
    X_train_final, y_train_final, undersampling_applied = handle_class_imbalance(X_train_subset, y_train_subset, target_class, control_class)

    # Specify the group (accession)
    accessions_train = X_train_final["accession"].values

    X_train_final = X_train_final.drop(columns=["accession"])
    
    # Create binary labels: 1 for the target class, 0 for the control class
    y_train_binary = (y_train_final == target_class).astype(int)
    pos = (y_train_binary["feature_type"] == 1).sum()
    neg = (y_train_binary["feature_type"] == 0).sum()
    pos_weight = neg / pos

    cv_with_groups = GroupKFold(n_splits=3)

    xgb_clf = XGBClassifier(
        scale_pos_weight=pos_weight,
        objective='binary:logistic',
        n_jobs=2, #2
        random_state=42
    )
    
    if comp:
        grid_search = GridSearchCV(
            estimator=xgb_clf,
            param_grid=param_grid,
            cv=cv_with_groups,
            n_jobs=5, 
            verbose=2,
            scoring='f1_weighted'
            )
    else:
        grid_search = RandomizedSearchCV(
            estimator=xgb_clf,
            param_distributions=param_grid,
            n_iter=n_size,
            cv=cv_with_groups,
            n_jobs=5, 
            verbose=2,
            scoring='f1_weighted',
            random_state=42
        )
    
    # Train the Random Forest with RandomizedSearchCV
    grid_search.fit(X_train_final, y_train_binary, groups=accessions_train)
    
    accessions_test = y_test_subset["accession"].values
    positions_test = y_test_subset["position"].values

    X_test_subset = X_test_subset.drop(columns=["accession"])
    y_test_subset = y_test_subset["feature_type"]
    y_test_binary = (y_test_subset == target_class).astype(int)

    # Get predictions on the test set
    y_pred = grid_search.best_estimator_.predict(X_test_subset)

    # Get the probabilities for the positive class (target class)
    y_pred_proba = grid_search.best_estimator_.predict_proba(X_test_subset)[:, 1]
    
    # Save the predictions
    predictions = pd.DataFrame({
        'Accession': accessions_test,
        'Position': positions_test,
        'True': y_test_binary, 
        'Prediction': y_pred,
        'Probability': y_pred_proba  # Add predicted probabilities
    })

    predictions.to_csv(os.path.join(pred_dir, f'pred_{target_class}.csv'), index=False)
    
    # Save the best model for the current class
    joblib.dump(grid_search.best_estimator_, os.path.join(model_dir, f'model_{target_class}.pkl'))

    return {
        'target_class': target_class,
        'best_params': grid_search.best_params_,
        'best_score': grid_search.best_score_,
        'undersampling': undersampling_applied
    }

#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#

# Prepare data
#---------------------------------------------------------------------------------------#

# Read the data
df = pd.read_pickle(file_name)
df = df.drop(columns=["category"])

df = df.sample(frac=1, random_state=42).reset_index(drop=True)

# Make result directories
pred_dir = os.path.join(res_dir, "predictions")
sum_dir = os.path.join(res_dir, "summary")
model_dir = os.path.join(res_dir, "models")
os.makedirs(pred_dir, exist_ok=True)
os.makedirs(sum_dir, exist_ok=True)
os.makedirs(model_dir, exist_ok=True)

# Split the dataset into train and test sets (you mentioned you have a column for this)
train_df = df[df["partition"] == "train"]
test_df = df[df['partition'] == 'test']

X_train = train_df.drop(columns=["position", "feature_type", "partition"])
y_train = train_df[["accession", "position", "feature_type"]]

X_test = test_df.drop(columns=["position", "feature_type", "partition"])
y_test = test_df[["accession", "position", "feature_type"]]

# Set the control group to be the unannotated group
control_class = 'unannotated'

# Define all the unique feature types
if small:
    classes = ['unannotated', 'domain', 'region_of_interest', 'transmembrane_region', 'topological_domain', 'disulfide_bond', 'zinc_finger_region', 'binding_site']
else:
    classes = y_train["feature_type"].unique()

# Run model training
#---------------------------------------------------------------------------------------#

# Initialize a dictionary to store results
results = {}

"""
# Train the different models in paralllel
all_results = Parallel(n_jobs=10, verbose=10)(
    delayed(train_model_for_class)(target_class, X_train, y_train, X_test, y_test, param_grid)
    for target_class in classes 
)
"""

all_results = []
for target_class in classes:
    result = train_model_for_class(target_class, X_train, y_train, X_test, y_test, param_grid)
    if result is not None:
        all_results.append(result)

# Organize the results
results = {res['target_class']: res for res in all_results if res is not None}

results_summary = pd.DataFrame.from_dict({
    class_id: {'Best Params': res['best_params'], 'Best Score': res['best_score'], 'Undersampling': res['undersampling']}
    for class_id, res in results.items()
}, orient='index')

results_summary.to_csv(os.path.join(sum_dir, 'model_comparison_results.csv'))

# Time 
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script took {elapsed_time/60:.2f} minutes to run.")