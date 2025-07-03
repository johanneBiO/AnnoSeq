##########################################################################################
# This script contains functions used by main.py to predict regions from pretrained 
# XGBoost model and extract and summarize corresponding ESM-2 embeddings.
##########################################################################################

import os
import pickle
import numpy as np
import pandas as pd
from joblib import load
from pathlib import Path
from scipy.ndimage import gaussian_filter1d 

#---------------------------------------------------------------------------------------#
# Functions
#---------------------------------------------------------------------------------------#

def load_xgb_models(model_folder):
    models = {}
    for file in os.listdir(model_folder):
        if file.endswith('.pkl'):
            model_name = file.replace('.pkl', '')
            model_path = os.path.join(model_folder, file)
            model = load(model_path)
            models[model_name] = model
    return models

def predict_and_smooth(model, X, smooth_window = 5, smooth_sigma=2):
    pred_proba = model.predict_proba(X)[:, 1]
    smooth_proba = gaussian_filter1d(pred_proba, radius = smooth_window, sigma=smooth_sigma)
    return smooth_proba

def extract_regions(probabilities, threshold=0.5):
    above = probabilities > threshold
    regions = []
    start = None
    for i, val in enumerate(above):
        if val and start is None:
            start = i
        elif not val and start is not None:
            regions.append((start, i - 1))
            start = None
    if start is not None:
        regions.append((start, len(probabilities) - 1))
    return regions

def process_data_with_models(test_df, model_folder, smooth_window = 5, smooth_sigma=2, threshold=0.5):
    
    models = load_xgb_models(model_folder)
    results = []

    feature_cols = [col for col in test_df.columns if col not in ['accession', 'position', 'feature_type', 'partition']]

    for accession, group in test_df.groupby('accession'):
        
        group = group.sort_values('position')
        X = group[feature_cols].values
        pos = group['position'].values

        for model_name, model in models.items():
            smooth_proba = predict_and_smooth(model, X, smooth_window = smooth_window, smooth_sigma=smooth_sigma)
            
            regions = extract_regions(smooth_proba, threshold=threshold)
            for start_idx, end_idx in regions:
                results.append({
                    'accession': accession,
                    'start_position': pos[start_idx],
                    'end_position': pos[end_idx],
                    'feature_type': model_name
                })
    
    return pd.DataFrame(results)

def summarize_region(embedding, start, end):
    region = embedding[start-1:end, :]  
    return region.mean(axis=0)