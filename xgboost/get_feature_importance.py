import os
import pandas as pd
from joblib import load  # assuming you use joblib to save/load models

def load_xgb_models(model_folder):
    models = {}
    for file in os.listdir(model_folder):
        if file.endswith('.pkl'):
            model_name = file.replace('.pkl', '')
            model_path = os.path.join(model_folder, file)
            model = load(model_path)
            models[model_name] = model
    return models

def save_feature_importance(models, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    
    for model_name, model in models.items():
        booster = model.get_booster()
        importance_dict = booster.get_score(importance_type='gain')
        
        # Convert to DataFrame
        importance_df = pd.DataFrame(list(importance_dict.items()), columns=['Feature', 'Importance'])
        importance_df = importance_df.sort_values(by='Importance', ascending=False)
        
        # Save CSV
        csv_path = os.path.join(output_folder, f'{model_name}_feature_importance.csv')
        importance_df.to_csv(csv_path, index=False)
        print(f'Saved feature importance for {model_name} to {csv_path}')

# Example usage:
model_folder = "/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/results/xgboost_headQuan90_colQuan_01000_100iter/models/"
output_folder = "/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/results/xgboost_headQuan90_colQuan_01000_100iter/feature_importance/"

models = load_xgb_models(model_folder)
save_feature_importance(models, output_folder)