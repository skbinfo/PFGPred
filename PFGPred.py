import pandas as pd
import numpy as np
import pickle
import xgboost as xgb
import tensorflow as tf
import joblib
import re
import warnings
import argparse
import os

warnings.filterwarnings('ignore')

def sanitize_feature_name(name):

    return re.sub(r'[\[\]<>\-,\(\):]', '_', str(name))

def run_prediction(input_path, model_type, confidence_cutoff, output_dir):
    print(f"--- Starting Prediction ---")
    print(f"Input file: {input_path}")
    print(f"Model type: {model_type}")
    print(f"Confidence cutoff: {confidence_cutoff}")
    print(f"Output directory: {output_dir}")

    if model_type == 'Arabidopsis':
        model_dir = 'Arab'
        prefix = 'a_'
    else: 
        model_dir = 'All_Plant'
        prefix = ''

    scaler_path = os.path.join(model_dir, f'{prefix}scaler.pkl')
    xgb_model_path = os.path.join(model_dir, f'{prefix}XGB.model')
    lstm_model_path = os.path.join(model_dir, f'{prefix}lstm_model.keras')
    rf_model_path = os.path.join(model_dir, f'{prefix}RF_model.pkl')
    meta_learner_path = os.path.join(model_dir, f'{prefix}meta_learner.pkl')
    splice_map_path = os.path.join(model_dir, f'{prefix}splice_site_mapping.pkl')
    encoders_path = os.path.join(model_dir, f'{prefix}binary_encoders.pkl') # Assuming you will create this file
    features_path = os.path.join(model_dir, f'{prefix}encoded_features.csv')
    
    try:
        scaler = joblib.load(scaler_path)
        xgb_model = xgb.Booster()
        xgb_model.load_model(xgb_model_path)
        lstm_model = tf.keras.models.load_model(lstm_model_path)
        rf_model = joblib.load(rf_model_path)
        meta_learner = joblib.load(meta_learner_path)
    except Exception as e:
        print(f"\n--- Error loading a model or scaler file: {e}. ---")
        print("--- Please ensure all model files are correctly placed in their directories. ---")
        exit()

    try:
        original_data = pd.read_csv(input_path)
        independent_data = original_data.copy()
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found. Exiting.")
        exit()

    columns_to_drop = [
        'fusion_pair', '5_geneid', '3_geneid', 'CDS_LEFT_ID', 'CDS_RIGHT_ID',
        '5_gene_go_term', '3_gene_go_term', 'shs_sequence', 'shs_length',
        'Chromosome1', 'Chromosome2'
    ]
    independent_data = independent_data.drop(columns=[col for col in columns_to_drop if col in independent_data.columns], errors='ignore')

    # Define feature types
    numerical_features = [
        'LeftBreakpoint', 'RightBreakpoint', 'FFPM', 'Total_Mapped_Reads', 'Left_Exon', 'Right_Exon',
        '5_gene_start', '5_gene_end', '5_gene_length', '3_gene_start', '3_gene_end', '3_gene_length',
        'exon_count5', 'exon_count3', 'Total_Count_(SC+RC)', 'alternate_junction_count',
        'Splice_Site'
    ]
    categorical_features = [
        'LeftStrand', 'RightStrand', 'Chromosome_Feature', 'Same_Strand', 'Reciprocal_Fusion', 
        'Splice_Pattern', 'Splice_Pattern_Class', '5_loc', '3_loc', 'alternate_junction'
    ]

    independent_data = independent_data.replace('NF', 'Unknown')
    for col in numerical_features:
        if col in independent_data.columns and col != 'Splice_Site':
            independent_data[col] = pd.to_numeric(independent_data[col], errors='coerce')
            independent_data[col] = independent_data[col].fillna(independent_data[col].median())
    for col in categorical_features:
        if col in independent_data.columns:
            independent_data[col] = independent_data[col].fillna('Unknown')

    try:
        with open(splice_map_path, 'rb') as f:
            splice_site_mapping = pickle.load(f)
        independent_data['Splice_Site'] = independent_data['Splice_Site'].map(splice_site_mapping).fillna(0)
    except FileNotFoundError:
        print(f"Warning: '{splice_map_path}' not found. Splice_Site encoding may be inconsistent.")
        unique_vals = independent_data['Splice_Site'].astype(str).unique()
        fallback_map = {val: i for i, val in enumerate(unique_vals)}
        independent_data['Splice_Site'] = independent_data['Splice_Site'].map(fallback_map).fillna(0)

    try:
        with open(encoders_path, 'rb') as f:
            binary_encoders = pickle.load(f)
    except FileNotFoundError:
        print(f"Warning: '{encoders_path}' not found. This file is crucial for correct encoding.")
        binary_encoders = {col: independent_data[col].unique().tolist() for col in categorical_features}

    encoded_data = independent_data.copy()
    for col in categorical_features:
        if col in encoded_data.columns:
            train_unique_values = binary_encoders.get(col, [])
            for value in train_unique_values:
                new_col_name = f"{col}_{sanitize_feature_name(value)}"
                encoded_data[new_col_name] = (encoded_data[col] == value).astype(int)
            encoded_data.drop(columns=[col], inplace=True)

    try:
        training_features_df = pd.read_csv(features_path)
        training_features = [col for col in training_features_df.columns if col != 'label']
        X_independent = encoded_data.reindex(columns=training_features, fill_value=0)
    except FileNotFoundError:
        print(f"Warning: '{features_path}' not found. Feature alignment may be incorrect.")
        X_independent = encoded_data

    X_independent = X_independent.fillna(0)

    X_independent_scaled = scaler.transform(X_independent)
    X_independent_lstm = X_independent_scaled.reshape((X_independent_scaled.shape[0], 1, X_independent_scaled.shape[1]))

    dtest = xgb.DMatrix(X_independent)
    xgb_predictions_prob = xgb_model.predict(dtest)
    lstm_predictions_prob = lstm_model.predict(X_independent_lstm, verbose=0).flatten()
    rf_predictions_prob = rf_model.predict_proba(X_independent_scaled)[:, 1]

    X_meta = np.hstack([
        xgb_predictions_prob.reshape(-1, 1),
        lstm_predictions_prob.reshape(-1, 1),
        rf_predictions_prob.reshape(-1, 1)
    ])
    ensemble_predictions_prob = meta_learner.predict_proba(X_meta)[:, 1]
    
    ensemble_pred_labels = (ensemble_predictions_prob >= confidence_cutoff).astype(int)

    predictions_df = pd.DataFrame({
        'Ensemble_Probability': ensemble_predictions_prob,
        'Predicted_Label': ensemble_pred_labels
    })
    
    results_df = pd.concat([original_data.reset_index(drop=True), predictions_df], axis=1)

    per_entry_output = os.path.join(output_dir, 'prediction_results.csv')
    results_df.to_csv(per_entry_output, index=False)
    print(f"Detailed prediction results saved to '{per_entry_output}'")
    
    positive_count = np.sum(ensemble_pred_labels)
    total_count = len(ensemble_pred_labels)
    summary_data = {
        'Metric': ['Positive Predictions (Fusion)', 'Negative Predictions (Non-Fusion)'],
        'Count': [positive_count, total_count - positive_count]
    }
    summary_df = pd.DataFrame(summary_data)
    summary_output = os.path.join(output_dir, 'prediction_summary.csv')
    summary_df.to_csv(summary_output, index=False)
    print(f"Prediction summary saved to '{summary_output}'")
    print("\n--- Prediction Complete ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="StackFG: Fusion Gene Prediction Tool")
    parser.add_argument("input_file", type=str, help="Path to the input CSV file for prediction.")
    parser.add_argument("--model", type=str, default="All_Plant", choices=["All_Plant", "Arabidopsis"],
                        help="The prediction model to use. Defaults to 'All_Plant'.")
    parser.add_argument("--cutoff", type=float, default=0.9,
                        help="Confidence cutoff for positive prediction (0.0 to 1.0). Defaults to 0.9.")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="The directory where result files will be saved.")
    
    args = parser.parse_args()

    if not 0.0 <= args.cutoff <= 1.0:
        print("Error: Confidence cutoff must be between 0.0 and 1.0.")
        exit()

    run_prediction(
        input_path=args.input_file,
        model_type=args.model,
        confidence_cutoff=args.cutoff,
        output_dir=args.output_dir
    )

