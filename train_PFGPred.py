import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, LSTM
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score,
                             roc_auc_score, matthews_corrcoef, confusion_matrix, roc_curve)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import xgboost as xgb
import joblib
import argparse
import os
import warnings
import pickle
import re
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR') # Suppress Keras/TF messages
warnings.filterwarnings('ignore')

def sanitize_feature_name(name):
    return re.sub(r'[\[\]<>\-,\(\):]', '_', str(name))

def specificity_score(y_true, y_pred):
    try:
        # Ensure predictions are binary (0 or 1)
        y_pred = np.round(y_pred).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        # Handle potential division by zero
        return tn / (tn + fp) if (tn + fp) > 0 else 0
    except ValueError:
        # Handle cases where confusion matrix doesn't have 4 values (e.g., only one class predicted)
        return 0

def validate_target_column(series, column_name, dataset_name):
    unique_values = series.unique()
    if not np.all(np.isin(unique_values, [0, 1])):
        print(f"Error: Target column '{column_name}' in the {dataset_name} data must contain only 0s and 1s.")
        print(f"Found unique values: {unique_values}")
        exit(1)
    print(f"--- {dataset_name.capitalize()} Data Sample Counts ---\n{series.value_counts().to_string()}\n")

def run_training(args):
    print("--- Starting PFGPred Training Process ---")

    train_df = pd.read_csv(args.train_file)
    print(f"Loaded training data: {train_df.shape[0]} rows.")
    validate_target_column(train_df[args.train_target_column], args.train_target_column, "training")

    columns_to_drop = ['fusion_pair', '5_geneid', '3_geneid', 'CDS_LEFT_ID', 'CDS_RIGHT_ID', '5_gene_go_term', '3_gene_go_term', 'shs_sequence', 'shs_length', 'Chromosome1', 'Chromosome2']
    X_train_full = train_df.drop(columns=[col for col in columns_to_drop if col in train_df.columns], errors='ignore')
    y_train_full = X_train_full.pop(args.train_target_column) # Extract target

    if args.test_file:
        test_df = pd.read_csv(args.test_file)
        print(f"Loaded test data: {test_df.shape[0]} rows.")
        validate_target_column(test_df[args.test_target_column], args.test_target_column, "testing")
        X_test_original = test_df.drop(columns=[col for col in columns_to_drop if col in test_df.columns], errors='ignore')
        y_test = X_test_original.pop(args.test_target_column) # Extract target
        X_train = X_train_full
        y_train = y_train_full
    else:
        split_ratio = args.split_ratio
        print(f"No test file provided. Splitting training data into {int((1-split_ratio)*100)}% training and {int(split_ratio*100)}% testing sets.")
        X_train, X_test_original, y_train, y_test = train_test_split(
            X_train_full, y_train_full, test_size=split_ratio, random_state=42, stratify=y_train_full
        )

    print("\n--- Preprocessing Data... ---")
    numerical_features = ['LeftBreakpoint', 'RightBreakpoint', 'FFPM', 'Total_Mapped_Reads', 'Left_Exon', 'Right_Exon', '5_gene_start', '5_gene_end', '5_gene_length', '3_gene_start', '3_gene_end', '3_gene_length', 'exon_count5', 'exon_count3', 'Total_Count_(SC+RC)', 'alternate_junction_count']
    categorical_features = ['LeftStrand', 'RightStrand', 'Chromosome_Feature', 'Same_Strand', 'Reciprocal_Fusion', 'Splice_Pattern', 'Splice_Pattern_Class', '5_loc', '3_loc', 'alternative_junction'] # Corrected spelling as per reference

    X_test = X_test_original.copy() # Make a copy to avoid modifying original split data

    for df in [X_train, X_test]:
        df.replace(['NF', '-', '.'], np.nan, inplace=True) # Replace before numeric conversion
        for col in categorical_features:
            if col in df.columns:
                df[col].fillna('Unknown', inplace=True)
        for col in numerical_features:
             if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        train_median = X_train[numerical_features].median() # Calculate median on cleaned training data
        df.fillna(train_median, inplace=True) # Apply median imputation to both train and test

    splice_site_mapping = {}
    if 'Splice_Site' in X_train.columns:
        X_train['Splice_Site'].fillna('Unknown', inplace=True) # Fill NaNs before creating map
        splice_site_values = X_train['Splice_Site'].unique()
        splice_site_mapping = {value: i for i, value in enumerate(splice_site_values)}
        if 'Unknown' in splice_site_mapping:
            unknown_val = splice_site_mapping.pop('Unknown')
            temp_map = {k: i + 1 for i, k in enumerate(splice_site_mapping.keys())}
            splice_site_mapping = {'Unknown': 0, **temp_map} # Ensure Unknown maps to 0

        X_train['Splice_Site'] = X_train['Splice_Site'].map(splice_site_mapping)
        if 'Splice_Site' in X_test.columns:
            X_test['Splice_Site'].fillna('Unknown', inplace=True)
            # Use mapping learned from training data, fill unknown test values with 0
            X_test['Splice_Site'] = X_test['Splice_Site'].map(splice_site_mapping).fillna(0)

    binary_encoders = {}
    X_train_encoded = X_train.copy()
    X_test_encoded = X_test.copy()

    for col in categorical_features:
        if col in X_train.columns:
            unique_values = X_train[col].unique()
            binary_encoders[col] = unique_values
            for value in unique_values:
                sanitized_value = sanitize_feature_name(value)
                X_train_encoded[f"{col}_{sanitized_value}"] = (X_train[col] == value).astype(int)
                X_test_encoded[f"{col}_{sanitized_value}"] = (X_test[col] == value).astype(int)
            # Drop original column from both
            X_train_encoded.drop(columns=[col], inplace=True)
            X_test_encoded.drop(columns=[col], inplace=True)
    final_feature_list = list(X_train_encoded.columns)
    X_test_encoded = X_test_encoded.reindex(columns=final_feature_list, fill_value=0)

    X_train_encoded.fillna(0, inplace=True)
    X_test_encoded.fillna(0, inplace=True)

    print("Preprocessing complete.")

    print("\n--- Saving preprocessing objects for future predictions... ---")
    # Fit scaler ONLY on training data
    scaler_to_save = StandardScaler().fit(X_train_encoded)
    joblib.dump(scaler_to_save, os.path.join(args.output_dir, 'scaler.pkl'))
    # Save the learned mappings and encoders
    with open(os.path.join(args.output_dir, 'splice_site_mapping.pkl'), 'wb') as f:
        pickle.dump(splice_site_mapping, f)
    with open(os.path.join(args.output_dir, 'binary_encoders.pkl'), 'wb') as f:
        pickle.dump(binary_encoders, f)
    # Save the list of final features used for training
    pd.DataFrame(columns=final_feature_list).to_csv(os.path.join(args.output_dir, 'encoded_features.csv'), index=False)
    print("Saved scaler.pkl, splice_site_mapping.pkl, binary_encoders.pkl, and encoded_features.csv")

    # Apply the fitted scaler to both train and test data
    X_train_scaled = scaler_to_save.transform(X_train_encoded)
    X_test_scaled = scaler_to_save.transform(X_test_encoded)

    X_train_lstm = np.reshape(X_train_scaled, (X_train_scaled.shape[0], 1, X_train_scaled.shape[1]))
    X_test_lstm = np.reshape(X_test_scaled, (X_test_scaled.shape[0], 1, X_test_scaled.shape[1]))

    print("\n--- Training Base Models... ---")
    # Prepare XGBoost data AFTER scaling is done and preprocessing objects saved
    dtrain = xgb.DMatrix(X_train_encoded, label=y_train) # XGBoost uses unscaled data
    dtest = xgb.DMatrix(X_test_encoded, label=y_test)   # XGBoost uses unscaled data

    xgb_params = {'objective': 'binary:logistic', 'eval_metric': 'logloss', 'seed': 42, 'eta': args.xgb_eta, 'max_depth': args.xgb_max_depth}
    xgb_model = xgb.train(xgb_params, dtrain, num_boost_round=100, verbose_eval=False)

    rf_model = RandomForestClassifier(random_state=42, n_estimators=args.rf_n_estimators, max_depth=args.rf_max_depth)
    rf_model.fit(X_train_scaled, y_train) # RF uses scaled data

    def create_lstm_model(input_shape):
        model = Sequential([LSTM(64, input_shape=input_shape), Dropout(0.4), Dense(16, activation='relu'), Dropout(0.4), Dense(1, activation='sigmoid')])
        model.compile(optimizer='adam', loss='binary_crossentropy')
        return model
    lstm_model = create_lstm_model(input_shape=(X_train_lstm.shape[1], X_train_lstm.shape[2]))
    lstm_model.fit(X_train_lstm, y_train, epochs=args.lstm_epochs, batch_size=args.lstm_batch_size, verbose=0) # LSTM uses scaled data
    print("Base models trained successfully.")

    print("\n--- Training PFGPred Meta-Learner... ---")
    # Generate meta-features using predictions on TRAINING data
    X_train_meta = np.column_stack((
        xgb_model.predict(dtrain),                 # Predict on unscaled train
        rf_model.predict_proba(X_train_scaled)[:, 1], # Predict on scaled train
        lstm_model.predict(X_train_lstm, verbose=0).flatten() # Predict on scaled train
    ))
    X_train_meta = np.nan_to_num(X_train_meta, nan=0.5) # Impute NaNs if any

    meta_learner = LogisticRegression(random_state=42, C=args.meta_c)
    meta_learner.fit(X_train_meta, y_train)
    print("Meta-learner trained successfully.")

    print("\n--- Evaluating PFGPred model on test data... ---")
    # Generate meta-features using predictions on TEST data
    X_test_meta = np.column_stack((
        xgb_model.predict(dtest),                 # Predict on unscaled test
        rf_model.predict_proba(X_test_scaled)[:, 1], # Predict on scaled test
        lstm_model.predict(X_test_lstm, verbose=0).flatten() # Predict on scaled test
    ))
    X_test_meta = np.nan_to_num(X_test_meta, nan=0.5) # Impute NaNs if any

    ensemble_preds_prob = meta_learner.predict_proba(X_test_meta)[:, 1]
    ensemble_preds_label = (ensemble_preds_prob > 0.5).astype(int)

    print("\n--- Generating Output Files... ---")
    metrics_df = pd.DataFrame({'Metric': ['Accuracy', 'Precision', 'Sensitivity', 'Specificity', 'F1_Score', 'MCC', 'AUC'], 'Value': [accuracy_score(y_test, ensemble_preds_label), precision_score(y_test, ensemble_preds_label, zero_division=0), recall_score(y_test, ensemble_preds_label, zero_division=0), specificity_score(y_test, ensemble_preds_label), f1_score(y_test, ensemble_preds_label, zero_division=0), matthews_corrcoef(y_test, ensemble_preds_label), roc_auc_score(y_test, ensemble_preds_prob)]})
    metrics_df.to_csv(os.path.join(args.output_dir, 'metrics.csv'), index=False)
    print(f"Saved metrics to {os.path.join(args.output_dir, 'metrics.csv')}")

    fpr, tpr, _ = roc_curve(y_test, ensemble_preds_prob)
    pd.DataFrame({'FPR': fpr, 'TPR': tpr}).to_csv(os.path.join(args.output_dir, 'roc_data.csv'), index=False)
    print(f"Saved ROC curve data to {os.path.join(args.output_dir, 'roc_data.csv')}")

    predictions_df = pd.DataFrame({'Ensemble_Probability': ensemble_preds_prob, 'Predicted_Label': ensemble_preds_label})
    # Use X_test_original which still has the original categorical values before encoding
    final_results_df = pd.concat([X_test_original.reset_index(drop=True), y_test.reset_index(drop=True), predictions_df], axis=1)
    final_results_df.to_csv(os.path.join(args.output_dir, 'test_predictions.csv'), index=False)
    print(f"Saved test predictions to {os.path.join(args.output_dir, 'test_predictions.csv')}")

    print("\n--- Training and Evaluation Complete ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PFGPred Model Training Script")
    parser.add_argument('--train-file', required=True)
    parser.add_argument('--train-target-column', required=True)
    parser.add_argument('--test-file', default=None)
    parser.add_argument('--test-target-column', default=None)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--split-ratio', type=float, default=0.2, help="Test set ratio if no test file provided.")

    parser.add_argument('--xgb-eta', type=float, default=0.01)
    parser.add_argument('--xgb-max-depth', type=int, default=6)
    parser.add_argument('--rf-n-estimators', type=int, default=100)
    parser.add_argument('--rf-max-depth', type=int, default=6)
    parser.add_argument('--lstm-epochs', type=int, default=35)
    parser.add_argument('--lstm-batch-size', type=int, default=32)
    parser.add_argument('--meta-c', type=float, default=0.1)

    args = parser.parse_args()

    if args.test_file and not args.test_target_column:
        parser.error("--test-target-column is required when --test-file is provided.")

    run_training(args)
