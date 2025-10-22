#!/bin/bash
# This script creates the required Conda environment for PFGPred.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Checking for Conda... ---"
if ! command -v conda &> /dev/null
then
    echo ""
    echo "--- ERROR: Conda could not be found. ---"
    echo "--- Please install Conda (Miniconda or Anaconda) first. ---"
    echo "--- Visit: https://docs.conda.io/projects/miniconda/en/latest/ ---"
    exit 1
fi

ENV_NAME="pfgpred_env" # Specific name for the user's environment

# Check if the environment already exists and remove it to ensure a clean install
if conda env list | grep -q "$ENV_NAME"; then
    echo ""
    echo "--- Found existing environment '$ENV_NAME'. Removing it for a clean setup... ---"
    conda env remove -n "$ENV_NAME"
fi

echo ""
echo "--- Creating Conda environment '$ENV_NAME' from environment.yml... ---"
echo "--- This might take several minutes depending on your internet connection. ---"

# Create the environment directly from the YAML file.
conda env create -f environment.yml -n "$ENV_NAME"

echo ""
echo "--- Environment '$ENV_NAME' created successfully! ---"
echo ""
echo "--- To use the tools, activate the environment first: ---"
echo "conda activate $ENV_NAME"
echo ""
echo "--- Example Usage: ---"
echo "python PFGPred.py your_input.csv --model All_Plant --output-dir ./results"
echo "python train_PFGPred.py --train-file train.csv --train-target-column label --output-dir ./training_output"
echo ""
echo "--- When finished, deactivate the environment with: ---"
echo "conda deactivate"
