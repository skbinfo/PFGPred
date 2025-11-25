# PFGPred: A stack ensemble classifier for the detection of fusion genes in plants from RNA-seq data
PFGPred (Plant Fusion Gene Predictor) is a comprehensive bioinformatics pipeline for the identification of fusion genes in plants using RNA-Seq data. 
It integrates conventional fusion detection tools, a feature extraction framework, and a stacked ensemble machine-learning model specifically trained on plant datasets to identify high-confidence fusion events with significantly reduced false positives.

# PFGPred operates in two major steps:
**1. Fusion Detection & Feature Extraction:** RNA-Seq data are processed using conventional fusion callers to detect candidate fusion transcripts. For each candidate, a set of features derived from RNA-Seq and genomic annotation was extracted.

**2.Classification of True vs. False Fusions:** Extracted features were evaluated using a plant-optimized stacked ensemble model (XGBoost, Random Forest, and LSTM with a logistic regression meta-classifier) to distinguish true fusion transcripts from false positives. This ensures the identification of biologically relevant, high-confidence fusion genes.

PFGPred thus provides an end-to-end, accurate, and plant-specific framework for fusion gene discovery without the need for validation on WGS data.

**Key features of PFGPred:**
-**Integrating RNA-Seq fusion predictions with machine learning**
-**Trained on fusion events validated on WGS datasets**
-**Supports any RNA-Seq fusion caller**
-**Allows species-specific retraining**

# Installation
1. **Clone the repository**:
```bash
git clone
cd PFGPred
chmod +x setup.sh
./setup.sh
conda activate PFGPred_env
```
