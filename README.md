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
git clone https://github.com/skbinfo/PFGPred.git
```
2. **Create Environment**:
```bash
cd PFGPred
chmod +x setup.sh
./setup.sh
```
3. **Activate Environment**:
```bash
conda activate PFGPred_env
```

# Running the Pipeline
1. **Prepare the Configuration File**
Edit config.yaml and provide paths to the required reference files:
```bash
Reference genome (.fasta)
Genome annotation (.gtf)
CTAT resource library (required for STAR-Fusion)
```

2. **How to Prepare a CTAT Genome Library**
PFGPred uses STAR-Fusion as its default fusion detection tool, which requires a CTAT Genome Resource Library built for your species.
You can prepare this library by following the official CTAT tutorial:

**🔗 CTAT Genome Lib Preparation Guide**
 https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki


3. **Specify Samples**
List the sample identifiers in one of the following files based on your sample type:
```
single.txt for single-end reads
paired.txt for paired-end reads
```
Each line should contain one sample ID (matching FASTQ file names).


4. **Create Output Directory**
```bash
mkdir Fusion_output
```


5. **Organize Input Files**
Place all input FASTQ files into the directory where ft.py is located.
Run Feature Extraction
```bash
python ft.py \
--start step1 \
--config config.yaml \
–gtf path/to/annotation.gtf
```
This step:
- ***Processes outputs from fusion callers (Default: STAR-Fusion)***
-***Extracts RNA-Seq based features for each fusion transcript***
-***Generates a feature table compatible with the model***


6. **Run the Prediction Model**
To classify fusions using the pretrained stacked ensemble model:
```bash

```
***Or upload the generated feature table to the PFGPred Web Server (http://223.31.159.15/PFGPred/predict.php) if the number of entries in the file is less than 10000.***


The output will provide:
1.**probability scores**


2.**classification labels**




3.**high-confidence filtered fusion list (Default threshold 0.9)**
