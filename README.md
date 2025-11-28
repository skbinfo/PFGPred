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
(pfgpred_env) [sk-202@localhost PFGPred]$ git clone https://github.com/skbinfo/PFGPred.git
```
2. **Create Environment**:
```bash
(base) [sk-202@localhost PFGPred]$ cd PFGPred
(base) [sk-202@localhost PFGPred]$ chmod +x setup_env.sh
(base) [sk-202@localhost PFGPred]$ ./setup_env.sh
```
3. **Activate Environment**:
```bash
conda activate pfgpred_env
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
(pfgpred_env) [sk-202@localhost PFGPred]$ mkdir Fusion_output
```


5. **Organize Input Files**

Place all input FASTQ files into the directory where ft.py is located.
Run Feature Extraction
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python ft.py --help
usage: ft.py [-h] [--start {step1,step2}] [--config CONFIG] [--merged MERGED] [--gtf GTF] [--outdir OUTDIR]

Fusion pipeline (resumable).

options:
  -h, --help            show this help message and exit
  --start {step1,step2}
                        Which step to start from. step1 runs everything (including STAR-Fusion). step2 runs GTF processing & downstream.
  --config CONFIG       Path to config.yaml (required for step1)
  --merged MERGED       Path to 1_initial_merged_predictions.tsv (required for step3)
  --gtf GTF             Path to genome GTF (required for step1/2)
  --outdir OUTDIR       Output directory



```

This step:

-**Processes outputs from fusion callers (Default: STAR-Fusion)**
  
-**Extracts RNA-Seq based features for each fusion transcript**

-**Generates a feature table compatible with the model**


6. **Run the Prediction Model**
To classify fusions using the pretrained stacked ensemble model:
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python PFGPred.py --help
usage: PFGPred.py [-h] [--model {All_Plant,Arabidopsis}] [--cutoff CUTOFF] --output-dir OUTPUT_DIR input_file

PFGPred: Plant Fusion Gene Predictor

positional arguments:
  input_file            Path to the input CSV file for prediction.

options:
  -h, --help            show this help message and exit
  --model {All_Plant,Arabidopsis}  The prediction model to use. Defaults to 'All_Plant'.                 
  --cutoff CUTOFF       Confidence cutoff for positive prediction (0.0 to 1.0). Defaults to 0.9.
  --output-dir OUTPUT_DIR  The directory where result files will be saved.

(pfgpred_env) [sk-202@localhost PFGPred]$ mkdir OutFile #Output Directory must be pre-existing

(pfgpred_env) [sk-202@localhost PFGPred]$ python PFGPred.py data.csv --model All_Plant --output-dir OutFile
To enable the following instructions: AVX512F, in other operations, rebuild TensorFlow with the appropriate compiler flags.
--- Starting Prediction ---
Input file: data.csv
Model type: All_Plant
Confidence cutoff: 0.9
Output directory: OutFile
Detailed prediction results saved to 'OutFile/prediction_results.csv'
Prediction summary saved to 'OutFile/prediction_summary.csv'

--- Prediction Complete ---
```
**Or upload the generated feature table to the PFGPred Web Server (http://223.31.159.15/PFGPred/predict.php) if the number of entries in the file is less than 10000.**


The output will provide:


1.**probability scores**


2.**classification labels**




3.**high-confidence filtered fusion list (Default threshold 0.9)**


# Kickstart Mode (Using Outputs from Any Fusion Caller)


PFGPred also supports a Kickstart Mode, allowing you to use fusion predictions generated by any RNA-Seq fusion detection tool. This mode is useful when STAR-Fusion is not preferred or when you already have fusion caller outputs from another workflow. In Kickstart mode, you simply provide the output of your fusion caller to the feature extraction module and then run the PFGPred model.

1. **Input Format**

Prepare the fusion caller output in the following format. If any field is missing or not available, use **. (dot)** or **NF** to indicate missing information. This ensures consistent parsing during feature extraction.

| Sample_ID | #FusionName | 5_geneid | 3_geneid | FFPM | CDS_LEFT_ID | CDS_RIGHT_ID | Splice_Pattern | Chromosome1 | LeftBreakpoint | LeftStrand | Chromosome2 | RightBreakpoint | RightStrand | Splice_Site | Total_Count_(SC+RC) | Splice_Pattern_Class |
| :-------- | :---------- | :------- | :------- | :---- | :----------- | :------------ | :-------------- | :----------- | :-------------- | :---------- | :----------- | :--------------- | :----------- | :----------- | :-------------------- | :--------------------- |
| ERR1588794 | AT2G47640->AT4G18230 | AT2G47640 | AT4G18230 | 2.1745 | . | . | . | 2 | 19537529 | + | 4 | 10081731 | - | GT-AG | 18 | Canonical |
| ERR1588794 | AT2G37600->AT5G02450 | AT2G37600 | AT5G02450 | 2.1746 | transcript:AT2G37600.1 | transcript:AT5G02450.1 | INFRAME | 2 | 15774593 | - | 5 | 534358 | + | GT-AG | 18 | Canonical |
| ERR1588794 | AT3G23050->AT2G22670 | AT3G23050 | AT2G22670 | 0.2416 | transcript:AT3G23050.1 | transcript:AT2G22670.1 | INFRAME | 3 | 8196210 | + | 2 | 9638275 | + | GT-AG | 2 | Canonical |
| SRR3458666 | AT3G53870->AT5G35530 | AT3G53870 | AT5G35530 | 1.5013 | transcript:AT3G53870.1 | transcript:AT5G35530.1 | INFRAME | 3 | 19952063 | + | 5 | 13710916 | - | AT-AC | 6 | Canonical |
| SRR3458666 | AT4G27040->AT4G27080 | AT4G27040 | AT4G27080 | 0.7507 | transcript:AT4G27040.1 | transcript:AT4G27080.1 | INFRAME | 4 | 13573491 | - | 4 | 13591228 | + | GT-AG | 3 | Canonical |
| SRR3458666 | AT2G24060->AT1G56350 | AT2G24060 | AT1G56350 | 0.7507 | transcript:AT2G24060.1 | transcript:AT1G56350.1 | INFRAME | 2 | 10230006 | + | 1 | 21097040 | + | AT-AG | 3 | Non-Canonical |
| SRR3458666 | AT1G60800->AT1G04880 | AT1G60800 | AT1G04880 | 0.7507 | . | . | . | 1 | 22385933 | - | 1 | 1378266 | - | GT-AG | 3 | Canonical |
| SRR3458666 | AT1G71830->AT2G13790 | AT1G71830 | AT2G13790 | 0.5004 | transcript:AT1G71830.1 | transcript:AT2G13790.1 | FRAMESHIFT | 1 | 27019380 | + | 2 | 5743213 | + | GT-AG | 2 | Canonical |



2. **To construct a PFGPred-compatible feature table from your fusion caller outputs, follow the steps below:**
   
-**Create an output directory named Fusion_output.** 
-**Copy your fusion caller output into the Fusion_output directory and rename the file to 1_initial_merged_predictions.tsv**
-**Run the feature extraction script:**
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python ft.py --help
usage: ft.py [-h] [--start {step1,step2}] [--config CONFIG] [--merged MERGED] [--gtf GTF] [--outdir OUTDIR]

Fusion pipeline (resumable).

options:
  -h, --help            show this help message and exit
  --start {step1,step2}
                        Which step to start from. step1 runs everything (including STAR-Fusion). step2 runs GTF processing & downstream.
  --config CONFIG       Path to config.yaml (required for step1)
  --merged MERGED       Path to 1_initial_merged_predictions.tsv (required for step3)
  --gtf GTF             Path to genome GTF (required for step1/2)
  --outdir OUTDIR       Output directory

(pfgpred_env) [sk-202@localhost PFGPred]$ python ft.py --start step2 --gtf Oryza_sativa.IRGSP-1.0.56.gtf
Parsing GTF (exons) from Oryza_sativa.IRGSP-1.0.56.gtf ...
Saved 2_filtered_gtf.tsv
Loading merged predictions from Fusion_output/1_initial_merged_predictions.tsv
Saved 3_merged_filled.tsv
Parsing GTF for exon coordinates from Oryza_sativa.IRGSP-1.0.56.gtf ...
Saved 4_final_merged_predictions.tsv
Saved 5_fusion_annotated.tsv
Saved 6_gene_info.tsv
Saved 7_fusion_with_gene_info.tsv
/home/sk-202/sakshi/try/PFGPred/ft.py:502: FutureWarning: DataFrameGroupBy.apply operated on the grouping columns. This behavior is deprecated, and in a future version of pandas the grouping columns will be excluded from the operation. Either pass `include_groups=False` to exclude the groupings or explicitly select the grouping columns after groupby to silence this warning.
  fusion_counts = df.groupby(["5_geneid","3_geneid"]).apply(
Saved 8_fusion_with_junction_counts.tsv and 9_fusion_counts_formatted.tsv
Saved 10_fusion_with_breakpoint_locations.tsv
Saved 11_final_fusion_annotated.tsv
Saved 12_annotated_with_reads.tsv
Saved fts_features.csv
Cleanup (temp_dir) done.

```

3. **Run the Prediction Model**

Then upload the generated feature table **(fts_features.csv)** to the **PFGPred Web Server (http://223.31.159.15/PFGPred/predict.php)** or use the following command to classify fusion transcripts using the PFGPred model:
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python PFGPred.py --help
usage: PFGPred.py [-h] [--model {All_Plant,Arabidopsis}] [--cutoff CUTOFF] --output-dir OUTPUT_DIR input_file

PFGPred: Plant Fusion Gene Predictor

positional arguments:
  input_file            Path to the input CSV file for prediction.

options:
  -h, --help            show this help message and exit
  --model {All_Plant,Arabidopsis}  The prediction model to use. Defaults to 'All_Plant'.                 
  --cutoff CUTOFF       Confidence cutoff for positive prediction (0.0 to 1.0). Defaults to 0.9.
  --output-dir OUTPUT_DIR  The directory where result files will be saved.

(pfgpred_env) [sk-202@localhost PFGPred]$ mkdir OutFile #Output Directory must be pre-existing

(pfgpred_env) [sk-202@localhost PFGPred]$ python PFGPred.py data.csv --model All_Plant --output-dir OutFile
To enable the following instructions: AVX512F, in other operations, rebuild TensorFlow with the appropriate compiler flags.
--- Starting Prediction ---
Input file: data.csv
Model type: All_Plant
Confidence cutoff: 0.9
Output directory: OutFile
Detailed prediction results saved to 'OutFile/prediction_results.csv'
Prediction summary saved to 'OutFile/prediction_summary.csv'
```

# Model Retraining Workflow

If you want to retrain PFGPred for a new plant species:

1. **Fusion Identification & Feature Extraction (for training)**
   
Run the standard extraction step:
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python ft.py --start step1 --config config.yaml --gtf Oryza_sativa.IRGSP-1.0.56.gtf
Saved read lengths to Fusion_output/reads.tsv
....................................................
ALL DONE: final files produced (11_final_fusion_annotated.tsv, 12_annotated_with_reads.tsv, fts_features.csv)
```
This step will generate the fts_features.csv file.


2. **Validation of Fusion Transcripts Using WGS**

To validate fusion breakpoints at the DNA level, run:

** Use the published WGS fusion pipeline:**
```bash
https://github.com/VolundurH/wgs_fusion_pipeline
```

This step produces a set of WGS-validated (true) fusion genes. This step will generate new_fusion_search_table.txt and fusion_summary_table_breakpoint_supporting_reads.txt files.


3. **Generate Positive and Negative Datasets**

To construct the training dataset, run:

```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python segregate.py --help
usage: segregate.py [-h] --summary SUMMARY --annot ANNOT --feat FEAT

Fusion processing pipeline (produce final merged CSV with labels)

options:
  -h, --help         show this help message and exit
  --summary SUMMARY  fusion_summary_table_breakpoint_supporting_reads.txt
  --annot ANNOT      new_fusion_search_table.txt
  --feat FEAT        fts_features.csv
(pfgpred_env) [sk-202@localhost PFGPred]$ python segregate.py  --summary fusion_summary_table_breakpoint_supporting_reads.txt --annot new_fusion_search_table.txt --feat Fusion_output/fts_features.csv
STEP 1: Removing duplicates...
 → unique_fusion_summary.tsv created.
STEP 2: Formatting to split chr:bpt...
 → formatted_fusion_summary.tsv created.
STEP 3: Adding annotations...
 → formatted_fusion_summary_annotated.tsv created.
STEP 4: Adding gene start/end...
 → merged_fusion_with_coords.tsv created.
STEP 5: Checking overlap...
 → overlapping_genes.csv and non_overlapping_genes.csv created.
STEP 6: Matching against feature table and creating final merged CSV with labels...
Cleaning up intermediary files...

 PIPELINE COMPLETE!
 → Final merged CSV with labels: final_with_labels.csv

```

This script generates RNA-Seq fusions supported by WGS as positive (post.tsv) and unsupported as negative (neg.tsv).

4. **Train Your Own Model**
```bash
(pfgpred_env) [sk-202@localhost PFGPred]$ python train_PFGPred.py --help

usage: train_PFGPred.py [-h] --train-file TRAIN_FILE --train-target-column TRAIN_TARGET_COLUMN [--test-file TEST_FILE] [--test-target-column TEST_TARGET_COLUMN]
                        --output-dir OUTPUT_DIR [--split-ratio SPLIT_RATIO] [--xgb-eta XGB_ETA] [--xgb-max-depth XGB_MAX_DEPTH] [--rf-n-estimators RF_N_ESTIMATORS]
                        [--rf-max-depth RF_MAX_DEPTH] [--lstm-epochs LSTM_EPOCHS] [--lstm-batch-size LSTM_BATCH_SIZE] [--meta-c META_C]

PFGPred Model Training Script

options:
  -h, --help            show this help message and exit
  --train-file TRAIN_FILE
  --train-target-column TRAIN_TARGET_COLUMN
  --test-file TEST_FILE
  --test-target-column TEST_TARGET_COLUMN
  --output-dir OUTPUT_DIR
  --split-ratio SPLIT_RATIO
                        Test set ratio if no test file provided.
  --xgb-eta XGB_ETA
  --xgb-max-depth XGB_MAX_DEPTH
  --rf-n-estimators RF_N_ESTIMATORS
  --rf-max-depth RF_MAX_DEPTH
  --lstm-epochs LSTM_EPOCHS
  --lstm-batch-size LSTM_BATCH_SIZE
  --meta-c META_C

(pfgpred_env) [sk-202@localhost PFGPred]$ python train_PFGPred.py --train-file final_with_labels.csv --train-target-column 'label' --split-ratio 0.2 --output-dir OutFile
--- Starting PFGPred Training Process ---
Loaded training data: 10000 rows.
--- Training Data Sample Counts ---
label
1    5000
0    5000

No test file provided. Splitting training data into 80% training and 20% testing sets.

--- Preprocessing Data... ---
Preprocessing complete.

--- Saving preprocessing objects for future predictions... ---
Saved scaler.pkl, splice_site_mapping.pkl, binary_encoders.pkl, and encoded_features.csv

--- Training Base Models... ---
Base models trained successfully.

--- Training PFGPred Meta-Learner... ---
Meta-learner trained successfully.

--- Evaluating PFGPred model on test data... ---

--- Generating Output Files... ---
Saved metrics to OutFile/metrics.csv
Saved ROC curve data to OutFile/roc_data.csv
Saved test predictions to OutFile/test_predictions.csv

```
Or Users can upload the positive and negative datasets to the web server ***(http://223.31.159.15/PFGPred/train.php)*** to train a custom species-specific model.
