import os
import requests
import json
import pandas as pd
import re
import shutil
import subprocess
import csv
import yaml

# Load config.yaml
with open("config.yaml") as f:
    config = yaml.safe_load(f)

fasta_file = config["fasta_file"]
gtf_file = config["gtf_file"] 
GENOME_LIB_DIR = config["genome_lib_dir"]

# Input files
single_file = "single.txt"
paired_file = "paired.txt"
output_file = "Fusion_output/reads.tsv"

# Ensure output directory exists
os.makedirs("Fusion_output", exist_ok=True)

# Function to extract read length using FASTQ format logic
def get_read_length(fq_file):
    try:
        with open(fq_file, "r") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:  # 2nd line of first FASTQ entry
                    return len(line.strip())
    except FileNotFoundError:
        return "File not found"

# ADDED: Helper function to find read files with flexible extensions
def find_read_file(patterns):
    """Checks a list of patterns and returns the first file that exists."""
    for pat in patterns:
        if os.path.isfile(pat):
            return pat
    return None

# Collect results
results = []

# Process single-end
with open(single_file, "r") as f:
    for line in f:
        sample = line.strip()
        # MODIFIED: Check for multiple file extensions
        patterns = [f"{sample}_val.fq", f"{sample}.fq", f"{sample}.fa", f"{sample}.fasta"]
        fq_file = find_read_file(patterns)
        
        if fq_file:
            read_len = get_read_length(fq_file)
        else:
            read_len = "File not found"
        results.append((read_len, sample))

# Process paired-end
with open(paired_file, "r") as f:
    for line in f:
        sample = line.strip()
        # MODIFIED: Check for multiple file extensions for R1
        patterns = [
            f"{sample}_1_val_1.fq", 
            f"{sample}_1.fq", f"{sample}_R1.fq",
            f"{sample}_1.fa", f"{sample}_R1.fa",
            f"{sample}_1.fasta", f"{sample}_R1.fasta"
        ]
        fq_file = find_read_file(patterns)
        
        if fq_file:
            read_len = get_read_length(fq_file)
        else:
            read_len = "File not found"
        results.append((read_len, sample))

# Write output
with open(output_file, "w") as out:
    out.write("Total_Mapped_Reads\tsample\n")
    for read_len, sample in results:
        out.write(f"{read_len}\t{sample}\n")

print(f"✅ Read lengths saved to {output_file}")


###### STAR FUSION  ####

# Load sample lists
with open("paired.txt") as f:
    paired_samples = {line.strip() for line in f if line.strip()}

with open("single.txt") as f:
    single_samples = {line.strip() for line in f if line.strip()}

# Output directories
PE_OUTPUT_DIR = "pe_output"
SE_OUTPUT_DIR = "se_output"
os.makedirs(PE_OUTPUT_DIR, exist_ok=True)
os.makedirs(SE_OUTPUT_DIR, exist_ok=True)

# Run Paired-end STAR-Fusion
for sample in paired_samples:
    # MODIFIED: Find R1 file
    r1_patterns = [
        f"{sample}_1_val_1.fq", 
        f"{sample}_1.fq", f"{sample}_R1.fq",
        f"{sample}_1.fa", f"{sample}_R1.fa",
        f"{sample}_1.fasta", f"{sample}_R1.fasta"
    ]
    r1 = find_read_file(r1_patterns)
    
    # MODIFIED: Find R2 file
    r2_patterns = [
        f"{sample}_2_val_2.fq", 
        f"{sample}_2.fq", f"{sample}_R2.fq",
        f"{sample}_2.fa", f"{sample}_R2.fa",
        f"{sample}_2.fasta", f"{sample}_R2.fasta"
    ]
    r2 = find_read_file(r2_patterns)

    if not r1 or not r2: # MODIFIED
        print(f"❌ Skipping {sample} (paired): missing file(s)")
        continue

    output_dir = os.path.join(PE_OUTPUT_DIR, sample)
    os.makedirs(output_dir, exist_ok=True)

    print(f"🔄 Running PE STAR-Fusion for sample: {sample}")
    subprocess.run([
        "STAR-Fusion",
        "--left_fq", r1,
        "--right_fq", r2,
        "--genome_lib_dir", GENOME_LIB_DIR,
        "-O", output_dir,
        "--CPU", "20",
        #"--FusionInspector", "validate",
        "--examine_coding_effect",
        #"--denovo_reconstruct"
    ])

# Run Single-end STAR-Fusion
for sample in single_samples:
    # MODIFIED: Find single-end file
    patterns = [f"{sample}_val.fq", f"{sample}.fq", f"{sample}.fa", f"{sample}.fasta"]
    fq = find_read_file(patterns)

    if not fq: # MODIFIED
        print(f"❌ Skipping {sample} (single): missing file")
        continue

    output_dir = os.path.join(SE_OUTPUT_DIR, sample)
    os.makedirs(output_dir, exist_ok=True)

    print(f"🔄 Running SE STAR-Fusion for sample: {sample}")
    subprocess.run([
        "STAR-Fusion",
        "--left_fq", fq,
        "--genome_lib_dir", GENOME_LIB_DIR,
        "--CPU", "20",
        "-O", output_dir,
        "--examine_coding_effect",
        #"--denovo_reconstruct"
    ])

print("STAR-Fusion completed for all samples.")


# --------------------------- PART 1: Process Star-Fusion outputs --------------------------- #
print("Process Star-Fusion outputs")

# Define source and destination directories
src_dirs = ['se_output', 'pe_output']
dest_dir = 'Fusion_output'

# Create the destination directory if it doesn't exist
os.makedirs(dest_dir, exist_ok=True)

# Loop through both source directories
for src in src_dirs:
    if not os.path.exists(src):
        print(f"Source directory not found: {src}")
        continue

    for subdir in os.listdir(src):
        src_path = os.path.join(src, subdir)
        dest_path = os.path.join(dest_dir, subdir)

        if os.path.isdir(src_path):
            print(f"Moving {src_path} → {dest_path}")
            shutil.move(src_path, dest_path)
        else:
            print(f"Skipping non-directory: {src_path}")

# Change working directory to Fusion_output/
fusion_dir = os.path.join(os.getcwd(), "Fusion_output")
os.chdir(fusion_dir)

print(f"Working directory set to: {os.getcwd()}")

# Process Star-Fusion files
current_dir = os.getcwd()
dataframes = []

for dir_name in os.listdir(current_dir):
    dir_path = os.path.join(current_dir, dir_name)
    if os.path.isdir(dir_path):
        input_file_path = os.path.join(dir_path, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
        output_file_path = os.path.join(dir_path, "updated_star-fusion.fusion_predictions.abridged.coding_effect.tsv")

        if os.path.exists(input_file_path):
            df = pd.read_csv(input_file_path, sep="\t", engine="python")
            df.insert(0, "Sample_ID", dir_name)
            df.to_csv(output_file_path, sep="\t", index=False)
            print(f"Updated file created: {output_file_path}")
            dataframes.append(df)
        else:
            print(f"File not found: {input_file_path}")

if dataframes:
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Clean and process the merged dataframe
    columns_to_remove = [
        "est_J", "est_S", "SpliceType", "LargeAnchorSupport", "LeftBreakEntropy",
        "RightBreakEntropy", "annots", "CDS_LEFT_RANGE", "CDS_RIGHT_RANGE", "FUSION_MODEL",
        "FUSION_CDS", "FUSION_TRANSL", "PFAM_LEFT", "PFAM_RIGHT"
    ]
    merged_df = merged_df.drop([col for col in columns_to_remove if col in merged_df.columns], axis=1, errors='ignore')
    
    # Process breakpoints and genes
    merged_df[['Chromosome1', 'LeftBreakpoint_Pos', 'LeftStrand']] = merged_df['LeftBreakpoint'].str.split(':', expand=True)
    merged_df[['Chromosome2', 'RightBreakpoint_Pos', 'RightStrand']] = merged_df['RightBreakpoint'].str.split(':', expand=True)
    merged_df = merged_df.drop(['LeftBreakpoint', 'RightBreakpoint'], axis=1)
    merged_df["LeftGene"] = merged_df["LeftGene"].str.extract(r"([^:^]+)$")
    merged_df["RightGene"] = merged_df["RightGene"].str.extract(r"([^:^]+)$")
    
    if '#FusionName' in merged_df.columns:
        merged_df['#FusionName'] = merged_df['#FusionName'].str.replace('--', '->', regex=False)
    
    merged_df['Splice_Site'] = merged_df['LeftBreakDinuc'] + '-' + merged_df['RightBreakDinuc']
    merged_df = merged_df.drop(['LeftBreakDinuc', 'RightBreakDinuc'], axis=1)
    
    # ADDED: Standardize fusion type labels before renaming column
    if 'PROT_FUSION_TYPE' in merged_df.columns:
        merged_df['PROT_FUSION_TYPE'] = merged_df['PROT_FUSION_TYPE'].replace({
            'INFRAME': 'InFrame',
            'FRAMESHIFT': 'FrameShift',
            '.': 'NF'
        })
    
    merged_df.rename(columns={
        'LeftGene': '5_geneid',
        'RightGene': '3_geneid',
        'PROT_FUSION_TYPE': 'Splice_Pattern',
        'RightBreakpoint_Pos': 'RightBreakpoint',
        'LeftBreakpoint_Pos': 'LeftBreakpoint'
    }, inplace=True)
    
    merged_df['Total_Count_(SC+RC)'] = merged_df['JunctionReadCount'] + merged_df['SpanningFragCount']
    merged_df = merged_df.drop(['JunctionReadCount', 'SpanningFragCount'], axis=1)
    
    # MODIFIED: Update splice pattern classification based on new rules
    canonical_sites = {'GT-AG', 'GC-AG', 'AT-AC'}
    merged_df['Splice_Pattern_Class'] = merged_df['Splice_Site'].apply(
        lambda x: 'CanonicalPattern' if x in canonical_sites else 'NonCanonicalPattern'
    )
    
    merged_output_path = os.path.join(current_dir, "1_initial_merged_predictions.tsv")
    merged_df.to_csv(merged_output_path, sep="\t", index=False)
    print(f"Initial merged file created: {merged_output_path}")

# --------------------------- PART 2: Process GTF to extract longest transcript per gene --------------------------- #
#gtf_file = "Arabidopsis_thaliana.TAIR10.51.gtf"
columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=columns, dtype=str)
gtf_df = gtf_df[gtf_df["feature"] == "exon"]

gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(r'transcript_id "([^"]+)"')
gtf_df["exon_number"] = gtf_df["attribute"].str.extract(r'exon_number "(\d+)"').astype(float)

max_exon_per_gene = gtf_df.groupby("gene_id")["exon_number"].max().reset_index()
result = gtf_df.merge(max_exon_per_gene, on=["gene_id", "exon_number"], how="inner")
result = result[["gene_id", "transcript_id", "exon_number"]].drop_duplicates()

output_file = "2_filtered_gtf.tsv"
result.to_csv(output_file, sep="\t", index=False, header=["Gene_ID", "Transcript_ID", "Exon_Count"])
print(f"Filtered GTF saved as {output_file}")

# --------------------------- PART 3: Fill missing transcript IDs --------------------------- #
fusion_df = pd.read_csv("1_initial_merged_predictions.tsv", sep="\t")
filtered_gtf = pd.read_csv("2_filtered_gtf.tsv", sep="\t")

gene_to_transcripts = filtered_gtf.groupby("Gene_ID")["Transcript_ID"].apply(list).to_dict()

def fill_transcript_ids(row):
    left_transcripts = gene_to_transcripts.get(row["5_geneid"], ["."]) if row["CDS_LEFT_ID"] == "." else [row["CDS_LEFT_ID"]]
    right_transcripts = gene_to_transcripts.get(row["3_geneid"], ["."]) if row["CDS_RIGHT_ID"] == "." else [row["CDS_RIGHT_ID"]]
    
    expanded_rows = []
    for left in left_transcripts:
        for right in right_transcripts:
            new_row = row.copy()
            new_row["CDS_LEFT_ID"] = left
            new_row["CDS_RIGHT_ID"] = right
            expanded_rows.append(new_row)
    return expanded_rows

expanded_data = fusion_df.apply(fill_transcript_ids, axis=1)
expanded_df = pd.DataFrame([item for sublist in expanded_data for item in sublist])

expanded_df["CDS_LEFT_ID"] = expanded_df["CDS_LEFT_ID"].str.replace("transcript:", "", regex=False)
expanded_df["CDS_RIGHT_ID"] = expanded_df["CDS_RIGHT_ID"].str.replace("transcript:", "", regex=False)

expanded_df.to_csv("3_merged_filled.tsv", sep="\t", index=False)
print("Temporary expanded DataFrame saved as '3_merged_filled.tsv'.")

# --------------------------- PART 4: Add exon information --------------------------- #
def parse_gtf(gtf_file):
    exon_map = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[2] == "exon":
                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                attr = cols[8]
                transcript_match = re.search(r'transcript_id "([^"]+)"', attr)
                exon_num_match = re.search(r'exon_number "([^"]+)"', attr)
                if transcript_match and exon_num_match:
                    transcript_id = transcript_match.group(1)
                    exon_number = exon_num_match.group(1)
                    exon_map.setdefault(transcript_id, []).append({
                        'start': start,
                        'end': end,
                        'exon_number': exon_number,
                        'chrom': chrom,
                        'strand': strand
                    })
    return exon_map

print("Parsing GTF file for exon information...")
exon_data = parse_gtf(gtf_file)
print("GTF parsing completed.\n")

def find_exon(transcript_id, breakpoint, chrom):
    exons = exon_data.get(transcript_id, [])
    for exon in exons:
        if exon['chrom'] == chrom and exon['start'] <= breakpoint <= exon['end']:
            return exon['exon_number']
    return "NF"

print("Adding exon information to fusion data...")
expanded_df['Left_Exon'] = expanded_df.apply(
    lambda row: find_exon(row['CDS_LEFT_ID'], int(row['LeftBreakpoint']), str(row['Chromosome1'])),
    axis=1
)

expanded_df['Right_Exon'] = expanded_df.apply(
    lambda row: find_exon(row['CDS_RIGHT_ID'], int(row['RightBreakpoint']), str(row['Chromosome2'])),
    axis=1
)

final_output_path = os.path.join(current_dir, "4_final_merged_predictions.tsv")
expanded_df.to_csv(final_output_path, sep="\t", index=False)
print(f"Final output saved as: {final_output_path}")

if os.path.exists("3_merged_filled.tsv"):
    os.remove("3_merged_filled.tsv")
    print("Temporary file '3_merged_filled.tsv' removed.")

print("All processing completed successfully.")

# --------------------------- PART 5: Annotate Fusion Features --------------------------- #
df = pd.read_csv("4_final_merged_predictions.tsv", sep="\t")

# Initialize new columns
chromosome_feature = []
strand_classification = []
reciprocal_fusion_detected = []
fusion_dict = {}

# Step 1: Determine chromosomal/strand relationships & collect fusions
for _, row in df.iterrows():
    chromosome_feature.append("Intrachromosomal" if row['Chromosome1'] == row['Chromosome2'] else "Interchromosomal")
    strand_classification.append("Yes" if row['LeftStrand'] == row['RightStrand'] else "No")
    
    sample_id = row['Sample_ID']
    fusion_gene = row['#FusionName']
    if sample_id not in fusion_dict:
        fusion_dict[sample_id] = set()
    fusion_dict[sample_id].add(fusion_gene)

# Step 2: Check for reciprocal fusion presence
for _, row in df.iterrows():
    sample_id = row['Sample_ID']
    fusion_gene = row['#FusionName']
    reverse_fusion = "->".join(fusion_gene.split("->")[::-1])
    reciprocal_fusion_detected.append("Yes" if reverse_fusion in fusion_dict.get(sample_id, set()) else "No")

# Add results to DataFrame
df['Chromosome_Feature'] = chromosome_feature
df['Same_Strand'] = strand_classification
df['Reciprocal_Fusion'] = reciprocal_fusion_detected

df.to_csv("5_fusion_annotated.tsv", sep="\t", index=False)
print("✅ Step 1: Annotated fusion file saved as '5_fusion_annotated.tsv'")

# --------------------------- PART 6: Extract Gene Info from GTF --------------------------- #
print("⏳ Extracting gene info from GTF...")

#gtf_file = "../Arabidopsis_thaliana.TAIR10.51.gtf"
output_file = "6_gene_info.tsv"

gene_data = []

with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attributes = parts
        attr_dict = {}
        for attr in attributes.split(';'):
            if attr.strip():
                key_value = attr.strip().replace('"', '').split(' ')
                if len(key_value) >= 2:
                    attr_dict[key_value[0]] = ' '.join(key_value[1:])
        gene_id = attr_dict.get('gene_id', 'NA')
        gene_name = attr_dict.get('gene_name', 'NA')
        transcript_id = attr_dict.get('transcript_id', 'NA')

        if feature == "gene":
            gene_data.append([chrom, start, end, strand, gene_id, gene_name])

# Save gene-level data
pd.DataFrame(gene_data, columns=['chromosome', 'start', 'end', 'strand', 'gene_id', 'gene_name']) \
    .to_csv(output_file, sep='\t', index=False)
print(f"✅ Step 2: Gene data extracted to {output_file}")

# --------------------------- PART 7: Merge Gene Info with Fusion Data --------------------------- #
print("⏳ Merging gene info with fusion data...")

temp_df = pd.read_csv("5_fusion_annotated.tsv", sep="\t")
temp_df[["5_geneid", "3_geneid"]] = temp_df["#FusionName"].str.split("->", expand=True)

gene_info_df = pd.read_csv("6_gene_info.tsv", sep="\t")
gene_info_df.rename(columns={"gene_id": "geneid", "start": "gene_start", "end": "gene_end"}, inplace=True)
gene_info_df["gene_length"] = gene_info_df["gene_end"].astype(int) - gene_info_df["gene_start"].astype(int) + 1
gene_info_df = gene_info_df[["geneid", "gene_start", "gene_end", "gene_length"]]

# Merge to get 5' gene info
temp_df = temp_df.merge(gene_info_df, left_on="5_geneid", right_on="geneid", how="left")
temp_df.rename(columns={"gene_start": "5_gene_start", "gene_end": "5_gene_end", "gene_length": "5_gene_length"}, inplace=True)
temp_df.drop(columns=["geneid"], inplace=True)

# 3' gene
temp_df = temp_df.merge(gene_info_df, left_on="3_geneid", right_on="geneid", how="left")
temp_df.rename(columns={"gene_start": "3_gene_start", "gene_end": "3_gene_end", "gene_length": "3_gene_length"}, inplace=True)
temp_df.drop(columns=["geneid"], inplace=True)

# Convert columns to integers where applicable
int_columns = ['5_gene_start', '5_gene_end', '5_gene_length', '3_gene_start', '3_gene_end', '3_gene_length']
temp_df[int_columns] = temp_df[int_columns].astype('Int64')

temp_df.to_csv("7_fusion_with_gene_info.tsv", sep="\t", index=False)
print("✅ Step 3: Gene start/end/length merged and saved as '7_fusion_with_gene_info.tsv'")

# --------------------------- PART 8: Count Junctions Per Fusion Pair --------------------------- #
print("⏳ Calculating junction counts...")

df = pd.read_csv("7_fusion_with_gene_info.tsv", sep="\t")
df["fusion_pair"] = df["5_geneid"].astype(str) + "->" + df["3_geneid"].astype(str)

fusion_counts = df.groupby(["5_geneid", "3_geneid"]).apply(
    lambda x: x[["LeftBreakpoint", "RightBreakpoint"]].drop_duplicates().shape[0]
).reset_index(name="alternate_junction_count")

df = df.merge(fusion_counts, how="left", on=["5_geneid", "3_geneid"])
df["alternative_junction"] = df["alternate_junction_count"].apply(lambda x: "Yes" if x > 1 else "No")

# Save temporary
df.to_csv("8_fusion_with_junction_counts.tsv", sep="\t", index=False)

# Optional formatted counts
fusion_counts["formatted"] = fusion_counts["alternate_junction_count"].astype(str) + " " + fusion_counts["5_geneid"].astype(str) + "->" + fusion_counts["3_geneid"].astype(str)
fusion_counts[["formatted"]].to_csv("9_fusion_counts_formatted.tsv", index=False, header=False)

print(f"✅ Step 4: Junction counts added and saved as '8_fusion_with_junction_counts.tsv'")

# --------------------------- PART 9: Breakpoint Location Classification --------------------------- #
print("⏳ Classifying breakpoint locations relative to exons...")

def extract_exons(gtf_file):
    exon_map = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = parts
            if feature != 'exon':
                continue
            attr_dict = {}
            for attr in attributes.split(';'):
                if attr.strip():
                    key_value = attr.strip().replace('"', '').split(' ')
                    if len(key_value) >= 2:
                        attr_dict[key_value[0]] = ' '.join(key_value[1:])
            transcript_id = attr_dict.get('transcript_id')
            gene_id = attr_dict.get('gene_id')
            if gene_id not in exon_map:
                exon_map[gene_id] = []
            exon_map[gene_id].append({"start": int(start), "end": int(end)})
    return exon_map

def classify_breakpoint_location(breakpoint, exons):
    for exon in exons:
        if breakpoint == exon["start"]:
            return "S"
        elif breakpoint == exon["end"]:
            return "E"
        elif exon["start"] < breakpoint < exon["end"]:
            return "M"
    return "O"

# Load exon map
exon_map = extract_exons(gtf_file)

# Initialize columns for output
left_exon_location = []
right_exon_location = []

for _, row in df.iterrows():
    left_bp = row['LeftBreakpoint']
    right_bp = row['RightBreakpoint']
    left_gene = row['5_geneid']
    right_gene = row['3_geneid']

    left_loc = classify_breakpoint_location(left_bp, exon_map.get(left_gene, []))
    right_loc = classify_breakpoint_location(right_bp, exon_map.get(right_gene, []))

    left_exon_location.append(left_loc)
    right_exon_location.append(right_loc)

df['5_loc'] = left_exon_location
df['3_loc'] = right_exon_location

df.to_csv("10_fusion_with_breakpoint_locations.tsv", sep='\t', index=False)

# --------------------------- PART 10: Exon Count Integration --------------------------- #
print("⏳ Incorporating exon count information...")

# Load the main data file
main_df = pd.read_csv("10_fusion_with_breakpoint_locations.tsv", sep="\t")

# Load the exon count file
exon_df = pd.read_csv("2_filtered_gtf.tsv", sep="\t", names=["Gene_ID", "Transcript_ID", "Exon_Count"])

# Aggregate exon counts for each Gene_ID
exon_counts = exon_df.groupby("Gene_ID")["Exon_Count"].unique().to_dict()

# Function to get exon count for a given gene_id
def get_exon_count(gene_id):
    return exon_counts.get(gene_id, [])

# Expand rows if exon counts are non-unique
data_expanded = []
for _, row in main_df.iterrows():
    exon5 = get_exon_count(row["5_geneid"])
    exon3 = get_exon_count(row["3_geneid"])

    # Handle cases where exon counts might be empty
    if not exon5: exon5 = [""]
    if not exon3: exon3 = [""]

    for e5 in exon5:
        for e3 in exon3:
            new_row = row.copy()
            new_row["exon_count5"] = int(float(e5)) if pd.notna(e5) and e5 != "" else ""
            new_row["exon_count3"] = int(float(e3)) if pd.notna(e3) and e3 != "" else ""
            data_expanded.append(new_row)

# Convert back to DataFrame
expanded_df = pd.DataFrame(data_expanded)

# Save the final output
expanded_df.to_csv("11_final_fusion_annotated.tsv", sep='\t', index=False)
print("✅ FINAL STEP: Exon count added. Output saved as '11_final_fusion_annotated.tsv'")



###total mapped reads. 
input_file = '11_final_fusion_annotated.tsv'
output_file = '12_annotated_with_reads.tsv'

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    fieldnames = reader.fieldnames + ['Total_Mapped_Reads']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    for row in reader:
        ffpm = row['FFPM']
        total_count = row['Total_Count_(SC+RC)']

        if ffpm not in ['NA', 'NF', ''] and total_count not in ['NA', 'NF', '']:
            try:
                ffpm_value = float(ffpm)
                count_value = int(total_count)
                if ffpm_value > 0:
                    total_reads = (count_value * 1_000_000) / ffpm_value
                    row['Total_Mapped_Reads'] = int(total_reads)
                else:
                    row['Total_Mapped_Reads'] = 'NF'
            except ValueError:
                row['Total_Mapped_Reads'] = 'NF'
        else:
            row['Total_Mapped_Reads'] = 'NF'

        writer.writerow(row)

###########################################################
# MODIFIED SECTION
###########################################################

###
input_file = "12_annotated_with_reads.tsv"
output_file = "fts_features.csv" # MODIFIED: Changed extension to .csv

# Read the input file (still reads the .tsv file)
df = pd.read_csv(input_file, sep="\t")

# Columns to remove
columns_to_remove = [
    "#FusionName",
    "3_geneid",
    "CDS_LEFT_ID",
    "CDS_RIGHT_ID"
]

# Drop the columns (only if they exist, to avoid errors)
df = df.drop(columns=[col for col in columns_to_remove if col in df.columns])

# Save to output file (MODIFIED: saves as .csv)
df.to_csv(output_file, sep=",", index=False)

print(f"Processed file saved as: {output_file}")

##Temp Dir
output = "fts_features.csv" # MODIFIED: Changed to .csv to protect it from cleanup
# Temporary directory
temp_dir = "temp_dir"
os.makedirs(temp_dir, exist_ok=True)

# Loop through all files in current directory
for f in os.listdir("."):
    # Match files that start with number(s) + underscore
    if re.match(r"^\d+_.*", f) and f != output:
        shutil.move(f, os.path.join(temp_dir, f))
        print(f"Moved {f} → {temp_dir}/")

print("Cleanup done.")
