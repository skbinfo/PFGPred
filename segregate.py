#!/usr/bin/env python3
import pandas as pd
import csv
import argparse
import os
import sys

# ----------------------------------------------------------
# ARGPARSE — Add command-line input support
# ----------------------------------------------------------
parser = argparse.ArgumentParser(description="Fusion processing pipeline (produce final merged CSV with labels)")

parser.add_argument("--summary", required=True,
                    help="fusion_summary_table_breakpoint_supporting_reads.txt")

parser.add_argument("--annot", required=True,
                    help="new_fusion_search_table.txt")

parser.add_argument("--feat", required=True,
                    help="fts_features.csv")

args = parser.parse_args()

summary_file = args.summary
annot_file = args.annot
feature_file = args.feat

# Helper to remove intermediate files safely
def safe_remove(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception as e:
        print(f"Warning: could not remove {path}: {e}", file=sys.stderr)

# ----------------------------------------------------------
# STEP 1 — Remove duplicates
# ----------------------------------------------------------
print("STEP 1: Removing duplicates...")

df = pd.read_csv(summary_file, sep='\t')
unique_df = df.drop_duplicates(
    subset=['threeprime_breakpoint_coordinate', 'fiveprime_breakpoint_coordinate']
)
unique_df = unique_df.sort_values(
    by=['threeprime_breakpoint_coordinate', 'fiveprime_breakpoint_coordinate']
)
unique_out = 'unique_fusion_summary.tsv'
unique_df.to_csv(unique_out, sep='\t', index=False)

print(f" → {unique_out} created.")

# ----------------------------------------------------------
# STEP 2 — Split chr:bpt into columns
# ----------------------------------------------------------
print("STEP 2: Formatting to split chr:bpt...")

formatted_out = 'formatted_fusion_summary.tsv'
with open(unique_out, 'r') as csvfile, \
     open(formatted_out, 'w', newline='') as outfile:

    reader = csv.DictReader(csvfile, delimiter='\t')
    fieldnames = ['fusion_id', 'softclip_ID', 'chr1', 'bpt1', 'chr2', 'bpt2', 'softclip_sequence']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    for row in reader:
        try:
            chr1, bpt1 = row['threeprime_breakpoint_coordinate'].split(':')
            chr2, bpt2 = row['fiveprime_breakpoint_coordinate'].split(':')
        except Exception:
            # If parsing fails, write empty coords and continue
            chr1, bpt1, chr2, bpt2 = '', '', '', ''

        writer.writerow({
            'fusion_id': row.get('fusion_id', ''),
            'softclip_ID': row.get('softclip_ID', ''),
            'chr1': chr1,
            'bpt1': bpt1,
            'chr2': chr2,
            'bpt2': bpt2,
            'softclip_sequence': row.get('softclip_sequence', '')
        })

print(f" → {formatted_out} created.")

# ----------------------------------------------------------
# STEP 3 — Add annotations
# ----------------------------------------------------------
print("STEP 3: Adding annotations...")

summary_df = pd.read_csv(formatted_out, sep="\t")
fusion_table = pd.read_csv(annot_file, sep="\t")

columns_to_add = [
    "fusion_id", "sample_id", "fiveprimepartner", "threeprimepartner",
    "fiveprime_chr", "threeprime_chr",
    "fiveprime_junction", "threeprime_junction",
    "fiveprime_strand", "threeprime_strand"
]

# Ensure columns exist in annotation table
missing_cols = [c for c in columns_to_add if c not in fusion_table.columns]
if missing_cols:
    print(f"ERROR: annotation file is missing required columns: {missing_cols}", file=sys.stderr)
    sys.exit(1)

fusion_table = fusion_table[columns_to_add]
merged_df = summary_df.merge(fusion_table, on="fusion_id", how="left")
annotated_out = "formatted_fusion_summary_annotated.tsv"
merged_df.to_csv(annotated_out, sep="\t", index=False)

print(f" → {annotated_out} created.")

# ----------------------------------------------------------
# STEP 4 — Add gene start / end based on ±5bp matching
# ----------------------------------------------------------
print("STEP 4: Adding gene start/end...")

df_filt = pd.read_csv(annotated_out, sep="\t", dtype=str)
df_annot = pd.read_csv(annot_file, sep="\t", dtype=str)

# Keep only rows with valid strand info
df_filt = df_filt[df_filt['fiveprime_strand'].isin(['+', '-']) &
                  df_filt['threeprime_strand'].isin(['+', '-'])].copy()

match_cols = [
    "fiveprimepartner", "threeprimepartner",
    "fiveprime_chr", "threeprime_chr",
    "fiveprime_strand", "threeprime_strand"
]

df_filt_match = df_filt[match_cols + ["fiveprime_junction", "threeprime_junction"]].copy()
df_annot_match = df_annot[match_cols + [
    "fiveprime_junction", "threeprime_junction",
    "fiveprime_search_start", "fiveprime_search_end",
    "threeprime_search_start", "threeprime_search_end"
]].copy()

# Normalize keys and coerce types with safe conversions
for col in match_cols:
    df_filt_match[col] = df_filt_match[col].astype(str).str.strip().str.lower()
    df_annot_match[col] = df_annot_match[col].astype(str).str.strip().str.lower()

def safe_int(x):
    try:
        return int(float(x))
    except Exception:
        return None

df_filt_match["fiveprime_junction"] = df_filt_match["fiveprime_junction"].apply(safe_int)
df_filt_match["threeprime_junction"] = df_filt_match["threeprime_junction"].apply(safe_int)
df_annot_match["fiveprime_junction"] = df_annot_match["fiveprime_junction"].apply(safe_int)
df_annot_match["threeprime_junction"] = df_annot_match["threeprime_junction"].apply(safe_int)

# Drop rows with missing junctions in either side for matching
df_filt_match = df_filt_match.dropna(subset=["fiveprime_junction", "threeprime_junction"]).copy()
df_annot_match = df_annot_match.dropna(subset=["fiveprime_junction", "threeprime_junction"]).copy()

df_filt_match["match_key"] = df_filt_match[match_cols].agg("__".join, axis=1)
df_annot_match["match_key"] = df_annot_match[match_cols].agg("__".join, axis=1)

# Build annot dictionary
annot_dict = {}
for _, row in df_annot_match.iterrows():
    annot_dict.setdefault(row["match_key"], []).append(row)

final_rows = []

# We will iterate over df_filt rows but need to map index to df_filt_match indices.
# Build mapping from df_filt index to df_filt_match row (by matching keys and junctions)
df_filt_match_indexed = df_filt_match.reset_index()

for _, filtrow in df_filt_match_indexed.iterrows():
    original_idx = filtrow['index']
    row = df_filt.loc[original_idx].to_dict()
    key = filtrow["match_key"]
    candidates = annot_dict.get(key, [])

    row_out = row.copy()
    row_out["5_gene_start"] = ""
    row_out["5_gene_end"] = ""
    row_out["3_gene_start"] = ""
    row_out["3_gene_end"] = ""

    for cand in candidates:
        cand_fj = cand["fiveprime_junction"]
        cand_tj = cand["threeprime_junction"]
        if cand_fj is None or cand_tj is None:
            continue
        if abs(cand_fj - int(filtrow["fiveprime_junction"])) <= 5 and \
           abs(cand_tj - int(filtrow["threeprime_junction"])) <= 5:
            row_out["5_gene_start"] = cand.get("fiveprime_search_start", "")
            row_out["5_gene_end"] = cand.get("fiveprime_search_end", "")
            row_out["3_gene_start"] = cand.get("threeprime_search_start", "")
            row_out["3_gene_end"] = cand.get("threeprime_search_end", "")
            break

    final_rows.append(row_out)

df_final = pd.DataFrame(final_rows)
merged_coords_out = "merged_fusion_with_coords.tsv"
df_final.to_csv(merged_coords_out, sep="\t", index=False)

print(f" → {merged_coords_out} created.")

# ----------------------------------------------------------
# STEP 5 — Separate overlapping vs non-overlapping
# ----------------------------------------------------------
print("STEP 5: Checking overlap...")

def overlap(s1, e1, s2, e2):
    try:
        s1f = float(s1)
        e1f = float(e1)
        s2f = float(s2)
        e2f = float(e2)
    except Exception:
        # If conversion fails, treat as non-overlapping (conservative)
        return False
    return not (e1f < s2f or s1f > e2f)

input_file = merged_coords_out
overlap_file = "overlapping_genes.csv"
non_overlap_file = "non_overlapping_genes.csv"

with open(input_file, 'r') as infile, \
     open(overlap_file, 'w', newline='') as overlap_outfile, \
     open(non_overlap_file, 'w', newline='') as non_overlap_outfile:

    reader = csv.DictReader(infile, delimiter='\t')
    header = reader.fieldnames

    overlap_writer = csv.writer(overlap_outfile)
    non_overlap_writer = csv.writer(non_overlap_outfile)

    overlap_writer.writerow(header)
    non_overlap_writer.writerow(header)

    for row in reader:
        s1 = row.get('5_gene_start', '')
        e1 = row.get('5_gene_end', '')
        s2 = row.get('3_gene_start', '')
        e2 = row.get('3_gene_end', '')

        if overlap(s1, e1, s2, e2):
            overlap_writer.writerow([row.get(h, '') for h in header])
        else:
            non_overlap_writer.writerow([row.get(h, '') for h in header])

print(f" → {overlap_file} and {non_overlap_file} created.")

# ----------------------------------------------------------
# STEP 6 — Compare with feature table and produce final merged CSV with labels
# ----------------------------------------------------------
print("STEP 6: Matching against feature table and creating final merged CSV with labels...")

validated_df = pd.read_csv(non_overlap_file)
temp_df = pd.read_csv(feature_file)

column_mapping = {
    "fiveprimepartner": "5_geneid",
    "threeprimepartner": "3_geneid",
    "fiveprime_chr": "Chromosome1",
    "fiveprime_junction": "LeftBreakpoint",
    "fiveprime_strand": "LeftStrand",
    "threeprime_chr": "Chromosome2",
    "threeprime_junction": "RightBreakpoint",
    "threeprime_strand": "RightStrand"
}

validated_renamed = validated_df.rename(columns=column_mapping)

compare_columns = [
    "5_geneid", "3_geneid",
    "LeftBreakpoint", "LeftStrand",
    "RightBreakpoint", "RightStrand"
]

# Coerce to string to allow reliable merging
for col in compare_columns:
    validated_renamed[col] = validated_renamed[col].astype(str)
    temp_df[col] = temp_df[col].astype(str)

# Merge so every row in temp_df is present and get indicator to mark matches
merged = temp_df.merge(validated_renamed[compare_columns],
                       on=compare_columns,
                       how='left',
                       indicator=True)

# Add 'label' column: 1 for matches (aligned to WGS validated data), 0 for non-matches
merged['label'] = merged['_merge'].apply(lambda x: 1 if x == 'both' else 0)

# Drop internal merge indicator and deduplicate
final_df = merged.drop(columns=['_merge']).drop_duplicates().reset_index(drop=True)

if "Sample_ID" in final_df.columns:
    final_df = final_df.drop(columns=["Sample_ID"])

final_out = "final_with_labels.csv"
final_df.to_csv(final_out, index=False)

# Count positives and negatives
pos_count = int(final_df['label'].sum())
neg_count = int((final_df['label'] == 0).sum())

# ----------------------------------------------------------
# Cleanup intermediary files
# ----------------------------------------------------------
print("Cleaning up intermediary files...")
intermediate_files = [
    unique_out,
    formatted_out,
    annotated_out,
    merged_coords_out,
    overlap_file,
    non_overlap_file,
    "post.tsv",
    "neg.tsv",
    "overlapping_genes.csv",
    "non_overlapping_genes.csv"
]
for f in intermediate_files:
    safe_remove(f)

print("\n PIPELINE COMPLETE!")
print(f" → Final merged CSV with labels: {final_out}")
print(f" → Positive (label=1): {pos_count} rows")
print(f" → Negative (label=0): {neg_count} rows")
