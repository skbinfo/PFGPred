#!/usr/bin/env python3
"""
fusion_pipeline_resume.py

Usage examples:
  # Resume from 1_initial_merged_predictions.tsv (PART 3 onward)
  python3 fusion_pipeline_resume.py --start step2 --merged 1_initial_merged_predictions.tsv --gtf path/to/annotation.gtf

  # Run full pipeline (requires config.yaml, single.txt, paired.txt, FASTQ files)
  python3 fusion_pipeline_resume.py --start step1 --config config.yaml
"""
import argparse
import os
import re
import shutil
import subprocess
import csv
import yaml
import pandas as pd
from typing import Dict, List, Any

# -------------------------
# ARGPARSE
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Fusion pipeline (resumable).")
    p.add_argument("--start", choices=["step1","step2"], default="step1",
                   help="Which step to start from. step1 runs everything (including STAR-Fusion). "
                        "step2 runs GTF processing & downstream.")
    p.add_argument("--config", help="Path to config.yaml (required for step1)")
    p.add_argument("--merged", help="Path to 1_initial_merged_predictions.tsv (required for step3)")
    p.add_argument("--gtf", help="Path to genome GTF (required for step1/2)")
    p.add_argument("--outdir", default="Fusion_output", help="Output directory")
    return p.parse_args()

args = parse_args()

# -------------------------
# SIMPLE HELPERS
# -------------------------
def safe_read_tsv(path, **kwargs):
    return pd.read_csv(path, sep=kwargs.pop("sep", "\t"), dtype=str, **kwargs)

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)

# -------------------------
# STEP 1: (Optional) STAR-Fusion and generate 1_initial_merged_predictions.tsv
# This implements the original behavior. If you already have merged file, start at step2.
# -------------------------
def step1_run_star_and_merge(config_path: str, out_root: str):
    """
    Runs read-length extraction, optionally runs STAR-Fusion if FASTQs present,
    collects outputs, and builds 1_initial_merged_predictions.tsv in out_root.
    NOTE: Running STAR-Fusion requires STAR-Fusion binary in PATH.
    """
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    fasta_file = cfg.get("fasta_file")
    gtf_file = cfg.get("gtf_file")
    GENOME_LIB_DIR = cfg.get("genome_lib_dir")
    single_file = cfg.get("single_list", "single.txt")
    paired_file = cfg.get("paired_list", "paired.txt")

    ensure_dir(out_root)
    reads_tsv = os.path.join(out_root, "reads.tsv")

    # helper functions
    def get_read_length(fq_file):
        try:
            with open(fq_file, "r") as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        return len(line.strip())
        except Exception:
            return "NF"

    def find_read_file(patterns):
        for pat in patterns:
            if os.path.isfile(pat):
                return pat
        return None

    # collect read lengths
    results = []
    if os.path.exists(single_file):
        with open(single_file) as fh:
            for line in fh:
                s = line.strip()
                if not s: 
                    continue
                patterns = [f"{s}_val.fq", f"{s}.fq", f"{s}.fa", f"{s}.fasta"]
                fq = find_read_file(patterns)
                results.append((get_read_length(fq) if fq else "NF", s))
    if os.path.exists(paired_file):
        with open(paired_file) as fh:
            for line in fh:
                s = line.strip()
                if not s: 
                    continue
                patterns = [f"{s}_1_val_1.fq", f"{s}_1.fq", f"{s}_R1.fq"]
                r1 = find_read_file(patterns)
                results.append((get_read_length(r1) if r1 else "NF", s))

    # write reads.tsv
    with open(reads_tsv, "w") as out:
        out.write("Total_Mapped_Reads\tsample\n")
        for rl, s in results:
            out.write(f"{rl}\t{s}\n")
    print(f"Saved read lengths to {reads_tsv}")

    # run STAR-Fusion if FASTQs found and genome lib present
    # This block mimics the original script but is optional and user must have STAR-Fusion installed.
    if GENOME_LIB_DIR and (os.path.exists(paired_file) or os.path.exists(single_file)):
        ensure_dir("pe_output"); ensure_dir("se_output")
        # paired
        if os.path.exists(paired_file):
            with open(paired_file) as fh:
                for line in fh:
                    s = line.strip()
                    if not s: continue
                    r1 = find_read_file([f"{s}_1_val_1.fq", f"{s}_1.fq", f"{s}_R1.fq"])
                    r2 = find_read_file([f"{s}_2_val_2.fq", f"{s}_2.fq", f"{s}_R2.fq"])
                    if not r1 or not r2:
                        print(f"Skipping paired sample {s}: missing R1/R2")
                        continue
                    outdir = os.path.join("pe_output", s)
                    ensure_dir(outdir)
                    print(f"Running STAR-Fusion for paired sample {s} → {outdir}")
                    subprocess.run([
                        "STAR-Fusion", "--left_fq", r1, "--right_fq", r2,
                        "--genome_lib_dir", GENOME_LIB_DIR, "-O", outdir, "--CPU", "4",
                        "--examine_coding_effect"
                    ])
        # single
        if os.path.exists(single_file):
            with open(single_file) as fh:
                for line in fh:
                    s = line.strip()
                    if not s: continue
                    fq = find_read_file([f"{s}_val.fq", f"{s}.fq"])
                    if not fq:
                        print(f"Skipping single sample {s}: missing fq")
                        continue
                    outdir = os.path.join("se_output", s)
                    ensure_dir(outdir)
                    print(f"Running STAR-Fusion for single sample {s} → {outdir}")
                    subprocess.run([
                        "STAR-Fusion", "--left_fq", fq, "--genome_lib_dir", GENOME_LIB_DIR,
                        "-O", outdir, "--CPU", "4", "--examine_coding_effect"
                    ])

    # collect STAR-Fusion outputs into Fusion_output and create merged_df as earlier
    src_dirs = ['se_output', 'pe_output']
    dest_dir = out_root
    ensure_dir(dest_dir)
    for src in src_dirs:
        if not os.path.exists(src):
            continue
        for sub in os.listdir(src):
            sp = os.path.join(src, sub)
            dp = os.path.join(dest_dir, sub)
            if os.path.isdir(sp):
                if os.path.exists(dp):
                    shutil.rmtree(dp)
                shutil.move(sp, dp)

    # run the merging of star-fusion outputs into 1_initial_merged_predictions.tsv
    fusion_root = os.path.abspath(out_root)
    dataframes = []
    for name in os.listdir(fusion_root):
        d = os.path.join(fusion_root, name)
        if os.path.isdir(d):
            tsv = os.path.join(d, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
            if os.path.exists(tsv):
                df = pd.read_csv(tsv, sep="\t", engine="python", dtype=str)
                df.insert(0, "Sample_ID", name)
                dataframes.append(df)
    if dataframes:
        merged_df = pd.concat(dataframes, ignore_index=True)

        # cleanup like original script
        to_drop = [
            "est_J","est_S","SpliceType","LargeAnchorSupport","LeftBreakEntropy",
            "RightBreakEntropy","annots","CDS_LEFT_RANGE","CDS_RIGHT_RANGE",
            "FUSION_MODEL","FUSION_CDS","FUSION_TRANSL","PFAM_LEFT","PFAM_RIGHT"
        ]
        merged_df = merged_df.drop(
            [c for c in to_drop if c in merged_df.columns],
            axis=1,
            errors='ignore'
        )

        # process breakpoints
        if 'LeftBreakpoint' in merged_df.columns and 'RightBreakpoint' in merged_df.columns:
            merged_df[['Chromosome1','LeftBreakpoint_Pos','LeftStrand']] = \
                merged_df['LeftBreakpoint'].astype(str).str.split(':', expand=True)
            merged_df[['Chromosome2','RightBreakpoint_Pos','RightStrand']] = \
                merged_df['RightBreakpoint'].astype(str).str.split(':', expand=True)
            merged_df = merged_df.drop(
                ['LeftBreakpoint','RightBreakpoint'],
                axis=1,
                errors='ignore'
            )

        # Left/Right gene extraction
        if 'LeftGene' in merged_df.columns:
            merged_df["LeftGene"] = merged_df["LeftGene"].astype(str).str.extract(r"([^:^]+)$")
        if 'RightGene' in merged_df.columns:
            merged_df["RightGene"] = merged_df["RightGene"].astype(str).str.extract(r"([^:^]+)$")

        if '#FusionName' in merged_df.columns:
            merged_df['#FusionName'] = merged_df['#FusionName'].astype(str).str.replace('--','->', regex=False)

        # splice site column
        if 'LeftBreakDinuc' in merged_df.columns and 'RightBreakDinuc' in merged_df.columns:
            merged_df['Splice_Site'] = (
                merged_df['LeftBreakDinuc'].astype(str) + '-' +
                merged_df['RightBreakDinuc'].astype(str)
            )
            merged_df = merged_df.drop(
                ['LeftBreakDinuc','RightBreakDinuc'],
                axis=1,
                errors='ignore'
            )

        # **Normalize Splice_Site to uppercase**
        if 'Splice_Site' in merged_df.columns:
            merged_df['Splice_Site'] = merged_df['Splice_Site'].astype(str).str.upper()

        # standardize fusion type
        if 'PROT_FUSION_TYPE' in merged_df.columns:
            merged_df['PROT_FUSION_TYPE'] = merged_df['PROT_FUSION_TYPE'].replace(
                {'INFRAME':'InFrame','FRAMESHIFT':'FrameShift','.':'NF'}
            )

        # rename
        renamed = {}
        if 'LeftGene' in merged_df.columns:
            renamed['LeftGene'] = '5_geneid'
        if 'RightGene' in merged_df.columns:
            renamed['RightGene'] = '3_geneid'
        if 'PROT_FUSION_TYPE' in merged_df.columns:
            renamed['PROT_FUSION_TYPE'] = 'Splice_Pattern'
        if 'RightBreakpoint_Pos' in merged_df.columns:
            renamed['RightBreakpoint_Pos'] = 'RightBreakpoint'
        if 'LeftBreakpoint_Pos' in merged_df.columns:
            renamed['LeftBreakpoint_Pos'] = 'LeftBreakpoint'
        merged_df = merged_df.rename(columns=renamed)

        # compute counts if present
        if 'JunctionReadCount' in merged_df.columns and 'SpanningFragCount' in merged_df.columns:
            merged_df['Total_Count_(SC+RC)'] = (
                merged_df['JunctionReadCount'].astype(float).fillna(0) +
                merged_df['SpanningFragCount'].astype(float).fillna(0)
            )
            merged_df = merged_df.drop(
                ['JunctionReadCount','SpanningFragCount'],
                axis=1,
                errors='ignore'
            )

        # splice pattern class (using normalized, upper-case Splice_Site)
        canonical_sites = {'GT-AG', 'GC-AG', 'AT-AC'}
        if 'Splice_Site' in merged_df.columns:
            merged_df['Splice_Pattern_Class'] = merged_df['Splice_Site'].apply(
                lambda x: 'CanonicalPattern' if x in canonical_sites else 'NonCanonicalPattern'
            )

        # save
        outpath = os.path.join(fusion_root, "1_initial_merged_predictions.tsv")
        merged_df.to_csv(outpath, sep="\t", index=False)
        print(f"Created {outpath}")
    else:
        print("No STAR-Fusion outputs found to merge.")

# -------------------------
# STEP 2: Process GTF to extract longest transcript per gene
# -------------------------
def step2_process_gtf(gtf_file: str):
    cols = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    print(f"Parsing GTF (exons) from {gtf_file} ...")
    gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=cols, dtype=str)
    gtf_df = gtf_df[gtf_df["feature"] == "exon"]
    gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(r'transcript_id "([^"]+)"')
    gtf_df["exon_number"] = gtf_df["attribute"].str.extract(r'exon_number "(\d+)"').astype(float)
    max_exon_per_gene = gtf_df.groupby("gene_id")["exon_number"].max().reset_index()
    result = gtf_df.merge(max_exon_per_gene, on=["gene_id","exon_number"], how="inner")
    result = result[["gene_id","transcript_id","exon_number"]].drop_duplicates()
    result.to_csv("2_filtered_gtf.tsv", sep="\t", index=False, header=["Gene_ID","Transcript_ID","Exon_Count"])
    print("Saved 2_filtered_gtf.tsv")

# -------------------------
# STEP 3: Fill missing transcript IDs (expects 1_initial_merged_predictions.tsv)
# -------------------------
def step3_fill_transcripts(merged_path: str):
    print(f"Loading merged predictions from {merged_path}")
    fusion_df = pd.read_csv(merged_path, sep="\t", dtype=str)
    filtered_gtf = pd.read_csv("2_filtered_gtf.tsv", sep="\t", dtype=str)
    gene_to_transcripts = filtered_gtf.groupby("Gene_ID")["Transcript_ID"].apply(list).to_dict()

    def fill_transcript_ids(row):
        left_transcripts = gene_to_transcripts.get(row.get("5_geneid"), ["."]) if str(row.get("CDS_LEFT_ID",".")) == "." else [row.get("CDS_LEFT_ID")]
        right_transcripts = gene_to_transcripts.get(row.get("3_geneid"), ["."]) if str(row.get("CDS_RIGHT_ID",".")) == "." else [row.get("CDS_RIGHT_ID")]
        out = []
        for L in left_transcripts:
            for R in right_transcripts:
                nr = row.copy()
                nr["CDS_LEFT_ID"] = L
                nr["CDS_RIGHT_ID"] = R
                out.append(nr)
        return out

    expanded = fusion_df.apply(fill_transcript_ids, axis=1)
    expanded_df = pd.DataFrame([x for sub in expanded for x in sub])
    # strip 'transcript:' prefix if present
    for col in ["CDS_LEFT_ID","CDS_RIGHT_ID"]:
        if col in expanded_df.columns:
            expanded_df[col] = expanded_df[col].astype(str).str.replace("transcript:","",regex=False)
    expanded_df.to_csv("3_merged_filled.tsv", sep="\t", index=False)
    print("Saved 3_merged_filled.tsv")

# -------------------------
# STEP 4: Add exon information (per-transcript exon coordinates)
# -------------------------
def step4_add_exons(gtf_file: str):
    def parse_gtf_to_exons(gtf):
        exon_map = {}
        with open(gtf) as fh:
            for ln in fh:
                if ln.startswith("#"): 
                    continue
                cols = ln.strip().split("\t")
                if len(cols) < 9: 
                    continue
                if cols[2] != "exon": 
                    continue
                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                attrs = cols[8]
                tidm = re.search(r'transcript_id "([^"]+)"', attrs)
                exn = re.search(r'exon_number "([^"]+)"', attrs)
                if not tidm or not exn: 
                    continue
                tid = tidm.group(1)
                exnum = exn.group(1)
                exon_map.setdefault(tid, []).append({
                    "start":start,
                    "end":end,
                    "exon_number":exnum,
                    "chrom":chrom,
                    "strand":strand
                })
        return exon_map

    print(f"Parsing GTF for exon coordinates from {gtf_file} ...")
    exon_map = parse_gtf_to_exons(gtf_file)
    df = pd.read_csv("3_merged_filled.tsv", sep="\t", dtype=str)

    def find_exon(tid, bp, chrom):
        if pd.isna(tid) or tid == ".":
            return "NF"
        exons = exon_map.get(tid, [])
        try:
            bp_int = int(bp)
        except Exception:
            return "NF"
        for e in exons:
            if e["chrom"] == str(chrom) and e["start"] <= bp_int <= e["end"]:
                return e["exon_number"]
        return "NF"

    df["Left_Exon"] = df.apply(
        lambda r: find_exon(r.get("CDS_LEFT_ID","."), r.get("LeftBreakpoint",0), r.get("Chromosome1","")),
        axis=1
    )
    df["Right_Exon"] = df.apply(
        lambda r: find_exon(r.get("CDS_RIGHT_ID","."), r.get("RightBreakpoint",0), r.get("Chromosome2","")),
        axis=1
    )
    df.to_csv("4_final_merged_predictions.tsv", sep="\t", index=False)
    print("Saved 4_final_merged_predictions.tsv")
    # remove temporary
    try:
        os.remove("3_merged_filled.tsv")
    except Exception:
        pass

# -------------------------
# STEP 5: Annotate Fusion Features (chromosome feature, reciprocal, same strand)
# -------------------------
def step5_annotate_features():
    df = pd.read_csv("4_final_merged_predictions.tsv", sep="\t", dtype=str)
    # chromosome feature
    df["Chromosome_Feature"] = df.apply(
        lambda r: "Intrachromosomal" if r.get("Chromosome1") == r.get("Chromosome2") else "Interchromosomal",
        axis=1
    )
    df["Same_Strand"] = df.apply(
        lambda r: "Yes" if r.get("LeftStrand") == r.get("RightStrand") else "No",
        axis=1
    )
    # reciprocal: collect per-sample set
    fusion_per_sample = {}
    for _, r in df.iterrows():
        s = r.get("Sample_ID","")
        fn = r.get("#FusionName","")
        fusion_per_sample.setdefault(s, set()).add(fn)
    df["Reciprocal_Fusion"] = df.apply(
        lambda r: "Yes" if "->".join(str(r.get("#FusionName","")).split("->")[::-1]) in fusion_per_sample.get(r.get("Sample_ID",""), set())
        else "No",
        axis=1
    )
    df.to_csv("5_fusion_annotated.tsv", sep="\t", index=False)
    print("Saved 5_fusion_annotated.tsv")

# -------------------------
# STEP 6: Extract gene-level info from GTF (gene start/end)
# -------------------------
def step6_extract_gene_info(gtf_file: str):
    gene_rows = []
    with open(gtf_file) as fh:
        for ln in fh:
            if ln.startswith("#"): 
                continue
            parts = ln.strip().split("\t")
            if len(parts) != 9: 
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene": 
                continue
            ad = {}
            for a in attrs.split(";"):
                if not a.strip(): 
                    continue
                kv = a.strip().replace('"','').split(" ")
                if len(kv) >= 2:
                    ad[kv[0]] = " ".join(kv[1:])
            gid = ad.get("gene_id","NA")
            gname = ad.get("gene_name","NA")
            gene_rows.append([chrom, int(start), int(end), strand, gid, gname])
    gene_df = pd.DataFrame(gene_rows, columns=["chromosome","start","end","strand","gene_id","gene_name"])
    gene_df.to_csv("6_gene_info.tsv", sep="\t", index=False)
    print("Saved 6_gene_info.tsv")

# -------------------------
# STEP 7: Merge gene info with fusion data (adds 5_gene_start, 3_gene_start, lengths)
# -------------------------
def step7_merge_gene_info():
    df = pd.read_csv("5_fusion_annotated.tsv", sep="\t", dtype=str)
    # ensure #FusionName exists and split
    if "#FusionName" not in df.columns:
        raise RuntimeError("Missing #FusionName column required to extract 5_geneid and 3_geneid.")
    df[["5_geneid","3_geneid"]] = df["#FusionName"].astype(str).str.split("->", expand=True)
    gene_df = pd.read_csv("6_gene_info.tsv", sep="\t", dtype=str)
    gene_df.rename(columns={"gene_id":"geneid","start":"gene_start","end":"gene_end"}, inplace=True)
    gene_df["gene_length"] = gene_df["gene_end"].astype(int) - gene_df["gene_start"].astype(int) + 1
    gene_df = gene_df[["geneid","gene_start","gene_end","gene_length"]]
    # merge 5'
    df = df.merge(
        gene_df,
        left_on="5_geneid",
        right_on="geneid",
        how="left"
    ).rename(
        columns={
            "gene_start":"5_gene_start",
            "gene_end":"5_gene_end",
            "gene_length":"5_gene_length"
        }
    ).drop(columns=["geneid"])
    # merge 3'
    df = df.merge(
        gene_df,
        left_on="3_geneid",
        right_on="geneid",
        how="left"
    ).rename(
        columns={
            "gene_start":"3_gene_start",
            "gene_end":"3_gene_end",
            "gene_length":"3_gene_length"
        }
    ).drop(columns=["geneid"])
    # convert to Int64 where possible
    for c in ['5_gene_start','5_gene_end','5_gene_length','3_gene_start','3_gene_end','3_gene_length']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce').astype('Int64')
    df.to_csv("7_fusion_with_gene_info.tsv", sep="\t", index=False)
    print("Saved 7_fusion_with_gene_info.tsv")

# -------------------------
# STEP 8: Count alternate junctions per fusion pair
# -------------------------
def step8_count_junctions():
    df = pd.read_csv("7_fusion_with_gene_info.tsv", sep="\t", dtype=str)
    df["fusion_pair"] = df["5_geneid"].astype(str) + "->" + df["3_geneid"].astype(str)
    fusion_counts = df.groupby(["5_geneid","3_geneid"]).apply(
        lambda x: x[["LeftBreakpoint","RightBreakpoint"]].drop_duplicates().shape[0]
    ).reset_index(name="alternate_junction_count")
    df = df.merge(fusion_counts, on=["5_geneid","3_geneid"], how="left")
    df["alternative_junction"] = df["alternate_junction_count"].apply(
        lambda x: "Yes" if int(x) > 1 else "No"
    )
    df.to_csv("8_fusion_with_junction_counts.tsv", sep="\t", index=False)
    fusion_counts["formatted"] = (
        fusion_counts["alternate_junction_count"].astype(str) + " " +
        fusion_counts["5_geneid"].astype(str) + "->" +
        fusion_counts["3_geneid"].astype(str)
    )
    fusion_counts[["formatted"]].to_csv("9_fusion_counts_formatted.tsv", index=False, header=False)
    print("Saved 8_fusion_with_junction_counts.tsv and 9_fusion_counts_formatted.tsv")

# -------------------------
# STEP 9: Classify breakpoint locations relative to exons
# -------------------------
def step9_classify_breakpoint_locations(gtf_file: str):
    # Build gene->exon lists (start,end)
    exon_map = {}
    with open(gtf_file) as fh:
        for ln in fh:
            if ln.startswith("#"): 
                continue
            parts = ln.strip().split("\t")
            if len(parts) != 9: 
                continue
            chrom,_,feature,start,end,_,strand,_,attrs = parts
            if feature != "exon": 
                continue
            ad = {}
            for a in attrs.split(";"):
                if not a.strip(): 
                    continue
                kv = a.strip().replace('"','').split(" ")
                if len(kv) >= 2:
                    ad[kv[0]] = " ".join(kv[1:])
            gid = ad.get("gene_id")
            if not gid: 
                continue
            exon_map.setdefault(gid, []).append({
                "start":int(start),
                "end":int(end)
            })
    def classify(bp, exons):
        for e in exons:
            if bp == e["start"]: 
                return "S"
            if bp == e["end"]: 
                return "E"
            if e["start"] < bp < e["end"]: 
                return "M"
        return "O"

    df = pd.read_csv("8_fusion_with_junction_counts.tsv", sep="\t", dtype=str)
    left_loc = []
    right_loc = []
    for _, r in df.iterrows():
        try:
            lb = int(r.get("LeftBreakpoint",0))
            rb = int(r.get("RightBreakpoint",0))
        except Exception:
            lb = rb = 0
        left_loc.append(classify(lb, exon_map.get(r.get("5_geneid",""), [])))
        right_loc.append(classify(rb, exon_map.get(r.get("3_geneid",""), [])))
    df["5_loc"] = left_loc
    df["3_loc"] = right_loc
    df.to_csv("10_fusion_with_breakpoint_locations.tsv", sep="\t", index=False)
    print("Saved 10_fusion_with_breakpoint_locations.tsv")

# -------------------------
# STEP 10: Integrate exon count info and finalize
# -------------------------
def step10_exon_counts_and_finalize():
    main_df = pd.read_csv("10_fusion_with_breakpoint_locations.tsv", sep="\t", dtype=str)
    exon_df = pd.read_csv("2_filtered_gtf.tsv", sep="\t", dtype=str, names=["Gene_ID","Transcript_ID","Exon_Count"])
    exon_counts = exon_df.groupby("Gene_ID")["Exon_Count"].unique().to_dict()
    rows = []
    for _, r in main_df.iterrows():
        exon5 = exon_counts.get(r.get("5_geneid"), [""])
        exon3 = exon_counts.get(r.get("3_geneid"), [""])
        if not exon5: 
            exon5 = [""]
        if not exon3: 
            exon3 = [""]
        for e5 in exon5:
            for e3 in exon3:
                nr = r.copy()
                try:
                    nr["exon_count5"] = int(float(e5)) if e5 != "" else ""
                except Exception:
                    nr["exon_count5"] = ""
                try:
                    nr["exon_count3"] = int(float(e3)) if e3 != "" else ""
                except Exception:
                    nr["exon_count3"] = ""
                rows.append(nr)
    out = pd.DataFrame(rows)
    out.to_csv("11_final_fusion_annotated.tsv", sep="\t", index=False)
    print("Saved 11_final_fusion_annotated.tsv")

    # compute Total_Mapped_Reads (FFPM->reads)
    in_f = "11_final_fusion_annotated.tsv"
    out_f = "12_annotated_with_reads.tsv"
    with open(in_f) as infile, open(out_f, "w", newline='') as outfile:
        r = csv.DictReader(infile, delimiter='\t')
        fieldnames = r.fieldnames + ["Total_Mapped_Reads"]
        w = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in r:
            ffpm = row.get("FFPM","")
            total_count = row.get("Total_Count_(SC+RC)","")
            tm = "NF"
            if ffpm not in ["","NA","NF"] and total_count not in ["","NA","NF"]:
                try:
                    ffpm_v = float(ffpm)
                    cnt = float(total_count)
                    if ffpm_v > 0:
                        tm = int((cnt * 1_000_000) / ffpm_v)
                except Exception:
                    tm = "NF"
            row["Total_Mapped_Reads"] = tm
            w.writerow(row)
    print("Saved 12_annotated_with_reads.tsv")

    # prepare fts_features.csv
    df = pd.read_csv("12_annotated_with_reads.tsv", sep="\t", dtype=str)

    # Ensure 3_geneid is present in the final CSV by splitting fusion_pair
    # Example: fusion_pair = "Os05g0352800->Os05g0333200"
    #          3_geneid   = "Os05g0333200"
    if "fusion_pair" in df.columns:
        df["3_geneid"] = df["fusion_pair"].astype(str).str.split("->").str[-1]

    # Drop columns we don't want in the final features CSV
    # (Note: we KEEP 3_geneid now)
    cols_to_remove = ["#FusionName", "CDS_LEFT_ID", "CDS_RIGHT_ID"]
    df = df.drop(
        columns=[c for c in cols_to_remove if c in df.columns],
        errors='ignore'
    )

    df.to_csv("Fusion_output/fts_features.csv", sep=",", index=False)
    print("Saved fts_features.csv")

    # move any numeric-prefixed files to temp_dir (original cleanup)
    temp_dir = "Fusion_output/temp_dir"
    ensure_dir(temp_dir)
    for fn in os.listdir("."):
        if re.match(r"^\d+_.*", fn) and fn != "fts_features.csv":
            shutil.move(fn, os.path.join(temp_dir, fn))
    print("Cleanup (temp_dir) done.")

# -------------------------
# MAIN orchestrator
# -------------------------
def main():
    start = args.start
    outdir = args.outdir
    ensure_dir(outdir)

    # if start==step1 -> run full STAR-Fusion & merging
    if start == "step1":
        if not args.config:
            raise SystemExit("Error: --config required for step1")
        step1_run_star_and_merge(args.config, outdir)
        # after step1 we still need to run subsequent steps, so fallthrough:
        start = "step2"

    if start == "step2":
        if not args.gtf:
            raise SystemExit("Error: --gtf required for step2")
        step2_process_gtf(args.gtf)
        start = "step3"

    if start == "step3":
        # For step3 we need the merged file
        merged_path = args.merged if args.merged else os.path.join(outdir, "1_initial_merged_predictions.tsv")
        if not os.path.exists(merged_path):
            raise SystemExit(f"Error: merged file not found at {merged_path}. Provide --merged or run step1.")
        if not args.gtf:
            raise SystemExit("Error: --gtf required for exon/gene processing")

        # run part 3..10 in order
        step3_fill_transcripts(merged_path)
        step4_add_exons(args.gtf)
        step5_annotate_features()
        step6_extract_gene_info(args.gtf)
        step7_merge_gene_info()
        step8_count_junctions()
        step9_classify_breakpoint_locations(args.gtf)
        step10_exon_counts_and_finalize()
        print("\nALL DONE: final files produced (11_final_fusion_annotated.tsv, 12_annotated_with_reads.tsv, fts_features.csv)")

if __name__ == "__main__":
    main()
