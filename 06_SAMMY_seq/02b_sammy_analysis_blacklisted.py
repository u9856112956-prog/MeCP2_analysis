#!/usr/bin/env python3
"""
SAMMY-seq Chromatin State Analysis Pipeline

This script processes SAMMY-seq data to identify heterochromatin and euchromatin regions
by comparing the S2S (accessible) and S3 (inaccessible) fractions.

Based on the methodology described in the SAMMY-seq paper where:
- Euchromatin: regions with log2(S2S/S3) > 0.1
- Heterochromatin: regions with log2(S2S/S3) < -0.1
"""

import os
import sys
import glob
import argparse
import subprocess
import re
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tempfile
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import pybedtools

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("sammy_seq_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('SAMMY-seq')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='SAMMY-seq Chromatin State Analysis')
    parser.add_argument('--workdir', required=True, help='Working directory containing processed BAM files')
    parser.add_argument('--outdir', help='Output directory for results', default=None)
    parser.add_argument('--bin-size', type=int, default=50000, help='Bin size for genome segmentation (default: 50000)')
    parser.add_argument('--eu-threshold', type=float, default=0.1, help='Log2 ratio threshold for euchromatin (default: 0.1)')
    parser.add_argument('--hetero-threshold', type=float, default=-0.1, help='Log2 ratio threshold for heterochromatin (default: -0.1)')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use (default: 8)')
    parser.add_argument('--blacklist', help='BED file with blacklisted regions')
    return parser.parse_args()

def find_bam_files(bam_dir):
    """Find and categorize BAM files by condition and fraction"""
    logger.info(f"Searching for BAM files in {bam_dir}")
    samples = {
        "Neu_GFP": {"S2S": [], "S3": []},
        "Neu_M2": {"S2S": [], "S3": []},
        "NSC_GFP": {"S2S": [], "S3": []},
        "NSC_M2": {"S2S": [], "S3": []}
    }
    
    bam_files = glob.glob(os.path.join(bam_dir, "*.filtered.bam"))
    logger.info(f"Found {len(bam_files)} BAM files")
    
    for bam_file in bam_files:
        basename = os.path.basename(bam_file)
        
        # Determine condition and fraction using regex
        if re.search(r'Neu[123]?_GFP', basename):
            condition = "Neu_GFP"
        elif re.search(r'Neu[123]?_M2', basename):
            condition = "Neu_M2"
        elif re.search(r'NSC[123]?_GFP', basename):
            condition = "NSC_GFP"
        elif re.search(r'NSC[123]?_M2', basename):
            condition = "NSC_M2"
        else:
            logger.warning(f"Skipping {basename} - could not determine condition")
            continue
        
        if "_S2S_" in basename:
            fraction = "S2S"
        elif "_S3_" in basename:
            fraction = "S3"
        else:
            logger.warning(f"Skipping {basename} - could not determine fraction")
            continue
        
        # Ensure BAM file is indexed
        if not os.path.exists(f"{bam_file}.bai"):
            logger.info(f"Indexing {basename}")
            subprocess.run(["samtools", "index", bam_file], check=True)
        
        samples[condition][fraction].append(bam_file)
        logger.info(f"Added {basename} to {condition} {fraction}")
    
    return samples

def merge_bam_files(samples, results_dir):
    """Merge BAM files for replicates of each condition and fraction"""
    logger.info("Merging BAM files for replicates")
    merged_dir = os.path.join(results_dir, "merged_bams")
    os.makedirs(merged_dir, exist_ok=True)
    
    merged_samples = {
        "Neu_GFP": {"S2S": None, "S3": None},
        "Neu_M2": {"S2S": None, "S3": None},
        "NSC_GFP": {"S2S": None, "S3": None},
        "NSC_M2": {"S2S": None, "S3": None}
    }
    
    for condition in samples:
        for fraction in samples[condition]:
            if not samples[condition][fraction]:
                logger.warning(f"No files found for {condition} {fraction}")
                continue
            
            output_file = os.path.join(merged_dir, f"{condition}_{fraction}_merged.bam")
            
            # If only one file, create a symlink
            if len(samples[condition][fraction]) == 1:
                logger.info(f"Only one file for {condition} {fraction}, creating symlink")
                if os.path.exists(output_file):
                    os.remove(output_file)
                os.symlink(samples[condition][fraction][0], output_file)
                if not os.path.exists(f"{output_file}.bai"):
                    subprocess.run(["samtools", "index", output_file], check=True)
                merged_samples[condition][fraction] = output_file
                continue
            
            # If merged file already exists, use it
            if os.path.exists(output_file):
                logger.info(f"Merged file already exists for {condition} {fraction}")
                if not os.path.exists(f"{output_file}.bai"):
                    subprocess.run(["samtools", "index", output_file], check=True)
                merged_samples[condition][fraction] = output_file
                continue
            
            # Merge files
            logger.info(f"Merging {len(samples[condition][fraction])} files for {condition} {fraction}")
            cmd = ["samtools", "merge", "-f", output_file] + samples[condition][fraction]
            subprocess.run(cmd, check=True)
            
            # Index merged file
            subprocess.run(["samtools", "index", output_file], check=True)
            
            merged_samples[condition][fraction] = output_file
    
    return merged_samples

def create_coverage_files(merged_samples, results_dir, bin_size, threads, blacklist=None):
    """Create coverage BigWig files for each condition and fraction"""
    logger.info("Creating coverage BigWig files for S2S and S3 fractions")
    
    coverage_dir = os.path.join(results_dir, "coverage")
    os.makedirs(coverage_dir, exist_ok=True)
    
    coverage_files = {
        "Neu_GFP": {"S2S": None, "S3": None},
        "Neu_M2": {"S2S": None, "S3": None},
        "NSC_GFP": {"S2S": None, "S3": None},
        "NSC_M2": {"S2S": None, "S3": None}
    }
    
    for condition in merged_samples:
        for fraction in ["S2S", "S3"]:
            if not merged_samples[condition][fraction]:
                logger.warning(f"Skipping {condition} {fraction} - missing BAM file")
                continue
            
            bam_file = merged_samples[condition][fraction]
            output_file = os.path.join(coverage_dir, f"{condition}_{fraction}_coverage.bw")
            coverage_files[condition][fraction] = output_file
            
            # Skip if file already exists
            if os.path.exists(output_file):
                logger.info(f"Coverage file already exists for {condition} {fraction}")
                continue
            
            logger.info(f"Creating coverage file for {condition} {fraction}")
            
            # Create command with proper parameters based on methods paper
            cmd = [
                "bamCoverage",
                "--bam", bam_file,
                "--outFileName", output_file,
                "--binSize", str(bin_size),  # 50kb binning as in the paper
                "--normalizeUsing", "RPKM",  # RPKM normalization as in the paper
                "--extendReads", "250",      # Extend reads to 250bp as in the paper
                "--numberOfProcessors", str(threads)
            ]
            
            if blacklist:
                cmd.extend(["--blackListFileName", blacklist])
            
            try:
                process = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Successfully created coverage file for {condition} {fraction}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error creating coverage file for {condition} {fraction}:\n{e.stderr}")
                coverage_files[condition][fraction] = None
    
    return coverage_files

def create_ratio_files(coverage_files, results_dir, bin_size, threads):
    """Create log2 ratio BigWig files (S2S/S3) for each condition"""
    logger.info("Creating log2 ratio files for S2S/S3")
    
    ratio_dir = os.path.join(results_dir, "ratio")
    os.makedirs(ratio_dir, exist_ok=True)
    
    ratio_files = {}
    
    for condition in coverage_files:
        if not coverage_files[condition]["S2S"] or not coverage_files[condition]["S3"]:
            logger.warning(f"Skipping {condition} - missing S2S or S3 coverage file")
            continue
        
        s2s_file = coverage_files[condition]["S2S"]
        s3_file = coverage_files[condition]["S3"]
        
        output_file = os.path.join(ratio_dir, f"{condition}_S2S_vs_S3_ratio.bw")
        ratio_files[condition] = output_file
        
        # Skip if file already exists
        if os.path.exists(output_file):
            logger.info(f"Ratio file already exists for {condition}")
            continue
        
        logger.info(f"Creating ratio file for {condition}")
        
        # Create command with proper parameters to mimic SPP package functionality
        cmd = [
            "bigwigCompare",
            "-b1", s2s_file,
            "-b2", s3_file,
            "--operation", "log2",
            "--pseudocount", "1",     # Add pseudocount to avoid log2(0) issues
            "--binSize", str(bin_size),
            "--numberOfProcessors", str(threads),
            "--outFileName", output_file
        ]
        
        try:
            process = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Successfully created ratio file for {condition}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating ratio file for {condition}:\n{e.stderr}")
            
            # Try alternative approach with different parameters
            logger.info(f"Retrying with alternative parameters for {condition}")
            cmd = [
                "bigwigCompare",
                "-b1", s2s_file,
                "-b2", s3_file,
                "--operation", "log2",
                "--pseudocount", "1",
                "--binSize", str(bin_size),
                "--scaleFactorsMethod", "None",
                "--numberOfProcessors", str(threads),
                "--outFileName", output_file
            ]
            
            try:
                process = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Successfully created ratio file with alternative parameters for {condition}")
            except subprocess.CalledProcessError as e:
                logger.error(f"All attempts failed for {condition}:\n{e.stderr}")
                ratio_files[condition] = None
    
    return ratio_files

def convert_bigwig_to_bedgraph(ratio_files, results_dir):
    """Convert BigWig files to BedGraph for easier processing"""
    logger.info("Converting ratio BigWig files to BedGraph")
    
    bedgraph_dir = os.path.join(results_dir, "bedgraph")
    os.makedirs(bedgraph_dir, exist_ok=True)
    
    bedgraph_files = {}
    
    for condition, bigwig_file in ratio_files.items():
        if not bigwig_file or not os.path.exists(bigwig_file):
            logger.warning(f"Skipping {condition} - missing ratio file")
            continue
        
        output_file = os.path.join(bedgraph_dir, f"{condition}_S2S_vs_S3_ratio.bedgraph")
        bedgraph_files[condition] = output_file
        
        # Skip if file already exists
        if os.path.exists(output_file):
            logger.info(f"BedGraph file already exists for {condition}")
            continue
        
        logger.info(f"Converting BigWig to BedGraph for {condition}")
        
        cmd = ["bigWigToBedGraph", bigwig_file, output_file]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"Successfully converted BigWig to BedGraph for {condition}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error converting BigWig to BedGraph for {condition}")
            bedgraph_files[condition] = None
    
    return bedgraph_files

def identify_chromatin_states(bedgraph_files, results_dir, eu_threshold, hetero_threshold):
    """
    Identify euchromatin and heterochromatin regions based on thresholds
    
    As stated in the methods paper:
    - Euchromatin: log2 ratio > 0.1 (S2S > S3)
    - Heterochromatin: log2 ratio < -0.1 (S2S < S3)
    """
    logger.info("Identifying chromatin states")
    
    chromatin_dir = os.path.join(results_dir, "chromatin_states")
    os.makedirs(chromatin_dir, exist_ok=True)
    
    chromatin_state_results = {}
    
    for condition, bedgraph_file in bedgraph_files.items():
        if not bedgraph_file or not os.path.exists(bedgraph_file):
            logger.warning(f"Skipping {condition} - missing bedgraph file")
            continue
        
        logger.info(f"Processing {condition}")
        
        try:
            # Read BedGraph file
            chromatin_states = pd.read_csv(
                bedgraph_file, 
                sep='\t', 
                header=None,
                names=['chrom', 'start', 'end', 'score']
            )
            
            # Filter out lines with "NA" or NaN scores
            chromatin_states = chromatin_states.dropna(subset=['score'])
            
            # Apply thresholds from the paper
            euchromatin = chromatin_states[chromatin_states['score'] > eu_threshold]
            heterochromatin = chromatin_states[chromatin_states['score'] < hetero_threshold]
            
            # Save results to BED files
            euchromatin_bed = os.path.join(chromatin_dir, f"{condition}_euchromatin.bed")
            heterochromatin_bed = os.path.join(chromatin_dir, f"{condition}_heterochromatin.bed")
            
            euchromatin[['chrom', 'start', 'end', 'score']].to_csv(
                euchromatin_bed, sep='\t', header=False, index=False
            )
            heterochromatin[['chrom', 'start', 'end', 'score']].to_csv(
                heterochromatin_bed, sep='\t', header=False, index=False
            )
            
            # Calculate statistics
            euchromatin_count = len(euchromatin)
            heterochromatin_count = len(heterochromatin)
            
            # Calculate total bases
            euchromatin_bp = int(euchromatin['end'].sum() - euchromatin['start'].sum())
            heterochromatin_bp = int(heterochromatin['end'].sum() - heterochromatin['start'].sum())
            total_bp = euchromatin_bp + heterochromatin_bp
            
            # Store results
            chromatin_state_results[condition] = {
                'euchromatin_count': euchromatin_count,
                'heterochromatin_count': heterochromatin_count,
                'euchromatin_bp': euchromatin_bp,
                'heterochromatin_bp': heterochromatin_bp,
                'euchromatin_ratio': euchromatin_bp / total_bp if total_bp > 0 else 0,
                'heterochromatin_ratio': heterochromatin_bp / total_bp if total_bp > 0 else 0,
                'euchromatin_bed': euchromatin_bed,
                'heterochromatin_bed': heterochromatin_bed
            }
            
            logger.info(f"Processed {condition}: {euchromatin_count} euchromatin regions, "
                        f"{heterochromatin_count} heterochromatin regions")
            
        except Exception as e:
            logger.error(f"Error processing {condition}: {e}")
            continue
    
    # Save summary to TSV
    summary_file = os.path.join(results_dir, "chromatin_state_summary.tsv")
    summary_df = pd.DataFrame.from_dict(chromatin_state_results, orient='index')
    summary_df.to_csv(summary_file, sep='\t')
    logger.info(f"Saved summary to {summary_file}")
    
    return chromatin_state_results

def compare_conditions(results_dir, condition1, condition2, state_type="heterochromatin"):
    """Compare chromatin states between two conditions"""
    logger.info(f"Comparing {state_type} regions between {condition1} and {condition2}")
    
    # Create temporary directory for intermediate files
    tmp_dir = os.path.join(results_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    file1 = os.path.join(results_dir, "chromatin_states", f"{condition1}_{state_type}.bed")
    file2 = os.path.join(results_dir, "chromatin_states", f"{condition2}_{state_type}.bed")
    
    if not os.path.exists(file1) or not os.path.exists(file2):
        logger.error(f"Missing BED files for comparison")
        return None
    
    # Find regions specific to condition1
    specific_to_1_file = os.path.join(tmp_dir, f"{state_type}_specific_to_{condition1}.bed")
    specific_to_1_cmd = f"bedtools subtract -a {file1} -b {file2} -A > {specific_to_1_file}"
    subprocess.run(specific_to_1_cmd, shell=True, check=True)
    
    # Find regions specific to condition2
    specific_to_2_file = os.path.join(tmp_dir, f"{state_type}_specific_to_{condition2}.bed")
    specific_to_2_cmd = f"bedtools subtract -a {file2} -b {file1} -A > {specific_to_2_file}"
    subprocess.run(specific_to_2_cmd, shell=True, check=True)
    
    # Find common regions
    common_file = os.path.join(tmp_dir, f"{state_type}_common_{condition1}_and_{condition2}.bed")
    common_cmd = f"bedtools intersect -a {file1} -b {file2} > {common_file}"
    subprocess.run(common_cmd, shell=True, check=True)
    
    # Count regions
    count_1 = sum(1 for _ in open(specific_to_1_file))
    count_2 = sum(1 for _ in open(specific_to_2_file))
    count_common = sum(1 for _ in open(common_file))
    
    # Calculate base pairs
    if count_1 > 0:
        specific_1_df = pd.read_csv(specific_to_1_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
        specific_1_bp = int(specific_1_df['end'].sum() - specific_1_df['start'].sum())
    else:
        specific_1_bp = 0
        
    if count_2 > 0:
        specific_2_df = pd.read_csv(specific_to_2_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
        specific_2_bp = int(specific_2_df['end'].sum() - specific_2_df['start'].sum())
    else:
        specific_2_bp = 0
        
    if count_common > 0:
        common_df = pd.read_csv(common_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
        common_bp = int(common_df['end'].sum() - common_df['start'].sum())
    else:
        common_bp = 0
    
    return {
        'specific_to_1_count': count_1,
        'specific_to_2_count': count_2,
        'common_count': count_common,
        'specific_to_1_bp': specific_1_bp,
        'specific_to_2_bp': specific_2_bp,
        'common_bp': common_bp
    }

def create_summary_plots(chromatin_state_results, comparisons, results_dir):
    """Create summary plots for chromatin states"""
    logger.info("Creating summary plots")
    
    if not chromatin_state_results:
        logger.error("No chromatin state results to plot")
        return
    
    plot_dir = os.path.join(results_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Create DataFrame from results
    plot_columns = [
        'euchromatin_count', 'heterochromatin_count', 
        'euchromatin_ratio', 'heterochromatin_ratio'
    ]
    results_df = pd.DataFrame.from_dict(chromatin_state_results, orient='index')
    results_df = results_df[plot_columns]
    
    # 2. Plot chromatin state proportions
    plt.figure(figsize=(14, 10))
    
    # A. Stacked bar chart showing proportion of heterochromatin vs euchromatin
    ax1 = plt.subplot(2, 2, 1)
    results_df[['heterochromatin_ratio', 'euchromatin_ratio']].plot(
        kind='bar', stacked=True, ax=ax1,
        color=['#1f77b4', '#aec7e8']
    )
    ax1.set_ylabel('Proportion of genome')
    ax1.set_title('Heterochromatin vs Euchromatin Proportion')
    ax1.legend(['Heterochromatin', 'Euchromatin'])
    ax1.set_ylim(0, 1)
    
    # B. Count of domains
    ax2 = plt.subplot(2, 2, 2)
    results_df[['heterochromatin_count', 'euchromatin_count']].plot(
        kind='bar', ax=ax2,
        color=['#1f77b4', '#aec7e8']
    )
    ax2.set_ylabel('Number of domains')
    ax2.set_title('Count of Heterochromatin and Euchromatin Domains')
    ax2.legend(['Heterochromatin', 'Euchromatin'])
    
    # If we have comparisons, plot those too
    if comparisons:
        # C. Comparison of heterochromatin regions
        if "Neu_GFP_vs_Neu_M2" in comparisons:
            comparison = comparisons["Neu_GFP_vs_Neu_M2"]["heterochromatin"]
            ax3 = plt.subplot(2, 2, 3)
            
            labels = ['Neu_GFP specific', 'Common', 'Neu_M2 specific']
            sizes = [comparison['specific_to_1_count'], 
                     comparison['common_count'], 
                     comparison['specific_to_2_count']]
            
            ax3.pie(sizes, labels=labels, autopct='%1.1f%%',
                   shadow=False, startangle=90,
                   colors=['#ff9999','#66b3ff','#99ff99'])
            ax3.axis('equal')
            ax3.set_title('Comparison of Heterochromatin Regions\nNeu_GFP vs Neu_M2')
        
        # D. Comparison of heterochromatin regions for NSC
        if "NSC_GFP_vs_NSC_M2" in comparisons:
            comparison = comparisons["NSC_GFP_vs_NSC_M2"]["heterochromatin"]
            ax4 = plt.subplot(2, 2, 4)
            
            labels = ['NSC_GFP specific', 'Common', 'NSC_M2 specific']
            sizes = [comparison['specific_to_1_count'], 
                     comparison['common_count'], 
                     comparison['specific_to_2_count']]
            
            ax4.pie(sizes, labels=labels, autopct='%1.1f%%',
                   shadow=False, startangle=90,
                   colors=['#ff9999','#66b3ff','#99ff99'])
            ax4.axis('equal')
            ax4.set_title('Comparison of Heterochromatin Regions\nNSC_GFP vs NSC_M2')
    
    plt.tight_layout()
    plot_file = os.path.join(plot_dir, "chromatin_states_summary.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()
    
    # 3. Create heatmap visualization for genome browser
    for condition in chromatin_state_results:
        try:
            # Read the ratio file
            bedgraph_file = os.path.join(results_dir, "bedgraph", f"{condition}_S2S_vs_S3_ratio.bedgraph")
            if not os.path.exists(bedgraph_file):
                continue
                
            df = pd.read_csv(bedgraph_file, sep='\t', header=None, 
                            names=['chrom', 'start', 'end', 'score'])
            
            # Create a custom colormap: blue for heterochromatin, red for euchromatin
            cmap = LinearSegmentedColormap.from_list(
                'custom_diverging',
                [(0, '#1f77b4'), (0.5, 'white'), (1, '#d62728')],
                N=256
            )
            
            # Clip values to a reasonable range
            df['score'] = np.clip(df['score'], -2, 2)
            
            # For each chromosome, create a separate heatmap
            chromosomes = df['chrom'].unique()
            for i, chrom in enumerate(chromosomes[:5]):  # Limit to first 5 chromosomes
                chrom_df = df[df['chrom'] == chrom]
                
                plt.figure(figsize=(10, 2))
                plt.title(f"{condition} - {chrom} - S2S/S3 Log2 Ratio")
                
                # Create a scatter plot with color mapped to score
                plt.scatter(chrom_df['start'], np.ones(len(chrom_df)), 
                           c=chrom_df['score'], cmap=cmap, 
                           s=10, marker='s')
                
                plt.colorbar(label='Log2 Ratio')
                plt.ylim(0.5, 1.5)
                plt.axis('off')
                
                # Save the figure
                chrom_plot_file = os.path.join(plot_dir, f"{condition}_{chrom}_heatmap.png")
                plt.savefig(chrom_plot_file, dpi=150, bbox_inches='tight')
                plt.close()
                
        except Exception as e:
            logger.error(f"Error creating heatmap for {condition}: {e}")
    
    logger.info(f"Saved summary plots to {plot_dir}")

def merge_overlapping_blacklist(blacklist_file, output_dir):
    """Merge overlapping regions in a blacklist file"""
    if not blacklist_file:
        return None
        
    logger.info(f"Processing blacklist file: {blacklist_file}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Output file path
    merged_blacklist = os.path.join(output_dir, "merged_blacklist.bed")
    
    try:
        # Use pybedtools to merge overlapping regions
        bed = pybedtools.BedTool(blacklist_file)
        merged = bed.sort().merge()
        merged.saveas(merged_blacklist)
        
        # Count original and merged regions
        original_count = len(bed)
        merged_count = len(merged)
        
        logger.info(f"Original blacklist had {original_count} regions")
        logger.info(f"After merging, blacklist has {merged_count} regions")
        logger.info(f"Merged blacklist saved to {merged_blacklist}")
        
        return merged_blacklist
    except Exception as e:
        logger.error(f"Error merging blacklist regions: {str(e)}")
        return blacklist_file  # Return original file if there's an error

def run_sammy_seq_analysis(args):
    """Run the full SAMMY-seq analysis pipeline"""
    logger.info("Starting SAMMY-seq chromatin state analysis")
    
    # Set up directories
    bam_dir = os.path.join(args.workdir, "results/bam")
    if not os.path.exists(bam_dir):
        logger.error(f"BAM directory not found: {bam_dir}")
        sys.exit(1)
    
    results_dir = args.outdir if args.outdir else os.path.join(args.workdir, "chromatin_states_results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Process blacklist file if provided
    processed_blacklist = None
    if args.blacklist:
        processed_blacklist = merge_overlapping_blacklist(args.blacklist, results_dir)
    
    # 1. Find and categorize BAM files
    samples = find_bam_files(bam_dir)
    
    # 2. Merge BAM files for replicates
    merged_samples = merge_bam_files(samples, results_dir)
    
    # 3. Create individual coverage files for S2S and S3 fractions
    coverage_files = create_coverage_files(
        merged_samples, 
        results_dir, 
        args.bin_size, 
        args.threads,
        processed_blacklist if processed_blacklist else args.blacklist
    )
    
    # 4. Create log2 ratio files (S2S/S3)
    ratio_files = create_ratio_files(
        coverage_files, 
        results_dir, 
        args.bin_size, 
        args.threads
    )
    
    # 5. Convert BigWig to BedGraph
    bedgraph_files = convert_bigwig_to_bedgraph(ratio_files, results_dir)
    
    # 6. Identify chromatin states
    chromatin_state_results = identify_chromatin_states(
        bedgraph_files, 
        results_dir, 
        args.eu_threshold, 
        args.hetero_threshold
    )
    
    # 7. Compare conditions if we have multiple conditions
    comparisons = {}
    if "Neu_GFP" in chromatin_state_results and "Neu_M2" in chromatin_state_results:
        comparisons["Neu_GFP_vs_Neu_M2"] = {
            "heterochromatin": compare_conditions(results_dir, "Neu_GFP", "Neu_M2", "heterochromatin"),
            "euchromatin": compare_conditions(results_dir, "Neu_GFP", "Neu_M2", "euchromatin")
        }

    if "NSC_GFP" in chromatin_state_results and "NSC_M2" in chromatin_state_results:
        comparisons["NSC_GFP_vs_NSC_M2"] = {
            "heterochromatin": compare_conditions(results_dir, "NSC_GFP", "NSC_M2", "heterochromatin"),
            "euchromatin": compare_conditions(results_dir, "NSC_GFP", "NSC_M2", "euchromatin")
        }
    
    if "Neu_GFP" in chromatin_state_results and "NSC_GFP" in chromatin_state_results:
        comparisons["Neu_GFP_vs_NSC_GFP"] = {
            "heterochromatin": compare_conditions(results_dir, "Neu_GFP", "NSC_GFP", "heterochromatin"),
            "euchromatin": compare_conditions(results_dir, "Neu_GFP", "NSC_GFP", "euchromatin")
        }
    
    # 8. Create summary plots
    create_summary_plots(chromatin_state_results, comparisons, results_dir)
    
    # Save comparison results to TSV
    comparison_file = os.path.join(results_dir, "condition_comparison_summary.tsv")
    comparison_results = {}
    
    for comp_name, comp_data in comparisons.items():
        for state_type, data in comp_data.items():
            if data:
                key = f"{comp_name}_{state_type}"
                comparison_results[key] = data
    
    if comparison_results:
        comparison_df = pd.DataFrame.from_dict(comparison_results, orient='index')
        comparison_df.to_csv(comparison_file, sep='\t')
        logger.info(f"Saved comparison results to {comparison_file}")
    
    logger.info("SAMMY-seq chromatin state analysis complete")
    logger.info(f"Results saved to {results_dir}")

if __name__ == "__main__":
    args = parse_arguments()
    run_sammy_seq_analysis(args)