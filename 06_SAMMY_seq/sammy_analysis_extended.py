#!/usr/bin/env python3
"""
Enhanced SAMMY-seq Chromatin State Analysis Pipeline

This script processes SAMMY-seq data to identify heterochromatin and euchromatin regions

Features:
- Processes all four SAMMY-seq fractions (S2S, S2L, S3, S4)
- Implements SPP-like methodology for ratio calculation
- Performs quantile normalization
- Conducts statistical testing with z-tests and multiple testing correction
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
import seaborn as sns
import tempfile
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.stats import zscore
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import eigh
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import quantile_transform

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

# Configuration - set CHROM_SIZES_FILE to your chromosome sizes file
CHROM_SIZES_FILE = ""  # e.g., "/path/to/mm10.chrom.sizes"

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Enhanced SAMMY-seq Chromatin State Analysis')
    parser.add_argument('--workdir', required=True, help='Working directory containing processed BAM files')
    parser.add_argument('--outdir', help='Output directory for results', default=None)
    parser.add_argument('--coverage-bin-size', type=int, default=50, 
                       help='Bin size for initial coverage calculation (default: 50bp as in methods)')
    parser.add_argument('--analysis-bin-size', type=int, default=50000, 
                       help='Bin size for genomic segmentation analysis (default: 50000bp as in methods)')
    parser.add_argument('--eu-threshold', type=float, default=0.1, 
                       help='Log2 ratio threshold for euchromatin (default: 0.1)')
    parser.add_argument('--hetero-threshold', type=float, default=-0.1, 
                       help='Log2 ratio threshold for heterochromatin (default: -0.1)')
    parser.add_argument('--threads', type=int, default=8, 
                       help='Number of threads to use (default: 8)')
    parser.add_argument('--blacklist', help='BED file with blacklisted regions')
    parser.add_argument('--blacklist-merge-distance', type=int, default=50000,
                       help='Distance to merge adjacent blacklisted regions (default: 50000)')
    parser.add_argument('--confidence-interval', type=float, default=0.99,
                       help='Confidence interval for z-test (default: 0.99)')
    parser.add_argument('--no-compartment', action='store_true',
                      help='Skip compartment analysis')
    parser.add_argument('--conditions', help='Comma-separated list of condition names to analyze')
    return parser.parse_args()

def find_bam_files(bam_dir):
    """Find and categorize BAM files by condition and fraction"""
    logger.info(f"Searching for BAM files in {bam_dir}")
    
    # Use a dictionary with dynamic keys for conditions
    samples = {}
    
    bam_files = glob.glob(os.path.join(bam_dir, "*.filtered.bam"))
    logger.info(f"Found {len(bam_files)} BAM files")
    
    for bam_file in bam_files:
        basename = os.path.basename(bam_file)
        
        # Determine condition using regex
        condition_match = re.search(r'(Neu[123]?_GFP|Neu[123]?_M2|NSC[123]?_GFP|NSC[123]?_M2|WT|LmnaKO|LBDKO)', basename)
        if condition_match:
            condition = condition_match.group(1)
        else:
            logger.warning(f"Skipping {basename} - could not determine condition")
            continue
        
        # Initialize condition in dictionary if not present
        if condition not in samples:
            samples[condition] = {"S2S": [], "S2L": [], "S3": [], "S4": []}
        
        # Determine fraction
        if "_S2S_" in basename:
            fraction = "S2S"
        elif "_S2L_" in basename:
            fraction = "S2L"
        elif "_S3_" in basename:
            fraction = "S3"
        elif "_S4_" in basename:
            fraction = "S4"
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
    
    merged_samples = {}
    
    for condition in samples:
        merged_samples[condition] = {"S2S": None, "S2L": None, "S3": None, "S4": None}
        
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
    logger.info(f"Creating coverage BigWig files for all fractions with bin size {bin_size}bp")
    
    coverage_dir = os.path.join(results_dir, "coverage")
    os.makedirs(coverage_dir, exist_ok=True)
    
    coverage_files = {}
    
    for condition in merged_samples:
        coverage_files[condition] = {"S2S": None, "S2L": None, "S3": None, "S4": None}
        
        for fraction in ["S2S", "S2L", "S3", "S4"]:
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
                "--binSize", str(bin_size),   # 50bp binning as in the methods
                "--normalizeUsing", "RPKM",   # RPKM normalization as in the methods
                "--extendReads", "250",       # Extend reads to 250bp as in the methods
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

def rebin_coverage(coverage_files, results_dir, original_bin_size, new_bin_size):
    """Rebin coverage files to a larger bin size for analysis"""
    logger.info(f"Rebinning coverage files from {original_bin_size}bp to {new_bin_size}bp")
    
    rebinned_dir = os.path.join(results_dir, "rebinned_coverage")
    os.makedirs(rebinned_dir, exist_ok=True)
    
    rebinned_files = {}
    
    for condition in coverage_files:
        rebinned_files[condition] = {}
        
        for fraction in ["S2S", "S2L", "S3", "S4"]:
            if not coverage_files[condition][fraction]:
                continue
                
            input_file = coverage_files[condition][fraction]
            output_file = os.path.join(rebinned_dir, f"{condition}_{fraction}_rebinned_{new_bin_size}.bw")
            rebinned_files[condition][fraction] = output_file
            
            # Skip if file already exists
            if os.path.exists(output_file):
                logger.info(f"Rebinned file already exists for {condition} {fraction}")
                continue
            
            logger.info(f"Rebinning {condition} {fraction} to {new_bin_size}bp")
            
            # Use bigwigAverage to rebin the coverage file
            cmd = [
                "bigwigAverage",
                "--bigwigs", input_file,
                "--outFileName", output_file,
                "--binSize", str(new_bin_size),
                "--numberOfProcessors", "4"
            ]
            
            # If bigwigAverage is not available, try alternative with multiBigwigSummary
            try:
                subprocess.run(cmd, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                logger.warning(f"bigwigAverage not available, using alternative method")
                
                # Create temporary bedGraph and rebin
                temp_bedgraph = os.path.join(rebinned_dir, f"temp_{condition}_{fraction}.bedgraph")
                subprocess.run(["bigWigToBedGraph", input_file, temp_bedgraph], check=True)
                
                # Read bedGraph
                df = pd.read_csv(temp_bedgraph, sep='\t', header=None, 
                                names=['chrom', 'start', 'end', 'score'])
                
                # Group by windows of new_bin_size
                df['bin'] = df['start'] // new_bin_size
                df['chrom_bin'] = df['chrom'] + '_' + df['bin'].astype(str)
                
                # Calculate average for each bin
                rebinned = df.groupby(['chrom', 'bin']).agg({
                    'score': 'mean',
                    'start': lambda x: (x.min() // new_bin_size) * new_bin_size,
                }).reset_index()
                
                rebinned['end'] = rebinned['start'] + new_bin_size
                
                # Save to temporary bedGraph
                rebinned[['chrom', 'start', 'end', 'score']].to_csv(
                    temp_bedgraph, sep='\t', header=False, index=False
                )
                
                # Convert to bigWig
                subprocess.run(["bedGraphToBigWig", temp_bedgraph,
                              CHROM_SIZES_FILE, output_file], check=True)
                
                # Clean up temporary file
                os.remove(temp_bedgraph)
    
    return rebinned_files

def apply_quantile_normalization(bedgraph_files, results_dir, blacklist=None, blacklist_merge_distance=50000):
    """Apply quantile normalization to ratio files as described in the methods"""
    logger.info("Applying quantile normalization to ratio files")
    
    normalized_dir = os.path.join(results_dir, "normalized")
    os.makedirs(normalized_dir, exist_ok=True)
    
    # Step 1: Load all bedgraph files
    all_data = {}
    
    for condition, bedgraph_file in bedgraph_files.items():
        if not bedgraph_file or not os.path.exists(bedgraph_file):
            logger.warning(f"Skipping {condition} - missing bedgraph file")
            continue
            
        logger.info(f"Loading data for {condition}")
        
        # Read BedGraph file
        df = pd.read_csv(
            bedgraph_file, 
            sep='\t', 
            header=None,
            names=['chrom', 'start', 'end', 'score']
        )
        
        # Filter out lines with "NA" or NaN scores
        df = df.dropna(subset=['score'])
        
        # Process blacklist if provided
        if blacklist:
            logger.info(f"Processing blacklist regions for {condition}")
            
            # Read blacklist file
            blacklist_df = pd.read_csv(
                blacklist, 
                sep='\t', 
                header=None,
                names=['chrom', 'start', 'end', 'name', 'score', 'strand']
            )
            
            # If blacklist_merge_distance is provided, merge adjacent blacklist regions
            if blacklist_merge_distance > 0:
                # Sort by chromosome and start position
                blacklist_df = blacklist_df.sort_values(['chrom', 'start'])
                
                # Initialize merged blacklist
                merged_blacklist = []
                current_chrom = None
                current_start = None
                current_end = None
                
                # Merge adjacent regions
                for _, row in blacklist_df.iterrows():
                    if current_chrom is None or current_chrom != row['chrom'] or row['start'] > current_end + blacklist_merge_distance:
                        if current_chrom is not None:
                            merged_blacklist.append({
                                'chrom': current_chrom,
                                'start': current_start,
                                'end': current_end
                            })
                        current_chrom = row['chrom']
                        current_start = row['start']
                        current_end = row['end']
                    else:
                        current_end = max(current_end, row['end'])
                
                # Add the last region
                if current_chrom is not None:
                    merged_blacklist.append({
                        'chrom': current_chrom,
                        'start': current_start,
                        'end': current_end
                    })
                
                # Create dataframe from merged blacklist
                blacklist_df = pd.DataFrame(merged_blacklist)
            
            # Filter out blacklisted regions
            for _, blacklist_row in blacklist_df.iterrows():
                df = df[~((df['chrom'] == blacklist_row['chrom']) & 
                         (df['end'] > blacklist_row['start']) & 
                         (df['start'] < blacklist_row['end']))]
        
        all_data[condition] = df
    
    # Step 2: Create a common coordinate system
    logger.info("Creating common coordinate system for all samples")
    
    # Find all unique genomic regions across all samples
    all_regions = pd.concat([df[['chrom', 'start', 'end']] for df in all_data.values()])
    all_regions = all_regions.drop_duplicates().sort_values(['chrom', 'start'])
    
    # Step 3: Create a matrix with scores for each region in each sample
    region_matrix = pd.DataFrame()
    region_matrix[['chrom', 'start', 'end']] = all_regions[['chrom', 'start', 'end']]
    
    for condition, df in all_data.items():
        # Create a dictionary for faster lookup
        score_dict = dict(zip(zip(df['chrom'], df['start'], df['end']), df['score']))
        
        # Add scores for this condition
        region_matrix[condition] = [
            score_dict.get((chrom, start, end), np.nan)
            for chrom, start, end in zip(region_matrix['chrom'], region_matrix['start'], region_matrix['end'])
        ]
    
    # Step 4: Apply quantile normalization to the score matrix
    logger.info("Applying quantile normalization")
    
    # Extract just the score columns
    score_columns = [col for col in region_matrix.columns if col not in ['chrom', 'start', 'end']]
    
    # Handle missing values by filling with median of each condition
    for col in score_columns:
        region_matrix[col] = region_matrix[col].fillna(region_matrix[col].median())
    
    # Apply quantile normalization
    normalized_scores = pd.DataFrame(
        quantile_transform(region_matrix[score_columns], axis=0, n_quantiles=1000),
        columns=score_columns
    )
    
    # Replace original scores with normalized scores
    for col in score_columns:
        region_matrix[col] = normalized_scores[col]
    
    # Step 5: Save normalized data
    logger.info("Saving normalized data")
    
    normalized_files = {}
    
    for condition in score_columns:
        # Create a dataframe with just the data for this condition
        norm_df = region_matrix[['chrom', 'start', 'end', condition]].rename(
            columns={condition: 'score'}
        )
        
        # Save to bedgraph
        normalized_file = os.path.join(normalized_dir, f"{condition}_normalized.bedgraph")
        norm_df.to_csv(normalized_file, sep='\t', header=False, index=False)
        
        # Convert to bigwig if needed
        normalized_bigwig = os.path.join(normalized_dir, f"{condition}_normalized.bw")
        
        try:
            subprocess.run([
                "bedGraphToBigWig",
                normalized_file,
                CHROM_SIZES_FILE,
                normalized_bigwig
            ], check=True)
            normalized_files[condition] = normalized_bigwig
        except subprocess.CalledProcessError:
            logger.warning(f"Error converting to BigWig for {condition}, using BedGraph instead")
            normalized_files[condition] = normalized_file
    
    return normalized_files

def identify_chromatin_states(normalized_files, results_dir, eu_threshold, hetero_threshold):
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
    
    for condition, normalized_file in normalized_files.items():
        if not normalized_file or not os.path.exists(normalized_file):
            logger.warning(f"Skipping {condition} - missing normalized file")
            continue
        
        logger.info(f"Processing {condition}")
        
        try:
            # Read normalized file
            is_bigwig = normalized_file.endswith('.bw')
            
            if is_bigwig:
                # Convert BigWig to BedGraph first
                temp_bedgraph = os.path.join(chromatin_dir, f"{condition}_temp.bedgraph")
                subprocess.run(["bigWigToBedGraph", normalized_file, temp_bedgraph], check=True)
                file_to_read = temp_bedgraph
            else:
                file_to_read = normalized_file
            
            # Read the file
            chromatin_states = pd.read_csv(
                file_to_read, 
                sep='\t', 
                header=None,
                names=['chrom', 'start', 'end', 'score']
            )
            
            # Clean up temp file if created
            if is_bigwig:
                os.remove(temp_bedgraph)
            
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

def perform_differential_analysis(normalized_files, results_dir, confidence_interval=0.99):
    """
    Perform differential enrichment analysis between conditions and implement statistical tests
    """
    logger.info("Performing differential enrichment analysis")
    
    diff_dir = os.path.join(results_dir, "differential")
    os.makedirs(diff_dir, exist_ok=True)
    
    # Load normalized data for all conditions
    condition_data = {}
    
    for condition, normalized_file in normalized_files.items():
        if not normalized_file or not os.path.exists(normalized_file):
            continue
            
        is_bigwig = normalized_file.endswith('.bw')
        
        if is_bigwig:
            # Convert BigWig to BedGraph first
            temp_bedgraph = os.path.join(diff_dir, f"{condition}_temp.bedgraph")
            subprocess.run(["bigWigToBedGraph", normalized_file, temp_bedgraph], check=True)
            file_to_read = temp_bedgraph
        else:
            file_to_read = normalized_file
        
        # Read the file
        df = pd.read_csv(
            file_to_read, 
            sep='\t', 
            header=None,
            names=['chrom', 'start', 'end', 'score']
        )
        
        # Clean up temp file if created
        if is_bigwig:
            os.remove(temp_bedgraph)
        
        # Store the data
        condition_data[condition] = df
    
    # Create a common coordinate system
    all_regions = pd.concat([df[['chrom', 'start', 'end']] for df in condition_data.values()])
    all_regions = all_regions.drop_duplicates().sort_values(['chrom', 'start'])
    
    # Create a matrix with scores for each region
    region_matrix = pd.DataFrame()
    region_matrix[['chrom', 'start', 'end']] = all_regions[['chrom', 'start', 'end']]
    
    for condition, df in condition_data.items():
        # Create a dictionary for faster lookup
        score_dict = dict(zip(zip(df['chrom'], df['start'], df['end']), df['score']))

        # Add scores for this condition
        region_matrix[condition] = [
            score_dict.get((chrom, start, end), np.nan)
            for chrom, start, end in zip(region_matrix['chrom'], region_matrix['start'], region_matrix['end'])
        ]

    # Get all condition pairs to compare
    conditions = [col for col in region_matrix.columns if col not in ['chrom', 'start', 'end']]

    # Define comparison pairs based on condition names
    comparison_pairs = []
    for cond1 in conditions:
        for cond2 in conditions:
            if cond1 != cond2:
                # Compare same cell type between GFP and M2
                if ('GFP' in cond1 and 'M2' in cond2 and
                    cond1.replace('GFP', '') == cond2.replace('M2', '')):
                    comparison_pairs.append((cond1, cond2))
                # Compare NSC vs Neu for same treatment
                elif (('NSC' in cond1 and 'Neu' in cond2) and
                      (('GFP' in cond1 and 'GFP' in cond2) or
                       ('M2' in cond1 and 'M2' in cond2))):
                    comparison_pairs.append((cond1, cond2))

    # Remove duplicates and sort
    comparison_pairs = list(set(comparison_pairs))

    differential_results = {}

    for cond1, cond2 in comparison_pairs:
        logger.info(f"Comparing {cond1} vs {cond2}")

        # Get scores for both conditions
        scores1 = region_matrix[cond1].values
        scores2 = region_matrix[cond2].values

        # Calculate the difference (cond2 - cond1)
        diff = scores2 - scores1

        # Remove NaN values for Z-score calculation
        valid_mask = ~np.isnan(diff)
        diff_valid = diff[valid_mask]

        if len(diff_valid) < 10:
            logger.warning(f"Too few valid regions for {cond1} vs {cond2}")
            continue

        # Calculate Z-scores of the differences
        z_scores = np.full(len(diff), np.nan)
        z_scores[valid_mask] = zscore(diff_valid)

        # Calculate two-tailed p-values from Z-scores
        p_values = np.full(len(diff), np.nan)
        p_values[valid_mask] = 2 * stats.norm.sf(np.abs(z_scores[valid_mask]))

        # Apply Benjamini-Hochberg FDR correction
        padj = np.full(len(diff), np.nan)
        valid_pvals = p_values[valid_mask]
        _, padj_valid, _, _ = multipletests(valid_pvals, method='fdr_bh')
        padj[valid_mask] = padj_valid

        # Create results dataframe
        results_df = pd.DataFrame({
            'chrom': region_matrix['chrom'],
            'start': region_matrix['start'],
            'end': region_matrix['end'],
            f'score_{cond1}': scores1,
            f'score_{cond2}': scores2,
            'diff': diff,
            'z_score': z_scores,
            'pvalue': p_values,
            'padj': padj
        })

        # Identify significant differential regions based on confidence interval
        alpha = 1 - confidence_interval
        sig_regions = results_df[results_df['padj'] < alpha].copy()

        # Separate into up and down regulated
        sig_up = sig_regions[sig_regions['diff'] > 0]
        sig_down = sig_regions[sig_regions['diff'] < 0]

        # Save results
        comparison_name = f"{cond1}_vs_{cond2}"

        # Full results TSV
        results_file = os.path.join(diff_dir, f"{comparison_name}_differential_results.tsv")
        results_df.to_csv(results_file, sep='\t', index=False)

        # Significant regions BED files
        if len(sig_up) > 0:
            sig_up_file = os.path.join(diff_dir, f"{comparison_name}_sig_up.bed")
            sig_up[['chrom', 'start', 'end', 'diff']].to_csv(
                sig_up_file, sep='\t', header=False, index=False
            )

        if len(sig_down) > 0:
            sig_down_file = os.path.join(diff_dir, f"{comparison_name}_sig_down.bed")
            sig_down[['chrom', 'start', 'end', 'diff']].to_csv(
                sig_down_file, sep='\t', header=False, index=False
            )

        # Store summary statistics
        differential_results[comparison_name] = {
            'total_regions': len(results_df),
            'valid_regions': int(valid_mask.sum()),
            'sig_regions': len(sig_regions),
            'sig_up': len(sig_up),
            'sig_down': len(sig_down),
            'sig_up_bp': int(sig_up['end'].sum() - sig_up['start'].sum()) if len(sig_up) > 0 else 0,
            'sig_down_bp': int(sig_down['end'].sum() - sig_down['start'].sum()) if len(sig_down) > 0 else 0,
        }

        logger.info(f"  Total regions: {differential_results[comparison_name]['total_regions']}")
        logger.info(f"  Significant regions: {differential_results[comparison_name]['sig_regions']}")
        logger.info(f"  Up in {cond2}: {differential_results[comparison_name]['sig_up']}")
        logger.info(f"  Down in {cond2}: {differential_results[comparison_name]['sig_down']}")

    # Save summary
    summary_file = os.path.join(diff_dir, "differential_analysis_summary.tsv")
    summary_df = pd.DataFrame.from_dict(differential_results, orient='index')
    summary_df.to_csv(summary_file, sep='\t')
    logger.info(f"Saved differential analysis summary to {summary_file}")

    return differential_results


def main():
    """Main function to run the SAMMY-seq analysis pipeline"""
    args = parse_arguments()

    # Set up output directory
    if args.outdir:
        results_dir = args.outdir
    else:
        results_dir = os.path.join(args.workdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    logger.info(f"Starting SAMMY-seq analysis pipeline")
    logger.info(f"Working directory: {args.workdir}")
    logger.info(f"Output directory: {results_dir}")

    # Find BAM files
    samples = find_bam_files(args.workdir)
    if not samples:
        logger.error("No BAM files found")
        return

    # Merge BAM files
    merged_samples = merge_bam_files(samples, results_dir)

    # Create coverage files
    coverage_files = create_coverage_files(
        merged_samples, results_dir,
        args.coverage_bin_size, args.threads,
        args.blacklist
    )

    # Create ratio files
    ratio_files = create_ratio_files(
        coverage_files, results_dir,
        args.coverage_bin_size, args.threads
    )

    # Convert to bedgraph
    bedgraph_files = convert_bigwig_to_bedgraph(ratio_files, results_dir)

    # Apply quantile normalization
    normalized_files = apply_quantile_normalization(
        bedgraph_files, results_dir,
        args.blacklist, args.blacklist_merge_distance
    )

    # Identify chromatin states
    chromatin_states = identify_chromatin_states(
        normalized_files, results_dir,
        args.eu_threshold, args.hetero_threshold
    )

    # Perform differential analysis
    differential_results = perform_differential_analysis(
        normalized_files, results_dir,
        args.confidence_interval
    )

    logger.info("SAMMY-seq analysis pipeline complete")


if __name__ == "__main__":
    main()