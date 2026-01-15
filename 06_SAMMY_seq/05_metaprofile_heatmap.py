#!/usr/bin/env python3
"""
SAMMY-seq Metaprofile Heatmap Generator

This script creates metaprofile heatmaps for gene sets, showing the log2 ratio of S2S/S3
across gene bodies (from -5kb to TSS to TES to +5kb).

The heatmap uses a blue-white-red color scheme where:
- Blue represents heterochromatin (log2(S2S/S3) < -0.1)
- White represents undefined regions (-0.1 <= log2(S2S/S3) <= 0.1)
- Red represents euchromatin (log2(S2S/S3) > 0.1)
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pyBigWig
import pybedtools
from tqdm import tqdm
import seaborn as sns

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("mecp2_metaprofile_heatmap.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('SAMMY-seq-Metaprofile')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='SAMMY-seq Metaprofile Heatmap Generator')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-dir', required=True, help='Output directory for heatmaps')
    parser.add_argument('--genome-gtf', required=True, help='GTF file with gene annotations')
    parser.add_argument('--gene-lists-dir', required=True, help='Directory containing gene list files')
    parser.add_argument('--upstream-distance', type=int, default=5000, help='Distance upstream of TSS (default: 5000)')
    parser.add_argument('--downstream-distance', type=int, default=5000, help='Distance downstream of TES (default: 5000)')
    parser.add_argument('--bin-size', type=int, default=100, help='Bin size for metaprofile (default: 100)')
    parser.add_argument('--gene-body-bins', type=int, default=100, help='Number of bins for gene body (default: 100)')
    return parser.parse_args()

def load_gene_list(file_path):
    """Load gene list from file"""
    logger.info(f"Loading gene list from {file_path}")
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        logger.error(f"Error loading gene list from {file_path}: {e}")
        return []

def get_gene_coordinates(genes, gtf_file):
    """Get gene coordinates (TSS, TES, strand) from GTF file"""
    logger.info(f"Getting coordinates for {len(genes)} genes from GTF file")
    
    # Create a BED file with gene names for filtering
    gene_list_bed = pybedtools.BedTool('\n'.join([f"placeholder\t0\t1\t{gene}" for gene in genes]), from_string=True)
    
    # Parse GTF file to extract gene coordinates
    gene_coords = {}
    try:
        gtf = pybedtools.BedTool(gtf_file)
        
        # Filter for gene features and extract relevant information
        gene_features = gtf.filter(lambda x: x[2] == 'gene')
        
        for feature in gene_features:
            attributes = dict(item.strip().split(' ') for item in feature.fields[8].split(';') if item.strip())
            
            # Extract gene name (removing quotes if present)
            gene_name = None
            for attr in ['gene_name', 'gene_id']:
                if attr in attributes:
                    gene_name = attributes[attr].strip('"\'')
                    break
            
            if gene_name and gene_name in genes:
                chrom = feature.chrom
                start = int(feature.start)  # BedTool uses 0-based coordinates
                end = int(feature.end)
                strand = feature.strand
                
                gene_coords[gene_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'tss': start if strand == '+' else end,
                    'tes': end if strand == '+' else start
                }
        
        logger.info(f"Found coordinates for {len(gene_coords)} out of {len(genes)} genes")
        return gene_coords
    
    except Exception as e:
        logger.error(f"Error getting gene coordinates: {e}")
        return {}

def compute_metaprofile(gene_coords, bigwig_file, upstream_distance, downstream_distance, bin_size, gene_body_bins):
    """Compute metaprofile for a set of genes from a BigWig file"""
    logger.info(f"Computing metaprofile from {bigwig_file}")
    
    try:
        bw = pyBigWig.open(bigwig_file)
        
        # Initialize matrices for storing values
        upstream_bins = int(upstream_distance / bin_size)
        downstream_bins = int(downstream_distance / bin_size)
        total_bins = upstream_bins + gene_body_bins + downstream_bins
        
        # Initialize matrix for all genes
        matrix = np.zeros((len(gene_coords), total_bins))
        matrix.fill(np.nan)  # Fill with NaN initially
        
        # Process each gene
        for i, (gene_name, coords) in enumerate(tqdm(gene_coords.items(), desc="Processing genes")):
            chrom = coords['chrom']
            tss = coords['tss']
            tes = coords['tes']
            strand = coords['strand']
            
            # Skip genes on chromosomes not in the BigWig file
            if chrom not in bw.chroms():
                continue
            
            # Calculate regions
            if strand == '+':
                upstream_start = max(0, tss - upstream_distance)
                upstream_end = tss
                gene_body_start = tss
                gene_body_end = tes
                downstream_start = tes
                downstream_end = min(bw.chroms()[chrom], tes + downstream_distance)
            else:  # strand == '-'
                upstream_start = max(0, tes - upstream_distance)
                upstream_end = tes
                gene_body_start = tes
                gene_body_end = tss
                downstream_start = tss
                downstream_end = min(bw.chroms()[chrom], tss + downstream_distance)
            
            # Get values for upstream region
            if upstream_end > upstream_start:
                try:
                    upstream_values = bw.stats(chrom, upstream_start, upstream_end, 
                                              nBins=upstream_bins, type="mean")
                    for j, val in enumerate(upstream_values):
                        if val is not None:
                            matrix[i, j] = val
                except:
                    pass
            
            # Get values for gene body
            if gene_body_end > gene_body_start:
                try:
                    gene_body_values = bw.stats(chrom, gene_body_start, gene_body_end, 
                                               nBins=gene_body_bins, type="mean")
                    for j, val in enumerate(gene_body_values):
                        if val is not None:
                            matrix[i, upstream_bins + j] = val
                except:
                    pass
            
            # Get values for downstream region
            if downstream_end > downstream_start:
                try:
                    downstream_values = bw.stats(chrom, downstream_start, downstream_end, 
                                                nBins=downstream_bins, type="mean")
                    for j, val in enumerate(downstream_values):
                        if val is not None:
                            matrix[i, upstream_bins + gene_body_bins + j] = val
                except:
                    pass
        
        bw.close()
        return matrix
    
    except Exception as e:
        logger.error(f"Error computing metaprofile: {e}")
        return None

def create_metaprofile_heatmap(matrices, conditions, gene_list_name, output_file, 
                              upstream_bins, gene_body_bins, downstream_bins):
    """Create a metaprofile heatmap for multiple conditions"""
    logger.info(f"Creating metaprofile heatmap for {gene_list_name}")
    
    try:
        # Create a custom colormap: blue for heterochromatin, white for undefined, red for euchromatin
        cmap = LinearSegmentedColormap.from_list(
            'custom_diverging',
            [(0, '#1f77b4'), (0.45, '#1f77b4'), (0.5, 'white'), (0.55, '#d62728'), (1, '#d62728')],
            N=256
        )
        
        # Set up the figure
        fig, ax = plt.subplots(figsize=(10, len(conditions) * 1.5))
        
        # Combine matrices for visualization
        combined_matrix = np.vstack([np.nanmean(matrix, axis=0).reshape(1, -1) for matrix in matrices])
        
        # Set vmin and vmax for consistent color scaling
        vmin = -0.1
        vmax = 0.1
        
        # Create heatmap
        sns.heatmap(combined_matrix, cmap=cmap, vmin=vmin, vmax=vmax, 
                   cbar=True, ax=ax, yticklabels=conditions)
        
        # Add vertical lines to mark TSS and TES
        ax.axvline(x=upstream_bins, color='black', linestyle='--', linewidth=1)
        ax.axvline(x=upstream_bins + gene_body_bins, color='black', linestyle='--', linewidth=1)
        
        # Set x-axis labels
        total_bins = upstream_bins + gene_body_bins + downstream_bins
        xticks = [0, upstream_bins, upstream_bins + gene_body_bins, total_bins]
        xticklabels = [f"-{upstream_bins/10}kb", "TSS", "TES", f"+{downstream_bins/10}kb"]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        
        # Set title and labels
        plt.title(f"Mecp2 repressed genes {gene_list_name} metaprofile")
        plt.xlabel("genes")
        
        # Save the figure
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved heatmap to {output_file}")
    
    except Exception as e:
        logger.error(f"Error creating metaprofile heatmap: {e}")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define gene lists to process
    gene_lists = {
        "NEU_target": "NEU_exo.txt",
        "NEU_non_target": "no_NEU_exo.txt",
        "NSC_target": "NSC_exo.txt",
        "NSC_non_target": "no_NSC_exo.txt"
    }
    
    # Define conditions to process
    conditions = ["NSC_GFP", "NSC_M2", "Neu_GFP", "Neu_M2"]
    
    # Calculate bin counts
    upstream_bins = int(args.upstream_distance / args.bin_size)
    downstream_bins = int(args.downstream_distance / args.bin_size)
    
    # Process each gene list
    for gene_list_name, gene_list_file in gene_lists.items():
        logger.info(f"Processing gene list: {gene_list_name}")
        
        # Load gene list
        gene_list_path = os.path.join(args.gene_lists_dir, gene_list_file)
        genes = load_gene_list(gene_list_path)
        
        if not genes:
            logger.warning(f"No genes found in {gene_list_file}, skipping")
            continue
        
        # Get gene coordinates
        gene_coords = get_gene_coordinates(genes, args.genome_gtf)
        
        if not gene_coords:
            logger.warning(f"No gene coordinates found for {gene_list_file}, skipping")
            continue
        
        # Compute metaprofiles for each condition
        matrices = []
        valid_conditions = []
        
        for condition in conditions:
            # Find the ratio BigWig file
            ratio_file = os.path.join(args.results_dir, "ratio", f"{condition}_S2S_vs_S3_ratio.bw")
            
            if not os.path.exists(ratio_file):
                logger.warning(f"Ratio file not found for {condition}, skipping")
                continue
            
            # Compute metaprofile
            matrix = compute_metaprofile(
                gene_coords, 
                ratio_file, 
                args.upstream_distance, 
                args.downstream_distance, 
                args.bin_size, 
                args.gene_body_bins
            )
            
            if matrix is not None:
                matrices.append(matrix)
                valid_conditions.append(condition)
        
        if not matrices:
            logger.warning(f"No valid metaprofiles computed for {gene_list_name}, skipping")
            continue
        
        # Create metaprofile heatmap
        output_file = os.path.join(args.output_dir, f"{gene_list_name}_metaprofile_heatmap.png")
        create_metaprofile_heatmap(
            matrices, 
            valid_conditions, 
            gene_list_name, 
            output_file, 
            upstream_bins, 
            args.gene_body_bins, 
            downstream_bins
        )
    
    # Create a combined heatmap for all gene lists
    logger.info("Creating combined metaprofile heatmap")
    
    # Process each condition
    for condition in conditions:
        # Find the ratio BigWig file
        ratio_file = os.path.join(args.results_dir, "ratio", f"{condition}_S2S_vs_S3_ratio.bw")
        
        if not os.path.exists(ratio_file):
            logger.warning(f"Ratio file not found for {condition}, skipping")
            continue
        
        matrices = []
        gene_list_names = []
        
        # Process each gene list
        for gene_list_name, gene_list_file in gene_lists.items():
            # Load gene list
            gene_list_path = os.path.join(args.gene_lists_dir, gene_list_file)
            genes = load_gene_list(gene_list_path)
            
            if not genes:
                continue
            
            # Get gene coordinates
            gene_coords = get_gene_coordinates(genes, args.genome_gtf)
            
            if not gene_coords:
                continue
            
            # Compute metaprofile
            matrix = compute_metaprofile(
                gene_coords, 
                ratio_file, 
                args.upstream_distance, 
                args.downstream_distance, 
                args.bin_size, 
                args.gene_body_bins
            )
            
            if matrix is not None:
                matrices.append(matrix)
                gene_list_names.append(gene_list_name)
        
        if not matrices:
            logger.warning(f"No valid metaprofiles computed for {condition}, skipping")
            continue
        
        # Create metaprofile heatmap
        output_file = os.path.join(args.output_dir, f"{condition}_combined_metaprofile_heatmap.png")
        create_metaprofile_heatmap(
            matrices, 
            gene_list_names, 
            condition, 
            output_file, 
            upstream_bins, 
            args.gene_body_bins, 
            downstream_bins
        )
    
    logger.info("Metaprofile heatmap generation complete")

if __name__ == "__main__":
    main() 