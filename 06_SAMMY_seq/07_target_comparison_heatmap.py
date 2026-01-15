#!/usr/bin/env python3
"""
SAMMY-seq Mecp2 Target Comparison Heatmap Generator

This script creates metaprofile heatmaps for Mecp2 target and non-target genes in NSC and NEU conditions,
showing the log2 ratio of S2S/S3 across gene bodies (from -5kb to TSS to TES to +5kb).

The heatmaps use a blue-white-red color scheme where:
- Blue represents heterochromatin (log2(S2S/S3) < 0)
- White represents undefined regions (log2(S2S/S3) = 0)
- Red represents euchromatin (log2(S2S/S3) > 0)

The script generates a 2x2 grid of heatmaps for:
1. NEU Mecp2 target genes (NEU_exo.txt)
2. NEU non-target genes (no_NEU_exo.txt)
3. NSC Mecp2 target genes (NSC_exo.txt)
4. NSC non-target genes (no_NSC_exo.txt)
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import gaussian_filter1d
import pickle
import hashlib
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("mecp2_target_comparison_heatmap.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('SAMMY-seq-Mecp2-Comparison')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='SAMMY-seq Mecp2 Target Comparison Heatmap Generator')
    parser.add_argument('--results-dir', required=True, help='Directory containing SAMMY-seq analysis results')
    parser.add_argument('--output-dir', required=True, help='Output directory for heatmaps')
    parser.add_argument('--genome-gtf', required=True, help='GTF file with gene annotations')
    parser.add_argument('--gene-lists-dir', required=True, help='Directory containing gene list files')
    parser.add_argument('--upstream-distance', type=int, default=5000, help='Distance upstream of TSS (default: 5000)')
    parser.add_argument('--downstream-distance', type=int, default=5000, help='Distance downstream of TES (default: 5000)')
    parser.add_argument('--bin-size', type=int, default=100, help='Bin size for metaprofile (default: 100)')
    parser.add_argument('--gene-body-bins', type=int, default=100, help='Number of bins for gene body (default: 100)')
    parser.add_argument('--force-recompute', action='store_true', help='Force recomputation of metaprofiles even if cache exists')
    parser.add_argument('--cache-dir', default=None, help='Directory to store cached results (default: output_dir/cache)')
    parser.add_argument('--plot-only', action='store_true', help='Only generate plots using cached data, skip computation')
    return parser.parse_args()

def generate_cache_key(args, gene_lists):
    """Generate a unique cache key based on input parameters"""
    # Create a string representation of the parameters that affect computation
    param_str = (
        f"{args.genome_gtf}_{args.gene_lists_dir}_{args.results_dir}_"
        f"{args.upstream_distance}_{args.downstream_distance}_"
        f"{args.bin_size}_{args.gene_body_bins}_"
        f"{'_'.join(sorted(gene_lists.keys()))}_"
        f"{'_'.join(sorted(gene_lists.values()))}"
    )
    
    # Generate a hash of the parameter string
    return hashlib.md5(param_str.encode()).hexdigest()

def get_cache_path(cache_dir, cache_key, suffix=""):
    """Get the path to a cache file"""
    return os.path.join(cache_dir, f"{cache_key}{suffix}.pkl")

def save_to_cache(cache_dir, cache_key, data, suffix=""):
    """Save data to cache"""
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = get_cache_path(cache_dir, cache_key, suffix)
    
    logger.info(f"Saving data to cache: {cache_path}")
    start_time = time.time()
    
    try:
        with open(cache_path, 'wb') as f:
            pickle.dump(data, f)
        logger.info(f"Cache saved successfully in {time.time() - start_time:.2f} seconds")
        return True
    except Exception as e:
        logger.error(f"Error saving to cache: {e}")
        return False

def load_from_cache(cache_dir, cache_key, suffix=""):
    """Load data from cache"""
    cache_path = get_cache_path(cache_dir, cache_key, suffix)
    
    if not os.path.exists(cache_path):
        logger.info(f"Cache file not found: {cache_path}")
        return None
    
    logger.info(f"Loading data from cache: {cache_path}")
    start_time = time.time()
    
    try:
        with open(cache_path, 'rb') as f:
            data = pickle.load(f)
        logger.info(f"Cache loaded successfully in {time.time() - start_time:.2f} seconds")
        return data
    except Exception as e:
        logger.error(f"Error loading from cache: {e}")
        return None

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
    """
    Get gene coordinates (TSS, TES, strand) from GTF file.
    
    This function parses a GTF file to extract the genomic coordinates 
    (chromosome, start, end, strand, TSS, TES) for a given list of genes.
    
    Args:
        genes (list): A list of gene names to extract coordinates for.
        gtf_file (str): Path to the GTF file containing gene annotations.
    
    Returns:
        dict: A dictionary where keys are gene names and values are dictionaries
              containing the gene's coordinates ('chrom', 'start', 'end', 'strand', 'tss', 'tes').
              Returns an empty dictionary if an error occurs or no coordinates are found.
    """
    logger.info(f"Getting coordinates for {len(genes)} genes from GTF file")
    
    # Initialize an empty dictionary to store gene coordinates
    gene_coords = {}
    
    try:
        # Use pybedtools to open and parse the GTF file
        gtf = pybedtools.BedTool(gtf_file)
        
        # Filter the GTF file to only include 'gene' features
        gene_features = gtf.filter(lambda x: x[2] == 'gene')
        
        # Iterate through each gene feature in the GTF file
        for feature in gene_features:
            # Extract attributes from the GTF feature's attribute string
            # The attribute string is split into key-value pairs, and a dictionary is created
            attributes = dict(item.strip().split(' ') for item in feature.fields[8].split(';') if item.strip())
            
            # Extract the gene name from the attributes
            # It tries 'gene_name' first, then 'gene_id' if 'gene_name' is not found
            gene_name = None
            for attr in ['gene_name', 'gene_id']:
                if attr in attributes:
                    gene_name = attributes[attr].strip('"\'')  # Remove quotes from the gene name
                    break
            
            # If the extracted gene name is in the list of genes we're interested in
            if gene_name and gene_name in genes:
                # Extract genomic information from the GTF feature
                chrom = feature.chrom  # Chromosome
                start = int(feature.start)  # Start coordinate (0-based)
                end = int(feature.end)  # End coordinate
                strand = feature.strand  # Strand (+ or -)
                
                # Determine TSS (Transcription Start Site) and TES (Transcription End Site)
                # based on the gene's strand
                if strand == '+':
                    tss = start  # TSS is the start coordinate for positive strand genes
                    tes = end  # TES is the end coordinate for positive strand genes
                else:  # strand == '-'
                    tss = end  # TSS is the end coordinate for negative strand genes
                    tes = start  # TES is the start coordinate for negative strand genes
                
                # Store the extracted gene coordinates in the dictionary
                gene_coords[gene_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'tss': tss,
                    'tes': tes
                }
        
        logger.info(f"Found coordinates for {len(gene_coords)} out of {len(genes)} genes")
        return gene_coords
    
    except Exception as e:
        # Log any errors that occur during the process
        logger.error(f"Error getting gene coordinates: {e}")
        return {}  # Return an empty dictionary if an error occurred

def compute_metaprofile(gene_coords, bigwig_file, upstream_distance, downstream_distance, bin_size, gene_body_bins):
    """
    Compute metaprofile for a set of genes from a BigWig file.
    
    This function calculates the average signal profile around a set of genes, 
    using data from a BigWig file. It divides the region around each gene 
    (upstream, gene body, and downstream) into bins and calculates the average 
    signal within each bin.
    
    Args:
        gene_coords (dict): A dictionary containing the genomic coordinates of each gene.
                             Keys are gene names, and values are dictionaries with keys 
                             'chrom', 'tss', 'tes', and 'strand'.
        bigwig_file (str): Path to the BigWig file containing the signal data.
        upstream_distance (int): Distance upstream of the TSS to consider (bp).
        downstream_distance (int): Distance downstream of the TES to consider (bp).
        bin_size (int): Size of each bin in the metaprofile (bp).
        gene_body_bins (int): Number of bins to use for the gene body.
    
    Returns:
        np.ndarray: A matrix where each row represents a gene and each column 
                      represents a bin in the metaprofile.  Returns None if an error occurs.
                      The matrix is filled with NaN values if no data is available for a bin.
    """
    logger.info(f"Computing metaprofile from {bigwig_file}")
    
    try:
        # Open the BigWig file for reading
        bw = pyBigWig.open(bigwig_file)
        
        # Calculate the number of bins for the upstream and downstream regions
        upstream_bins = int(upstream_distance / bin_size)
        downstream_bins = int(downstream_distance / bin_size)
        
        # Calculate the total number of bins in the metaprofile
        total_bins = upstream_bins + gene_body_bins + downstream_bins
        
        # Initialize a matrix to store the metaprofile data for all genes
        matrix = np.zeros((len(gene_coords), total_bins))
        matrix.fill(np.nan)  # Initialize all values to NaN
        
        # Iterate through each gene and its coordinates
        for i, (gene_name, coords) in enumerate(tqdm(gene_coords.items(), desc="Processing genes")):
            chrom = coords['chrom']  # Chromosome
            tss = coords['tss']  # Transcription start site
            tes = coords['tes']  # Transcription end site
            strand = coords['strand']  # Strand
            
            # Skip genes located on chromosomes not present in the BigWig file
            if chrom not in bw.chroms():
                continue
            
            # Calculate the start and end coordinates for the upstream, gene body, and downstream regions
            # based on the gene's strand
            if strand == '+':
                upstream_start = max(0, tss - upstream_distance)
                upstream_end = tss
                gene_body_start = tss
                gene_body_end = tes
                downstream_start = tes
                downstream_end = min(bw.chroms()[chrom], tes + downstream_distance)
            else:  # strand == '-'
                upstream_start = max(0, tes)  # Upstream is after TES on the - strand
                upstream_end = min(bw.chroms()[chrom], tes + upstream_distance)
                gene_body_start = tes
                gene_body_end = tss
                downstream_start = max(0, tss - downstream_distance)
                downstream_end = tss
            
            # Retrieve the signal values for the upstream region
            if upstream_end > upstream_start:
                try:
                    upstream_values = bw.stats(chrom, upstream_start, upstream_end, 
                                              nBins=upstream_bins, type="mean")
                    
                    # Reverse the order of values for genes on the negative strand to maintain 5' to 3' direction
                    if strand == '-':
                        upstream_values = upstream_values[::-1] if upstream_values else upstream_values
                    
                    # Store the signal values in the matrix
                    for j, val in enumerate(upstream_values):
                        if val is not None:
                            matrix[i, j] = val
                except Exception as e:
                    logger.debug(f"Error processing upstream region for {gene_name}: {e}")
            
            # Retrieve the signal values for the gene body
            if gene_body_end > gene_body_start:
                try:
                    gene_body_values = bw.stats(chrom, gene_body_start, gene_body_end, 
                                               nBins=gene_body_bins, type="mean")
                    
                    # Reverse the order of values for genes on the negative strand to maintain 5' to 3' direction
                    if strand == '-':
                        gene_body_values = gene_body_values[::-1] if gene_body_values else gene_body_values
                    
                    # Store the signal values in the matrix
                    for j, val in enumerate(gene_body_values):
                        if val is not None:
                            matrix[i, upstream_bins + j] = val
                except Exception as e:
                    logger.debug(f"Error processing gene body for {gene_name}: {e}")
            
            # Retrieve the signal values for the downstream region
            if downstream_end > downstream_start:
                try:
                    downstream_values = bw.stats(chrom, downstream_start, downstream_end, 
                                                nBins=downstream_bins, type="mean")
                    
                    # Reverse the order of values for genes on the negative strand to maintain 5' to 3' direction
                    if strand == '-':
                        downstream_values = downstream_values[::-1] if downstream_values else downstream_values
                    
                    # Store the signal values in the matrix
                    for j, val in enumerate(downstream_values):
                        if val is not None:
                            matrix[i, upstream_bins + gene_body_bins + j] = val
                except Exception as e:
                    logger.debug(f"Error processing downstream region for {gene_name}: {e}")
        
        # Close the BigWig file
        bw.close()
        
        # Return the matrix containing the metaprofile data
        return matrix
    
    except Exception as e:
        # Log any errors that occur during the process
        logger.error(f"Error computing metaprofile: {e}")
        return None  # Return None if an error occurred

def create_comparison_heatmap(gene_sets, condition_matrices, output_file, 
                             upstream_bins, gene_body_bins, downstream_bins,
                             bin_size, upstream_distance, downstream_distance):
    """
    Create a 2x2 grid of heatmaps comparing Mecp2 target and non-target genes.
    
    This function generates a figure containing four heatmaps, each representing
    the average metaprofile of a specific gene set (NEU target, NEU non-target,
    NSC target, NSC non-target) across a region spanning upstream, gene body,
    and downstream of the genes. It also generates corresponding line plots
    for each gene set to visualize the average profiles more clearly.
    
    Args:
        gene_sets (list): A list of gene set names (e.g., 'NEU_target', 'NSC_non_target').
        condition_matrices (dict): A dictionary where keys are gene set names and
                                   values are dictionaries of condition-specific matrices.
                                   Each matrix contains the metaprofile data for the genes
                                   in that set under a specific condition.
        output_file (str): The path to save the generated heatmap figure.
        upstream_bins (int): The number of bins representing the upstream region.
        gene_body_bins (int): The number of bins representing the gene body.
        downstream_bins (int): The number of bins representing the downstream region.
        bin_size (int): Size of each bin in the metaprofile (bp).
        upstream_distance (int): Distance upstream of the TSS to consider (bp).
        downstream_distance (int): Distance downstream of the TES to consider (bp).
    """
    logger.info(f"Creating comparison heatmap")
    
    try:
        # Create a custom colormap with more contrast: blue for heterochromatin, white for undefined, red for euchromatin
        cmap = LinearSegmentedColormap.from_list(
            'custom_diverging',
            [(0, '#0000FF'), (0.45, '#8080FF'), (0.5, 'white'), (0.55, '#FF8080'), (1, '#FF0000')],
            N=256
        )
        
        # Set up the figure with 2x2 grid
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Set vmin and vmax for consistent color scaling across all heatmaps
        # Use tighter bounds to enhance contrast
        vmin = -0.05
        vmax = 0.05
        
        # Define the order of gene sets to display
        gene_set_order = ['NEU_target', 'NEU_non_target', 'NSC_target', 'NSC_non_target']
        
        # Define display names for gene sets
        gene_set_display = {
            'NEU_target': 'NEU Mecp2 target genes',
            'NEU_non_target': 'NEU non-target genes',
            'NSC_target': 'NSC Mecp2 target genes',
            'NSC_non_target': 'NSC non-target genes'
        }
        
        # Define the positions for each gene set in the 2x2 grid
        positions = {
            'NEU_target': (0, 0),  # Row 0, Column 0
            'NEU_non_target': (0, 1), # Row 0, Column 1
            'NSC_target': (1, 0),  # Row 1, Column 0
            'NSC_non_target': (1, 1)  # Row 1, Column 1
        }
        
        # Iterate through each gene set to create its corresponding heatmap
        for gene_set in gene_set_order:
            row, col = positions[gene_set]  # Get the row and column for the current gene set
            
            # Get the matrices for this gene set
            matrices = condition_matrices.get(gene_set, {})
            
            # If no matrices are found for this gene set, log a warning and skip to the next gene set
            if not matrices:
                logger.warning(f"No matrices found for {gene_set}, skipping")
                continue
            
            # Initialize variables to store the averaged matrix data and the conditions with valid data
            matrix_data = None
            conditions_with_data = []
            
            # Iterate through each condition to compute the average profile
            for condition in ['NSC_GFP', 'NSC_M2', 'Neu_GFP', 'Neu_M2']:  # Ensure consistent order
                # Check if the condition exists in the matrices and if the data is not None
                if condition in matrices and matrices[condition] is not None:
                    # Average across genes for each position
                    avg_profile = np.nanmean(matrices[condition], axis=0)
                    
                    # Apply smoothing to reduce noise and highlight patterns
                    avg_profile = gaussian_filter1d(avg_profile, sigma=2)
                    
                    # If matrix_data is None, initialize it with the first average profile
                    if matrix_data is None:
                        matrix_data = avg_profile.reshape(1, -1)
                    # Otherwise, stack the current average profile onto the existing matrix_data
                    else:
                        matrix_data = np.vstack([matrix_data, avg_profile.reshape(1, -1)])
                    
                    # Append the condition to the list of conditions with valid data
                    conditions_with_data.append(condition)
            
            # If no valid data is found for the gene set, log a warning and skip to the next gene set
            if matrix_data is None or matrix_data.size == 0:
                logger.warning(f"No valid data for {gene_set}, skipping")
                continue
            
            # Create heatmap
            ax = axes[row, col]  # Get the axis object for the current gene set
            
            # Create the heatmap with enhanced parameters
            sns.heatmap(matrix_data, cmap=cmap, vmin=vmin, vmax=vmax, 
                       cbar=True, ax=ax, yticklabels=conditions_with_data,
                       xticklabels=False, rasterized=True)
            
            # Add vertical lines to mark TSS and TES
            ax.axvline(x=upstream_bins, color='black', linestyle='--', linewidth=1)
            ax.axvline(x=upstream_bins + gene_body_bins, color='black', linestyle='--', linewidth=1)
            
            # Set title with larger font
            ax.set_title(gene_set_display[gene_set], fontsize=14, fontweight='bold')
            
            # Set x-axis labels
            total_bins = upstream_bins + gene_body_bins + downstream_bins
            xticks = [0, upstream_bins, upstream_bins + gene_body_bins, total_bins]
            xticklabels = [f"-{upstream_bins/10}kb", "TSS", "TES", f"+{downstream_bins/10}kb"]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels, fontsize=12)
            
            # Add colorbar label
            cbar = ax.collections[0].colorbar
            cbar.set_label('log2(S2S/S3)', fontsize=12)
            
            # Add a line plot below the heatmap to better visualize differences
            divider = make_axes_locatable(ax)
            ax_line = divider.append_axes("bottom", size="30%", pad=0.5)
            
            # Plot each condition as a line
            for i, condition in enumerate(conditions_with_data):
                ax_line.plot(matrix_data[i], label=condition, linewidth=2)
            
            # Add vertical lines to mark TSS and TES in line plot
            ax_line.axvline(x=upstream_bins, color='black', linestyle='--', linewidth=1)
            ax_line.axvline(x=upstream_bins + gene_body_bins, color='black', linestyle='--', linewidth=1)
            
            # Set y-limits for line plot to enhance visibility of differences
            y_min = np.min(matrix_data) - 0.01
            y_max = np.max(matrix_data) + 0.01
            ax_line.set_ylim(y_min, y_max)
            
            # Add legend to line plot
            ax_line.legend(loc='upper right', fontsize=10)
            
            # Set x-axis labels for line plot
            ax_line.set_xticks(xticks)
            ax_line.set_xticklabels(xticklabels, fontsize=12)
            
            # Set y-axis label for line plot
            ax_line.set_ylabel('log2(S2S/S3)', fontsize=12)
        
        # Set overall title
        plt.suptitle("Mecp2 repressed genes SAMMY metaprofile", fontsize=18, fontweight='bold')
        
        # Add a color legend and bin information to the figure
        fig.subplots_adjust(bottom=0.2)  # Make room for the legend at the bottom
        
        # Create a new axis for the color legend
        color_legend_ax = fig.add_axes([0.15, 0.08, 0.7, 0.05])
        
        # Create a gradient for the color legend
        gradient = np.linspace(vmin, vmax, 256).reshape(1, -1)
        color_legend_ax.imshow(gradient, aspect='auto', cmap=cmap)
        color_legend_ax.set_xticks([0, 128, 255])
        color_legend_ax.set_xticklabels([f"{vmin} (Heterochromatin)", "0 (Undefined)", f"{vmax} (Euchromatin)"])
        color_legend_ax.set_yticks([])
        color_legend_ax.set_title("Color Scale: log2(S2S/S3)", fontsize=12)
        
        # Add bin information
        bin_info_text = (
            f"Bin Information:\n"
            f"Upstream region: {upstream_bins} bins (each {bin_size}bp) = {upstream_distance}bp\n"
            f"Gene body: {gene_body_bins} bins (scaled to gene length)\n"
            f"Downstream region: {downstream_bins} bins (each {bin_size}bp) = {downstream_distance}bp"
        )
        fig.text(0.5, 0.02, bin_info_text, ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
        
        # Adjust layout
        plt.tight_layout(rect=[0, 0.2, 1, 0.96])  # Make room for suptitle and legend
        
        # Save the figure
        plt.savefig(output_file, dpi=300)
        
        # Create a second figure with individual line plots for clearer comparison
        plt.figure(figsize=(16, 12))
        
        # Iterate through each gene set to create individual line plots
        for i, gene_set in enumerate(gene_set_order):
            plt.subplot(2, 2, i+1)  # Create a subplot for each gene set
            
            matrices = condition_matrices.get(gene_set, {})
            if not matrices:
                continue
                
            # Plot the average profile for each condition
            for condition in ['NSC_GFP', 'NSC_M2', 'Neu_GFP', 'Neu_M2']:
                if condition in matrices and matrices[condition] is not None:
                    avg_profile = np.nanmean(matrices[condition], axis=0)
                    avg_profile = gaussian_filter1d(avg_profile, sigma=2)  # Apply smoothing
                    plt.plot(avg_profile, label=condition, linewidth=2)
            
            # Add vertical lines to mark TSS and TES
            plt.axvline(x=upstream_bins, color='black', linestyle='--', linewidth=1)
            plt.axvline(x=upstream_bins + gene_body_bins, color='black', linestyle='--', linewidth=1)
            
            # Set title
            plt.title(gene_set_display[gene_set], fontsize=14, fontweight='bold')
            
            # Set x-axis labels
            plt.xticks(xticks, xticklabels, fontsize=12)
            
            # Set y-axis label
            plt.ylabel('log2(S2S/S3)', fontsize=12)
            
            # Add legend
            plt.legend(loc='upper right', fontsize=10)
            
            # Add grid for better readability
            plt.grid(True, linestyle='--', alpha=0.3)
        
        plt.suptitle("Mecp2 repressed genes SAMMY metaprofile - Line Plots", fontsize=18, fontweight='bold')
        
        # Add color interpretation and bin information to the line plot figure
        plt.subplots_adjust(bottom=0.2)  # Make room for the legend
        
        # Add color interpretation text
        color_info_text = (
            "Color Interpretation:\n"
            "Blue (log2(S2S/S3) < 0): Heterochromatin\n"
            "White (log2(S2S/S3) = 0): Undefined regions\n"
            "Red (log2(S2S/S3) > 0): Euchromatin"
        )
        plt.figtext(0.25, 0.08, color_info_text, ha='left', fontsize=12, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
        
        # Add bin information
        bin_info_text = (
            f"Bin Information:\n"
            f"Upstream region: {upstream_bins} bins (each {bin_size}bp) = {upstream_distance}bp\n"
            f"Gene body: {gene_body_bins} bins (scaled to gene length)\n"
            f"Downstream region: {downstream_bins} bins (each {bin_size}bp) = {downstream_distance}bp"
        )
        plt.figtext(0.75, 0.08, bin_info_text, ha='right', fontsize=12, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
        
        plt.tight_layout(rect=[0, 0.2, 1, 0.96])
        
        # Save the line plot figure
        line_plot_file = os.path.splitext(output_file)[0] + "_line_plots.png"
        plt.savefig(line_plot_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved comparison heatmap to {output_file}")
        logger.info(f"Saved line plots to {line_plot_file}")
    
    except Exception as e:
        logger.error(f"Error creating comparison heatmap: {e}")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up cache directory
    if args.cache_dir is None:
        args.cache_dir = os.path.join(args.output_dir, "cache")
    os.makedirs(args.cache_dir, exist_ok=True)
    
    # Define gene lists to process
    gene_lists = {
        'NEU_target': 'NEU_exo.txt',
        'NEU_non_target': 'not_NEU_exo.txt',
        'NSC_target': 'NSC_exo.txt',
        'NSC_non_target': 'not_NSC_exo.txt'
    }
    
    # Define conditions to process
    conditions = ['NSC_GFP', 'NSC_M2', 'Neu_GFP', 'Neu_M2']
    
    # Calculate bin counts
    upstream_bins = int(args.upstream_distance / args.bin_size)
    downstream_bins = int(args.downstream_distance / args.bin_size)
    
    # Generate a unique cache key based on input parameters
    cache_key = generate_cache_key(args, gene_lists)
    logger.info(f"Cache key: {cache_key}")
    
    # Initialize variables for gene coordinates and matrices
    gene_coords = {}
    condition_matrices = {}
    
    # Check if we should use cached data
    use_cached_data = not args.force_recompute
    
    if use_cached_data:
        # Try to load gene coordinates from cache
        gene_coords = load_from_cache(args.cache_dir, cache_key, "_gene_coords")
        
        # Try to load matrices from cache
        condition_matrices = load_from_cache(args.cache_dir, cache_key, "_matrices")
        
        # If either cache is missing, we need to recompute
        if gene_coords is None or condition_matrices is None:
            use_cached_data = False
            logger.info("Cache incomplete or missing, will recompute")
    
    # Compute data if needed
    if not use_cached_data and not args.plot_only:
        logger.info("Computing metaprofiles from scratch")
        
        # Dictionary to store gene coordinates for each gene set
        gene_coords = {}
        
        # Dictionary to store matrices for each gene set and condition
        condition_matrices = {}
        
        # Process each gene list
        for gene_set, gene_list_file in gene_lists.items():
            logger.info(f"Processing gene list: {gene_set}")
            
            # Load gene list
            gene_list_path = os.path.join(args.gene_lists_dir, gene_list_file)
            genes = load_gene_list(gene_list_path)
            
            if not genes:
                logger.warning(f"No genes found in {gene_list_file}, skipping")
                continue
            
            # Get gene coordinates
            gene_coords[gene_set] = get_gene_coordinates(genes, args.genome_gtf)
            
            if not gene_coords[gene_set]:
                logger.warning(f"No gene coordinates found for {gene_list_file}, skipping")
                continue
            
            # Initialize matrices dictionary for this gene set
            condition_matrices[gene_set] = {}
            
            # Compute metaprofiles for each condition
            for condition in conditions:
                # Find the ratio BigWig file
                ratio_file = os.path.join(args.results_dir, "ratio", f"{condition}_S2S_vs_S3_ratio.bw")
                
                if not os.path.exists(ratio_file):
                    logger.warning(f"Ratio file not found for {condition}, skipping")
                    continue
                
                # Compute metaprofile
                matrix = compute_metaprofile(
                    gene_coords[gene_set], 
                    ratio_file, 
                    args.upstream_distance, 
                    args.downstream_distance, 
                    args.bin_size, 
                    args.gene_body_bins
                )
                
                if matrix is not None:
                    # Apply data normalization to enhance differences
                    # Center the data around zero
                    matrix = matrix - np.nanmean(matrix)
                    
                    condition_matrices[gene_set][condition] = matrix
        
        # Save results to cache
        save_to_cache(args.cache_dir, cache_key, gene_coords, "_gene_coords")
        save_to_cache(args.cache_dir, cache_key, condition_matrices, "_matrices")
    elif args.plot_only and (gene_coords is None or condition_matrices is None):
        logger.error("Cannot generate plots: cached data not found and --plot-only specified")
        sys.exit(1)
    
    # Create comparison heatmap
    output_file = os.path.join(args.output_dir, "mecp2_target_comparison_heatmap.png")
    create_comparison_heatmap(
        gene_lists, 
        condition_matrices, 
        output_file, 
        upstream_bins, 
        args.gene_body_bins, 
        downstream_bins,
        args.bin_size,
        args.upstream_distance,
        args.downstream_distance
    )
    
    logger.info("Mecp2 target comparison heatmap generation complete")

if __name__ == "__main__":
    main() 