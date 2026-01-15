#!/usr/bin/env python3
"""
Associate Shifted Regions with Promoters

This script takes the BED files containing shifted chromatin regions and associates them
with gene promoters from the mm10 genome using the provided GTF annotation file.
"""

import os
import sys
import pandas as pd
import logging
from collections import defaultdict
import argparse
import glob
import matplotlib.pyplot as plt

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger()

def parse_gtf(gtf_file):
    """
    Parse GTF file and extract gene information.
    
    Returns a dictionary with gene coordinates and information.
    """
    logger.info(f"Parsing GTF file: {gtf_file}")
    genes = {}
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                
                # Parse attributes
                attr_dict = {}
                for attr in fields[8].split(';'):
                    if not attr.strip():
                        continue
                    try:
                        key, value = attr.strip().split(' ', 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        continue
                
                gene_id = attr_dict.get('gene_id', '')
                gene_name = attr_dict.get('gene_name', gene_id)
                gene_type = attr_dict.get('gene_type', '')
                
                genes[gene_id] = {
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'gene_name': gene_name,
                    'gene_type': gene_type
                }
        
        logger.info(f"Successfully parsed {len(genes)} genes from GTF file")
        return genes
    
    except Exception as e:
        logger.error(f"Error parsing GTF file: {e}")
        sys.exit(1)

def find_promoter_overlaps(region, genes, promoter_size=5000, include_gene_body=False):
    """
    Find genes whose promoters overlap with a given region.
    
    Parameters:
    -----------
    region : dict
        Dictionary containing region information (chrom, start, end)
    genes : dict
        Dictionary of gene information
    promoter_size : int
        Number of base pairs upstream of TSS to consider as promoter region
    include_gene_body : bool
        Whether to include gene body in the promoter definition
    
    Returns:
    --------
    list
        List of overlapping gene IDs and their overlap type
    """
    overlapping = []
    
    for gene_id, gene in genes.items():
        if gene['chrom'] != region['chrom']:
            continue
        
        # Define promoter region based on strand
        if gene['strand'] == '+':
            promoter_start = gene['start'] - promoter_size
            promoter_end = gene['start'] + (1 if not include_gene_body else gene['end'] - gene['start'])
        else:  # '-' strand
            promoter_start = gene['end'] - (1 if not include_gene_body else gene['end'] - gene['start'])
            promoter_end = gene['end'] + promoter_size
        
        # Check for overlap with promoter
        if not (region['end'] < promoter_start or region['start'] > promoter_end):
            # Determine overlap type
            if gene['strand'] == '+':
                if region['start'] >= gene['start'] and include_gene_body:
                    overlap_type = 'gene_body'
                else:
                    overlap_type = 'promoter'
            else:  # '-' strand
                if region['end'] <= gene['end'] and include_gene_body:
                    overlap_type = 'gene_body'
                else:
                    overlap_type = 'promoter'
            
            overlapping.append({
                'gene_id': gene_id,
                'gene_name': gene['gene_name'],
                'gene_type': gene['gene_type'],
                'overlap_type': overlap_type,
                'strand': gene['strand'],
                'tss': gene['start'] if gene['strand'] == '+' else gene['end'],
                'distance_to_tss': min(
                    abs(region['start'] - (gene['start'] if gene['strand'] == '+' else gene['end'])),
                    abs(region['end'] - (gene['start'] if gene['strand'] == '+' else gene['end']))
                )
            })
    
    return overlapping

def process_bed_file(bed_file, genes, output_file, promoter_size, include_gene_body):
    """
    Process a BED file and find genes whose promoters overlap with each region.
    """
    logger.info(f"Processing BED file: {bed_file}")
    
    try:
        # Read BED file
        regions_df = pd.read_csv(bed_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'ratio_change'])
        
        if len(regions_df) == 0:
            logger.warning(f"No regions found in {bed_file}")
            return
        
        logger.info(f"Found {len(regions_df)} regions in {bed_file}")
        
        # Process each region
        results = []
        for _, region in regions_df.iterrows():
            overlapping = find_promoter_overlaps(region, genes, promoter_size, include_gene_body)
            
            for gene in overlapping:
                results.append({
                    'chrom': region['chrom'],
                    'start': region['start'],
                    'end': region['end'],
                    'ratio_change': region['ratio_change'],
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['gene_name'],
                    'gene_type': gene['gene_type'],
                    'overlap_type': gene['overlap_type'],
                    'strand': gene['strand'],
                    'tss': gene['tss'],
                    'distance_to_tss': gene['distance_to_tss']
                })
        
        # Create results DataFrame and save
        if results:
            results_df = pd.DataFrame(results)
            
            # Ensure output directory exists
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            results_df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Found {len(results)} gene promoter associations")
            logger.info(f"Results saved to: {output_file}")
            
            # Create summary statistics
            summary_file = output_file.replace('.tsv', '_summary.tsv')
            gene_type_summary = results_df.groupby('gene_type').size().reset_index(name='count')
            gene_type_summary = gene_type_summary.sort_values('count', ascending=False)
            gene_type_summary.to_csv(summary_file, sep='\t', index=False)
            logger.info(f"Summary statistics saved to: {summary_file}")
            
            # Create gene list file (just gene names)
            gene_list_file = output_file.replace('.tsv', '_gene_list.txt')
            unique_genes = results_df['gene_name'].unique()
            with open(gene_list_file, 'w') as f:
                for gene in unique_genes:
                    f.write(f"{gene}\n")
            logger.info(f"Gene list saved to: {gene_list_file} ({len(unique_genes)} unique genes)")
            
            # Create distance to TSS histogram
            if 'distance_to_tss' in results_df.columns:
                plt.figure(figsize=(10, 6))
                plt.hist(results_df['distance_to_tss'], bins=50)
                plt.xlabel('Distance to TSS (bp)')
                plt.ylabel('Number of Regions')
                plt.title(f'Distance to TSS Distribution - {os.path.basename(bed_file)}')
                hist_file = output_file.replace('.tsv', '_tss_distance.png')
                plt.savefig(hist_file, dpi=300, bbox_inches='tight')
                plt.close()
                logger.info(f"TSS distance histogram saved to: {hist_file}")
            
        else:
            logger.warning("No gene promoter associations found")
            
    except Exception as e:
        logger.error(f"Error processing BED file: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return

def main():
    parser = argparse.ArgumentParser(description='Associate shifted regions with gene promoters')
    parser.add_argument('--results-dir', required=True, help='Directory containing shifted regions BED files')
    parser.add_argument('--gtf-file', required=True, 
                       help='Path to GTF annotation file')
    parser.add_argument('--promoter-size', type=int, default=5000,
                       help='Number of base pairs upstream of TSS to consider as promoter region')
    parser.add_argument('--include-gene-body', action='store_true',
                       help='Include gene body in promoter definition')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # Check if results directory exists
    if not os.path.exists(args.results_dir):
        logger.error(f"Results directory does not exist: {args.results_dir}")
        sys.exit(1)
    
    # Check if GTF file exists
    if not os.path.exists(args.gtf_file):
        logger.error(f"GTF file does not exist: {args.gtf_file}")
        sys.exit(1)
    
    logger.info(f"Starting gene promoter association analysis")
    logger.info(f"Results directory: {args.results_dir}")
    logger.info(f"GTF file: {args.gtf_file}")
    logger.info(f"Promoter size: {args.promoter_size} bp upstream of TSS")
    logger.info(f"Include gene body: {args.include_gene_body}")
    
    # Parse GTF file
    genes = parse_gtf(args.gtf_file)
    
    # Create output directory
    output_dir = os.path.join(args.results_dir, 'promoter_associations')
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Find all BED files in the shifted regions directory
    bed_files = glob.glob(os.path.join(args.results_dir, '*.bed'))
    
    if not bed_files:
        logger.error(f"No BED files found in {args.results_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(bed_files)} BED files to process")
    
    # Process each BED file in the shifted regions directory
    for bed_file in bed_files:
        filename = os.path.basename(bed_file)
        output_file = os.path.join(output_dir, filename.replace('.bed', '_promoters.tsv'))
        process_bed_file(bed_file, genes, output_file, args.promoter_size, args.include_gene_body)
    
    logger.info("Analysis complete")

if __name__ == "__main__":
    main() 