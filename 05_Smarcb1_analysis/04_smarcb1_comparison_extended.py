#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import os
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_gene_lists(base_dir: str) -> dict:
    """Load the three gene lists with their peak regions."""
    gene_lists = {}
    for category in ['up', 'down', 'no_deg']:
        file_path = os.path.join(base_dir, f'extended_{category}.csv')
        df = pd.read_csv(file_path)
        # Filter out rows without coordinates
        df = df.dropna(subset=['seqnames', 'start', 'end'])
        # Convert coordinates to integers
        df['start'] = df['start'].astype(float).astype(int)
        df['end'] = df['end'].astype(float).astype(int)
        gene_lists[category] = df
        logger.info(f"Loaded {len(df)} regions for {category}")
    return gene_lists

def get_smarcb1_signal(bw_file: str, region: pd.Series, window_size: int = 500) -> float:
    """Get average SMARCB1 signal for a region."""
    with pyBigWig.open(bw_file) as bw:
        try:
            chrom = str(region['seqnames'])
            center = (int(region['start']) + int(region['end'])) // 2
            window_start = max(0, center - window_size)
            window_end = center + window_size
            
            # Verify coordinates are valid
            if window_start >= window_end:
                logger.warning(f"Invalid coordinates for region {chrom}:{window_start}-{window_end}")
                return 0.0
                
            values = bw.stats(chrom, window_start, window_end, type="mean")
            if values is None or values[0] is None:
                return 0.0
                
            return values[0]
            
        except Exception as e:
            logger.warning(f"Error processing region {chrom}: {str(e)}")
            return 0.0

def analyze_smarcb1_levels(gene_lists: dict, smarcb1_dir: str) -> pd.DataFrame:
    """Analyze SMARCB1 levels for each category."""
    results = []
    
    # Get file paths
    bm_file = os.path.join(smarcb1_dir, "BM3_RPKM.bw")
    bg_files = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]
    
    for category, df in gene_lists.items():
        logger.info(f"Processing {category} regions...")
        
        for idx, region in df.iterrows():
            if idx % 100 == 0:
                logger.info(f"  Processing region {idx}/{len(df)}")
            
            # Get signals
            bm_signal = get_smarcb1_signal(bm_file, region)
            bg_signals = [get_smarcb1_signal(bg_file, region) for bg_file in bg_files]
            
            # Calculate metrics
            avg_bg = np.mean(bg_signals) if bg_signals else 0
            fold_change = bm_signal / avg_bg if avg_bg > 0 else 0
            is_higher = bm_signal > max(bg_signals) if bg_signals else False
            
            results.append({
                'category': category,
                'gene': region.get('Gene', 'Unknown'),
                'chromosome': region['seqnames'],
                'start': region['start'],
                'end': region['end'],
                'bm3_signal': bm_signal,
                'avg_bg_signal': avg_bg,
                'fold_change': fold_change,
                'is_higher_than_bg': is_higher
            })
    
    return pd.DataFrame(results)

def create_summary_plots(results_df: pd.DataFrame, output_dir: str):
    """Create summary plots for the analysis."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save genes with higher SMARCB1 in BM3 for each category
    for category in ['up', 'down', 'no_deg']:
        category_df = results_df[results_df['category'] == category]
        higher_genes = category_df[category_df['is_higher_than_bg']].copy()
        
        # Sort by fold change for better interpretation
        higher_genes = higher_genes.sort_values('fold_change', ascending=False)
        
        # Save to CSV
        output_file = os.path.join(output_dir, f'{category}_higher_smarcb1_genes.csv')
        higher_genes.to_csv(output_file, index=False)
        
        # Log statistics
        total_genes = len(category_df)
        higher_count = len(higher_genes)
        percentage = (higher_count / total_genes * 100) if total_genes > 0 else 0
        logger.info(f"\n{category.upper()} genes statistics:")
        logger.info(f"Total genes: {total_genes}")
        logger.info(f"Genes with higher SMARCB1 in BM3: {higher_count} ({percentage:.2f}%)")
    
    # Original plotting code
    # 1. Fold change distribution
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=results_df, x='category', y='fold_change')
    plt.title('SMARCB1 Signal Fold Change (BM3/BG) by Category')
    plt.yscale('log')
    plt.ylabel('Fold Change (log scale)')
    plt.savefig(os.path.join(output_dir, 'fold_change_distribution.png'))
    plt.close()
    
    # 2. Percentage of regions with higher signal in BM3
    summary = results_df.groupby('category').agg({
        'is_higher_than_bg': ['count', 'sum']
    })
    summary.columns = ['total', 'higher_in_bm3']
    summary['percentage'] = (summary['higher_in_bm3'] / summary['total']) * 100
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x=summary.index, y='percentage', data=summary.reset_index())
    plt.title('Percentage of Regions with Higher SMARCB1 in BM3')
    plt.ylabel('Percentage (%)')
    plt.savefig(os.path.join(output_dir, 'percentage_higher_in_bm3.png'))
    plt.close()
    
    # Save summary statistics
    summary.to_csv(os.path.join(output_dir, 'summary_statistics.csv'))
    
    return summary

def main():
    # Set paths - configure these for your data
    base_dir = ""  # e.g., "/path/to/cutandtag/analysis"
    data_dir = os.path.join(base_dir, "data")
    smarcb1_dir = ""  # e.g., "/path/to/smarcb1/bigwig"
    output_dir = os.path.join(base_dir, "results/smarcb1_comparison_2")
    
    # Load gene lists
    gene_lists = load_gene_lists(data_dir)
    
    # Analyze SMARCB1 levels
    results_df = analyze_smarcb1_levels(gene_lists, smarcb1_dir)
    
    # Save detailed results
    os.makedirs(output_dir, exist_ok=True)
    results_df.to_csv(os.path.join(output_dir, 'detailed_results.csv'), index=False)
    
    # Create and save summary plots
    summary = create_summary_plots(results_df, output_dir)
    
    # Print summary
    logger.info("\nAnalysis Summary:")
    logger.info(summary.to_string())

if __name__ == "__main__":
    main() 