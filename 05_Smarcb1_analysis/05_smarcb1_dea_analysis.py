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
import argparse
import gtfparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_dea_genes(dea_file: str) -> pd.DataFrame:
    """Load DEA results and return gene list."""
    df = pd.read_csv(dea_file)
    logger.info(f"Loaded {len(df)} genes from DEA file")
    return df

def load_gene_annotations(gtf_file: str) -> pd.DataFrame:
    """Load gene annotations from GTF file."""
    logger.info("Loading gene annotations...")
    df = gtfparse.read_gtf(gtf_file)
    # Convert to pandas and filter for genes
    df = df.to_pandas()
    genes_df = df[df['feature'] == 'gene'].copy()
    genes_df['tss'] = genes_df.apply(
        lambda x: x['start'] if x['strand'] == '+' else x['end'], 
        axis=1
    )
    return genes_df

def get_smarcb1_signal(bw_file: str, chrom: str, tss: int, window_size: int = 500) -> float:
    """Get average SMARCB1 signal around TSS."""
    with pyBigWig.open(bw_file) as bw:
        try:
            window_start = max(0, tss - window_size)
            window_end = tss + window_size
            
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

def analyze_smarcb1_levels(genes_df: pd.DataFrame, annotations_df: pd.DataFrame, smarcb1_dir: str) -> pd.DataFrame:
    """Analyze SMARCB1 levels for DEA genes."""
    results = []
    
    # Get file paths
    bm_file = os.path.join(smarcb1_dir, "BM3_RPKM.bw")
    bg_files = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]
    
    # Merge gene info with annotations
    merged_df = pd.merge(
        genes_df, 
        annotations_df[['gene_name', 'seqname', 'tss', 'strand']], 
        left_on='gene', 
        right_on='gene_name',
        how='inner'
    )
    
    logger.info(f"Processing {len(merged_df)} genes...")
    
    for idx, gene in merged_df.iterrows():
        if idx % 100 == 0:
            logger.info(f"  Processing gene {idx}/{len(merged_df)}")
        
        # Get signals
        bm_signal = get_smarcb1_signal(bm_file, gene['seqname'], gene['tss'])
        bg_signals = [get_smarcb1_signal(bg_file, gene['seqname'], gene['tss']) for bg_file in bg_files]
        
        # Calculate metrics
        avg_bg = np.mean(bg_signals) if bg_signals else 0
        fold_change = bm_signal / avg_bg if avg_bg > 0 else 0
        is_higher = bm_signal > max(bg_signals) if bg_signals else False
        
        results.append({
            'gene': gene['gene'],
            'chromosome': gene['seqname'],
            'tss': gene['tss'],
            'strand': gene['strand'],
            'log2FoldChange': gene['log2FoldChange'],
            'padj': gene['padj'],
            'bm3_signal': bm_signal,
            'avg_bg_signal': avg_bg,
            'fold_change': fold_change,
            'is_higher_than_bg': is_higher
        })
    
    return pd.DataFrame(results)

def create_summary_plots(results_df: pd.DataFrame, output_dir: str):
    """Create summary plots for the analysis."""
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Distribution of SMARCB1 signals
    plt.figure(figsize=(10, 6))
    sns.histplot(data=results_df, x='fold_change', bins=50)
    plt.title('Distribution of SMARCB1 Fold Change at TSS')
    plt.xlabel('Fold Change (BM3/avg BG)')
    plt.savefig(os.path.join(output_dir, 'smarcb1_fold_change_dist.png'))
    plt.close()
    
    # 2. Correlation with DEA
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=results_df, x='log2FoldChange', y='fold_change')
    plt.title('SMARCB1 Enrichment vs Gene Expression Change')
    plt.xlabel('log2 Fold Change (Expression)')
    plt.ylabel('Fold Change (SMARCB1)')
    plt.savefig(os.path.join(output_dir, 'smarcb1_vs_expression.png'))
    plt.close()
    
    # 3. Save results to CSV
    results_df.to_csv(os.path.join(output_dir, 'smarcb1_tss_analysis.csv'), index=False)
    
    # Calculate some statistics
    high_smarcb1 = results_df[results_df['is_higher_than_bg']]
    logger.info(f"Found {len(high_smarcb1)} genes with higher SMARCB1 in BM3")
    
    # Correlation test
    corr, pval = stats.pearsonr(results_df['log2FoldChange'], results_df['fold_change'])
    logger.info(f"Correlation between expression and SMARCB1: r={corr:.3f}, p={pval:.3e}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze SMARCB1 binding levels at DEA gene TSS regions'
    )
    parser.add_argument('--dea-file', type=str, required=True,
                        help='Path to DEA results CSV file (e.g., DEA_NSC.csv)')
    parser.add_argument('--gtf-file', type=str, required=True,
                        help='Path to GTF annotation file')
    parser.add_argument('--smarcb1-dir', type=str, required=True,
                        help='Directory containing SMARCB1 bigWig files (BM3_RPKM.bw, BG1-3_RPKM.bw)')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Output directory for results')
    args = parser.parse_args()

    # Load data
    dea_genes = load_dea_genes(args.dea_file)
    gene_annotations = load_gene_annotations(args.gtf_file)

    # Analyze SMARCB1 levels
    results = analyze_smarcb1_levels(dea_genes, gene_annotations, args.smarcb1_dir)

    # Create plots and save results
    create_summary_plots(results, args.output_dir)

if __name__ == "__main__":
    main()
