import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
from typing import Dict, List, Tuple
from pathlib import Path
import scipy.stats as stats

# Set up more detailed logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class RNASeqIntegrator:
    def __init__(self, enrichment_file: str, rna_seq_file: str):
        """
        Initialize with pre-calculated enrichment data and RNA-seq results
        
        Args:
            enrichment_file: Path to the enrichment results from analyze_mecp2_cpg_enrichment.py
            rna_seq_file: Path to RNA-seq differential expression results
        """
        # Load pre-calculated enrichment data
        self.enrichment_df = pd.read_csv(enrichment_file)
        logger.info(f"Loaded {len(self.enrichment_df)} CpG regions with enrichment data")
        
        # Load RNA-seq data
        self.rna_seq_df = pd.read_csv(rna_seq_file)
        logger.info(f"Loaded {len(self.rna_seq_df)} genes with RNA-seq data")
        
        # Verify enrichment data format
        required_cols = ['chr', 'start', 'end', 'enrichment', 'significant', 
                        'exo_signal', 'endo_signal', 'binding_type']
        missing_cols = [col for col in required_cols if col not in self.enrichment_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in enrichment data: {missing_cols}")
        
        # Filter for significant MeCP2 enrichment
        self.significant_regions = self.enrichment_df[
            (self.enrichment_df['significant']) & 
            (self.enrichment_df['binding_type'].isin(['both', 'exo_only']))
        ]
        logger.info(f"Found {len(self.significant_regions)} significantly enriched CpG regions")
        
        # Log binding type distribution
        binding_counts = self.enrichment_df['binding_type'].value_counts()
        logger.info("\nBinding type distribution:")
        for binding_type, count in binding_counts.items():
            logger.info(f"{binding_type}: {count}")
        
        # Verify RNA-seq data format
        required_cols = ['gene', 'baseMean', 'log2FoldChange', 'padj']
        missing_cols = [col for col in required_cols if col not in self.rna_seq_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in RNA-seq data: {missing_cols}")
        
    def integrate_data(self, gene_annotations: pd.DataFrame) -> pd.DataFrame:
        """Map significantly enriched CpG regions to genes and integrate with RNA-seq data"""
        # First, find the nearest gene for each CpG region
        integrated_data = []
        
        for idx, region in self.significant_regions.iterrows():
            # Find nearest gene
            nearby_genes = gene_annotations[
                (gene_annotations['chr'] == region['chr']) &
                ((gene_annotations['start'] - 2000 <= region['end']) &
                 (gene_annotations['end'] + 2000 >= region['start']))
            ]
            
            if nearby_genes.empty:
                continue
            
            # Get nearest gene only
            distances = nearby_genes.apply(lambda x: min(
                abs(region['start'] - x['start']),
                abs(region['end'] - x['end'])
            ), axis=1)
            nearest_gene = nearby_genes.iloc[distances.argmin()]
            
            gene_data = {
                'gene': nearest_gene['gene_name'],
                'chr': region['chr'],
                'cpg_start': region['start'],
                'cpg_end': region['end'],
                'mecp2_enrichment': region['enrichment'],
                'exo_signal': region['exo_signal'],
                'endo_signal': region['endo_signal'],
                'binding_type': region['binding_type'],
                'distance_to_gene': distances.min()
            }
            
            # Add RNA-seq data if available
            rna_data = self.rna_seq_df[self.rna_seq_df['gene'] == nearest_gene['gene_name']]
            if not rna_data.empty:
                expr = rna_data.iloc[0]
                gene_data.update({
                    'baseMean': expr['baseMean'],
                    'log2FoldChange': expr['log2FoldChange'],
                    'padj': expr['padj']
                })
            
            integrated_data.append(gene_data)
        
        result_df = pd.DataFrame(integrated_data)
        
        # Remove duplicate genes, keeping the closest CpG region
        result_df = result_df.sort_values('distance_to_gene').groupby('gene').first().reset_index()
        
        return result_df
    
    def categorize_genes(self, integrated_results: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Categorize genes based on differential expression status"""
        categories = {}
        
        # Ensure we have expression data
        mask = ~integrated_results['log2FoldChange'].isna()
        expr_data = integrated_results[mask].copy()
        
        # Define categories based on fold change and significance
        categories['up-regulated'] = expr_data[
            (expr_data['log2FoldChange'] > 0.5) &  # log2FC > 0.5
            (expr_data['padj'] < 0.05)             # significant
        ]
        
        categories['down-regulated'] = expr_data[
            (expr_data['log2FoldChange'] < -0.5) & # log2FC < -0.5
            (expr_data['padj'] < 0.05)             # significant
        ]
        
        categories['non-deregulated'] = expr_data[
            ((expr_data['log2FoldChange'].abs() <= 0.5) |  # small fold change
             (expr_data['padj'] >= 0.05))                  # not significant
        ]
        
        # Log category sizes
        for category, df in categories.items():
            logger.info(f"{category}: {len(df)} genes")
            logger.info(f"Mean enrichment: {df['mecp2_enrichment'].mean():.2f}")
        
        return categories
    
    def plot_enrichment_by_category(self, integrated_results: pd.DataFrame, output_dir: str):
        """Create visualization of enrichment distribution by gene category"""
        # Categorize genes
        category_results = self.categorize_genes(integrated_results)
        
        # Create enrichment distribution plot by category
        plt.figure(figsize=(12, 7))  # Increased figure size for better readability
        
        # Define colors and labels for each category
        category_specs = {
            'up-regulated': {
                'color': 'orange',
                'label': 'up-regulated'
            },
            'down-regulated': {
                'color': 'green',
                'label': 'down-regulated'
            },
            'non-deregulated': {
                'color': 'blue',
                'label': 'non-deregulated'
            }
        }
        
        # Plot density for each category
        for category, df in category_results.items():
            if not df.empty:
                # Calculate the 95th percentile for x-axis limit
                upper_limit = df['mecp2_enrichment'].quantile(0.95)
                
                # Plot the KDE with cleaned data
                sns.kdeplot(
                    data=df[df['mecp2_enrichment'] <= upper_limit]['mecp2_enrichment'],
                    label=category_specs[category]['label'],
                    color=category_specs[category]['color'],
                    linewidth=2  # Thicker lines for better visibility
                )
        
        plt.xlabel('MeCP2 Enrichment (Exo/Endo)', fontsize=12)
        plt.ylabel('Density', fontsize=12)
        plt.title('MeCP2 Enrichment Distribution by Gene Category', fontsize=14, pad=20)
        
        # Add legend with category sizes
        legend_labels = [
            f"{cat} (n={len(df)})" 
            for cat, df in category_results.items()
        ]
        plt.legend(legend_labels, fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        
        # Save plot
        plt.savefig(
            os.path.join(output_dir, 'enrichment_distribution_by_category.pdf'),
            bbox_inches='tight',
            dpi=300
        )
        plt.close()
    
    def plot_integration_summary(self, integrated_results: pd.DataFrame, output_dir: str):
        """Create summary plots of the integrated analysis"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Original plots
        # Plot 1: MeCP2 Enrichment Distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=integrated_results, x='mecp2_enrichment', bins=50)
        plt.title('Distribution of MeCP2 Enrichment at CpG Islands')
        plt.xlabel('MeCP2 Enrichment (Exo/Endo)')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'mecp2_enrichment_distribution.pdf'))
        plt.close()
        
        # Plot 2: MeCP2 Enrichment vs Expression Change
        plt.figure(figsize=(10, 6))
        mask = ~integrated_results['log2FoldChange'].isna()
        sns.scatterplot(
            data=integrated_results[mask],
            x='log2FoldChange',
            y='mecp2_enrichment',
            alpha=0.6
        )
        plt.title('MeCP2 Enrichment vs Gene Expression Change')
        plt.xlabel('log2 Fold Change (RNA-seq)')
        plt.ylabel('MeCP2 Enrichment (Exo/Endo)')
        plt.savefig(os.path.join(output_dir, 'mecp2_vs_expression.pdf'))
        plt.close()
        
        # New plot: Enrichment by category
        # self.plot_enrichment_by_category(integrated_results, output_dir)

def load_gtf(gtf_file: str) -> pd.DataFrame:
    """Load and parse GTF file, extracting only gene records"""
    logger.info(f"Loading GTF file: {gtf_file}")
    
    # Read GTF file with correct format
    df = pd.read_csv(gtf_file, sep='\t', comment='#',
                     names=['chr', 'source', 'feature', 'start', 'end', 
                           'score', 'strand', 'frame', 'attributes'])
    
    # Filter for gene entries only
    genes = df[df['feature'] == 'gene'].copy()
    logger.info(f"Found {len(genes)} gene records")
    
    # Extract gene_name from attributes
    def extract_gene_name(attr_str):
        for attr in attr_str.split('; '):
            if attr.startswith('gene_name'):
                return attr.split('"')[1]
        return None
    
    genes['gene_name'] = genes['attributes'].apply(extract_gene_name)
    
    # Select required columns
    result = genes[['chr', 'start', 'end', 'gene_name']].copy()
    logger.info(f"Processed gene annotations shape: {result.shape}")
    logger.info("Sample of processed annotations:")
    logger.info(result.head().to_string())
    
    return result

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Integrate MeCP2 enrichment with RNA-seq data')
    parser.add_argument('--enrichment-file', type=str, required=True,
                       help='Path to MeCP2 enrichment results')
    parser.add_argument('--rna-seq-file', type=str, required=True,
                       help='Path to RNA-seq results')
    parser.add_argument('--gene-annotations', type=str, required=True,
                       help='Path to gene annotations file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    args = parser.parse_args()
    
    try:
        # Initialize integrator
        integrator = RNASeqIntegrator(args.enrichment_file, args.rna_seq_file)
        
        # Load gene annotations using custom GTF parser
        logger.info("Loading gene annotations...")
        gene_annotations = load_gtf(args.gene_annotations)
        
        # Integrate data
        logger.info("Integrating MeCP2 enrichment with gene annotations and RNA-seq data...")
        integrated_results = integrator.integrate_data(gene_annotations)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Save full results
        integrated_results.to_csv(
            os.path.join(args.output_dir, 'mecp2_enriched_genes.csv'),
            index=False
        )
        
        # Save MeCP2-enriched genes only (without RNA-seq data)
        mecp2_enriched = integrated_results[['gene', 'chr', 'cpg_start', 'cpg_end', 
                                           'mecp2_enrichment', 'exo_signal', 'endo_signal']]
        mecp2_enriched = mecp2_enriched.drop_duplicates()
        mecp2_enriched.to_csv(
            os.path.join(args.output_dir, 'mecp2_enriched_genes_only.csv'),
            index=False
        )
        
        # Create visualizations
        integrator.plot_integration_summary(integrated_results, args.output_dir)
        
        # Generate summary statistics
        logger.info("\nAnalysis Summary:")
        logger.info(f"Total enriched CpG regions: {len(integrator.significant_regions)}")
        logger.info(f"Total genes near enriched regions: {len(mecp2_enriched)}")
        logger.info(f"Genes with RNA-seq data: {integrated_results['log2FoldChange'].notna().sum()}")
        
    except Exception as e:
        logger.error(f"Error during integration: {str(e)}")
        raise

if __name__ == "__main__":
    main() 