#!/usr/bin/env python3
"""
Integrate Gene Body MeCP2 Binding with RNA-seq Data

This script integrates the enhanced gene body MeCP2 binding analysis results
with RNA-seq differential expression data to understand the relationship between
MeCP2 binding patterns in gene bodies and gene expression changes.

Features:
- Integration of position-specific MeCP2 binding with gene expression
- Analysis of binding vs expression relationships
- Identification of genes with coordinated binding and expression changes
- Statistical analysis of binding-expression correlations
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneBodyRNASeqIntegrator:
    """Integrates gene body MeCP2 binding results with RNA-seq data"""
    
    def __init__(self, output_dir: str, cell_type: str):
        self.output_dir = Path(output_dir)
        self.cell_type = cell_type
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def load_gene_body_results(self, gene_body_file: str) -> pd.DataFrame:
        """Load gene body analysis results"""
        
        logger.info(f"Loading gene body results from: {gene_body_file}")
        gene_body_df = pd.read_csv(gene_body_file)
        
        logger.info(f"Loaded gene body results for {len(gene_body_df)} genes")
        
        return gene_body_df
    
    def load_rnaseq_data(self, rnaseq_file: str, gene_id_col: str = 'gene_id', 
                        gene_name_col: str = 'gene_name') -> pd.DataFrame:
        """
        Load RNA-seq differential expression data
        
        Expected columns:
        - gene_id or gene_name: Gene identifier
        - log2FoldChange: Log2 fold change
        - padj: Adjusted p-value
        - baseMean: Base mean expression
        """
        
        logger.info(f"Loading RNA-seq data from: {rnaseq_file}")
        rnaseq_df = pd.read_csv(rnaseq_file)
        
        # Standardize column names
        if 'gene' in rnaseq_df.columns:
            rnaseq_df['gene_identifier'] = rnaseq_df['gene']
        elif gene_id_col in rnaseq_df.columns:
            rnaseq_df['gene_identifier'] = rnaseq_df[gene_id_col]
        elif gene_name_col in rnaseq_df.columns:
            rnaseq_df['gene_identifier'] = rnaseq_df[gene_name_col]
        else:
            raise ValueError(f"Neither 'gene', {gene_id_col} nor {gene_name_col} found in RNA-seq data")
        
        # Clean gene identifiers (remove version numbers if present)
        rnaseq_df['gene_identifier'] = rnaseq_df['gene_identifier'].str.replace(r'\.\d+$', '', regex=True)
        
        # Ensure required columns exist
        required_cols = ['log2FoldChange', 'padj']
        missing_cols = [col for col in required_cols if col not in rnaseq_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in RNA-seq data: {missing_cols}")
        
        # Fill missing values
        rnaseq_df['padj'] = rnaseq_df['padj'].fillna(1.0)
        rnaseq_df['log2FoldChange'] = rnaseq_df['log2FoldChange'].fillna(0.0)
        
        if 'baseMean' not in rnaseq_df.columns:
            rnaseq_df['baseMean'] = 1.0  # Default value
        
        if gene_id_col in rnaseq_df.columns:
            rnaseq_df['gene_id_clean_source'] = rnaseq_df[gene_id_col].astype(str).str.replace(r'\.\d+$', '', regex=True)
        else:
            rnaseq_df['gene_id_clean_source'] = rnaseq_df['gene_identifier']
        
        if gene_name_col in rnaseq_df.columns:
            rnaseq_df['gene_name_clean'] = rnaseq_df[gene_name_col].astype(str).str.replace(r'\.\d+$', '', regex=True)
        elif 'gene' in rnaseq_df.columns:
            rnaseq_df['gene_name_clean'] = rnaseq_df['gene'].astype(str).str.replace(r'\.\d+$', '', regex=True)
        else:
            rnaseq_df['gene_name_clean'] = rnaseq_df['gene_identifier']
        
        logger.info(f"Loaded RNA-seq data for {len(rnaseq_df)} genes")
        
        return rnaseq_df
    
    def integrate_datasets(self, 
                         gene_body_df: pd.DataFrame, 
                         rnaseq_df: pd.DataFrame,
                         match_by: str = 'gene_id') -> pd.DataFrame:
        """
        Integrate gene body and RNA-seq datasets
        
        Args:
            match_by: Column to use for matching ('gene_id' or 'gene_name')
        """
        
        # Prepare gene identifiers for matching
        if match_by == 'gene_id':
            gene_body_df = gene_body_df.copy()
            gene_body_df['gene_id_clean'] = gene_body_df['gene_id'].astype(str).str.replace(r'\.\d+$', '', regex=True)
            gene_body_df['gene_name_clean'] = gene_body_df['gene_name'].astype(str).str.replace(r'\.\d+$', '', regex=True)
            merge_col_left = 'gene_id_clean'
        elif match_by == 'gene_name':  # match by gene name
            gene_body_df = gene_body_df.copy()
            gene_body_df['gene_name_clean'] = gene_body_df['gene_name'].astype(str).str.replace(r'\.\d+$', '', regex=True)
            merge_col_left = 'gene_name_clean'
        else:
            raise ValueError(f"Invalid match_by argument: {match_by}. Must be 'gene_id' or 'gene_name'.")
        
        # Merge datasets
        logger.info(f"Integrating datasets by {match_by}...")
        integrated_df = pd.merge(
            gene_body_df, rnaseq_df, 
            left_on=merge_col_left, right_on='gene_identifier',
            how='inner'
        )
        
        if match_by == 'gene_id' and 'gene_name_clean' in gene_body_df.columns and 'gene_name_clean' in rnaseq_df.columns:
            matched_ids = set(integrated_df.get('gene_id_clean', []))
            remaining = gene_body_df[~gene_body_df['gene_id_clean'].isin(matched_ids)]
            if not remaining.empty:
                fallback = pd.merge(
                    remaining,
                    rnaseq_df,
                    left_on='gene_name_clean',
                    right_on='gene_name_clean',
                    how='inner'
                )
                if not fallback.empty:
                    logging.info("Added %d fallback matches using gene names", len(fallback))
                    integrated_df = pd.concat([integrated_df, fallback], ignore_index=True)
        
        integrated_df = integrated_df.drop_duplicates(subset=['gene_id', 'gene_name', 'gene_identifier'])
        integrated_df = integrated_df.drop(columns=[col for col in ['gene_id_clean', 'gene_name_clean', 'gene_id_clean_source'] if col in integrated_df.columns], errors='ignore')
        logger.info(f"Successfully integrated {len(integrated_df)} genes")
        
        # Categorize genes by expression change
        integrated_df['expression_category'] = 'unchanged'
        integrated_df.loc[
            (integrated_df['log2FoldChange'] > 0.5) & (integrated_df['padj'] < 0.05),
            'expression_category'
        ] = 'upregulated'
        integrated_df.loc[
            (integrated_df['log2FoldChange'] < -0.5) & (integrated_df['padj'] < 0.05),
            'expression_category'
        ] = 'downregulated'
        
        # Add expression significance
        integrated_df['expression_significant'] = integrated_df['padj'] < 0.05
        
        logger.info("Expression categorization:")
        logger.info(integrated_df['expression_category'].value_counts())
        
        return integrated_df
    
    def analyze_binding_expression_correlations(self, integrated_df: pd.DataFrame) -> Dict:
        """Analyze correlations between MeCP2 binding and gene expression changes"""
        
        logger.info("Analyzing binding-expression correlations...")
        
        results = {}
        
        # Define regions to analyze
        regions = ['promoter', 'gene_body_5prime', 'gene_body_middle', 'gene_body_3prime', 'gene_body_full']
        region_names = ['Promoter', '5\' Gene Body', 'Middle Gene Body', '3\' Gene Body', 'Full Gene Body']
        
        for region, region_name in zip(regions, region_names):
            exo_col = f'{region}_exo_signal'
            endo_col = f'{region}_endo_signal'
            enrich_col = f'{region}_enrichment'
            
            if all(col in integrated_df.columns for col in [exo_col, endo_col, enrich_col]):
                # Filter for genes with measurable signals
                valid_data = integrated_df[
                    (integrated_df[exo_col] > 0) & (integrated_df[endo_col] > 0)
                ].copy()
                
                if len(valid_data) > 10:
                    # Correlations with log2FoldChange
                    try:
                        # Exogenous signal vs expression
                        corr_exo, p_exo = stats.pearsonr(
                            valid_data[exo_col], valid_data['log2FoldChange']
                        )
                        
                        # Endogenous signal vs expression
                        corr_endo, p_endo = stats.pearsonr(
                            valid_data[endo_col], valid_data['log2FoldChange']
                        )
                        
                        # Enrichment vs expression
                        enrichments = valid_data[enrich_col].replace([np.inf, -np.inf], np.nan)
                        enrichments = enrichments[(enrichments > 0) & (enrichments < 100)]
                        
                        if len(enrichments) > 10:
                            corr_enrich, p_enrich = stats.pearsonr(
                                enrichments, valid_data.loc[enrichments.index, 'log2FoldChange']
                            )
                        else:
                            corr_enrich, p_enrich = np.nan, np.nan
                        
                        results[region] = {
                            'region_name': region_name,
                            'n_genes': len(valid_data),
                            'exo_signal_correlation': corr_exo,
                            'exo_signal_pvalue': p_exo,
                            'endo_signal_correlation': corr_endo,
                            'endo_signal_pvalue': p_endo,
                            'enrichment_correlation': corr_enrich,
                            'enrichment_pvalue': p_enrich
                        }
                        
                    except Exception as e:
                        logger.warning(f"Error analyzing correlations for {region}: {e}")
                        
        return results
    
    def analyze_expression_categories(self, integrated_df: pd.DataFrame) -> Dict:
        """Analyze MeCP2 binding patterns across expression categories"""
        
        logger.info("Analyzing binding patterns by expression category...")
        
        results = {}
        
        # Define regions to analyze
        regions = ['promoter', 'gene_body_5prime', 'gene_body_middle', 'gene_body_3prime', 'gene_body_full']
        region_names = ['Promoter', '5\' Gene Body', 'Middle Gene Body', '3\' Gene Body', 'Full Gene Body']
        
        expression_categories = ['upregulated', 'downregulated', 'unchanged']
        
        for region, region_name in zip(regions, region_names):
            exo_col = f'{region}_exo_signal'
            endo_col = f'{region}_endo_signal'
            enrich_col = f'{region}_enrichment'
            
            if all(col in integrated_df.columns for col in [exo_col, endo_col]):
                region_results = {'region_name': region_name}
                
                for category in expression_categories:
                    category_data = integrated_df[
                        (integrated_df['expression_category'] == category) &
                        (integrated_df[exo_col] > 0) & (integrated_df[endo_col] > 0)
                    ]
                    
                    if len(category_data) > 5:
                        region_results[f'{category}_n_genes'] = len(category_data)
                        region_results[f'{category}_exo_median'] = category_data[exo_col].median()
                        region_results[f'{category}_endo_median'] = category_data[endo_col].median()
                        
                        if enrich_col in category_data.columns:
                            enrichments = category_data[enrich_col].replace([np.inf, -np.inf], np.nan)
                            enrichments = enrichments[(enrichments > 0) & (enrichments < 100)]
                            region_results[f'{category}_enrichment_median'] = enrichments.median()
                    else:
                        region_results[f'{category}_n_genes'] = len(category_data)
                        region_results[f'{category}_exo_median'] = np.nan
                        region_results[f'{category}_endo_median'] = np.nan
                        region_results[f'{category}_enrichment_median'] = np.nan
                
                results[region] = region_results
        
        return results
    
    def create_integration_visualizations(self, integrated_df: pd.DataFrame, correlation_results: Dict):
        """Create comprehensive visualizations for the integrated analysis"""
        
        # 1. Correlation matrix heatmap
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        correlation_data = []
        regions = []
        
        for region, results in correlation_results.items():
            if 'region_name' in results:
                regions.append(results['region_name'])
                correlation_data.append([
                    results.get('exo_signal_correlation', 0),
                    results.get('endo_signal_correlation', 0),
                    results.get('enrichment_correlation', 0)
                ])
        
        if correlation_data:
            corr_df = pd.DataFrame(
                correlation_data,
                index=regions,
                columns=['Exogenous Signal', 'Endogenous Signal', 'Enrichment']
            )
            
            sns.heatmap(corr_df, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                       vmin=-0.5, vmax=0.5, ax=ax)
            ax.set_title(f'MeCP2 Binding vs Gene Expression Correlations - {self.cell_type}')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f'binding_expression_correlations_{self.cell_type}.png',
                   dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / f'binding_expression_correlations_{self.cell_type}.pdf',
                   bbox_inches='tight')
        plt.close()
        
        # 2. Scatter plots of binding vs expression
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'MeCP2 Binding vs Gene Expression - {self.cell_type}', fontsize=16)
        
        regions = ['promoter', 'gene_body_5prime', 'gene_body_middle', 'gene_body_3prime', 'gene_body_full']
        region_labels = ['Promoter', '5\' Gene Body', 'Middle Gene Body', '3\' Gene Body', 'Full Gene Body']
        
        for i, (region, label) in enumerate(zip(regions, region_labels)):
            if i >= 6:
                break
                
            row, col = i // 3, i % 3
            ax = axes[row, col]
            
            enrich_col = f'{region}_enrichment'
            if enrich_col in integrated_df.columns:
                # Filter data for visualization
                plot_data = integrated_df[
                    (integrated_df[enrich_col] > 0) & 
                    (integrated_df[enrich_col] < 10) &
                    (integrated_df['log2FoldChange'].abs() < 5)
                ].copy()
                
                if len(plot_data) > 10:
                    # Color by expression category
                    colors = {'upregulated': 'red', 'downregulated': 'blue', 'unchanged': 'gray'}
                    
                    for category in ['upregulated', 'downregulated', 'unchanged']:
                        cat_data = plot_data[plot_data['expression_category'] == category]
                        if len(cat_data) > 0:
                            ax.scatter(cat_data[enrich_col], cat_data['log2FoldChange'],
                                     alpha=0.6, color=colors[category], label=category, s=20)
                    
                    ax.set_xlabel('MeCP2 Enrichment (Exo/Endo)')
                    ax.set_ylabel('Log2 Fold Change (Expression)')
                    ax.set_title(label)
                    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                    ax.axvline(x=1, color='black', linestyle='-', alpha=0.3)
                    
                    if i == 0:  # Add legend to first subplot
                        ax.legend()
        
        # Remove empty subplot
        if len(regions) < 6:
            axes[1, 2].remove()
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f'binding_vs_expression_scatter_{self.cell_type}.png',
                   dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / f'binding_vs_expression_scatter_{self.cell_type}.pdf',
                   bbox_inches='tight')
        plt.close()
        
        # 3. Expression category comparison
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'MeCP2 Binding by Expression Category - {self.cell_type}', fontsize=14)
        
        regions = ['promoter', 'gene_body_full']
        conditions = ['exo', 'endo']
        
        for i, (region, condition) in enumerate([(r, c) for r in regions for c in conditions]):
            ax = axes[i // 2, i % 2]
            
            signal_col = f'{region}_{condition}_signal'
            if signal_col in integrated_df.columns:
                # Prepare data for boxplot
                plot_data = []
                for category in ['upregulated', 'downregulated', 'unchanged']:
                    cat_data = integrated_df[
                        (integrated_df['expression_category'] == category) &
                        (integrated_df[signal_col] > 0)
                    ]
                    
                    for val in cat_data[signal_col]:
                        plot_data.append({
                            'Signal': val,
                            'Expression_Category': category.title()
                        })
                
                if plot_data:
                    plot_df = pd.DataFrame(plot_data)
                    sns.boxplot(data=plot_df, x='Expression_Category', y='Signal', ax=ax)
                    ax.set_title(f'{region.replace("_", " ").title()} - {condition.title()}')
                    ax.set_ylabel('MeCP2 Signal')
                    plt.setp(ax.get_xticklabels(), rotation=45)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f'binding_by_expression_category_{self.cell_type}.png',
                   dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / f'binding_by_expression_category_{self.cell_type}.pdf',
                   bbox_inches='tight')
        plt.close()
        
    def identify_coordinated_genes(self, integrated_df: pd.DataFrame, 
                                 enrichment_threshold: float = 1.5,
                                 expression_threshold: float = 0.5,
                                 pvalue_threshold: float = 0.05) -> Dict[str, pd.DataFrame]:
        """Identify genes with coordinated binding and expression changes"""
        
        logger.info("Identifying genes with coordinated binding and expression changes...")
        
        # Define coordination categories
        coordinated_genes = {}
        
        # High binding + upregulated
        high_binding_upregulated = integrated_df[
            (integrated_df['gene_body_full_enrichment'] > enrichment_threshold) &
            (integrated_df['log2FoldChange'] > expression_threshold) &
            (integrated_df['padj'] < pvalue_threshold)
        ].copy()
        
        # Low binding + downregulated
        low_binding_downregulated = integrated_df[
            (integrated_df['gene_body_full_enrichment'] < (1/enrichment_threshold)) &
            (integrated_df['log2FoldChange'] < -expression_threshold) &
            (integrated_df['padj'] < pvalue_threshold)
        ].copy()
        
        # High binding + downregulated (unexpected)
        high_binding_downregulated = integrated_df[
            (integrated_df['gene_body_full_enrichment'] > enrichment_threshold) &
            (integrated_df['log2FoldChange'] < -expression_threshold) &
            (integrated_df['padj'] < pvalue_threshold)
        ].copy()
        
        # Low binding + upregulated (unexpected)
        low_binding_upregulated = integrated_df[
            (integrated_df['gene_body_full_enrichment'] < (1/enrichment_threshold)) &
            (integrated_df['log2FoldChange'] > expression_threshold) &
            (integrated_df['padj'] < pvalue_threshold)
        ].copy()
        
        coordinated_genes['high_binding_upregulated'] = high_binding_upregulated
        coordinated_genes['low_binding_downregulated'] = low_binding_downregulated
        coordinated_genes['high_binding_downregulated'] = high_binding_downregulated
        coordinated_genes['low_binding_upregulated'] = low_binding_upregulated
        
        # Log results
        logger.info("Coordinated gene categories:")
        for category, genes_df in coordinated_genes.items():
            logger.info(f"  {category}: {len(genes_df)} genes")
        
        return coordinated_genes
    
    def save_integration_results(self, 
                               integrated_df: pd.DataFrame,
                               correlation_results: Dict,
                               category_results: Dict,
                               coordinated_genes: Dict[str, pd.DataFrame]):
        """Save all integration results"""
        
        # Save integrated dataset
        integrated_df.to_csv(
            self.output_dir / f'integrated_gene_body_rnaseq_{self.cell_type}.csv',
            index=False
        )
        
        # Save correlation results
        if correlation_results:
            corr_df = pd.DataFrame(correlation_results).T
            corr_df.to_csv(
                self.output_dir / f'binding_expression_correlations_{self.cell_type}.csv'
            )
        
        # Save category results
        if category_results:
            cat_df = pd.DataFrame(category_results).T
            cat_df.to_csv(
                self.output_dir / f'binding_by_expression_category_{self.cell_type}.csv'
            )
        
        # Save coordinated genes
        for category, genes_df in coordinated_genes.items():
            if not genes_df.empty:
                genes_df.to_csv(
                    self.output_dir / f'coordinated_genes_{category}_{self.cell_type}.csv',
                    index=False
                )
        
        # Create summary report
        with open(self.output_dir / f'integration_summary_{self.cell_type}.txt', 'w') as f:
            f.write(f"Gene Body MeCP2 Binding - RNA-seq Integration Summary\n")
            f.write(f"Cell Type: {self.cell_type}\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Dataset Overview:\n")
            f.write(f"- Total integrated genes: {len(integrated_df)}\n")
            f.write(f"- Upregulated genes: {(integrated_df['expression_category'] == 'upregulated').sum()}\n")
            f.write(f"- Downregulated genes: {(integrated_df['expression_category'] == 'downregulated').sum()}\n")
            f.write(f"- Unchanged genes: {(integrated_df['expression_category'] == 'unchanged').sum()}\n\n")
            
            if correlation_results:
                f.write("Binding-Expression Correlations:\n")
                f.write("-" * 35 + "\n")
                for region, results in correlation_results.items():
                    if 'region_name' in results:
                        f.write(f"{results['region_name']}:\n")
                        f.write(f"  - Exogenous signal correlation: {results.get('exo_signal_correlation', 'N/A'):.3f}\n")
                        f.write(f"  - Endogenous signal correlation: {results.get('endo_signal_correlation', 'N/A'):.3f}\n")
                        f.write(f"  - Enrichment correlation: {results.get('enrichment_correlation', 'N/A'):.3f}\n\n")
            
            f.write("Coordinated Gene Categories:\n")
            f.write("-" * 30 + "\n")
            for category, genes_df in coordinated_genes.items():
                f.write(f"  {category}: {len(genes_df)} genes\n")
        
        logger.info(f"Integration results saved to {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Integrate gene body MeCP2 binding with RNA-seq data'
    )
    parser.add_argument('--gene-body-results', required=True, 
                       help='CSV file with combined gene body analysis results')
    parser.add_argument('--rnaseq-data', required=True,
                       help='CSV file with RNA-seq differential expression results')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for integration results')
    parser.add_argument('--cell-type', required=True, choices=['NSC', 'Neu'],
                       help='Cell type being analyzed')
    parser.add_argument('--match-by', choices=['gene_id', 'gene_name'], default='gene_id',
                       help='Column to use for matching datasets')
    parser.add_argument('--gene-id-col', default='gene_id',
                       help='Gene ID column name in RNA-seq data')
    parser.add_argument('--gene-name-col', default='gene_name',
                       help='Gene name column name in RNA-seq data')
    
    args = parser.parse_args()
    
    # Initialize integrator
    integrator = GeneBodyRNASeqIntegrator(args.output_dir, args.cell_type)
    
    # Load datasets
    logger.info("Loading datasets...")
    gene_body_df = integrator.load_gene_body_results(args.gene_body_results)
    rnaseq_df = integrator.load_rnaseq_data(args.rnaseq_data, args.gene_id_col, args.gene_name_col)
    
    # Integrate datasets
    logger.info("Integrating datasets...")
    integrated_df = integrator.integrate_datasets(gene_body_df, rnaseq_df, args.match_by)
    
    # Perform analyses
    logger.info("Analyzing correlations...")
    correlation_results = integrator.analyze_binding_expression_correlations(integrated_df)
    
    logger.info("Analyzing expression categories...")
    category_results = integrator.analyze_expression_categories(integrated_df)
    
    logger.info("Identifying coordinated genes...")
    coordinated_genes = integrator.identify_coordinated_genes(integrated_df)
    
    # Create visualizations
    logger.info("Creating visualizations...")
    integrator.create_integration_visualizations(integrated_df, correlation_results)
    
    # Save results
    logger.info("Saving results...")
    integrator.save_integration_results(
        integrated_df, correlation_results, category_results, coordinated_genes
    )
    
    logger.info("Integration analysis complete!")

if __name__ == "__main__":
    main()