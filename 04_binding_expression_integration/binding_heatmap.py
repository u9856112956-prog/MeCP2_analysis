#!/usr/bin/env python3
"""
MeCP2 Binding vs Expression Heatmap Visualization

This script creates comprehensive heatmap visualizations showing the relationship
between MeCP2 binding differences and gene expression changes across different
genomic regions and gene expression categories.

Creates several types of heatmaps:
1. Correlation matrix heatmap (regions vs expression categories)
2. Binding signal heatmap for top differentially expressed genes
3. Combined binding ratio and expression fold change heatmap
4. Statistical significance overlay heatmaps

Usage:
    python binding_expression_heatmap_visualization.py \
        --nsc-results /path/to/binding_expression_correlation_nsc.csv \
        --neu-results /path/to/binding_expression_correlation_neu.csv \
        --nsc-integrated /path/to/integrated_gene_body_rnaseq_NSC.csv \
        --neu-integrated /path/to/integrated_gene_body_rnaseq_Neu.csv \
        --output-dir /path/to/output
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr, spearmanr

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Constants
REGION_ORDER = ["promoter", "gene_body_5prime", "gene_body_middle", "gene_body_3prime", "gene_body_full"]
REGION_LABELS = {
    "promoter": "Promoter",
    "gene_body_5prime": "5' Gene Body",
    "gene_body_middle": "Middle Gene Body",
    "gene_body_3prime": "3' Gene Body",
    "gene_body_full": "Full Gene Body"
}

SUBSET_MAP = {
    'all_dea': 'all',
    'all': 'all',
    'padj_lt_0.05': 'padj_lt_0.05',
    'padj_gte_0.05': 'padj_gte_0.05',
    'upregulated': 'upregulated',
    'downregulated': 'downregulated'
}
DEFAULT_SUBSET_ORDER = ['all', 'upregulated', 'downregulated', 'padj_lt_0.05', 'padj_gte_0.05']

def standardize_correlation_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Harmonize column names and subset labels across correlation result files."""
    df = df.copy()

    if 'correlation_log2_ratio' not in df.columns:
        for candidate in ('correlation_log2_ratio_abs', 'correlation_log2_ratio_signed'):
            if candidate in df.columns:
                df['correlation_log2_ratio'] = df[candidate]
                break

    if 'pvalue_log2_ratio' not in df.columns:
        for candidate in ('pvalue_log2_ratio_abs', 'pvalue_log2_ratio_signed'):
            if candidate in df.columns:
                df['pvalue_log2_ratio'] = df[candidate]
                break

    if 'subset' in df.columns:
        df['subset'] = df['subset'].map(lambda x: SUBSET_MAP.get(x, x))

    return df

def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

def load_correlation_results(nsc_path: Path, neu_path: Optional[Path] = None) -> pd.DataFrame:
    """Load and combine correlation results from NSC and optionally NEU."""
    logging.info("Loading correlation results")

    nsc_df = standardize_correlation_dataframe(pd.read_csv(nsc_path))
    nsc_df['cell_type'] = 'NSC'

    combined_frames = [nsc_df]

    if neu_path and neu_path.exists():
        neu_df = standardize_correlation_dataframe(pd.read_csv(neu_path))
        neu_df['cell_type'] = 'Neu'
        combined_frames.append(neu_df)
    else:
        logging.warning("NEU results not found, using NSC only")

    combined_df = pd.concat(combined_frames, ignore_index=True)
    logging.info("Loaded correlation results for %d combinations", len(combined_df))
    return combined_df

def load_integrated_data(nsc_path: Path, neu_path: Optional[Path] = None) -> pd.DataFrame:
    """Load and combine integrated gene body + RNA-seq data."""
    logging.info("Loading integrated datasets")

    nsc_df = pd.read_csv(nsc_path)
    nsc_df['cell_type'] = 'NSC'

    if neu_path and neu_path.exists():
        neu_df = pd.read_csv(neu_path)
        neu_df['cell_type'] = 'Neu'
        combined_df = pd.concat([nsc_df, neu_df], ignore_index=True)
    else:
        logging.warning("NEU integrated data not found, using NSC only")
        combined_df = nsc_df

    logging.info("Loaded integrated data for %d genes across cell types", len(combined_df))
    return combined_df

def create_correlation_matrix_heatmap(results_df: pd.DataFrame, output_dir: Path) -> None:
    """Create correlation coefficient heatmap across regions and subsets."""
    logging.info("Creating correlation matrix heatmap")

    available_subsets = [s for s in DEFAULT_SUBSET_ORDER if s in results_df['subset'].unique()]
    if not available_subsets:
        available_subsets = sorted(results_df['subset'].unique())

    pivot_data = results_df.pivot_table(
        index=['cell_type', 'region'],
        columns='subset',
        values='correlation_log2_ratio',
        aggfunc='first'
    ).reindex(columns=available_subsets)

    pval_pivot = results_df.pivot_table(
        index=['cell_type', 'region'],
        columns='subset',
        values='pvalue_log2_ratio',
        aggfunc='first'
    ).reindex(columns=available_subsets)


    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = sns.diverging_palette(250, 10, as_cmap=True, center='light')

    sns.heatmap(
        pivot_data,
        annot=True,
        fmt='.3f',
        cmap=cmap,
        center=0,
        cbar_kws={'label': 'Spearman Correlation (r)'},
        ax=ax,
        mask=pivot_data.isna(),
        linewidths=0.5
    )

    for (i, j), value in np.ndenumerate(pivot_data.values):
        if np.isnan(value):
            continue
        p_value = np.nan
        if pval_pivot is not None:
            try:
                p_value = pval_pivot.values[i, j]
            except IndexError:
                p_value = np.nan
        if not np.isnan(p_value) and p_value <= 0.05:
            ax.text(j + 0.7, i + 0.3, '*', fontsize=12, fontweight='bold', color='white')

    ax.set_title('MeCP2 Binding-Expression Correlations\n(* p < 0.05)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Gene Subset', fontsize=12)
    ax.set_ylabel('Cell Type & Genomic Region', fontsize=12)

    plt.tight_layout()
    output_path = output_dir / "correlation_matrix_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    logging.info("Saved correlation matrix heatmap to %s", output_path)
def create_top_genes_binding_heatmap(integrated_df: pd.DataFrame, output_dir: Path,
                                   n_genes: int = 50) -> None:
    """Create heatmap showing binding signals for top DE genes."""
    logging.info("Creating top genes binding heatmap")

    # Get top upregulated and downregulated genes
    top_up = integrated_df.nlargest(n_genes//2, 'log2FoldChange')
    top_down = integrated_df.nsmallest(n_genes//2, 'log2FoldChange')
    top_genes = pd.concat([top_up, top_down])

    # Prepare binding data matrix
    binding_cols = [f"{region}_{condition}_signal"
                   for region in REGION_ORDER
                   for condition in ['exo', 'endo']]

    # Filter for available columns
    available_cols = [col for col in binding_cols if col in top_genes.columns]
    binding_data = top_genes[available_cols + ['gene_name', 'log2FoldChange', 'cell_type']].copy()

    # Log transform binding signals
    for col in available_cols:
        binding_data[col] = np.log2(binding_data[col] + 1)

    # Create hierarchical clustering
    binding_matrix = binding_data[available_cols].fillna(0)
    linkage_matrix = linkage(binding_matrix, method='ward')

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10),
                                  gridspec_kw={'width_ratios': [4, 1]})

    # Plot dendrogram
    dendro = dendrogram(linkage_matrix, ax=ax2, orientation='right',
                       labels=binding_data['gene_name'].values, leaf_font_size=6)
    ax2.set_title('Gene Clustering')

    # Reorder data according to clustering
    gene_order = dendro['leaves']
    binding_matrix_ordered = binding_matrix.iloc[gene_order]

    # Plot heatmap
    im = ax1.imshow(binding_matrix_ordered.T, aspect='auto', cmap='viridis',
                   interpolation='nearest')

    # Customize heatmap
    ax1.set_xlabel('Genes (ordered by clustering)')
    ax1.set_ylabel('Genomic Region & Condition')
    ax1.set_title(f'MeCP2 Binding Signals for Top {n_genes} DE Genes\n(log2 transformed)')

    # Set y-tick labels
    region_labels = []
    for col in available_cols:
        parts = col.replace('_signal', '').split('_')
        region = '_'.join(parts[:-1])
        condition = parts[-1]
        region_labels.append(f"{REGION_LABELS.get(region, region)} ({condition.title()})")

    ax1.set_yticks(range(len(available_cols)))
    ax1.set_yticklabels(region_labels, fontsize=8)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('log2(Signal + 1)')

    plt.tight_layout()

    # Save
    output_path = output_dir / f"top_{n_genes}_genes_binding_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    logging.info("Saved top genes binding heatmap to %s", output_path)

def create_binding_ratio_expression_heatmap(integrated_df: pd.DataFrame, output_dir: Path,
                                          n_genes: int = 100) -> None:
    """Create dual heatmap showing binding ratios and expression changes."""
    logging.info("Creating binding ratio vs expression heatmap")

    # Select genes with strongest expression changes
    top_genes = integrated_df.nlargest(n_genes//2, 'log2FoldChange').copy()
    bottom_genes = integrated_df.nsmallest(n_genes//2, 'log2FoldChange').copy()
    selected_genes = pd.concat([top_genes, bottom_genes])

    # Calculate binding ratios for each region
    ratio_data = []
    for region in REGION_ORDER:
        exo_col = f"{region}_exo_signal"
        endo_col = f"{region}_endo_signal"

        if exo_col in selected_genes.columns and endo_col in selected_genes.columns:
            ratio = np.log2((selected_genes[endo_col] + 1e-6) / (selected_genes[exo_col] + 1e-6))
            ratio_data.append(ratio.values)

    if not ratio_data:
        logging.warning("No binding ratio data available")
        return

    # Create data matrices
    binding_ratios = np.array(ratio_data).T  # genes x regions
    expression_changes = selected_genes['log2FoldChange'].values.reshape(-1, 1)

    # Sort by expression change
    sort_idx = np.argsort(expression_changes.flatten())
    binding_ratios_sorted = binding_ratios[sort_idx]
    expression_sorted = expression_changes[sort_idx]
    gene_names_sorted = selected_genes.iloc[sort_idx]['gene_name'].values

    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 10),
                                       gridspec_kw={'width_ratios': [4, 1, 0.5]})

    # Plot binding ratios heatmap
    norm1 = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    im1 = ax1.imshow(binding_ratios_sorted.T, aspect='auto',
                    cmap='RdBu_r', norm=norm1, interpolation='nearest')

    ax1.set_title('MeCP2 Binding Ratio (log2 Endo/Exo)')
    ax1.set_xlabel('Genes (sorted by expression change)')
    ax1.set_ylabel('Genomic Region')
    ax1.set_yticks(range(len(REGION_ORDER)))
    ax1.set_yticklabels([REGION_LABELS[r] for r in REGION_ORDER])

    # Plot expression changes
    norm2 = TwoSlopeNorm(vmin=expression_sorted.min(), vcenter=0, vmax=expression_sorted.max())
    im2 = ax2.imshow(expression_sorted.T, aspect='auto',
                    cmap='RdYlBu_r', norm=norm2, interpolation='nearest')

    ax2.set_title('Expression\nlog2FC')
    ax2.set_xlabel('Genes')
    ax2.set_yticks([])

    # Add colorbars
    cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
    cbar1.set_label('log2(MeCP2 Endo/Exo)')

    cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
    cbar2.set_label('log2(Expression Endo/Exo)')

    # Add gene annotations for extreme cases
    ax3.axis('off')
    ax3.text(0.1, 0.9, 'Top Upregulated:', fontweight='bold', transform=ax3.transAxes)
    for i, gene in enumerate(gene_names_sorted[-5:]):
        ax3.text(0.1, 0.85-i*0.04, gene, fontsize=8, transform=ax3.transAxes)

    ax3.text(0.1, 0.6, 'Top Downregulated:', fontweight='bold', transform=ax3.transAxes)
    for i, gene in enumerate(gene_names_sorted[:5]):
        ax3.text(0.1, 0.55-i*0.04, gene, fontsize=8, transform=ax3.transAxes)

    plt.tight_layout()

    # Save
    output_path = output_dir / "binding_ratio_expression_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    logging.info("Saved binding ratio vs expression heatmap to %s", output_path)

def create_cell_type_comparison_heatmap(results_df: pd.DataFrame, output_dir: Path) -> None:
    """Create heatmap comparing correlations between cell types."""
    if 'NSC' not in results_df['cell_type'].values or 'Neu' not in results_df['cell_type'].values:
        logging.warning("Both cell types not available for comparison")
        return

    logging.info("Creating cell type comparison heatmap")

    nsc_data = results_df[results_df['cell_type'] == 'NSC'].copy()
    neu_data = results_df[results_df['cell_type'] == 'Neu'].copy()

    subset_candidates = [s for s in DEFAULT_SUBSET_ORDER if s in results_df['subset'].unique()]
    if not subset_candidates:
        subset_candidates = sorted(results_df['subset'].unique())

    comparison_data = []
    for region in REGION_ORDER:
        for subset in subset_candidates:
            nsc_row = nsc_data[(nsc_data['region'] == region) & (nsc_data['subset'] == subset)]
            neu_row = neu_data[(neu_data['region'] == region) & (neu_data['subset'] == subset)]

            if nsc_row.empty or neu_row.empty:
                continue

            nsc_corr = nsc_row['correlation_log2_ratio'].iloc[0]
            neu_corr = neu_row['correlation_log2_ratio'].iloc[0]
            nsc_pval = nsc_row.get('pvalue_log2_ratio', pd.Series([1.0])).iloc[0]
            neu_pval = neu_row.get('pvalue_log2_ratio', pd.Series([1.0])).iloc[0]

            comparison_data.append({
                'region': region,
                'subset': subset,
                'nsc_correlation': nsc_corr,
                'neu_correlation': neu_corr,
                'nsc_pvalue': nsc_pval,
                'neu_pvalue': neu_pval,
                'difference': nsc_corr - neu_corr
            })

    comparison_df = pd.DataFrame(comparison_data)

    if comparison_df.empty:
        logging.warning("No comparable data between cell types")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    nsc_pivot = comparison_df.pivot(index='region', columns='subset', values='nsc_correlation').reindex(
        index=REGION_ORDER, columns=subset_candidates
    )
    neu_pivot = comparison_df.pivot(index='region', columns='subset', values='neu_correlation').reindex(
        index=REGION_ORDER, columns=subset_candidates
    )
    diff_pivot = comparison_df.pivot(index='region', columns='subset', values='difference').reindex(
        index=REGION_ORDER, columns=subset_candidates
    )

    sns.heatmap(nsc_pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                ax=axes[0,0], cbar_kws={'label': 'Correlation'})
    axes[0,0].set_title('NSC Correlations')

    sns.heatmap(neu_pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                ax=axes[0,1], cbar_kws={'label': 'Correlation'})
    axes[0,1].set_title('Neuron Correlations')

    sns.heatmap(diff_pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                ax=axes[1,0], cbar_kws={'label': 'Difference (NSC - Neu)'})
    axes[1,0].set_title('Correlation Differences')

    sig_matrix = np.zeros_like(nsc_pivot.fillna(0).values)
    for i, region in enumerate(nsc_pivot.index):
        for j, subset in enumerate(nsc_pivot.columns):
            nsc_sig = comparison_df[(comparison_df['region'] == region) &
                                   (comparison_df['subset'] == subset)]['nsc_pvalue'].values
            neu_sig = comparison_df[(comparison_df['region'] == region) &
                                   (comparison_df['subset'] == subset)]['neu_pvalue'].values
            if len(nsc_sig) == 0 or len(neu_sig) == 0:
                sig_matrix[i, j] = 0
                continue
            if nsc_sig[0] < 0.05 and neu_sig[0] < 0.05:
                sig_matrix[i, j] = 3
            elif nsc_sig[0] < 0.05:
                sig_matrix[i, j] = 1
            elif neu_sig[0] < 0.05:
                sig_matrix[i, j] = 2

    sig_df = pd.DataFrame(sig_matrix, index=nsc_pivot.index, columns=nsc_pivot.columns)
    colors = ['white', 'lightblue', 'lightcoral', 'purple']
    cmap_sig = LinearSegmentedColormap.from_list('sig', colors, N=4)

    sns.heatmap(sig_df, annot=False, cmap=cmap_sig,
                ax=axes[1,1], cbar_kws={'label': 'Significance'})
    axes[1,1].set_title('Significance Pattern\n(0=Neither, 1=NSC, 2=Neu, 3=Both)')

    plt.tight_layout()

    output_path = output_dir / "cell_type_comparison_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    logging.info("Saved cell type comparison heatmap to %s", output_path)
def create_comprehensive_summary_plot(results_df: pd.DataFrame, integrated_df: pd.DataFrame,
                                    output_dir: Path) -> None:
    """Create a comprehensive summary visualization."""
    logging.info("Creating comprehensive summary plot")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # 1. Correlation strength by region (all genes)
    all_genes_data = results_df[results_df['subset'] == 'all']
    if not all_genes_data.empty:
        pivot_corr = all_genes_data.pivot(index='region', columns='cell_type',
                                         values='correlation_log2_ratio')
        sns.heatmap(pivot_corr, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                   ax=axes[0,0], cbar_kws={'label': 'Correlation'})
        axes[0,0].set_title('Correlation Strength by Region')

    # 2. P-value significance
    if not all_genes_data.empty:
        pivot_pval = all_genes_data.pivot(index='region', columns='cell_type',
                                         values='pvalue_log2_ratio')
        # Convert to -log10 for better visualization
        pivot_pval_log = -np.log10(pivot_pval)
        sns.heatmap(pivot_pval_log, annot=True, fmt='.1f', cmap='Reds',
                   ax=axes[0,1], cbar_kws={'label': '-log10(p-value)'})
        axes[0,1].set_title('Statistical Significance')
        axes[0,1].axhline(y=0, color='black', linewidth=2)  # Significance threshold line

    # 3. Effect sizes (R-squared)
    if not all_genes_data.empty and 'r_squared' in all_genes_data.columns:
        pivot_r2 = all_genes_data.pivot(index='region', columns='cell_type',
                                       values='r_squared')
        sns.heatmap(pivot_r2, annot=True, fmt='.4f', cmap='Greens',
                   ax=axes[0,2], cbar_kws={'label': 'R²'})
        axes[0,2].set_title('Effect Size (R²)')

    # 4. Gene count distribution
    subset_data = results_df[results_df['subset'].isin(['upregulated', 'downregulated'])]
    if not subset_data.empty:
        subset_pivot = subset_data.pivot_table(index=['cell_type', 'subset'],
                                              columns='region', values='n_genes')
        sns.heatmap(subset_pivot, annot=True, fmt='d', cmap='Blues',
                   ax=axes[1,0], cbar_kws={'label': 'Number of genes'})
        axes[1,0].set_title('Gene Count by Direction')

    # 5. Binding signal distribution
    if not integrated_df.empty:
        # Calculate mean binding signals across regions
        binding_means = []
        for region in REGION_ORDER:
            exo_col = f"{region}_exo_signal"
            endo_col = f"{region}_endo_signal"
            if exo_col in integrated_df.columns and endo_col in integrated_df.columns:
                mean_exo = integrated_df.groupby('cell_type')[exo_col].mean()
                mean_endo = integrated_df.groupby('cell_type')[endo_col].mean()
                for cell_type in mean_exo.index:
                    binding_means.append({
                        'region': region,
                        'cell_type': cell_type,
                        'condition': 'Exo',
                        'signal': mean_exo[cell_type]
                    })
                    binding_means.append({
                        'region': region,
                        'cell_type': cell_type,
                        'condition': 'Endo',
                        'signal': mean_endo[cell_type]
                    })

        if binding_means:
            binding_df = pd.DataFrame(binding_means)
            binding_pivot = binding_df.pivot_table(
                index=['cell_type', 'condition'], columns='region', values='signal'
            )
            sns.heatmap(np.log2(binding_pivot + 1), annot=True, fmt='.1f', cmap='viridis',
                       ax=axes[1,1], cbar_kws={'label': 'log2(Signal + 1)'})
            axes[1,1].set_title('Mean Binding Signals')

    # 6. Summary statistics
    axes[1,2].axis('off')

    # Calculate summary statistics
    stats_text = "Summary Statistics:\n\n"

    promoter_data = results_df[(results_df['region'] == 'promoter') &
                              (results_df['subset'] == 'all')]
    if not promoter_data.empty:
        for cell_type in promoter_data['cell_type'].unique():
            ct_data = promoter_data[promoter_data['cell_type'] == cell_type]
            if not ct_data.empty:
                corr = ct_data['correlation_log2_ratio'].iloc[0]
                pval = ct_data['pvalue_log2_ratio'].iloc[0]
                n_genes = ct_data['n_genes'].iloc[0]
                stats_text += f"{cell_type} Promoter:\n"
                stats_text += f"  r = {corr:.3f}\n"
                stats_text += f"  p = {pval:.2e}\n"
                stats_text += f"  n = {n_genes:,}\n\n"

    # Add interpretation
    stats_text += "Key Findings:\n"
    sig_regions = results_df[(results_df['pvalue_log2_ratio'] < 0.05) &
                            (results_df['subset'] == 'all')]['region'].unique()
    stats_text += f"• {len(sig_regions)} regions show significant correlations\n"

    if 'promoter' in sig_regions:
        stats_text += "• Promoter binding correlates with expression\n"

    negative_corrs = results_df[(results_df['correlation_log2_ratio'] < 0) &
                               (results_df['pvalue_log2_ratio'] < 0.05) &
                               (results_df['subset'] == 'all')]
    if not negative_corrs.empty:
        stats_text += "• Negative correlations suggest MeCP2 repression\n"

    axes[1,2].text(0.05, 0.95, stats_text, transform=axes[1,2].transAxes,
                  fontsize=10, verticalalignment='top', fontfamily='monospace')

    plt.tight_layout()

    # Save
    output_path = output_dir / "comprehensive_summary_heatmap.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()

    logging.info("Saved comprehensive summary heatmap to %s", output_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create heatmap visualizations for MeCP2 binding-expression correlations"
    )
    parser.add_argument(
        "--nsc-results",
        type=Path,
        required=True,
        help="CSV file with NSC correlation results"
    )
    parser.add_argument(
        "--neu-results",
        type=Path,
        help="CSV file with NEU correlation results (optional)"
    )
    parser.add_argument(
        "--nsc-integrated",
        type=Path,
        required=True,
        help="CSV file with NSC integrated gene body + RNA-seq data"
    )
    parser.add_argument(
        "--neu-integrated",
        type=Path,
        help="CSV file with NEU integrated gene body + RNA-seq data (optional)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to save heatmap visualizations"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    return parser.parse_args()

def main() -> None:
    args = parse_args()
    configure_logging(args.verbose)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    results_df = load_correlation_results(args.nsc_results, args.neu_results)
    integrated_df = load_integrated_data(args.nsc_integrated, args.neu_integrated)

    # Create visualizations
    create_correlation_matrix_heatmap(results_df, args.output_dir)
    create_top_genes_binding_heatmap(integrated_df, args.output_dir)
    create_binding_ratio_expression_heatmap(integrated_df, args.output_dir)

    # Create cell type comparison if both available
    if len(results_df['cell_type'].unique()) > 1:
        create_cell_type_comparison_heatmap(results_df, args.output_dir)

    create_comprehensive_summary_plot(results_df, integrated_df, args.output_dir)

    logging.info("All heatmap visualizations completed successfully")

if __name__ == "__main__":
    main()




