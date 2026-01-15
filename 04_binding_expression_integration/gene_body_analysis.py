#!/usr/bin/env python3
"""
Enhanced Gene Body MeCP2 Binding Analysis

This script provides comprehensive analysis of MeCP2 binding patterns in gene bodies,
including position-specific analysis, length normalization, and detailed genomic context.

Key Features:
- Position-specific analysis (5', middle, 3' regions)
- Length-normalized signal quantification
- Exon vs intron binding analysis
- Promoter vs gene body comparison
- RNA-seq integration capabilities
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pyBigWig
import os
import logging
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from pathlib import Path
from glob import glob
import argparse
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class GeneBodyConfig:
    """Configuration for enhanced gene body analysis"""
    # Signal thresholds
    min_signal_threshold: float = 0.1
    min_fold_change: float = 1.2
    max_qvalue: float = 0.05
    
    # Gene length filters
    min_gene_length: int = 1000  # Minimum gene length for analysis
    max_gene_length: int = 1000000  # Maximum gene length
    
    # Position-specific analysis
    region_5prime_fraction: float = 0.2  # First 20% of gene
    region_3prime_fraction: float = 0.2  # Last 20% of gene
    # Middle fraction is calculated as 1 - (5prime + 3prime)
    
    # Promoter definition
    promoter_upstream: int = 2000
    promoter_downstream: int = 500
    
    # Analysis windows
    bin_size: int = 100  # Size of bins for detailed analysis
    smoothing_window: int = 5  # Window for signal smoothing

class EnhancedGeneBodyAnalyzer:
    """Enhanced analyzer for MeCP2 binding in gene bodies"""

    def __init__(self, config: GeneBodyConfig):
        self.config = config
        self._bigwig_cache = {}

    def _get_bigwig_handle(self, path: str):
        """Return a cached pyBigWig handle for the provided path"""
        handle = self._bigwig_cache.get(path)
        if handle is None:
            try:
                handle = pyBigWig.open(path)
                self._bigwig_cache[path] = handle
            except Exception as exc:
                logger.warning(f"Failed to open bigWig file {path}: {exc}")
                raise
        return handle

    def close_bigwig_handles(self):
        """Close any cached bigWig handles"""
        for path, handle in list(self._bigwig_cache.items()):
            try:
                handle.close()
            except Exception as exc:
                logger.debug(f"Error closing bigWig handle for {path}: {exc}")
        self._bigwig_cache.clear()

    def load_gtf_annotations(self, gtf_file: str, gene_types: List[str] = None) -> pd.DataFrame:
        """
        Load comprehensive gene annotations from GTF file
        
        Returns DataFrame with genes, exons, and promoter information
        """
        if gene_types is None:
            gene_types = ['protein_coding', 'lncRNA']
            
        logger.info(f"Loading GTF annotations from: {gtf_file}")
        
        genes = []
        exons = []
        
        def parse_attributes(attr_string):
            """Parse GTF attribute string"""
            attributes = {}
            for item in attr_string.split(';'):
                if item.strip():
                    parts = item.strip().split(' ', 1)
                    if len(parts) == 2:
                        key = parts[0]
                        value = parts[1].strip('"')
                        attributes[key] = value
            return attributes
        
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                    
                if line_num % 100000 == 0:
                    logger.info(f"Processed {line_num} GTF lines...")
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom, source, feature_type, start, end, score, strand, phase, attributes_str = fields
                start, end = int(start) - 1, int(end)  # Convert to 0-based
                attributes = parse_attributes(attributes_str)
                
                gene_type = attributes.get('gene_type', 'unknown')
                if gene_types and gene_type not in gene_types:
                    continue
                
                if feature_type == 'gene':
                    gene_length = end - start
                    if self.config.min_gene_length <= gene_length <= self.config.max_gene_length:
                        # Calculate promoter region
                        if strand == '+':
                            promoter_start = max(0, start - self.config.promoter_upstream)
                            promoter_end = start + self.config.promoter_downstream
                        else:
                            promoter_start = max(0, end - self.config.promoter_downstream)
                            promoter_end = end + self.config.promoter_upstream
                        
                        # Calculate strand-aware position-specific regions
                        region_5_len = max(1, int(round(gene_length * self.config.region_5prime_fraction)))
                        region_3_len = max(1, int(round(gene_length * self.config.region_3prime_fraction)))
                        max_alloc = max(1, gene_length - 2)
                        if region_5_len + region_3_len >= gene_length:
                            scale = max_alloc / float(region_5_len + region_3_len)
                            region_5_len = max(1, int(round(region_5_len * scale)))
                            region_3_len = max(1, int(round(region_3_len * scale)))

                        forward_boundary = min(end, start + region_5_len)
                        reverse_boundary = max(start, end - region_3_len)

                        if strand == '+':
                            region_5prime_start, region_5prime_end = start, forward_boundary
                            region_3prime_start, region_3prime_end = reverse_boundary, end
                        else:
                            region_5prime_start, region_5prime_end = reverse_boundary, end
                            region_3prime_start, region_3prime_end = start, forward_boundary

                        region_middle_start = max(start, min(region_5prime_end, region_3prime_end))
                        region_middle_end = min(end, max(region_5prime_start, region_3prime_start))
                        if region_middle_end < region_middle_start:
                            midpoint = start + gene_length // 2
                            region_middle_start, region_middle_end = sorted((region_middle_start, midpoint))

                        region_5prime_start = max(start, min(region_5prime_start, region_5prime_end))
                        region_5prime_end = min(end, max(region_5prime_start, region_5prime_end))
                        region_3prime_start = max(start, min(region_3prime_start, region_3prime_end))
                        region_3prime_end = min(end, max(region_3prime_start, region_3prime_end))
                        region_middle_start, region_middle_end = sorted((max(start, region_middle_start), min(end, region_middle_end)))

                        genes.append({
                            'chr': chrom,
                            'start': start,
                            'end': end,
                            'length': gene_length,
                            'strand': strand,
                            'gene_id': attributes.get('gene_id', 'unknown'),
                            'gene_name': attributes.get('gene_name', attributes.get('gene_id', 'unknown')),
                            'gene_type': gene_type,
                            
                            # Promoter region
                            'promoter_start': promoter_start,
                            'promoter_end': promoter_end,
                            
                            # Position-specific regions
                            'region_5prime_start': region_5prime_start,
                            'region_5prime_end': region_5prime_end,
                            'region_middle_start': region_middle_start,
                            'region_middle_end': region_middle_end,
                            'region_3prime_start': region_3prime_start,
                            'region_3prime_end': region_3prime_end,
                            
                            # Region lengths
                            'region_5prime_length': max(0, region_5prime_end - region_5prime_start),
                            'region_middle_length': max(0, region_middle_end - region_middle_start),
                            'region_3prime_length': max(0, region_3prime_end - region_3prime_start),
                        })
                
                elif feature_type == 'exon':
                    exons.append({
                        'chr': chrom,
                        'start': start,
                        'end': end,
                        'length': end - start,
                        'strand': strand,
                        'gene_id': attributes.get('gene_id', 'unknown'),
                        'transcript_id': attributes.get('transcript_id', 'unknown'),
                        'exon_number': attributes.get('exon_number', 'unknown')
                    })
        
        genes_df = pd.DataFrame(genes)
        exons_df = pd.DataFrame(exons)
        
        logger.info(f"Loaded {len(genes_df)} genes and {len(exons_df)} exons")
        
        return genes_df, exons_df
    
    def calculate_signal_from_bigwig(self,                                    bigwig_files: Dict[str, str],
                                   chrom: str,
                                   start: int,
                                   end: int,
                                   normalize_by_length: bool = True) -> Tuple[float, List[float], Dict]:
        """Calculate signal statistics from bigWig files with cached handles."""
        region_length = max(0, end - start)
        if region_length == 0 or start >= end:
            zero_stats = {
                'mean_normalized_signal': 0.0,
                'std_normalized_signal': 0.0,
                'median_normalized_signal': 0.0,
                'mean_raw_signal': 0.0,
                'n_replicates': len(bigwig_files),
                'n_nonzero_replicates': 0,
                'region_length': 0
            }
            return 0.0, [0.0] * len(bigwig_files), zero_stats

        signals: List[float] = []
        raw_means: List[float] = []
        for sample, bw_path in bigwig_files.items():
            try:
                bw = self._get_bigwig_handle(bw_path)
                signal_mean = bw.stats(chrom, start, end, type="mean")[0]
                signal_sum = bw.stats(chrom, start, end, type="sum")[0]
            except Exception as exc:
                logger.warning(f"Error processing {sample} ({bw_path}): {exc}")
                signal_mean, signal_sum = None, None

            if signal_mean is None or signal_sum is None:
                normalized_signal = 0.0
                mean_value = 0.0
            else:
                mean_value = float(signal_mean)
                if normalize_by_length:
                    normalized_signal = (float(signal_sum) * 1000.0) / region_length
                else:
                    normalized_signal = mean_value

            signals.append(normalized_signal)
            raw_means.append(mean_value)

        detailed_stats = {
            'mean_normalized_signal': float(np.mean(signals)) if signals else 0.0,
            'std_normalized_signal': float(np.std(signals)) if signals else 0.0,
            'median_normalized_signal': float(np.median(signals)) if signals else 0.0,
            'mean_raw_signal': float(np.mean(raw_means)) if raw_means else 0.0,
            'n_replicates': len(signals),
            'n_nonzero_replicates': sum(1 for s in signals if s > 0),
            'region_length': region_length
        }

        mean_signal = detailed_stats['mean_normalized_signal'] if signals else 0.0
        return mean_signal, signals, detailed_stats
    
    def analyze_gene_body_positions(self,
                                  genes_df: pd.DataFrame,
                                  exo_bigwigs: Dict[str, str],
                                  endo_bigwigs: Dict[str, str],
                                  chunk_id: int = 0,
                                  total_chunks: int = 1) -> pd.DataFrame:
        """
        Analyze MeCP2 binding across different gene body positions
        
        Analyzes:
        1. Promoter regions
        2. Gene body 5' region
        3. Gene body middle region  
        4. Gene body 3' region
        5. Full gene body
        """
        
        # Calculate chunk boundaries
        chunk_size = math.ceil(len(genes_df) / total_chunks) if total_chunks else len(genes_df)
        start_idx = chunk_id * chunk_size
        end_idx = min(len(genes_df), start_idx + chunk_size)
        
        gene_chunk = genes_df.iloc[start_idx:end_idx].copy()
        logger.info(f"Processing chunk {chunk_id}: {len(gene_chunk)} genes")
        
        results = []
        
        for idx, gene in gene_chunk.iterrows():
            if idx % 100 == 0:
                logger.info(f"Processed {idx - start_idx} genes in chunk {chunk_id}")
            
            # Define regions to analyze
            regions = {
                'promoter': (gene['promoter_start'], gene['promoter_end']),
                'gene_body_full': (gene['start'], gene['end']),
                'gene_body_5prime': (gene['region_5prime_start'], gene['region_5prime_end']),
                'gene_body_middle': (gene['region_middle_start'], gene['region_middle_end']),
                'gene_body_3prime': (gene['region_3prime_start'], gene['region_3prime_end'])
            }
            
            gene_result = {
                'gene_id': gene['gene_id'],
                'gene_name': gene['gene_name'],
                'chr': gene['chr'],
                'start': gene['start'],
                'end': gene['end'],
                'length': gene['length'],
                'strand': gene['strand'],
                'gene_type': gene['gene_type']
            }
            
            # Analyze each region
            for region_name, (region_start, region_end) in regions.items():
                if region_start >= region_end:
                    continue
                
                # Calculate exogenous signals
                exo_mean, exo_signals, exo_stats = self.calculate_signal_from_bigwig(
                    exo_bigwigs, gene['chr'], region_start, region_end
                )
                
                # Calculate endogenous signals
                endo_mean, endo_signals, endo_stats = self.calculate_signal_from_bigwig(
                    endo_bigwigs, gene['chr'], region_start, region_end
                )
                
                # Statistical comparison
                if len([s for s in exo_signals if s > 0]) >= 1 and len([s for s in endo_signals if s > 0]) >= 1:
                    try:
                        _, p_value = stats.mannwhitneyu(
                            [s for s in exo_signals if s > 0],
                            [s for s in endo_signals if s > 0],
                            alternative='two-sided'
                        )
                    except:
                        p_value = 1.0
                else:
                    p_value = 1.0
                
                # Calculate enrichment
                if endo_mean > 0:
                    enrichment = exo_mean / endo_mean
                else:
                    enrichment = float('inf') if exo_mean > 0 else 0.0
                
                # Store region-specific results
                gene_result.update({
                    f'{region_name}_exo_signal': exo_mean,
                    f'{region_name}_endo_signal': endo_mean,
                    f'{region_name}_enrichment': enrichment,
                    f'{region_name}_pvalue': p_value,
                    f'{region_name}_exo_n_nonzero': exo_stats['n_nonzero_replicates'],
                    f'{region_name}_endo_n_nonzero': endo_stats['n_nonzero_replicates'],
                    f'{region_name}_exo_std': exo_stats['std_normalized_signal'],
                    f'{region_name}_endo_std': endo_stats['std_normalized_signal'],
                    f'{region_name}_length': region_end - region_start,
                    f'{region_name}_start': region_start,
                    f'{region_name}_end': region_end
                })
            
            # Calculate promoter vs gene body ratios
            if (gene_result['promoter_exo_signal'] > 0 and 
                gene_result['gene_body_full_exo_signal'] > 0):
                gene_result['promoter_vs_genebody_ratio_exo'] = (
                    gene_result['promoter_exo_signal'] / gene_result['gene_body_full_exo_signal']
                )
            else:
                gene_result['promoter_vs_genebody_ratio_exo'] = 0.0
            
            if (gene_result['promoter_endo_signal'] > 0 and 
                gene_result['gene_body_full_endo_signal'] > 0):
                gene_result['promoter_vs_genebody_ratio_endo'] = (
                    gene_result['promoter_endo_signal'] / gene_result['gene_body_full_endo_signal']
                )
            else:
                gene_result['promoter_vs_genebody_ratio_endo'] = 0.0
            
            results.append(gene_result)
        
        results_df = pd.DataFrame(results)
        logger.info(f"Completed analysis for chunk {chunk_id}: {len(results_df)} genes")
        
        return results_df

    def analyze_exon_vs_intron_binding(self,
                                     genes_df: pd.DataFrame,
                                     exons_df: pd.DataFrame,
                                     exo_bigwigs: Dict[str, str],
                                     endo_bigwigs: Dict[str, str],
                                     chunk_id: int = 0,
                                     total_chunks: int = 1) -> pd.DataFrame:
        """Analyze MeCP2 binding in exonic vs intronic regions within gene bodies"""

        chunk_size = math.ceil(len(genes_df) / total_chunks) if total_chunks else len(genes_df)
        start_idx = chunk_id * chunk_size
        end_idx = len(genes_df) if chunk_id == total_chunks - 1 else min(len(genes_df), start_idx + chunk_size)

        gene_chunk = genes_df.iloc[start_idx:end_idx].copy()
        logger.info(f"Processing exon/intron analysis for chunk {chunk_id}: {len(gene_chunk)} genes")

        results = []

        for _, gene in gene_chunk.iterrows():
            gene_exons = exons_df[exons_df['gene_id'] == gene['gene_id']].copy()
            if gene_exons.empty:
                continue

            merged_exons: List[Tuple[int, int]] = []
            for start, end in gene_exons[['start', 'end']].sort_values('start').itertuples(index=False):
                if not merged_exons or start > merged_exons[-1][1]:
                    merged_exons.append([start, end])
                else:
                    merged_exons[-1][1] = max(merged_exons[-1][1], end)

            total_exonic_length = sum(end - start for start, end in merged_exons)
            total_intronic_length = max(0, gene['length'] - total_exonic_length)
            if total_exonic_length == 0 or total_intronic_length == 0:
                continue

            exonic_sum_exo = 0.0
            exonic_sum_endo = 0.0
            for start, end in merged_exons:
                exon_len = max(0, end - start)
                if exon_len == 0:
                    continue

                exo_mean, _, _ = self.calculate_signal_from_bigwig(
                    exo_bigwigs, gene['chr'], start, end, normalize_by_length=False
                )
                endo_mean, _, _ = self.calculate_signal_from_bigwig(
                    endo_bigwigs, gene['chr'], start, end, normalize_by_length=False
                )
                exonic_sum_exo += max(0.0, exo_mean) * exon_len
                exonic_sum_endo += max(0.0, endo_mean) * exon_len

            weighted_exonic_signal_exo = exonic_sum_exo / total_exonic_length if total_exonic_length else 0.0
            weighted_exonic_signal_endo = exonic_sum_endo / total_exonic_length if total_exonic_length else 0.0

            total_exo_mean, _, _ = self.calculate_signal_from_bigwig(
                exo_bigwigs, gene['chr'], gene['start'], gene['end'], normalize_by_length=False
            )
            total_endo_mean, _, _ = self.calculate_signal_from_bigwig(
                endo_bigwigs, gene['chr'], gene['start'], gene['end'], normalize_by_length=False
            )

            total_exo_sum = max(0.0, total_exo_mean) * gene['length']
            total_endo_sum = max(0.0, total_endo_mean) * gene['length']

            intronic_sum_exo = max(0.0, total_exo_sum - exonic_sum_exo)
            intronic_sum_endo = max(0.0, total_endo_sum - exonic_sum_endo)

            intronic_signal_exo = intronic_sum_exo / total_intronic_length if total_intronic_length else 0.0
            intronic_signal_endo = intronic_sum_endo / total_intronic_length if total_intronic_length else 0.0

            exon_intron_ratio_exo = (
                weighted_exonic_signal_exo / intronic_signal_exo if intronic_signal_exo > 0 else float('inf')
            )
            exon_intron_ratio_endo = (
                weighted_exonic_signal_endo / intronic_signal_endo if intronic_signal_endo > 0 else float('inf')
            )

            results.append({
                'gene_id': gene['gene_id'],
                'gene_name': gene['gene_name'],
                'chr': gene['chr'],
                'start': gene['start'],
                'end': gene['end'],
                'length': gene['length'],
                'strand': gene['strand'],
                'gene_type': gene['gene_type'],
                'n_exons': len(merged_exons),
                'total_exonic_length': total_exonic_length,
                'total_intronic_length': total_intronic_length,
                'exonic_fraction': total_exonic_length / gene['length'] if gene['length'] else 0,
                'exonic_signal_exo': weighted_exonic_signal_exo,
                'exonic_signal_endo': weighted_exonic_signal_endo,
                'intronic_signal_exo': intronic_signal_exo,
                'intronic_signal_endo': intronic_signal_endo,
                'exon_intron_ratio_exo': exon_intron_ratio_exo,
                'exon_intron_ratio_endo': exon_intron_ratio_endo,
                'total_gene_signal_exo': max(0.0, total_exo_mean),
                'total_gene_signal_endo': max(0.0, total_endo_mean)
            })

        results_df = pd.DataFrame(results)
        logger.info(f"Completed exon/intron analysis for chunk {chunk_id}: {len(results_df)} genes")

        return results_df
    
    def create_comprehensive_visualizations(self, 
                                          position_results: pd.DataFrame,
                                          exon_intron_results: pd.DataFrame,
                                          output_dir: str,
                                          cell_type: str):
        """Create comprehensive visualizations for gene body analysis"""
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Position-specific signal comparison
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'MeCP2 Binding Across Gene Body Positions - {cell_type} Cells', fontsize=16)
        
        regions = ['promoter', 'gene_body_5prime', 'gene_body_middle', 'gene_body_3prime', 'gene_body_full']
        region_labels = ['Promoter', '5\' Gene Body', 'Middle Gene Body', '3\' Gene Body', 'Full Gene Body']
        
        for i, (region, label) in enumerate(zip(regions, region_labels)):
            if i >= 6:  # We only have 6 subplots
                break
            
            row, col = i // 3, i % 3
            ax = axes[row, col]
            
            exo_col = f'{region}_exo_signal'
            endo_col = f'{region}_endo_signal'
            
            if exo_col in position_results.columns and endo_col in position_results.columns:
                # Remove zeros for better visualization
                data_subset = position_results[
                    (position_results[exo_col] > 0) | (position_results[endo_col] > 0)
                ].copy()
                
                if len(data_subset) > 0:
                    ax.scatter(data_subset[endo_col], data_subset[exo_col], 
                             alpha=0.6, s=20)
                    ax.set_xlabel('Endogenous MeCP2 Signal')
                    ax.set_ylabel('Exogenous MeCP2 Signal')
                    ax.set_title(f'{label}')
                    
                    # Add diagonal line
                    max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
                    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
        
        # Remove empty subplot
        if len(regions) < 6:
            axes[1, 2].remove()
        
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f'position_specific_signals_{cell_type}.png', 
                   dpi=300, bbox_inches='tight')
        plt.savefig(Path(output_dir) / f'position_specific_signals_{cell_type}.pdf', 
                   bbox_inches='tight')
        plt.close()
        
        # 2. Enrichment comparison across positions
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        enrichment_data = []
        for region, label in zip(regions, region_labels):
            enrich_col = f'{region}_enrichment'
            if enrich_col in position_results.columns:
                # Filter infinite values and extreme outliers
                enrichments = position_results[enrich_col].replace([np.inf, -np.inf], np.nan)
                enrichments = enrichments[(enrichments > 0) & (enrichments < 100)]
                
                for val in enrichments:
                    if not np.isnan(val):
                        enrichment_data.append({
                            'Region': label,
                            'Enrichment': val
                        })
        
        if enrichment_data:
            enrich_df = pd.DataFrame(enrichment_data)
            sns.boxplot(data=enrich_df, x='Region', y='Enrichment', ax=ax)
            ax.set_title(f'MeCP2 Enrichment (Exo/Endo) Across Gene Body Regions - {cell_type}')
            ax.set_ylabel('Enrichment Ratio (Exogenous/Endogenous)')
            plt.xticks(rotation=45)
            ax.axhline(y=1, color='red', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f'enrichment_by_position_{cell_type}.png', 
                   dpi=300, bbox_inches='tight')
        plt.savefig(Path(output_dir) / f'enrichment_by_position_{cell_type}.pdf', 
                   bbox_inches='tight')
        plt.close()
        
        # 3. Promoter vs Gene Body comparison
        if not position_results.empty:
            fig, axes = plt.subplots(1, 2, figsize=(15, 6))
            fig.suptitle(f'Promoter vs Gene Body MeCP2 Binding - {cell_type}', fontsize=14)
            
            # Exogenous
            valid_data = position_results[
                (position_results['promoter_exo_signal'] > 0) & 
                (position_results['gene_body_full_exo_signal'] > 0)
            ]
            
            if len(valid_data) > 0:
                axes[0].scatter(valid_data['gene_body_full_exo_signal'], 
                              valid_data['promoter_exo_signal'],
                              alpha=0.6, color='blue')
                axes[0].set_xlabel('Gene Body Exogenous Signal')
                axes[0].set_ylabel('Promoter Exogenous Signal')
                axes[0].set_title('Exogenous MeCP2')
                
                # Add diagonal line
                max_val = max(axes[0].get_xlim()[1], axes[0].get_ylim()[1])
                axes[0].plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
            
            # Endogenous
            valid_data_endo = position_results[
                (position_results['promoter_endo_signal'] > 0) & 
                (position_results['gene_body_full_endo_signal'] > 0)
            ]
            
            if len(valid_data_endo) > 0:
                axes[1].scatter(valid_data_endo['gene_body_full_endo_signal'], 
                              valid_data_endo['promoter_endo_signal'],
                              alpha=0.6, color='red')
                axes[1].set_xlabel('Gene Body Endogenous Signal')
                axes[1].set_ylabel('Promoter Endogenous Signal')
                axes[1].set_title('Endogenous MeCP2')
                
                # Add diagonal line
                max_val = max(axes[1].get_xlim()[1], axes[1].get_ylim()[1])
                axes[1].plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
            
            plt.tight_layout()
            plt.savefig(Path(output_dir) / f'promoter_vs_genebody_{cell_type}.png', 
                       dpi=300, bbox_inches='tight')
            plt.savefig(Path(output_dir) / f'promoter_vs_genebody_{cell_type}.pdf', 
                       bbox_inches='tight')
            plt.close()
        
        # 4. Exon vs Intron analysis (if available)
        if not exon_intron_results.empty:
            fig, axes = plt.subplots(1, 2, figsize=(15, 6))
            fig.suptitle(f'Exonic vs Intronic MeCP2 Binding - {cell_type}', fontsize=14)
            
            # Exogenous
            valid_exo = exon_intron_results[
                (exon_intron_results['exonic_signal_exo'] > 0) & 
                (exon_intron_results['intronic_signal_exo'] > 0)
            ]
            
            if len(valid_exo) > 0:
                axes[0].scatter(valid_exo['intronic_signal_exo'], 
                              valid_exo['exonic_signal_exo'],
                              alpha=0.6, color='blue')
                axes[0].set_xlabel('Intronic Exogenous Signal')
                axes[0].set_ylabel('Exonic Exogenous Signal')
                axes[0].set_title('Exogenous MeCP2')
                
                # Add diagonal line
                max_val = max(axes[0].get_xlim()[1], axes[0].get_ylim()[1])
                axes[0].plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
            
            # Endogenous
            valid_endo = exon_intron_results[
                (exon_intron_results['exonic_signal_endo'] > 0) & 
                (exon_intron_results['intronic_signal_endo'] > 0)
            ]
            
            if len(valid_endo) > 0:
                axes[1].scatter(valid_endo['intronic_signal_endo'], 
                              valid_endo['exonic_signal_endo'],
                              alpha=0.6, color='red')
                axes[1].set_xlabel('Intronic Endogenous Signal')
                axes[1].set_ylabel('Exonic Endogenous Signal')
                axes[1].set_title('Endogenous MeCP2')
                
                # Add diagonal line
                max_val = max(axes[1].get_xlim()[1], axes[1].get_ylim()[1])
                axes[1].plot([0, max_val], [0, max_val], 'r--', alpha=0.5)
            
            plt.tight_layout()
            plt.savefig(Path(output_dir) / f'exon_vs_intron_{cell_type}.png', 
                       dpi=300, bbox_inches='tight')
            plt.savefig(Path(output_dir) / f'exon_vs_intron_{cell_type}.pdf', 
                       bbox_inches='tight')
            plt.close()
        
        logger.info(f"Visualizations saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Enhanced Gene Body MeCP2 Binding Analysis'
    )
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--bigwig-dir', required=True, help='Directory containing bigWig files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--cell-type', required=True, choices=['NSC', 'Neu'], 
                       help='Cell type to analyze')
    parser.add_argument('--chunk-id', type=int, default=0, help='Chunk ID for parallel processing')
    parser.add_argument('--total-chunks', type=int, default=1, help='Total number of chunks')
    parser.add_argument('--analysis-type', choices=['position', 'exon-intron', 'both'], 
                       default='both', help='Type of analysis to perform')
    parser.add_argument('--gene-types', nargs='+', default=['protein_coding', 'lncRNA'],
                       help='Gene types to analyze')
    parser.add_argument('--min-gene-length', type=int, default=1000,
                       help='Minimum gene length for analysis')
    parser.add_argument('--max-gene-length', type=int, default=1000000,
                       help='Maximum gene length for analysis')
    
    args = parser.parse_args()
    
    # Create configuration
    config = GeneBodyConfig(
        min_gene_length=args.min_gene_length,
        max_gene_length=args.max_gene_length
    )
    
    # Initialize analyzer
    analyzer = EnhancedGeneBodyAnalyzer(config)
    
    # Load annotations
    logger.info("Loading GTF annotations...")
    genes_df, exons_df = analyzer.load_gtf_annotations(args.gtf, args.gene_types)
    
    # Get bigWig files
    bigwig_dir = Path(args.bigwig_dir)
    
    exo_bigwigs = {
        f.stem: str(f)
        for f in bigwig_dir.glob("*.bw")
        if f.stem.startswith(f'{args.cell_type}v') or f.stem.startswith(f'{args.cell_type}V')
    }
    
    endo_bigwigs = {
        f.stem: str(f)
        for f in bigwig_dir.glob("*.bw")
        if f.stem.startswith(f'{args.cell_type}M')
    }
    
    logger.info(f"Found {len(exo_bigwigs)} exogenous and {len(endo_bigwigs)} endogenous bigWig files")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    position_results = pd.DataFrame()
    exon_intron_results = pd.DataFrame()

    try:
        if args.analysis_type in ['position', 'both']:
            logger.info("Performing position-specific analysis...")
            position_results = analyzer.analyze_gene_body_positions(
                genes_df, exo_bigwigs, endo_bigwigs, args.chunk_id, args.total_chunks
            )

            position_file = Path(args.output_dir) / f'position_analysis_chunk_{args.chunk_id}.csv'
            position_results.to_csv(position_file, index=False)
            logger.info(f"Position analysis results saved to {position_file}")

        if args.analysis_type in ['exon-intron', 'both']:
            logger.info("Performing exon vs intron analysis...")
            exon_intron_results = analyzer.analyze_exon_vs_intron_binding(
                genes_df, exons_df, exo_bigwigs, endo_bigwigs, args.chunk_id, args.total_chunks
            )

            exon_intron_file = Path(args.output_dir) / f'exon_intron_analysis_chunk_{args.chunk_id}.csv'
            exon_intron_results.to_csv(exon_intron_file, index=False)
            logger.info(f"Exon/intron analysis results saved to {exon_intron_file}")

        if args.total_chunks == 1 and not position_results.empty:
            logger.info("Creating comprehensive visualizations...")
            analyzer.create_comprehensive_visualizations(
                position_results, exon_intron_results, args.output_dir, args.cell_type
            )
    finally:
        analyzer.close_bigwig_handles()


    logger.info("Enhanced gene body analysis complete!")


if __name__ == "__main__":
    main()







