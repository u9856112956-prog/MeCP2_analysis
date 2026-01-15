#!/usr/bin/env python3
"""
CpG Enrichment Analysis using BigWig Signal

This script analyzes MeCP2 enrichment at CpG islands by extracting signal
intensities from BigWig files (RPKM-normalized) using pyBigWig.

Signal extraction method:
- Mean signal intensities extracted from bigWig files for each replicate
- Requires signal in at least 2 replicates for robust quantification
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pyBigWig
import os
import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
from glob import glob

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CPG_WINDOW = 0
NUM_REPS = 2

@dataclass
class EnrichmentConfig:
    """Configuration for enrichment analysis"""
    min_signal_threshold: float = 0.0
    min_fold_change: float = 1.0
    max_qvalue: float = 0.05

class MeCP2CpGAnalyzer:
    def __init__(self, config: EnrichmentConfig):
        self.config = config

    def load_cpg_islands(self, cpg_file: str) -> pd.DataFrame:
        """Load CpG islands data"""
        try:
            cpg = pd.read_csv(cpg_file, sep='\t',
                            names=['chr', 'start', 'end', 'name', 'score', 'strand'])
            logger.info(f"Loaded {len(cpg)} CpG islands")
            return cpg
        except Exception as e:
            logger.error(f"Error loading CpG islands: {str(e)}")
            raise

    def get_signal_from_bigwig(
        self,
        bigwig_files: Dict[str, str],
        chrom: str,
        start: int,
        end: int,
        method: str = "mean"
    ) -> Tuple[float, List[float]]:
        """
        Calculate signal from multiple bigWig replicates

        Returns:
        - Tuple of (mean signal across replicates, list of individual signals)
        """
        signals = []

        for sample, bw_path in bigwig_files.items():
            with pyBigWig.open(bw_path) as bw:
                try:
                    signal = bw.stats(
                        chrom,
                        start,
                        end,
                        type=method
                    )[0]

                    signals.append(signal if signal is not None else 0.0)
                except RuntimeError:
                    signals.append(0.0)

        return np.mean(signals) if signals else 0.0, signals

    def load_peak_files(self, peak_dir: str) -> List[pd.DataFrame]:
        """Load individual peak files from directory"""
        peak_files = glob(f"{peak_dir}/*_peaks.broadPeak")
        peaks = []
        for peak_file in peak_files:
            df = pd.read_csv(peak_file, sep='\t',
                            names=['chr', 'start', 'end', 'name', 'score',
                                  'strand', 'signalValue', 'pValue',
                                  'qValue', 'peak'])
            peaks.append(df)
        return peaks

    def calculate_enrichment_chunk(self,
                                 chunk_id: int,
                                 total_chunks: int,
                                 exo_bigwigs: Dict[str, str],
                                 endo_bigwigs: Dict[str, str],
                                 cpg_islands: pd.DataFrame,
                                 exo_peaks_list: List[pd.DataFrame],
                                 endo_peaks_list: List[pd.DataFrame]) -> pd.DataFrame:
        """Calculate enrichment using bigWig signal across peak-defined regions"""

        # Calculate chunk boundaries
        chunk_size = len(cpg_islands) // total_chunks
        start_idx = chunk_id * chunk_size
        end_idx = start_idx + chunk_size if chunk_id < total_chunks - 1 else len(cpg_islands)

        cpg_chunk = cpg_islands.iloc[start_idx:end_idx]

        enrichment_data = []
        for _, cpg in cpg_chunk.iterrows():
            # Find overlapping peaks in each replicate
            exo_overlaps_by_rep = []
            for exo_peaks in exo_peaks_list:
                overlaps = exo_peaks[
                    (exo_peaks['chr'] == cpg['chr']) &
                    (exo_peaks['start'] <= cpg['end'] + CPG_WINDOW) &
                    (exo_peaks['end'] >= cpg['start'] - CPG_WINDOW)
                ]
                exo_overlaps_by_rep.append(overlaps)

            endo_overlaps_by_rep = []
            for endo_peaks in endo_peaks_list:
                overlaps = endo_peaks[
                    (endo_peaks['chr'] == cpg['chr']) &
                    (endo_peaks['start'] <= cpg['end'] + CPG_WINDOW) &
                    (endo_peaks['end'] >= cpg['start'] - CPG_WINDOW)
                ]
                endo_overlaps_by_rep.append(overlaps)

            # Check if we have peaks in at least NUM_REPS replicates
            exo_rep_with_peaks = sum(1 for overlaps in exo_overlaps_by_rep if not overlaps.empty)
            endo_rep_with_peaks = sum(1 for overlaps in endo_overlaps_by_rep if not overlaps.empty)

            # Skip if neither condition has at least NUM_REPS replicates with peaks
            if exo_rep_with_peaks < NUM_REPS and endo_rep_with_peaks < NUM_REPS:
                continue

            # Get the outer bounds of all overlapping peaks
            all_starts = []
            all_ends = []

            for overlaps in exo_overlaps_by_rep + endo_overlaps_by_rep:
                if not overlaps.empty:
                    all_starts.extend(overlaps['start'].tolist())
                    all_ends.extend(overlaps['end'].tolist())

            # Use the outermost boundaries of overlapping peaks
            region_start = min(all_starts)
            region_end = max(all_ends)

            # Get signals using the peak-defined boundaries
            exo_signal, exo_replicates = self.get_signal_from_bigwig(
                exo_bigwigs,
                cpg['chr'],
                region_start,
                region_end
            )

            endo_signal, endo_replicates = self.get_signal_from_bigwig(
                endo_bigwigs,
                cpg['chr'],
                region_start,
                region_end
            )

            # Check for signal in at least NUM_REPS replicate
            exo_has_signal = sum(1 for s in exo_replicates if s > 0) >= NUM_REPS
            endo_has_signal = sum(1 for s in endo_replicates if s > 0) >= NUM_REPS

            if not (exo_has_signal or endo_has_signal):
                continue

            # Statistical test between replicates
            if exo_replicates and endo_replicates:
                nonzero_exo = [s for s in exo_replicates if s > 0]
                nonzero_endo = [s for s in endo_replicates if s > 0]

                if len(nonzero_exo) >= 1 and len(nonzero_endo) >= 1:
                    _, pvalue = stats.mannwhitneyu(
                        nonzero_exo,
                        nonzero_endo,
                        alternative='two-sided'
                    )
                else:
                    pvalue = 1.0
            else:
                pvalue = 1.0

            # Calculate enrichment
            if endo_signal > 0:
                enrichment = exo_signal / endo_signal
            else:
                enrichment = float('inf') if exo_signal > 0 else 0.0

            # Determine binding type
            if exo_has_signal and endo_has_signal:
                binding_type = 'both'
            elif exo_has_signal:
                binding_type = 'exo_only'
            elif endo_has_signal:
                binding_type = 'endo_only'
            else:
                continue

            # Determine binding type by peaks
            if exo_rep_with_peaks >= NUM_REPS and endo_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'both'
            elif exo_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'exo_only'
            elif endo_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'endo_only'
            else:
                binding_type_by_peaks = 'none'

            # Store results with replicate information
            enrichment_data.append({
                'chr': cpg['chr'],
                'start': cpg['start'],
                'end': cpg['end'],
                'exo_signal': exo_signal,
                'endo_signal': endo_signal,
                'enrichment': enrichment,
                'pvalue': pvalue,
                'binding_type': binding_type,
                'binding_type_by_peaks': binding_type_by_peaks,
                'significant': (
                    (enrichment >= self.config.min_fold_change if endo_signal > 0 else exo_signal > 0) and
                    exo_signal > self.config.min_signal_threshold and
                    pvalue < self.config.max_qvalue
                ),
                'exo_replicates_with_signal': sum(1 for s in exo_replicates if s > 0),
                'endo_replicates_with_signal': sum(1 for s in endo_replicates if s > 0),
                'exo_replicate_signals': ','.join(map(str, exo_replicates)),
                'endo_replicate_signals': ','.join(map(str, endo_replicates)),
                'region_length': region_end - region_start,
                'cpg_length': cpg['end'] - cpg['start'],
                'cpg_score': cpg['score'],
                'cpg_name': cpg['name'],
                'exo_replicates_with_peaks': exo_rep_with_peaks,
                'endo_replicates_with_peaks': endo_rep_with_peaks,
                'exo_peak_scores_by_rep': ';'.join(','.join(map(str, rep['signalValue'].tolist()))
                                                  for rep in exo_overlaps_by_rep if not rep.empty),
                'endo_peak_scores_by_rep': ';'.join(','.join(map(str, rep['signalValue'].tolist()))
                                                   for rep in endo_overlaps_by_rep if not rep.empty),
                'region_start': region_start,
                'region_end': region_end
            })

        return pd.DataFrame(enrichment_data)

    def plot_enrichment_results(self,
                              enrichment_df: pd.DataFrame,
                              output_dir: str):
        """Create visualization of enrichment results"""
        os.makedirs(output_dir, exist_ok=True)

        # Plot 1: Enrichment distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=enrichment_df, x='enrichment', bins=50)
        plt.title('Distribution of MeCP2 Enrichment at CpG Islands')
        plt.xlabel('Enrichment (Exo/Endo)')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'enrichment_distribution.pdf'))
        plt.close()

        # Plot 2: Signal comparison
        plt.figure(figsize=(10, 6))
        plt.scatter(enrichment_df['endo_signal'],
                   enrichment_df['exo_signal'],
                   alpha=0.5)
        plt.xlabel('Endogenous Signal')
        plt.ylabel('Exogenous Signal')
        plt.title('Exogenous vs Endogenous MeCP2 Signal at CpG Islands')
        plt.savefig(os.path.join(output_dir, 'signal_comparison.pdf'))
        plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyze MeCP2 enrichment at CpG islands using bigWig files')
    parser.add_argument('--bigwig-dir', type=str, required=True,
                       help='Directory containing bigWig files')
    parser.add_argument('--cpg-islands', type=str, required=True,
                       help='Path to CpG islands file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    parser.add_argument('--chunk-id', type=int, required=True,
                       help='Chunk ID for parallel processing')
    parser.add_argument('--total-chunks', type=int, required=True,
                       help='Total number of chunks')
    parser.add_argument('--exo-peaks', type=str, required=True,
                       help='Directory containing exogenous peak files')
    parser.add_argument('--endo-peaks', type=str, required=True,
                       help='Directory containing endogenous peak files')
    parser.add_argument('--cell-type', type=str, required=True, choices=['Neu', 'NSC'],
                       help='Cell type to analyze (Neu or NSC)')
    args = parser.parse_args()

    # Initialize analyzer
    config = EnrichmentConfig()
    analyzer = MeCP2CpGAnalyzer(config)

    # Load CpG islands
    logger.info(f"Loading CpG islands from: {args.cpg_islands}")
    cpg_islands = analyzer.load_cpg_islands(args.cpg_islands)

    # Get bigWig files
    bigwig_dir = Path(args.bigwig_dir)

    # Get peak files from directories
    exo_peak_files = list(Path(args.exo_peaks).glob("*.broadPeak"))
    endo_peak_files = list(Path(args.endo_peaks).glob("*.broadPeak"))

    logger.info(f"Found {len(exo_peak_files)} exogenous and {len(endo_peak_files)} endogenous peak files")

    # Separate exo and endo bigWig files
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

    logger.info(f"Found {len(exo_bigwigs)} exogenous and {len(endo_bigwigs)} endogenous bigWig files for {args.cell_type}")

    # Load peak files
    exo_peaks_list = analyzer.load_peak_files(args.exo_peaks)
    endo_peaks_list = analyzer.load_peak_files(args.endo_peaks)

    logger.info(f"Loaded {len(exo_peaks_list)} exogenous and {len(endo_peaks_list)} endogenous peak files for {args.cell_type}")

    # Calculate enrichment with peaks
    chunk_results = analyzer.calculate_enrichment_chunk(
        args.chunk_id,
        args.total_chunks,
        exo_bigwigs,
        endo_bigwigs,
        cpg_islands,
        exo_peaks_list,
        endo_peaks_list
    )

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Save chunk results
    chunk_file = os.path.join(args.output_dir, f'chunk_{args.chunk_id}.csv')
    chunk_results.to_csv(chunk_file, index=False)

    logger.info(f"Chunk {args.chunk_id} complete. Results saved to {chunk_file}")

if __name__ == "__main__":
    main()
