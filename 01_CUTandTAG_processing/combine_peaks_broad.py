import pandas as pd
import numpy as np
from typing import List, Tuple
import os
import argparse
from intervaltree import IntervalTree, Interval

def load_peak_file(file_path: str) -> pd.DataFrame:
    """Load broadPeak format file."""
    columns = ['chr', 'start', 'end', 'name', 'score', 
              'strand', 'signalValue', 'pValue', 'qValue']
    
    df = pd.read_csv(file_path, sep='\t', names=columns)
    return df

def create_interval_tree(peaks_df: pd.DataFrame) -> dict:
    """Create interval trees for each chromosome from peaks."""
    trees = {}
    for _, peak in peaks_df.iterrows():
        chrom = peak['chr']
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom].add(Interval(peak['start'], peak['end'], peak))
    return trees

def find_overlapping_peaks(trees: dict, peak: pd.Series, overlap_threshold: float = 0.25, 
                          max_gap: int = 2000) -> List[Interval]:
    """Find overlapping peaks in interval trees.
    
    Args:
        trees: Dictionary of interval trees
        peak: Peak to find overlaps for
        overlap_threshold: Minimum overlap fraction required (lowered for low depth data)
        max_gap: Maximum distance between peaks to consider them as potentially mergeable
    
    Returns:
        List of overlapping intervals
    """
    chrom = peak['chr']
    if chrom not in trees:
        return []
    
    # Look for overlaps in an extended region
    overlaps = trees[chrom].overlap(peak['start'] - max_gap, peak['end'] + max_gap)
    matching_intervals = []
    
    for interval in overlaps:
        # Check if peaks overlap or are within max_gap distance
        if (min(peak['end'] + max_gap, interval.end) - 
            max(peak['start'] - max_gap, interval.begin) > -max_gap):
            matching_intervals.append(interval)
            
    return matching_intervals

def combine_replicates(peak_files: List[str], min_replicates: int = 2, 
                      overlap_threshold: float = 0.25, max_gap: int = 2000) -> pd.DataFrame:
    """Combine peaks from replicates, requiring peaks to be present in at least min_replicates."""
    # Load all peak files
    peak_dfs = [load_peak_file(f) for f in peak_files]
    
    # Create interval trees for all replicates except the first
    trees = [create_interval_tree(df) for df in peak_dfs[1:]]
    
    combined_peaks = []
    processed_peaks = set()  # Track which peaks we've already processed
    
    # Start with peaks from first replicate
    for idx, peak in peak_dfs[0].iterrows():
        if idx in processed_peaks:
            continue
            
        # Find all overlapping peaks from other replicates
        all_overlapping_peaks = [peak]
        peaks_to_check = [(0, idx, peak)]  # (replicate_idx, peak_idx, peak)
        
        while peaks_to_check:
            rep_idx, peak_idx, current_peak = peaks_to_check.pop(0)
            
            # Look for overlaps in all other replicates
            for other_rep_idx, tree in enumerate(trees):
                overlaps = find_overlapping_peaks(tree, current_peak, 
                                                overlap_threshold, max_gap)
                
                for overlap in overlaps:
                    other_peak = overlap.data
                    other_peak_idx = peak_dfs[other_rep_idx + 1].index[
                        (peak_dfs[other_rep_idx + 1]['chr'] == other_peak['chr']) &
                        (peak_dfs[other_rep_idx + 1]['start'] == other_peak['start']) &
                        (peak_dfs[other_rep_idx + 1]['end'] == other_peak['end'])
                    ][0]
                    
                    if other_peak_idx not in processed_peaks:
                        all_overlapping_peaks.append(other_peak)
                        processed_peaks.add(other_peak_idx)
                        peaks_to_check.append((other_rep_idx + 1, other_peak_idx, other_peak))
        
        # Count unique replicates represented
        unique_replicates = len(set(p['name'].split('_')[0] for p in all_overlapping_peaks))
        
        # If peak is found in enough replicates, add it to results
        if unique_replicates >= min_replicates:
            # Merge overlapping peaks
            merged_peak = {
                'chr': peak['chr'],
                'start': min(p['start'] for p in all_overlapping_peaks),
                'end': max(p['end'] for p in all_overlapping_peaks),
                'name': f"Peak_{len(combined_peaks)+1}",
                'score': max(p['score'] for p in all_overlapping_peaks),
                'strand': '.',
                'signalValue': max(p['signalValue'] for p in all_overlapping_peaks),
                'pValue': min(p['pValue'] for p in all_overlapping_peaks),
                'qValue': min(p['qValue'] for p in all_overlapping_peaks),
                'peak': 0
            }
            combined_peaks.append(merged_peak)
    
    return pd.DataFrame(combined_peaks)

def main():
    parser = argparse.ArgumentParser(description='Combine replicate peaks.')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing peak files')
    parser.add_argument('--condition', required=True, choices=['endo', 'exo'], 
                       help='Condition to process (endo or exo)')
    parser.add_argument('--min-replicates', type=int, default=2,
                       help='Minimum number of replicates required (default: 2)')
    args = parser.parse_args()
    
    if args.condition == 'exo':
        pattern = 'NSCv'
    else:
        pattern = 'NSCM'
    
    # Find all replicate files
    peak_files = [
        os.path.join(args.peaks_dir, f) for f in os.listdir(args.peaks_dir)
        if f.startswith(pattern) and f.endswith('filtered.broadPeak')
    ]
    
    if len(peak_files) < 3:
        print(f"Error: Found only {len(peak_files)} replicate files")
        return
    
    # Combine peaks with more lenient parameters
    combined_peaks = combine_replicates(
        peak_files, 
        min_replicates=args.min_replicates,
        overlap_threshold=0.25,  # More lenient overlap requirement
        max_gap=500  # Larger gap allowed between peaks
    )
    
    # Save combined peaks
    output_file = os.path.join(args.peaks_dir, f"NPCs_{args.condition}_combined.broadPeak")
    combined_peaks.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Combined peaks saved to: {output_file}")
    print(f"Total peaks in combined file: {len(combined_peaks)}")

if __name__ == "__main__":
    main() 