import pandas as pd
import numpy as np
from typing import List, Tuple
import os
import argparse
from intervaltree import IntervalTree, Interval

def load_peak_file(file_path: str) -> pd.DataFrame:
    """Load narrowPeak format file."""
    columns = ['chr', 'start', 'end', 'name', 'score', 
              'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    return pd.read_csv(file_path, sep='\t', names=columns)

def create_interval_tree(peaks_df: pd.DataFrame) -> dict:
    """Create interval trees for each chromosome from peaks."""
    trees = {}
    for _, peak in peaks_df.iterrows():
        chrom = peak['chr']
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom].add(Interval(peak['start'], peak['end'], peak))
    return trees

def find_overlapping_peaks(trees: dict, peak: pd.Series, overlap_threshold: float = 0.25) -> Interval:
    """Find overlapping peaks in interval trees."""
    chrom = peak['chr']
    if chrom not in trees:
        return None
    
    overlaps = trees[chrom].overlap(peak['start'], peak['end'])
    for interval in overlaps:
        overlap_length = min(peak['end'], interval.end) - max(peak['start'], interval.begin)
        peak_length = peak['end'] - peak['start']
        interval_length = interval.end - interval.begin
        
        if (overlap_length / peak_length >= overlap_threshold and 
            overlap_length / interval_length >= overlap_threshold):
            return interval
    return None

def combine_replicates(peak_files: List[str], min_replicates: int = 2) -> pd.DataFrame:
    """
    Combine peaks from replicates, requiring peaks to be present in at least min_replicates.
    Optimized for low-depth Cut&Tag data.
    """
    # Load all peak files
    peak_dfs = [load_peak_file(f) for f in peak_files]
    
    # Create interval trees for all replicates except the first
    trees = [create_interval_tree(df) for df in peak_dfs[1:]]
    
    combined_peaks = []
    
    # Start with peaks from first replicate
    for _, peak in peak_dfs[0].iterrows():
        # Count how many other replicates have this peak
        overlap_count = 1  # Count itself
        overlapping_peaks = [peak]
        
        for tree in trees:
            overlap = find_overlapping_peaks(tree, peak, overlap_threshold=0.25)  # More lenient overlap
            if overlap:
                overlap_count += 1
                overlapping_peaks.append(overlap.data)
        
        # If peak is found in enough replicates, add it to results
        if overlap_count >= min_replicates:
            # For low-depth data, take the strongest peak's properties
            strongest_peak = max(overlapping_peaks, key=lambda x: x['signalValue'])
            
            merged_peak = {
                'chr': peak['chr'],
                'start': min(p['start'] for p in overlapping_peaks),
                'end': max(p['end'] for p in overlapping_peaks),
                'name': f"Peak_{len(combined_peaks)+1}",
                'score': strongest_peak['score'],  # Take strongest peak's score
                'strand': '.',
                'signalValue': strongest_peak['signalValue'],  # Take strongest peak's signal
                'pValue': strongest_peak['pValue'],  # Take strongest peak's p-value
                'qValue': strongest_peak['qValue'],  # Take strongest peak's q-value
                'peak': strongest_peak['peak']  # Take strongest peak's summit
            }
            combined_peaks.append(merged_peak)
    
    return pd.DataFrame(combined_peaks)

def main():
    parser = argparse.ArgumentParser(description='Combine replicate peaks.')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing peak files')
    parser.add_argument('--condition', required=True, choices=['endo', 'exo'], 
                       help='Condition to process (endo or exo)')
    args = parser.parse_args()
    
    if args.condition == 'exo':
        pattern = 'NSCv'
    else:
        pattern = 'NSCM'
    
    # Find all replicate files
    peak_files = [
        os.path.join(args.peaks_dir, f) for f in os.listdir(args.peaks_dir)
        if f.startswith(pattern) and f.endswith('filtered.narrowPeak')
    ]
    
    if len(peak_files) < 3:
        print(f"Error: Found only {len(peak_files)} replicate files")
        return
    
    # Combine peaks
    combined_peaks = combine_replicates(peak_files, min_replicates=2)
    
    # Save combined peaks
    output_file = os.path.join(args.peaks_dir, f"NPCs_{args.condition}_combined.narrowPeak")
    combined_peaks.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Combined peaks saved to: {output_file}")
    print(f"Total peaks in combined file: {len(combined_peaks)}")

if __name__ == "__main__":
    main() 