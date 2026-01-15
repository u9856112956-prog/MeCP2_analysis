import os
import pandas as pd
import numpy as np
from typing import Dict, List
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def load_gene_lists(base_dir: str) -> Dict[str, pd.DataFrame]:
    """Load the three gene lists."""
    categories = ['up', 'down', 'no_deg']
    gene_lists = {}
    for category in categories:
        file_path = os.path.join(base_dir, "data", f'extended_{category}.csv')
        df = pd.read_csv(file_path)
        
        # Convert coordinates to integers and handle NaN values
        if 'start' in df.columns and 'end' in df.columns:
            df['start'] = pd.to_numeric(df['start'], errors='coerce').fillna(0).astype(int)
            df['end'] = pd.to_numeric(df['end'], errors='coerce').fillna(0).astype(int)
        
        # Remove rows with invalid coordinates
        df = df[
            (df['seqnames'].notna()) & 
            (df['start'] > 0) & 
            (df['end'] > 0)
        ]
        
        # Ensure chromosome names are strings and properly formatted
        df['seqnames'] = df['seqnames'].astype(str)
        # Add 'chr' prefix if not present
        df['seqnames'] = df['seqnames'].apply(lambda x: f"chr{x}" if not x.startswith('chr') else x)
        
        gene_lists[category] = df
        print(f"Loaded {len(df)} valid regions for {category}")
    return gene_lists

def get_smarcb1_signal(bw_file: str, region: pd.Series, window_size: int = 500) -> float:
    """Get average SMARCB1 signal for a region."""
    with pyBigWig.open(bw_file) as bw:
        try:
            chrom = str(region['seqnames'])
            start = int(region['start'])
            end = int(region['end'])
            center = (start + end) // 2
            window_start = max(0, center - window_size)
            window_end = center + window_size
            
            # Verify coordinates are valid
            if window_start >= window_end:
                print(f"Invalid coordinates for region {chrom}:{window_start}-{window_end}")
                return 0.0
                
            values = bw.values(chrom, window_start, window_end)
            if values is None:
                return 0.0
                
            values = [v for v in values if v is not None and not np.isnan(v)]
            return np.mean(values) if values else 0.0
            
        except Exception as e:
            print(f"Error processing region {chrom}:{start}-{end}: {str(e)}")
            return 0.0

def compare_smarcb1_levels(gene_lists: Dict[str, pd.DataFrame],
                          smarcb1_dir: str,
                          output_dir: str):
    """Compare SMARCB1 levels between BM3 and BGs for each gene list."""
    
    # Setup files
    bm_file = os.path.join(smarcb1_dir, "BM3_RPKM.bw")
    bg_files = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]
    
    # Verify files exist
    all_files = [bm_file] + bg_files
    for f in all_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"BigWig file not found: {f}")
    
    results = []
    
    for category, df in gene_lists.items():
        print(f"Processing {category} genes ({len(df)} regions)...")
        
        for idx, region in df.iterrows():
            if idx % 100 == 0:
                print(f"  Processing gene {idx}/{len(df)}")
            
            # Skip invalid regions
            if pd.isna(region['start']) or pd.isna(region['end']) or pd.isna(region['seqnames']):
                continue
                
            # Get signals
            bm_signal = get_smarcb1_signal(bm_file, region)
            bg_signals = [get_smarcb1_signal(bg_file, region) for bg_file in bg_files]
            
            # Calculate metrics
            avg_bg = np.mean(bg_signals) if bg_signals else 0
            is_higher = bm_signal > max(bg_signals) if bg_signals else False
            
            results.append({
                'category': category,
                'gene_id': region.get('gene_id', 'Unknown'),
                'gene_name': region.get('gene_name', 'Unknown'),
                'chromosome': region['seqnames'],
                'start': region['start'],
                'end': region['end'],
                'bm3_signal': bm_signal,
                'avg_bg_signal': avg_bg,
                'fold_change': bm_signal / avg_bg if avg_bg > 0 else 0,
                'is_higher_than_all_bg': is_higher
            })
    
    results_df = pd.DataFrame(results)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save results
    results_path = os.path.join(output_dir, 'smarcb1_comparison_results.csv')
    results_df.to_csv(results_path, index=False)
    print(f"Saved results to {results_path}")
    
    # Plot results only if we have valid data
    if results_df['fold_change'].max() > 0:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=results_df, x='category', y='fold_change')
        plt.title('SMARCB1 Signal Fold Change (BM3/BG) by Category')
        plt.yscale('log')
        plot_path = os.path.join(output_dir, 'smarcb1_comparison_boxplot.png')
        plt.savefig(plot_path)
        plt.close()
        print(f"Saved plot to {plot_path}")
    else:
        print("Warning: No valid fold changes found for plotting")

def main():
    # Set paths - configure these for your data
    base_dir = ""  # e.g., "/path/to/cutandtag/analysis"
    smarcb1_dir = ""  # e.g., "/path/to/smarcb1/bigwig"
    output_dir = os.path.join(base_dir, "results/smarcb1_comparison")
    
    # Load gene lists
    gene_lists = load_gene_lists(base_dir)
    
    # Run analysis
    compare_smarcb1_levels(gene_lists, smarcb1_dir, output_dir)

if __name__ == "__main__":
    main()