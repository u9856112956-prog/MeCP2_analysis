#!/usr/bin/env python3
"""
Analysis of chromatin states for specific genes in NSC and Neu samples (GFP condition)
This script visualizes the chromatin states of genes from a provided list
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import logging
import subprocess
import tempfile

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("gene_chromatin_state_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('Gene-Chromatin-Analysis')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Gene Chromatin State Analysis')
    parser.add_argument('--gene-list', required=True, help='Path to gene list file')
    parser.add_argument('--results-dir', required=True, help='Directory containing chromatin state results')
    parser.add_argument('--output-dir', required=True, help='Output directory for analysis results')
    parser.add_argument('--genome-gtf', required=True,
                        help='Path to genome GTF file for gene coordinates (e.g., Mus_musculus.GRCm38.102.gtf)')
    parser.add_argument('--promoter-size', type=int, default=2000, 
                        help='Size of promoter region in base pairs upstream of transcription start site (default: 2000)')
    
    args = parser.parse_args()
    
    # Ensure output_dir is an absolute path or relative to current directory
    if not os.path.isabs(args.output_dir):
        args.output_dir = os.path.join(os.getcwd(), args.output_dir)
    
    return args

def load_gene_list(file_path):
    """Load gene list from file"""
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        logger.error(f"Error loading gene list: {e}")
        sys.exit(1)

def load_chromatin_states(results_dir, condition, state_type):
    """Load chromatin state data for a specific condition and state type"""
    # Try multiple possible paths for the chromatin state files
    path = os.path.join(results_dir, "chromatin_states", f"{condition}_{state_type}.bed")
    
    if not os.path.exists(path):
        logger.error(f"Chromatin state file not found for {condition}_{state_type}.")
        return pd.DataFrame()
    
    try:
        # The 4th column is the score, not a state name
        states = pd.read_csv(path, sep='\t',
                            names=['chromosome', 'start', 'end', 'score'])
        # Add a state type column
        states['state_type'] = state_type
        states['condition'] = condition
        logger.info(f"Loaded {len(states)} {state_type} regions for {condition} from {path}")
        return states
    except Exception as e:
        logger.error(f"Error loading chromatin states from {path}: {e}")
        return pd.DataFrame()

def get_gene_promoter_coordinates(genes, gtf_file, promoter_size=2000):
    """Extract promoter coordinates for genes from GTF file
    
    Promoter is defined as a region upstream of the transcription start site (TSS)
    Default promoter size is 2kb upstream of TSS
    """
    logger.info(f"Extracting promoter coordinates for {len(genes)} genes from GTF file")
    
    # Check if GTF file exists
    if not os.path.exists(gtf_file):
        logger.error(f"GTF file not found: {gtf_file}")
        return {}
    
    # Use grep to extract all gene entries from GTF file once
    all_genes_cmd = f"grep -w 'gene' {gtf_file} > /tmp/all_genes.gtf"
    try:
        subprocess.run(all_genes_cmd, shell=True, check=True)
        logger.info(f"Extracted all gene entries from GTF file")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error extracting genes from GTF: {e}")
        return {}
    
    gene_promoters = {}
    found_count = 0
    
    # Process each gene
    for gene in genes:
        # Try different patterns to match gene names
        patterns = [
            f"grep 'gene_name \"{gene}\"' /tmp/all_genes.gtf",
            f"grep 'gene_name \"{gene};\"' /tmp/all_genes.gtf",
            f"grep -w 'gene_name' /tmp/all_genes.gtf | grep -w '{gene}'"
        ]
        
        for cmd in patterns:
            try:
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                if result.stdout.strip():
                    # Parse the GTF line
                    fields = result.stdout.strip().split('\t')
                    if len(fields) >= 8:  # Need at least 8 fields to get strand information
                        chromosome = fields[0]
                        feature_type = fields[2]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]  # '+' or '-'
                        
                        # Define promoter region based on strand
                        if strand == '+':
                            # For genes on + strand, promoter is upstream of start
                            promoter_start = max(1, start - promoter_size)  # Ensure start is not negative
                            promoter_end = start
                        else:  # strand == '-'
                            # For genes on - strand, promoter is upstream of end
                            promoter_start = end
                            promoter_end = end + promoter_size
                        
                        gene_promoters[gene] = {
                            'chromosome': chromosome,
                            'start': promoter_start,
                            'end': promoter_end,
                            'strand': strand,
                            'gene_start': start,
                            'gene_end': end
                        }
                        found_count += 1
                        break  # Found the gene, no need to try other patterns
            except subprocess.CalledProcessError:
                continue  # Try next pattern
    
    # Clean up temporary file
    subprocess.run("rm -f /tmp/all_genes.gtf", shell=True)
    
    logger.info(f"Found promoter coordinates for {found_count} out of {len(genes)} genes")
    return gene_promoters

def assign_chromatin_states(gene_coords, eu_states, hetero_states):
    """Assign chromatin states to genes based on overlap with chromatin regions
    
    Returns both categorical state and numerical values for signal strength
    """
    gene_states = {}
    
    for gene, coords in gene_coords.items():
        chrom = coords['chromosome']
        gene_start = coords['start']
        gene_end = coords['end']
        
        # Find overlapping euchromatin regions
        eu_overlaps = eu_states[
            (eu_states['chromosome'] == chrom) & 
            (eu_states['end'] >= gene_start) & 
            (eu_states['start'] <= gene_end)
        ]
        
        # Find overlapping heterochromatin regions
        hetero_overlaps = hetero_states[
            (hetero_states['chromosome'] == chrom) & 
            (hetero_states['end'] >= gene_start) & 
            (hetero_states['start'] <= gene_end)
        ]
        
        # Initialize numerical values
        eu_overlap_len = 0
        hetero_overlap_len = 0
        eu_signal = 0
        hetero_signal = 0
        eu_count = len(eu_overlaps)
        hetero_count = len(hetero_overlaps)
        
        # Calculate overlap lengths and signal strengths
        if not eu_overlaps.empty:
            for _, region in eu_overlaps.iterrows():
                overlap_start = max(region['start'], gene_start)
                overlap_end = min(region['end'], gene_end)
                overlap_len = overlap_end - overlap_start
                eu_overlap_len += overlap_len
                # Add signal strength (score column)
                if 'score' in region:
                    eu_signal += float(region['score']) * overlap_len  # Weight by overlap length
            
            # Calculate average signal
            if eu_overlap_len > 0:
                eu_signal = eu_signal / eu_overlap_len
        
        if not hetero_overlaps.empty:
            for _, region in hetero_overlaps.iterrows():
                overlap_start = max(region['start'], gene_start)
                overlap_end = min(region['end'], gene_end)
                overlap_len = overlap_end - overlap_start
                hetero_overlap_len += overlap_len
                # Add signal strength (score column)
                if 'score' in region:
                    hetero_signal += float(region['score']) * overlap_len  # Weight by overlap length
            
            # Calculate average signal
            if hetero_overlap_len > 0:
                hetero_signal = hetero_signal / hetero_overlap_len
        
        # Store numerical data
        state_data = {
            'eu_overlap_length': eu_overlap_len,
            'hetero_overlap_length': hetero_overlap_len,
            'eu_signal': eu_signal,
            'hetero_signal': hetero_signal,
            'eu_count': eu_count,
            'hetero_count': hetero_count,
            'total_promoter_length': gene_end - gene_start
        }
        
        # Determine the predominant state
        if not eu_overlaps.empty and not hetero_overlaps.empty:
            # Assign state based on which has more overlap
            if eu_overlap_len > hetero_overlap_len:
                state = "Euchromatin"
            else:
                state = "Heterochromatin"
        elif not eu_overlaps.empty:
            state = "Euchromatin"
        elif not hetero_overlaps.empty:
            state = "Heterochromatin"
        else:
            state = "Unknown"
        
        # Store both categorical state and numerical data
        state_data['state'] = state
        gene_states[gene] = state_data
    
    return gene_states

def analyze_gene_chromatin_states_gfp_vs_m2(genes, results_dir, gtf_file, promoter_size=2000):
    """Analyze chromatin states for gene promoters comparing GFP vs M2 conditions
    
    This function focuses on the promoter regions of genes rather than the entire gene body.
    Promoters are defined as regions upstream of the transcription start site (TSS).
    """
    # Get gene promoter coordinates
    gene_promoters = get_gene_promoter_coordinates(genes, gtf_file, promoter_size)
    
    # Load chromatin states for NSC_GFP
    nsc_gfp_eu = load_chromatin_states(results_dir, "NSC_GFP", "euchromatin")
    nsc_gfp_hetero = load_chromatin_states(results_dir, "NSC_GFP", "heterochromatin")
    
    # Load chromatin states for NSC_M2
    nsc_m2_eu = load_chromatin_states(results_dir, "NSC_M2", "euchromatin")
    nsc_m2_hetero = load_chromatin_states(results_dir, "NSC_M2", "heterochromatin")
    
    # Assign chromatin states for gene promoters
    nsc_gfp_states = assign_chromatin_states(gene_promoters, nsc_gfp_eu, nsc_gfp_hetero)
    nsc_m2_states = assign_chromatin_states(gene_promoters, nsc_m2_eu, nsc_m2_hetero)
    
    # Combine results
    combined_states = {}
    for gene in genes:
        if gene in gene_promoters:
            # Get state data or create default if not found
            gfp_data = nsc_gfp_states.get(gene, {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0})
            m2_data = nsc_m2_states.get(gene, {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0})
            
            # Create combined entry with both categorical and numerical data
            combined_states[gene] = {
                "GFP": gfp_data,
                "M2": m2_data
            }
        else:
            # Default data for genes not found in GTF
            default_data = {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0}
            combined_states[gene] = {
                "GFP": default_data,
                "M2": default_data
            }
    
    logger.info(f"Completed GFP vs M2 chromatin state analysis for promoters of {len(combined_states)} genes")
    return combined_states

def analyze_gene_chromatin_states_nsc_vs_neu(genes, results_dir, gtf_file, promoter_size=2000):
    """Analyze chromatin states for gene promoters comparing NSC vs Neu (GFP condition)
    
    This function focuses on the promoter regions of genes rather than the entire gene body.
    Promoters are defined as regions upstream of the transcription start site (TSS).
    """
    # Get gene promoter coordinates
    gene_promoters = get_gene_promoter_coordinates(genes, gtf_file, promoter_size)
    
    # Load chromatin states for NSC_GFP
    nsc_eu = load_chromatin_states(results_dir, "NSC_GFP", "euchromatin")
    nsc_hetero = load_chromatin_states(results_dir, "NSC_GFP", "heterochromatin")
    
    # Load chromatin states for Neu_GFP
    neu_eu = load_chromatin_states(results_dir, "Neu_GFP", "euchromatin")
    neu_hetero = load_chromatin_states(results_dir, "Neu_GFP", "heterochromatin")
    
    # Assign chromatin states for gene promoters
    nsc_states = assign_chromatin_states(gene_promoters, nsc_eu, nsc_hetero)
    neu_states = assign_chromatin_states(gene_promoters, neu_eu, neu_hetero)
    
    # Combine results
    combined_states = {}
    for gene in genes:
        if gene in gene_promoters:
            # Get state data or create default if not found
            nsc_data = nsc_states.get(gene, {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0})
            neu_data = neu_states.get(gene, {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0})
            
            # Create combined entry with both categorical and numerical data
            combined_states[gene] = {
                "NSC": nsc_data,
                "Neu": neu_data
            }
        else:
            # Default data for genes not found in GTF
            default_data = {"state": "Unknown", "eu_signal": 0, "hetero_signal": 0}
            combined_states[gene] = {
                "NSC": default_data,
                "Neu": default_data
            }
    
    logger.info(f"Completed NSC vs Neu chromatin state analysis for promoters of {len(combined_states)} genes")
    return combined_states

def plot_chromatin_states(gene_states, output_dir, gene_list_name):
    """Create visualization of chromatin states for the genes"""
    # Prepare data for plotting
    data = []
    
    # Determine which comparison we're dealing with based on the keys in the first gene's states
    first_gene = next(iter(gene_states))
    comparison_type = "unknown"
    if 'NSC' in gene_states[first_gene] and 'Neu' in gene_states[first_gene]:
        comparison_type = "nsc_vs_neu"
        col1, col2 = 'NSC', 'Neu'
        title1, title2 = 'NSC (GFP)', 'Neu (GFP)'
        transition_title = 'NSC to Neu'
    elif 'GFP' in gene_states[first_gene] and 'M2' in gene_states[first_gene]:
        comparison_type = "gfp_vs_m2"
        col1, col2 = 'GFP', 'M2'
        title1, title2 = 'NSC (GFP)', 'NSC (M2)'
        transition_title = 'GFP to M2'
    else:
        logger.error("Unknown comparison type in gene states")
        return
    
    logger.info(f"Plotting {comparison_type} comparison")
    
    # Format data for DataFrame
    for gene, states in gene_states.items():
        # Extract categorical states
        state1 = states[col1]['state'] if isinstance(states[col1], dict) else 'Unknown'
        state2 = states[col2]['state'] if isinstance(states[col2], dict) else 'Unknown'
        
        # Extract numerical values
        eu_signal1 = states[col1].get('eu_signal', 0) if isinstance(states[col1], dict) else 0
        hetero_signal1 = states[col1].get('hetero_signal', 0) if isinstance(states[col1], dict) else 0
        eu_signal2 = states[col2].get('eu_signal', 0) if isinstance(states[col2], dict) else 0
        hetero_signal2 = states[col2].get('hetero_signal', 0) if isinstance(states[col2], dict) else 0
        
        data.append({
            'Gene': gene,
            # Categorical states
            f'{col1}_state': state1,
            f'{col2}_state': state2,
            # Numerical values
            f'{col1}_eu_signal': eu_signal1,
            f'{col1}_hetero_signal': hetero_signal1,
            f'{col2}_eu_signal': eu_signal2,
            f'{col2}_hetero_signal': hetero_signal2,
            # Signal ratios
            f'{col1}_eu_hetero_ratio': eu_signal1 / hetero_signal1 if hetero_signal1 > 0 else 0,
            f'{col2}_eu_hetero_ratio': eu_signal2 / hetero_signal2 if hetero_signal2 > 0 else 0
        })
    df = pd.DataFrame(data)
    
    # Save the data to CSV
    df.to_csv(os.path.join(output_dir, 'gene_chromatin_states.csv'), index=False)
    logger.info(f"Saved gene chromatin states to CSV")
    
    # Create state distribution plot
    plt.figure(figsize=(14, 7))
    
    # Define consistent color mapping
    color_map = {
        'Unknown': '#66b3ff',      # Blue for Unknown
        'Euchromatin': '#99ff99',  # Green for Euchromatin
        'Heterochromatin': '#ff9999', # Red for Heterochromatin
        'Mixed': '#ffcc99'         # Orange for Mixed
    }
    
    # Plot for first condition
    col1_counts = df[f'{col1}_state'].value_counts()
    plt.subplot(1, 2, 1)
    
    # Map colors to the actual states present in the data
    col1_colors = [color_map[state] for state in col1_counts.index]
    plt.pie(col1_counts.values, labels=col1_counts.index, autopct='%1.1f%%', colors=col1_colors)
    plt.title(f'Chromatin States Distribution in {title1}')
    
    # Plot for second condition
    col2_counts = df[f'{col2}_state'].value_counts()
    plt.subplot(1, 2, 2)
    
    # Map colors to the actual states present in the data
    col2_colors = [color_map[state] for state in col2_counts.index]
    plt.pie(col2_counts.values, labels=col2_counts.index, autopct='%1.1f%%', colors=col2_colors)
    plt.title(f'Chromatin States Distribution in {title2}')
    
    plt.suptitle(f'Chromatin State Distribution for {gene_list_name} Genes', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_states_distribution.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin states distribution plot")
    
    # Create state transition plot
    transitions = {}
    for gene, states in gene_states.items():
        state1 = states[col1]['state'] if isinstance(states[col1], dict) else 'Unknown'
        state2 = states[col2]['state'] if isinstance(states[col2], dict) else 'Unknown'
        transition = f"{state1} → {state2}"
        transitions[transition] = transitions.get(transition, 0) + 1
    
    # Sort transitions to group them meaningfully
    sorted_transitions = {}
    # First add transitions that maintain the same state
    for state in ['Euchromatin', 'Heterochromatin', 'Mixed', 'Unknown']:
        key = f"{state} → {state}"
        if key in transitions:
            sorted_transitions[key] = transitions[key]
    
    # Then add transitions from open to closed chromatin
    key = "Euchromatin → Heterochromatin"
    if key in transitions:
        sorted_transitions[key] = transitions[key]
    
    # Then add transitions from closed to open chromatin
    key = "Heterochromatin → Euchromatin"
    if key in transitions:
        sorted_transitions[key] = transitions[key]
    
    # Add all other transitions
    for key, value in transitions.items():
        if key not in sorted_transitions:
            sorted_transitions[key] = value
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(sorted_transitions.keys(), sorted_transitions.values(), color='skyblue')
    plt.xlabel(f'Chromatin State Transition ({title1} → {title2})')
    plt.ylabel('Number of Genes')
    plt.title(f'Chromatin State Transitions for {gene_list_name} Genes')
    plt.xticks(rotation=45, ha='right')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                 f'{height}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_state_transitions.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin state transitions plot")
    
    # Create a focused plot on chromatin state changes (open to closed and vice versa)
    plt.figure(figsize=(10, 6))
    
    # Define the key transitions we're interested in
    key_transitions = {
        'Euchromatin → Heterochromatin': 'Open to Closed',
        'Heterochromatin → Euchromatin': 'Closed to Open',
        'Euchromatin → Euchromatin': 'Remained Open',
        'Heterochromatin → Heterochromatin': 'Remained Closed'
    }
    
    # Define consistent colors for transitions
    transition_color_map = {
        'Open to Closed': '#ff9999',      # Red for closing (like heterochromatin)
        'Closed to Open': '#99ff99',      # Green for opening (like euchromatin)
        'Remained Open': '#99cc99',       # Darker green for remained open
        'Remained Closed': '#cc9999'      # Darker red for remained closed
    }
    
    # Extract counts for these transitions
    transition_counts = []
    transition_labels = []
    transition_colors = []
    
    for key, label in key_transitions.items():
        if key in transitions:
            transition_counts.append(transitions[key])
            transition_labels.append(label)
            transition_colors.append(transition_color_map[label])
        else:
            transition_counts.append(0)
            transition_labels.append(label)
            transition_colors.append(transition_color_map[label])
    
    # Create the bar plot
    bars = plt.bar(transition_labels, transition_counts, color=transition_colors)
    plt.xlabel('Chromatin State Change')
    plt.ylabel('Number of Genes')
    plt.title(f'Key Chromatin State Changes from {title1} to {title2}')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:  # Only add text if the bar has height
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{height}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'key_chromatin_state_changes.png'), dpi=300)
    plt.close()
    logger.info(f"Created key chromatin state changes plot")
    
    # Create a pie chart showing the proportion of genes that changed state vs. remained the same
    changed_state = 0
    same_state = 0
    unknown_state = 0
    
    for gene, states in gene_states.items():
        state1 = states[col1]['state'] if isinstance(states[col1], dict) else 'Unknown'
        state2 = states[col2]['state'] if isinstance(states[col2], dict) else 'Unknown'
        
        if state1 == 'Unknown' or state2 == 'Unknown':
            unknown_state += 1
        elif state1 == state2:
            same_state += 1
        else:
            changed_state += 1
    
    # Use consistent colors - blue for unknown, red for changed (like heterochromatin), green for maintained (like euchromatin)
    change_color_map = {
        'Changed State': '#ff9999',    # Red for changes
        'Maintained State': '#99ff99', # Green for maintained
        'Unknown': '#66b3ff'           # Blue for unknown
    }
    
    plt.figure(figsize=(8, 8))
    plt.pie([changed_state, same_state, unknown_state], 
            labels=['Changed State', 'Maintained State', 'Unknown'],
            autopct='%1.1f%%', colors=[change_color_map['Changed State'], 
                                      change_color_map['Maintained State'], 
                                      change_color_map['Unknown']])
    plt.title(f'Proportion of Genes Changing Chromatin State from {title1} to {title2}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_state_change_proportion.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin state change proportion plot")
    
    # Create numerical signal plots
    logger.info(f"Creating numerical signal plots")
    
    # 1. Boxplot of euchromatin and heterochromatin signal strengths
    plt.figure(figsize=(14, 7))
    
    # First condition
    plt.subplot(1, 2, 1)
    boxplot_data1 = [
        df[f'{col1}_eu_signal'].dropna(),
        df[f'{col1}_hetero_signal'].dropna()
    ]
    box1 = plt.boxplot(boxplot_data1, patch_artist=True, labels=['Euchromatin', 'Heterochromatin'])
    for patch, color in zip(box1['boxes'], ['#99ff99', '#ff9999']):
        patch.set_facecolor(color)
    plt.title(f'Chromatin Signal Strength in {title1}')
    plt.ylabel('Signal Strength')
    
    # Second condition
    plt.subplot(1, 2, 2)
    boxplot_data2 = [
        df[f'{col2}_eu_signal'].dropna(),
        df[f'{col2}_hetero_signal'].dropna()
    ]
    box2 = plt.boxplot(boxplot_data2, patch_artist=True, labels=['Euchromatin', 'Heterochromatin'])
    for patch, color in zip(box2['boxes'], ['#99ff99', '#ff9999']):
        patch.set_facecolor(color)
    plt.title(f'Chromatin Signal Strength in {title2}')
    plt.ylabel('Signal Strength')
    
    plt.suptitle(f'Chromatin Signal Strength Distribution for {gene_list_name} Genes', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_signal_strength.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin signal strength boxplot")
    
    # 2. Scatter plot comparing signal strengths between conditions
    plt.figure(figsize=(10, 8))
    
    # Euchromatin signal comparison
    plt.subplot(1, 2, 1)
    plt.scatter(df[f'{col1}_eu_signal'], df[f'{col2}_eu_signal'], alpha=0.6, color='#99ff99')
    max_val = max(df[f'{col1}_eu_signal'].max(), df[f'{col2}_eu_signal'].max()) * 1.1
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)  # Diagonal line
    plt.xlabel(f'{title1} Euchromatin Signal')
    plt.ylabel(f'{title2} Euchromatin Signal')
    plt.title('Euchromatin Signal Comparison')
    
    # Heterochromatin signal comparison
    plt.subplot(1, 2, 2)
    plt.scatter(df[f'{col1}_hetero_signal'], df[f'{col2}_hetero_signal'], alpha=0.6, color='#ff9999')
    max_val = max(df[f'{col1}_hetero_signal'].max(), df[f'{col2}_hetero_signal'].max()) * 1.1
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)  # Diagonal line
    plt.xlabel(f'{title1} Heterochromatin Signal')
    plt.ylabel(f'{title2} Heterochromatin Signal')
    plt.title('Heterochromatin Signal Comparison')
    
    plt.suptitle(f'Chromatin Signal Comparison Between {title1} and {title2}', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'chromatin_signal_comparison.png'), dpi=300)
    plt.close()
    logger.info(f"Created chromatin signal comparison scatter plot")
    
    # 3. Signal ratio plot (Euchromatin/Heterochromatin ratio)
    plt.figure(figsize=(14, 7))
    
    # Filter out zeros and extreme values for better visualization
    ratio1 = df[f'{col1}_eu_hetero_ratio']
    ratio1 = ratio1[(ratio1 > 0) & (ratio1 < ratio1.quantile(0.95))]
    
    ratio2 = df[f'{col2}_eu_hetero_ratio']
    ratio2 = ratio2[(ratio2 > 0) & (ratio2 < ratio2.quantile(0.95))]
    
    plt.subplot(1, 2, 1)
    plt.hist(ratio1, bins=20, alpha=0.7, color='#99ccff')
    plt.axvline(x=1, color='k', linestyle='--', alpha=0.5)  # Line at ratio=1 (equal eu/hetero)
    plt.xlabel('Euchromatin/Heterochromatin Ratio')
    plt.ylabel('Number of Genes')
    plt.title(f'Signal Ratio in {title1}')
    
    plt.subplot(1, 2, 2)
    plt.hist(ratio2, bins=20, alpha=0.7, color='#99ccff')
    plt.axvline(x=1, color='k', linestyle='--', alpha=0.5)  # Line at ratio=1 (equal eu/hetero)
    plt.xlabel('Euchromatin/Heterochromatin Ratio')
    plt.ylabel('Number of Genes')
    plt.title(f'Signal Ratio in {title2}')
    
    plt.suptitle(f'Euchromatin/Heterochromatin Signal Ratio for {gene_list_name} Genes', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'eu_hetero_signal_ratio.png'), dpi=300)
    plt.close()
    logger.info(f"Created euchromatin/heterochromatin signal ratio histogram")

def main():
    # Parse arguments
    args = parse_arguments()
    
    # Extract gene list name from file path
    gene_list_name = os.path.basename(args.gene_list).split('.')[0]
    
    # Create output directory if it doesn't exist
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Analysis output will be saved to {args.output_dir}")
    except PermissionError:
        # If permission denied, try to use a directory in the user's home or current directory
        user_home = os.path.expanduser("~")
        fallback_dir = os.path.join(user_home, "chromatin_analysis_output")
        try:
            os.makedirs(fallback_dir, exist_ok=True)
            args.output_dir = fallback_dir
            logger.warning(f"Permission denied for original output directory. Using {fallback_dir} instead.")
        except PermissionError:
            # If still permission denied, use current directory
            fallback_dir = os.path.join(os.getcwd(), "chromatin_analysis_output")
            os.makedirs(fallback_dir, exist_ok=True)
            args.output_dir = fallback_dir
            logger.warning(f"Permission denied for home directory. Using {fallback_dir} instead.")
    
    logger.info(f"Analyzing promoter regions with size: {args.promoter_size} bp")
    
    # Load gene list
    genes = load_gene_list(args.gene_list)
    logger.info(f"Loaded {len(genes)} genes from {args.gene_list}")
    
    # Create subdirectories for each comparison
    gfp_vs_m2_dir = os.path.join(args.output_dir, 'gfp_vs_m2')
    nsc_vs_neu_dir = os.path.join(args.output_dir, 'nsc_vs_neu')
    os.makedirs(gfp_vs_m2_dir, exist_ok=True)
    os.makedirs(nsc_vs_neu_dir, exist_ok=True)
    
    # Analyze GFP vs M2 comparison (in NSCs)
    logger.info("Analyzing GFP vs M2 comparison in NSCs...")
    gfp_vs_m2_states = analyze_gene_chromatin_states_gfp_vs_m2(genes, args.results_dir, args.genome_gtf, args.promoter_size)
    plot_chromatin_states(gfp_vs_m2_states, gfp_vs_m2_dir, f"{gene_list_name}_GFP_vs_M2")
    
    # Analyze NSC vs Neu comparison (in GFP condition)
    logger.info("Analyzing NSC vs Neu comparison in GFP condition...")
    nsc_vs_neu_states = analyze_gene_chromatin_states_nsc_vs_neu(genes, args.results_dir, args.genome_gtf, args.promoter_size)
    plot_chromatin_states(nsc_vs_neu_states, nsc_vs_neu_dir, f"{gene_list_name}_NSC_vs_Neu")
    
    logger.info(f"Analysis complete. Results saved in {args.output_dir}")

if __name__ == "__main__":
    main()
