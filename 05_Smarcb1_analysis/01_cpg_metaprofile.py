#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg') # Use non-interactive backend BEFORE pyplot import

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import argparse
import logging
import time
from scipy.ndimage import gaussian_filter1d # Import Gaussian filter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# --- Default Parameters --- 
# DEFAULT_WINDOW_SIZE = 2000 # Default half-window size
DEFAULT_WINDOW_SIZE = 5000 # Default half-window size (+/- 5kb)
DEFAULT_BIN_SIZE = 50

# Global variable for worker process
bw_reader = None

def init_worker(bigwig_path):
    """Initializer function for each worker process."""
    global bw_reader
    try:
        bw_reader = pyBigWig.open(bigwig_path)
    except Exception as e:
        bw_reader = None

def extract_signal_for_region(args):
    """Worker function to extract signal for a single region."""
    # bigwig_path is no longer needed here, worker uses the global bw_reader
    chrom, center, window_size, n_bins = args

    if bw_reader is None:
        return np.full(n_bins, np.nan)

    start = max(0, center - window_size)
    end = center + window_size

    try:
        if chrom not in bw_reader.chroms():
            return np.full(n_bins, np.nan)

        chrom_len = bw_reader.chroms(chrom)
        if start >= chrom_len:
            return np.full(n_bins, np.nan)

        end = min(end, chrom_len)
        if start >= end:
            return np.full(n_bins, np.nan)

        values = bw_reader.stats(chrom, start, end, type="mean", nBins=n_bins)
        signal = [v if v is not None else np.nan for v in values]
        return np.array(signal)

    except Exception:
        return np.full(n_bins, np.nan)

def calculate_metaprofile(cpg_regions: pd.DataFrame,
                          bigwig_path: str,
                          window_size: int = DEFAULT_WINDOW_SIZE,
                          bin_size: int = DEFAULT_BIN_SIZE,
                          n_processes: int = 4) -> np.ndarray:
    """
    Calculate the meta-profile signal around the center of CpG regions.

    Args:
        cpg_regions: DataFrame with 'chrom', 'start', 'end' for CpG islands.
        bigwig_path: Path to the input BigWig file.
        window_size: Half-size of the window around the center (e.g., 2000 for +/- 2kb).
        bin_size: Size of bins for averaging signal.
        n_processes: Number of parallel processes to use.

    Returns:
        A numpy array representing the mean signal profile across all regions.
    """
    start_time = time.time()
    n_bins = window_size * 2 // bin_size
    logger.info(f"Calculating meta-profile for {len(cpg_regions)} CpG islands.")
    logger.info(f"BigWig file: {os.path.basename(bigwig_path)}")
    logger.info(f"Window: +/- {window_size} bp, Bin size: {bin_size} bp ({n_bins} bins)")

    # Prepare arguments for parallel processing
    tasks = []
    for _, region in cpg_regions.iterrows():
        center = (region['start'] + region['end']) // 2
        # Remove bigwig_path from task tuple
        tasks.append((region['chrom'], center, window_size, n_bins))

    # Process in parallel
    logger.info(f"Starting signal extraction with {n_processes} processes...")
    results = []
    # Use initializer to open BigWig file once per worker
    with ProcessPoolExecutor(max_workers=n_processes, initializer=init_worker, initargs=(bigwig_path,)) as executor:
        # Using tqdm to show progress
        results = list(tqdm(executor.map(extract_signal_for_region, tasks), total=len(tasks), desc="Processing regions"))

    # --- Important: Close the file handles in the worker processes --- #
    # Relying on process termination is standard practice.

    # Filter out potential None results if any error occurred during task submission itself (unlikely with map)
    results = [r for r in results if r is not None]

    if not results:
        logger.error("No valid signals could be extracted. Cannot compute meta-profile.")
        return np.full(n_bins, np.nan) # Return NaNs if no results

    # Stack results into a matrix (regions x bins)
    signal_matrix = np.vstack(results)

    # Calculate the mean profile, ignoring NaNs
    logger.info("Aggregating results and calculating mean profile...")
    # Calculate number of valid (non-NaN) regions per bin
    valid_counts = np.sum(~np.isnan(signal_matrix), axis=0)
    min_valid_regions = int(len(cpg_regions) * 0.1) # Require at least 10% of regions to have data for a bin

    mean_profile = np.nanmean(signal_matrix, axis=0)

    # Set bins with too few valid regions to NaN
    mean_profile[valid_counts < min_valid_regions] = np.nan

    # Report bins with insufficient data
    insufficient_bins = np.where(valid_counts < min_valid_regions)[0]
    if len(insufficient_bins) > 0:
        logger.warning(f"Bins {insufficient_bins} have data from less than {min_valid_regions} regions and were set to NaN.")

    total_regions_processed = signal_matrix.shape[0]
    regions_with_any_nan = np.sum(np.isnan(signal_matrix).any(axis=1))
    logger.info(f"Processed {total_regions_processed} regions. {regions_with_any_nan} regions had at least one NaN value (likely due to edge effects or missing chromosomes).")

    end_time = time.time()
    logger.info(f"Meta-profile calculation completed in {end_time - start_time:.2f} seconds.")

    return mean_profile

def plot_metaprofile(profile: np.ndarray,
                     window_size: int, # Actual window size used (e.g., 5000 for +/- 2500)
                     output_file: str,
                     title: str = "Methylation Signal at CpG Islands",
                     ylabel: str = "Average Signal"):
    """
    Plot the calculated meta-profile.
    """
    n_bins = len(profile)
    if n_bins == 0 or np.all(np.isnan(profile)):
        logger.error(f"Profile data is empty or all NaN for {output_file}. Cannot generate plot.")
        return

    # Recalculate half_window based on the actual window_size used for data generation
    half_window = window_size # The arg is the half-window size
    x = np.linspace(-half_window, half_window, n_bins)

    plt.figure(figsize=(8, 6))
    plt.plot(x, profile, color='red', linewidth=2)

    # Add vertical dashed line at the center
    plt.axvline(x=0, color='grey', linestyle='--', alpha=0.7)

    # Improve axis formatting using half_window
    plt.xticks([-half_window, 0, half_window], [f'-{half_window/1000:.1f}kb', 'CpG Center', f'+{half_window/1000:.1f}kb'])
    plt.xlabel("Distance from CpG Island Center (bp)")
    plt.ylabel(ylabel)
    plt.title(title)

    # Remove top and right spines for cleaner look
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add grid
    plt.grid(axis='y', linestyle=':', alpha=0.7)

    plt.tight_layout()

    try:
        plt.savefig(output_file, dpi=300)
        logger.info(f"Meta-profile plot saved to: {output_file}")
    except Exception as e:
        logger.error(f"Failed to save plot to {output_file}: {str(e)}")
    finally:
        plt.close() # Ensure plot is closed

def load_cpg_islands(bed_file: str) -> pd.DataFrame:
    """Load CpG islands from a BED file."""
    try:
        # Assuming standard BED like format: chr, start, end, ...
        # Use first 3 columns, ignore header if present
        regions = pd.read_csv(
            bed_file,
            sep='\t',
            header=None,
            usecols=[0, 1, 2],
            names=['chrom', 'start', 'end'],
            engine='python' # Handles potential parsing issues better
        )

        # Basic cleaning
        regions = regions.dropna()
        regions['start'] = regions['start'].astype(int)
        regions['end'] = regions['end'].astype(int)
        regions['chrom'] = regions['chrom'].astype(str)

        # Remove potentially invalid regions
        initial_count = len(regions)
        regions = regions[regions['end'] > regions['start']]
        regions = regions[regions['start'] >= 0]

        if len(regions) < initial_count:
            logger.warning(f"Removed {initial_count - len(regions)} invalid regions from {bed_file}")

        logger.info(f"Loaded {len(regions)} valid CpG regions from {bed_file}")
        return regions

    except FileNotFoundError:
        logger.error(f"Error: CpG BED file not found at {bed_file}")
        raise
    except Exception as e:
        logger.error(f"Error loading or parsing CpG BED file {bed_file}: {str(e)}")
        raise

# --- Comparison Plotting Function (Generalized) ---
def plot_comparison(input_npy_file1, input_npy_file2, label1, label2, output_plot_path, window_size):
    """
    Loads two specified smoothed metaprofile data (.npy) files and plots a comparison.

    Args:
        input_npy_file1 (str): Path to the first smoothed metaprofile .npy file.
        input_npy_file2 (str): Path to the second smoothed metaprofile .npy file.
        label1 (str): Label for the first profile in the legend.
        label2 (str): Label for the second profile in the legend.
        output_plot_path (str): Path to save the output comparison plot (.png).
        window_size (int): The half-window size used when generating the data.
    """
    logger.info(f"Attempting to generate comparison plot: {output_plot_path}")
    logger.info(f"  - File 1 ({label1}): {os.path.basename(input_npy_file1)}")
    logger.info(f"  - File 2 ({label2}): {os.path.basename(input_npy_file2)}")
    try:
        # Load profile data directly from .npy files
        profile1 = np.load(input_npy_file1)
        profile2 = np.load(input_npy_file2)
        logger.info(f"Loaded comparison data: {label1} ({len(profile1)} points), {label2} ({len(profile2)} points)")

        # Check if profiles have the same length (implies same binning)
        if len(profile1) != len(profile2):
            logger.error(f"Error: Profile lengths differ for comparison ({label1}: {len(profile1)}, {label2}: {len(profile2)}). Skipping comparison plot.")
            return # Skip plotting this comparison

        num_points = len(profile1)
        if num_points == 0:
            logger.warning(f"Warning: Profiles for comparison have zero length. Skipping plot.")
            return

        # Generate the x-axis bins based on window size and number of points
        half_window = window_size
        bins = np.linspace(-half_window, half_window, num_points)

        plt.figure(figsize=(10, 6))
        # Use specific colors as requested if needed, otherwise use defaults or specify here
        plt.plot(bins, profile1, color='lightcoral', label=label1) # Example: Using lightcoral for first curve (e.g., Neurons)
        plt.plot(bins, profile2, color='lightblue', label=label2) # Example: Using lightblue for second curve (e.g., NSCs)

        # --- Adjust x-axis ticks --- 
        tick_positions = [bins[0], 0, bins[-1]]
        tick_labels = [f'{-half_window/1000:.1f} kb', 'Center', f'{half_window/1000:.1f} kb']
        center_index = np.argmin(np.abs(bins))
        actual_tick_indices = [0, center_index, len(bins) - 1]
        actual_tick_positions = bins[actual_tick_indices]
        plt.xticks(actual_tick_positions, tick_labels)
        # --- End x-axis tick adjustment ---

        plt.xlabel("Position relative to CpG Target Center")
        plt.ylabel("Average Signal (Smoothed)") # Indicate smoothed
        plt.title(f"Smoothed Metaprofile Comparison") # Generalized title
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig(output_plot_path)
        logger.info(f"Saved comparison plot to {output_plot_path}")

    except FileNotFoundError as e:
        logger.error(f"Error: Input file not found for comparison plot - {e}. Skipping.")
    except Exception as e:
        logger.error(f"An unexpected error occurred while plotting comparison: {e}")
    finally:
        plt.close() # Ensure plot is closed
# --- End Comparison Plotting Function ---

def main():
    parser = argparse.ArgumentParser(description="Generate methylation meta-profile centered on CpG islands OR plot target comparisons.")

    # --- Arguments for Profile Calculation & Individual Plotting --- 
    calc_group = parser.add_argument_group('Profile Calculation & Plotting')
    calc_group.add_argument("--cpg-bed", help="Path to the CpG islands BED file (chrom, start, end).")
    calc_group.add_argument("--bigwig", help="Path to the methylation signal BigWig file.")
    calc_group.add_argument("--output-dir", help="Directory to save the output plot and data.")
    calc_group.add_argument("--window-size", type=int, default=DEFAULT_WINDOW_SIZE, help=f"Window size around the CpG center (half-size, default: {DEFAULT_WINDOW_SIZE} for +/- {DEFAULT_WINDOW_SIZE/1000}kb).")
    calc_group.add_argument("--bin-size", type=int, default=DEFAULT_BIN_SIZE, help=f"Bin size for signal calculation (default: {DEFAULT_BIN_SIZE} bp).")
    calc_group.add_argument("--threads", type=int, default=4, help="Number of parallel processes to use (default: 4).")
    calc_group.add_argument("--output-prefix", default=None, help="Optional prefix for output files (e.g., sample name). Required for calculation mode.")
    calc_group.add_argument("--force", action='store_true', help="Force recalculation even if cache files exist.")

    # --- Arguments for Comparison Plotting --- 
    compare_group = parser.add_argument_group('Comparison Plotting (Mutually Exclusive with Calculation)')
    compare_group.add_argument("--plot-comparison", action="store_true", help="If set, run in comparison plotting mode.")
    compare_group.add_argument("--input-npy-file1", help="Path to the first smoothed .npy file for comparison.")
    compare_group.add_argument("--input-npy-file2", help="Path to the second smoothed .npy file for comparison.")
    compare_group.add_argument("--label1", help="Label for the first curve in the comparison plot legend.")
    compare_group.add_argument("--label2", help="Label for the second curve in the comparison plot legend.")
    compare_group.add_argument("--comparison-output-file", help="Full path for the output comparison plot PNG file.")

    args = parser.parse_args()

    # --- Mode Selection --- 
    if args.plot_comparison:
        # --- Comparison Plotting Mode --- 
        logger.info("Running in Comparison Plotting Mode.")
        # Validate required arguments for comparison mode
        required_comp_args = [
            args.input_npy_file1,
            args.input_npy_file2,
            args.label1,
            args.label2,
            args.comparison_output_file
        ]
        if not all(required_comp_args):
            parser.error("Missing required arguments for --plot-comparison mode. Need --input-npy-file1, --input-npy-file2, --label1, --label2, and --comparison-output-file.")
            return

        # Construct input/output paths using smoothed suffixes - NO LONGER NEEDED
        # nsc_file = os.path.join(args.comparison_input_dir, f"{args.comparison_bw_label}{args.nsc_target_smoothed_suffix}")
        # neuron_file = os.path.join(args.comparison_input_dir, f"{args.comparison_bw_label}{args.neuron_target_smoothed_suffix}")
        
        # Use the window size argument (assuming it's passed correctly from the calculation step or is the default)
        comp_window_size = args.window_size # Get window size used for data generation
        
        # Make sure comparison output directory exists (extract from output file path)
        output_dir = os.path.dirname(args.comparison_output_file)
        if output_dir:
             os.makedirs(output_dir, exist_ok=True)
        else:
             logger.warning("No output directory specified in --comparison-output-file path. Saving to current directory.")

        # output_file = os.path.join(args.comparison_output_dir, f"{args.comparison_bw_label}_NSC_vs_Neuron_target_comparison.png") # Use direct path now

        plot_comparison(
            input_npy_file1=args.input_npy_file1,
            input_npy_file2=args.input_npy_file2,
            label1=args.label1,
            label2=args.label2,
            output_plot_path=args.comparison_output_file,
            window_size=comp_window_size # Pass the window size
        )
        logger.info("Comparison plotting finished.")
        return # Exit after comparison plot generation

    else:
        # --- Profile Calculation & Individual Plotting Mode --- 
        logger.info("Running in Profile Calculation/Plotting Mode.")
        # Validate required arguments for calculation mode
        required_calc_args = [args.cpg_bed, args.bigwig, args.output_dir, args.output_prefix]
        if not all(required_calc_args):
            parser.error("Missing required arguments for calculation mode. Need --cpg-bed, --bigwig, --output-dir, and --output-prefix.")
            return

        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)

        # --- Caching Logic --- 
        bw_basename = os.path.basename(args.bigwig).replace('.bw', '').replace('.bigWig', '')
        # Use the provided prefix directly now
        prefix = f"{args.output_prefix}_" 
        plot_filename = os.path.join(args.output_dir, f"{prefix}cpg_metaprofile.png")
        data_filename = os.path.join(args.output_dir, f"{prefix}cpg_metaprofile.npy")
        smoothed_plot_filename = os.path.join(args.output_dir, f"{prefix}cpg_metaprofile_smoothed.png")

        profile = None
        # Check if the data file already exists
        if not args.force and os.path.exists(data_filename):
            try:
                logger.info(f"Found cached meta-profile data: {data_filename}. Loading...")
                profile = np.load(data_filename)
                logger.info("Cached data loaded successfully.")
                # Check if both plots also exist
                plots_exist = os.path.exists(plot_filename) and os.path.exists(smoothed_plot_filename)
                if plots_exist:
                     logger.info(f"Found cached plot files: {plot_filename} and {smoothed_plot_filename}. Nothing more to do for this prefix.")
                     return # Exit early if data and both plots exist
                else:
                     logger.info(f"Cached plot file(s) missing. Will generate plot(s) from cached data.")
                     # Proceed to plotting step below

            except Exception as e:
                logger.warning(f"Failed to load cached data from {data_filename}: {str(e)}. Recalculating...")
                profile = None # Ensure profile is None if loading failed

        # If profile wasn't loaded from cache (or --force is used), calculate it
        if profile is None:
            if args.force:
                logger.info(f"--force flag detected for {args.output_prefix}. Recalculating and overwriting any existing cache.")
                # Optionally, explicitly delete old files first
                for f in [data_filename, plot_filename, smoothed_plot_filename]:
                     try:
                         if os.path.exists(f): os.remove(f)
                     except OSError as e:
                         logger.warning(f"Could not delete existing cache file {f}: {e}")
            else:
                 logger.info(f"Cached data for {args.output_prefix} not found or failed to load. Proceeding with calculation...")

            # Load CpG regions
            try:
                cpg_regions = load_cpg_islands(args.cpg_bed)
            except Exception:
                return # Exit if loading fails

            if cpg_regions.empty:
                logger.error(f"No valid CpG regions loaded from {args.cpg_bed}. Exiting for prefix {args.output_prefix}.")
                return

            # Calculate meta-profile
            profile = calculate_metaprofile(
                cpg_regions,
                args.bigwig,
                window_size=args.window_size,
                bin_size=args.bin_size,
                n_processes=args.threads
            )

            if np.all(np.isnan(profile)):
                 logger.error(f"Meta-profile calculation for {args.output_prefix} resulted in all NaN values. Cannot proceed.")
                 return

            # Save profile data only if it was calculated
            logger.info(f"Attempting to save profile data to: {data_filename}")
            try:
                np.save(data_filename, profile)
                logger.info(f"Meta-profile data saved successfully.")
            except Exception as e:
                logger.error(f"Failed to save profile data to {data_filename}: {str(e)}")
                # Decide if you want to stop here or still attempt plotting

        # --- Plotting Section (Individual Plots) ---
        if profile is not None and not np.all(np.isnan(profile)):
            # 1. Handle Original Plot
            if args.force or not os.path.exists(plot_filename):
                 logger.info(f"Generating original plot: {plot_filename}")
                 plot_title = f"Signal from {bw_basename} at {args.output_prefix} regions"
                 plot_metaprofile(
                     profile,
                     args.window_size, # Pass the half-window size
                     plot_filename,
                     title=plot_title
                 )
            else:
                 logger.info(f"Original plot file {plot_filename} already exists and --force not specified. Skipping.")

            # 2. Handle Smoothed Plot
            if args.force or not os.path.exists(smoothed_plot_filename):
                try:
                    smoothing_sigma = 3.0 # Adjust as needed
                    if np.isnan(profile).any():
                         logger.warning(f"Profile for {args.output_prefix} contains NaNs, cannot generate smoothed plot or save smoothed data.")
                    else:
                        logger.info(f"Generating smoothed profile for {args.output_prefix} (sigma={smoothing_sigma})...")
                        smoothed_profile = gaussian_filter1d(profile, sigma=smoothing_sigma)

                        # --- Save Smoothed Profile Data --- 
                        smoothed_data_filename = os.path.join(args.output_dir, f"{prefix}cpg_metaprofile_smoothed.npy")
                        logger.info(f"Attempting to save smoothed profile data to: {smoothed_data_filename}")
                        try:
                            np.save(smoothed_data_filename, smoothed_profile)
                            logger.info(f"Smoothed profile data saved successfully.")
                        except Exception as e:
                            logger.error(f"Failed to save smoothed profile data to {smoothed_data_filename}: {str(e)}")
                        # --- End Save Smoothed Profile Data ---

                        # --- Plot Smoothed Profile --- 
                        logger.info(f"Generating smoothed plot: {smoothed_plot_filename}...")
                        smoothed_plot_title = f"Smoothed Signal ({bw_basename} at {args.output_prefix}, Sigma={smoothing_sigma})"
                        plot_metaprofile(
                             smoothed_profile,
                             args.window_size, # Pass the half-window size
                             smoothed_plot_filename,
                             title=smoothed_plot_title
                        )
                except ImportError:
                     logger.warning("SciPy not found. Cannot generate smoothed plot. Please install SciPy (`pip install scipy`).")
                except Exception as e:
                     logger.error(f"Error generating smoothed plot for {args.output_prefix}: {str(e)}")
            else:
                 logger.info(f"Smoothed plot file {smoothed_plot_filename} already exists and --force not specified. Skipping.")

        elif profile is None:
             logger.error(f"Profile calculation failed or was skipped for {args.output_prefix}, cannot plot.")

        logger.info(f"Individual processing finished for prefix: {args.output_prefix}")

if __name__ == "__main__":
    main() 