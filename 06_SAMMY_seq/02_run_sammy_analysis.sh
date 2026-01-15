#!/bin/bash
#SBATCH --job-name=run_sammy_seq_analysis
#SBATCH --output=logs/run_sammy_seq_analysis.out
#SBATCH --error=logs/run_sammy_seq_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --partition=workq

# SAMMY-seq analysis pipeline for identifying heterochromatin and euchromatin regions
# Based on comparison of S2S (accessible) and S3 (inaccessible) chromatin fractions

# Set WORKDIR to your working directory containing processed BAM files
WORKDIR=""  # e.g., "/path/to/sammy_seq/analysis"
OUTDIR="${WORKDIR}/results/chromatin_states"

cd ${WORKDIR}
mkdir -p logs
mkdir -p ${OUTDIR}

# Activate conda environment
conda activate snakemake  # Adjust to your environment name

# Run the Python script with all necessary arguments
python sammy_seq_analysis.py \
  --workdir ${WORKDIR} \
  --outdir ${OUTDIR} \
  --bin-size 50000 \
  --eu-threshold 0.1 \
  --hetero-threshold -0.1 \
  --threads 16 \