#!/bin/bash
#SBATCH --job-name=0_run_sammy_seq_analysis_blacklisted
#SBATCH --output=logs/0_run_sammy_seq_analysis_blacklisted.out
#SBATCH --error=logs/0_run_sammy_seq_analysis_blacklisted.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --partition=workq

# SAMMY-seq analysis pipeline with ENCODE blacklist filtering
# Identifies heterochromatin and euchromatin regions from S2S/S3 fractions

# Set BASE_DIR to your SAMMY-seq data directory
BASE_DIR=""  # e.g., "/path/to/sammy_seq"
WORKDIR="${BASE_DIR}/ratio_analysis"
OUTDIR="${WORKDIR}/results/chromatin_states_blacklisted"
BLACKLIST="${WORKDIR}/ENCFF547MET.bed"

cd ${WORKDIR}
mkdir -p logs
mkdir -p ${OUTDIR}

# Activate conda environment
conda activate snakemake  # Adjust to your environment name

python 0_sammy_seq_analysis_blacklisted.py \
  --workdir ${WORKDIR} \
  --outdir ${OUTDIR} \
  --bin-size 50000 \
  --eu-threshold 0.1 \
  --hetero-threshold -0.1 \
  --threads 16 \
  --blacklist ${BLACKLIST}