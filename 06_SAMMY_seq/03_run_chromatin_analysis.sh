#!/bin/bash
#SBATCH --job-name=1_run_gene_chromatin_state_analysis
#SBATCH --output=logs/1_run_gene_chromatin_state_analysis.out
#SBATCH --error=logs/1_run_gene_chromatin_state_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --partition=workq

# Chromatin state analysis for specific genes in NSC and Neu samples
# Analyzes euchromatin/heterochromatin states at gene promoter regions

# Set BASE_DIR to your SAMMY-seq data directory
BASE_DIR=""  # e.g., "/path/to/sammy_seq"
WORKDIR="${BASE_DIR}/ratio_analysis"

cd ${WORKDIR}
mkdir -p logs

# Activate conda environment
conda activate snakemake  # Adjust to your environment name

# Configure paths - adjust these to your data
GENE_LIST="${BASE_DIR}/gene_lists/gene_list.txt"  # Your gene list file
RESULTS_DIR="${WORKDIR}/results/chromatin_states"
OUTPUT_DIR="${WORKDIR}/results/gene_chromatin_analysis"
GENOME_GTF="${BASE_DIR}/DATA/gencode.vM25.basic.annotation.gtf"  # Gene annotation
PROMOTER_SIZE=2000

mkdir -p ${OUTPUT_DIR}

python 1_gene_chromatin_state_analysis.py \
  --gene-list ${GENE_LIST} \
  --results-dir ${RESULTS_DIR} \
  --output-dir ${OUTPUT_DIR} \
  --genome-gtf ${GENOME_GTF} \
  --promoter-size ${PROMOTER_SIZE}