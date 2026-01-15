#!/bin/bash
#SBATCH --job-name=2_run_gene_chromatin_state_comparison_nsc
#SBATCH --output=logs/2_run_gene_chromatin_state_comparison_nsc.out
#SBATCH --error=logs/2_run_gene_chromatin_state_comparison_nsc.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=workq

# Chromatin state comparison between NSC and Neu samples (GFP condition)
# Classifies gene promoters as euchromatin (open) or heterochromatin (closed)

# Set BASE_DIR to your SAMMY-seq data directory
BASE_DIR=""  # e.g., "/path/to/sammy_seq"
WORKDIR="${BASE_DIR}/ratio_analysis"

cd ${WORKDIR}
mkdir -p logs

# Activate conda environment
conda activate snakemake  # Adjust to your environment name

# Configure paths - adjust these to your data
GENE_LIST="${BASE_DIR}/gene_lists/nsc_gene_list.txt"  # Your gene list
RESULTS_DIR="${WORKDIR}/results/chromatin_states"
OUTPUT_DIR="${WORKDIR}/results/nsc_comparison"
GENOME_GTF="${BASE_DIR}/DATA/gencode.vM25.basic.annotation.gtf"  # Gene annotation
PROMOTER_SIZE=2000

mkdir -p ${OUTPUT_DIR}

python 2_gene_chromatin_state_comparison.py \
    --gene-list ${GENE_LIST} \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --promoter-size ${PROMOTER_SIZE}

echo "Analysis completed."
