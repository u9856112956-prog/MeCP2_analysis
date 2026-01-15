#!/bin/bash
#SBATCH --job-name=3_run_mecp2_metaprofile_heatmap
#SBATCH --output=logs/3_run_mecp2_metaprofile_heatmap.out
#SBATCH --error=logs/3_run_mecp2_metaprofile_heatmap.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --partition=workq

# Generates metaprofile heatmaps for MeCP2 target and non-target genes
# Shows log2(S2S/S3) ratio across gene bodies with chromatin state coloring

# Set BASE_DIR to your SAMMY-seq data directory
BASE_DIR=""  # e.g., "/path/to/sammy_seq"
WORKDIR="${BASE_DIR}/ratio_analysis"

cd ${WORKDIR}
mkdir -p logs

# Activate conda environment
conda activate snakemake  # Adjust to your environment name

# Configure paths - adjust these to your data
RESULTS_DIR="${WORKDIR}/results/chromatin_states"
OUTPUT_DIR="${WORKDIR}/results/metaprofile_heatmaps"
GENOME_GTF="${BASE_DIR}/DATA/gencode.vM25.basic.annotation.gtf"  # Gene annotation
GENE_LISTS_DIR="${BASE_DIR}/gene_lists"  # Directory with gene list files
UPSTREAM_DISTANCE=5000
DOWNSTREAM_DISTANCE=5000
BIN_SIZE=100
GENE_BODY_BINS=100

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run the analysis
python 3_create_mecp2_metaprofile_heatmap.py \
    --results-dir ${RESULTS_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --genome-gtf ${GENOME_GTF} \
    --gene-lists-dir ${GENE_LISTS_DIR} \
    --upstream-distance ${UPSTREAM_DISTANCE} \
    --downstream-distance ${DOWNSTREAM_DISTANCE} \
    --bin-size ${BIN_SIZE} \
    --gene-body-bins ${GENE_BODY_BINS}

echo "Metaprofile heatmap generation completed." 