#!/bin/bash
#SBATCH --job-name=03_featureCounts
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --error="./logs/03_featurecounts.err"
#SBATCH --output="./logs/03_featurecounts.out"

# RNA-seq Pipeline Step 3: Gene Quantification with featureCounts
# This script counts reads mapping to genes using featureCounts

# Set up conda environment with featureCounts
conda activate rnaseq  # Adjust to your environment name

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Set PROJECT_DIR to your project directory
PROJECT_DIR=""  # e.g., "/path/to/project"
cd ${PROJECT_DIR}

echo "Starting featureCounts Quantification"
echo "Start time: $(date)"

# Define paths - set GTF to your annotation file
BAM_DIR="results/alignment"
GTF=""  # e.g., "/path/to/annotation.gtf"
OUTPUT_DIR="results/counts"
THREADS=16

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Get all BAM files
BAM_FILES=(${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

echo "Found ${#BAM_FILES[@]} BAM files to process:"
printf '%s\n' "${BAM_FILES[@]}"
echo ""

# Run featureCounts
echo "Running featureCounts..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -g gene_id \
    -t exon \
    -s 0 \
    -a ${GTF} \
    -o ${OUTPUT_DIR}/gene_counts.txt \
    ${BAM_FILES[@]}

echo ""
echo "Creating simplified count matrix..."
# Extract only sample names and counts (remove first column which is Geneid, and columns 2-6 which are Chr, Start, End, Strand, Length)
# Keep Geneid and count columns only
cut -f1,7- ${OUTPUT_DIR}/gene_counts.txt > ${OUTPUT_DIR}/gene_counts_matrix.txt

# Clean up sample names in header (remove path and suffix)
sed -i '1 s|results/alignment/||g; 1 s|_Aligned.sortedByCoord.out.bam||g' ${OUTPUT_DIR}/gene_counts_matrix.txt

echo ""
echo "Gene quantification complete!"
echo "Output files:"
echo "  - ${OUTPUT_DIR}/gene_counts.txt (full output)"
echo "  - ${OUTPUT_DIR}/gene_counts.txt.summary (summary statistics)"
echo "  - ${OUTPUT_DIR}/gene_counts_matrix.txt (simplified count matrix)"
echo ""
echo "Summary statistics:"
cat ${OUTPUT_DIR}/gene_counts.txt.summary
echo ""
echo "End time: $(date)"
echo "=========================================="
