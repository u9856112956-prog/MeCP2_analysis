#!/bin/bash
#SBATCH --job-name=metaprofiles
#SBATCH --mem=64GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL

# Activate conda environment with deepTools
conda activate snakemake  # Adjust to your environment name

# Set WORKDIR to your CUT&TAG analysis directory
WORKDIR=""  # e.g., "/path/to/cutandtag/analysis"
cd $WORKDIR

# Create output directory
mkdir -p results/metaprofiles

# Define file paths
GTF="/path/to/mm10.gtf"  # You'll need to specify the path to your mm10 GTF file
BLACKLIST="/path/to/mm10-blacklist.bed"  # Optional: specify path to blacklist regions

# Define sample groups
NEU_MECP2_SAMPLES=(
    "results/bigwig/NeuM2.bw"
    "results/bigwig/NeuM3.bw"
)

NEU_VECTOR_SAMPLES=(
    "results/bigwig/NeuV1.bw"
    "results/bigwig/NeuV2.bw"
    "results/bigwig/NeuV3.bw"
)

NSC_MECP2_SAMPLES=(
    "results/bigwig/NSCM1.bw"
    "results/bigwig/NSCM2.bw"
    "results/bigwig/NSCM3.bw"
)

NSC_VECTOR_SAMPLES=(
    "results/bigwig/NSCv1.bw"
    "results/bigwig/NSCv2.bw"
    "results/bigwig/NSCv3.bw"
)

# Step 1: Create regions file for analysis
echo "Creating regions file..."
computeMatrix scale-regions \
    -S "${NEU_MECP2_SAMPLES[@]}" "${NEU_VECTOR_SAMPLES[@]}" \
       "${NSC_MECP2_SAMPLES[@]}" "${NSC_VECTOR_SAMPLES[@]}" \
    -R $GTF \
    --beforeRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 5000 \
    --skipZeros \
    --numberOfProcessors 16 \
    --transcriptID gene \
    --metagene \
    -o results/metaprofiles/matrix.gz \
    --outFileSortedRegions results/metaprofiles/regions.bed

# Step 2: Create profile plot
echo "Generating profile plot..."
plotProfile \
    -m results/metaprofiles/matrix.gz \
    -out results/metaprofiles/metaprofile.pdf \
    --perGroup \
    --colors "#FF0000" "#0000FF" "#00FF00" "#000000" \
    --plotTitle "Mecp2 binding profile around genes" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "Genes" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

# Step 3: Create heatmap
echo "Generating heatmap..."
plotHeatmap \
    -m results/metaprofiles/matrix.gz \
    -out results/metaprofiles/heatmap.pdf \
    --colorMap RdYlBu \
    --whatToShow 'plot, heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --zMin 0 \
    --zMax 10

# Step 4: Generate TSS-centered profile
echo "Generating TSS-centered profile..."
computeMatrix reference-point \
    -S "${NEU_MECP2_SAMPLES[@]}" "${NEU_VECTOR_SAMPLES[@]}" \
       "${NSC_MECP2_SAMPLES[@]}" "${NSC_VECTOR_SAMPLES[@]}" \
    -R $GTF \
    --referencePoint TSS \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --skipZeros \
    --numberOfProcessors 16 \
    -o results/metaprofiles/tss_matrix.gz

plotProfile \
    -m results/metaprofiles/tss_matrix.gz \
    -out results/metaprofiles/tss_profile.pdf \
    --perGroup \
    --colors "#FF0000" "#0000FF" "#00FF00" "#000000" \
    --plotTitle "Mecp2 binding profile around TSS" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "TSS" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

# Step 5: Generate statistics
echo "Generating statistics..."
plotProfileWithStderr \
    -m results/metaprofiles/matrix.gz \
    -out results/metaprofiles/profile_with_stderr.pdf \
    --perGroup \
    --plotTitle "Mecp2 binding profile with standard error" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "Genes" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

echo "Analysis complete. Check results in results/metaprofiles/" 