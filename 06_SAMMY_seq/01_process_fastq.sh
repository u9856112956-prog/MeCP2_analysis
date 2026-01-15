#!/bin/bash
#SBATCH --job-name=process_samples
#SBATCH --output=logs/process_samples_%a.out
#SBATCH --error=logs/process_samples_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=workq
#SBATCH --array=1-34%10  # Update based on the number of files

# Error handling
set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools (BWA, samtools, picard, bamCoverage)
conda activate snakemake  # Adjust to your environment name

# Set WORKDIR to your working directory
WORKDIR=""  # e.g., "/path/to/sammy_seq/analysis"
cd $WORKDIR

# Directory setup
FASTQ_DIR="${WORKDIR}/../DATA/Sammy_Seq_fastq"
OUTPUT_DIR="${WORKDIR}/results"
# Need to set up BWA index for mm10 genome
BWA_GENOME_DIR="${WORKDIR}/mm10"
BLACKLIST_FILE="${WORKDIR}/mm10-blacklist.v2.bed"
THREADS=16  # Using all CPUs allocated to the job

# Create a directory for trimmed reads
mkdir -p $OUTPUT_DIR/trimmed
mkdir -p $OUTPUT_DIR/bam/temp
mkdir -p logs

# Download blacklist file if it doesn't exist
if [ ! -f "$BLACKLIST_FILE" ]; then
  echo "Downloading mm10 blacklist file..."
  mkdir -p $(dirname $BLACKLIST_FILE)
  wget -O $BLACKLIST_FILE https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz
  gunzip -f $BLACKLIST_FILE.gz
fi

# Check if the FASTQ directory exists
if [ ! -d "$FASTQ_DIR" ]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist"
  exit 1
fi
echo "Using FASTQ directory: $FASTQ_DIR"

# Create output directories
mkdir -p $OUTPUT_DIR/bam
mkdir -p logs

# Check if the filtered_sample_list.txt exists
if [ ! -f "filtered_sample_list.txt" ]; then
  echo "Error: filtered_sample_list.txt not found"
  echo "Please run check_number_of_files.sh first to generate the file list"
  exit 1
fi

# Count how many files we have
FILE_COUNT=$(wc -l < filtered_sample_list.txt)
echo "Found $FILE_COUNT files in the filtered sample list"

# Display the first few files to verify
echo -e "\nFirst 10 files to process:"
head -n 10 filtered_sample_list.txt

# Check if the array size matches the number of files
if [[ "$SLURM_ARRAY_TASK_MAX" ]]; then
  ARRAY_SIZE=$((SLURM_ARRAY_TASK_MAX))
  if [ "$ARRAY_SIZE" != "$FILE_COUNT" ]; then
    echo "Warning: Array size ($ARRAY_SIZE) does not match the number of files ($FILE_COUNT)"
    echo "For optimal processing, you should cancel this job and resubmit with:"
    echo "sbatch --array=1-${FILE_COUNT}%10 $0"
  fi
else
  echo "Note: Running in non-array mode or array information not available"
fi

# Get the FASTQ file for this array job
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filtered_sample_list.txt)
if [ -z "$FASTQ" ]; then
  echo "Error: No file found for task ID ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
FASTQ="$FASTQ_DIR/$FASTQ"

if [ ! -f "$FASTQ" ]; then
    echo "Error: FASTQ file $FASTQ not found"
    exit 1
fi

# Get base name for output files
BASE=$(basename $FASTQ _R1_001.fastq.gz)
echo "Processing $BASE..."

# Create a directory to store trimmed adapter files if it doesn't exist
ADAPTER_DIR="$WORKDIR/adapters"
mkdir -p $ADAPTER_DIR
ADAPTER_FILE="$ADAPTER_DIR/TruSeq3-SE-2.fa"

# Create adapter file if it doesn't exist
if [ ! -f "$ADAPTER_FILE" ]; then
  echo "Creating TruSeq3-SE-2.fa adapter file..."
  cat > "$ADAPTER_FILE" << 'EOL'
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
EOL
fi

# Step 1: Trim reads using Trimmomatic as specified in the paper
TRIMMED_FASTQ="$OUTPUT_DIR/trimmed/${BASE}.trimmed.fastq.gz"
if [ -f "$TRIMMED_FASTQ" ]; then
  echo "Trimmed file already exists, skipping trimming step..."
else
  echo "Trimming reads with Trimmomatic..."
  trimmomatic SE $FASTQ $TRIMMED_FASTQ \
    ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
fi

# Step 2: Align with BWA as specified in the paper
echo "Aligning with BWA..."
# Check if BWA index exists and is valid, if not create it
BWA_INDEX_FILES=("$BWA_GENOME_DIR/mm10.amb" "$BWA_GENOME_DIR/mm10.ann" "$BWA_GENOME_DIR/mm10.bwt" "$BWA_GENOME_DIR/mm10.pac" "$BWA_GENOME_DIR/mm10.sa")
INDEX_VALID=true

# Check if all index files exist and have non-zero size
for index_file in "${BWA_INDEX_FILES[@]}"; do
  if [ ! -s "$index_file" ]; then
    INDEX_VALID=false
    echo "BWA index file $index_file is missing or empty."
    break
  fi
done

# Even if all files exist, perform a more thorough validation by attempting to use the index
if [ "$INDEX_VALID" = true ]; then
  echo "Testing BWA index validity..."
  # Create a small test FASTQ with a few reads
  TEST_FASTQ="${BWA_GENOME_DIR}/test_reads.fastq"
  echo -e "@test_read\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII" > "$TEST_FASTQ"
  
  # Try to align the test read - redirect stderr to capture any errors
  TEST_OUTPUT=$(bwa mem "$BWA_GENOME_DIR/mm10" "$TEST_FASTQ" 2>&1 || true)
  
  # Check if the output contains error messages related to index files
  if echo "$TEST_OUTPUT" | grep -q -E "(Can't|Cannot|failed|error|Unexpected end of file)"; then
    INDEX_VALID=false
    echo "BWA index validation failed with error:"
    echo "$TEST_OUTPUT" | grep -E "(Can't|Cannot|failed|error|Unexpected end of file)"
  else
    echo "BWA index validation successful."
  fi
  
  # Clean up test file
  rm -f "$TEST_FASTQ"
fi

if [ "$INDEX_VALID" = false ]; then
  echo "BWA index is invalid or incomplete. Removing old index files and recreating..."
  # Remove any existing index files that might be corrupted
  rm -f ${BWA_INDEX_FILES[@]}
  # Also remove any potential temporary files that might have been created during a failed indexing
  rm -f "$BWA_GENOME_DIR/mm10."*
  mkdir -p $BWA_GENOME_DIR
  
  # Check if reference genome exists
  if [ ! -f "$WORKDIR/references/mm10.fa" ]; then
    echo "Error: Reference genome file $WORKDIR/references/mm10.fa not found"
    exit 1
  fi
  
  # Check available disk space before indexing
  GENOME_SIZE=$(du -b "$WORKDIR/references/mm10.fa" | cut -f1)
  REQUIRED_SPACE=$((GENOME_SIZE * 5)) # BWA index typically requires ~5x the genome size
  AVAILABLE_SPACE=$(df -B1 "$BWA_GENOME_DIR" | awk 'NR==2 {print $4}')
  
  if [ "$AVAILABLE_SPACE" -lt "$REQUIRED_SPACE" ]; then
    echo "Error: Not enough disk space for BWA indexing. Need at least $((REQUIRED_SPACE/1024/1024/1024)) GB, but only $((AVAILABLE_SPACE/1024/1024/1024)) GB available."
    exit 1
  fi
  
  echo "Creating BWA index. This may take a while..."
  bwa index -p $BWA_GENOME_DIR/mm10 "$WORKDIR/references/mm10.fa"
  
  # Verify index was created successfully
  for index_file in "${BWA_INDEX_FILES[@]}"; do
    if [ ! -s "$index_file" ]; then
      echo "Error: Failed to create BWA index file $index_file"
      exit 1
    fi
  done
else
  echo "BWA index is valid and complete."
fi

# Check if BAM file already exists
if [ -f "$OUTPUT_DIR/bam/${BASE}.bam" ]; then
  echo "BAM file already exists, skipping alignment step..."
else
  # Align with BWA using parameters from the paper
  echo "Performing BWA alignment..."
  bwa mem -t $THREADS $BWA_GENOME_DIR/mm10 $OUTPUT_DIR/trimmed/${BASE}.trimmed.fastq.gz > $OUTPUT_DIR/bam/${BASE}.sam

  # Convert SAM to BAM
  echo "Converting SAM to BAM..."
  samtools view -bS $OUTPUT_DIR/bam/${BASE}.sam > $OUTPUT_DIR/bam/${BASE}.bam
  rm $OUTPUT_DIR/bam/${BASE}.sam
fi

# Sort BAM if it exists but hasn't been sorted yet
if [ -f "$OUTPUT_DIR/bam/${BASE}.bam" ]; then
  # Check if BAM is already sorted
  if ! samtools view -H "$OUTPUT_DIR/bam/${BASE}.bam" | grep -q 'SO:coordinate'; then
    echo "Sorting BAM..."
    samtools sort -@ $THREADS $OUTPUT_DIR/bam/${BASE}.bam -o $OUTPUT_DIR/bam/${BASE}.sorted.bam
    mv $OUTPUT_DIR/bam/${BASE}.sorted.bam $OUTPUT_DIR/bam/${BASE}.bam
  else
    echo "BAM file already sorted, skipping sorting step..."
  fi
fi

# Step 3: Mark and remove PCR duplicates using Picard
if [ -f "$OUTPUT_DIR/bam/${BASE}.dedup.bam" ]; then
  echo "Deduplicated BAM already exists, skipping deduplication step..."
else
  echo "Marking and removing PCR duplicates with Picard..."
  picard MarkDuplicates \
    I=$OUTPUT_DIR/bam/${BASE}.bam \
    O=$OUTPUT_DIR/bam/${BASE}.dedup.bam \
    M=$OUTPUT_DIR/bam/${BASE}.metrics.txt \
    REMOVE_DUPLICATES=true \
    TMP_DIR=$OUTPUT_DIR/bam/temp
fi

# Step 4: Filter by mapping quality (q>=1) as in the paper
if [ -f "$OUTPUT_DIR/bam/${BASE}.filtered.bam" ]; then
  echo "Filtered BAM already exists, skipping filtering step..."
else
  echo "Filtering by mapping quality..."
  samtools view -q 1 -b $OUTPUT_DIR/bam/${BASE}.dedup.bam > $OUTPUT_DIR/bam/${BASE}.filtered.bam
fi

# Index BAM if needed
if [ ! -f "$OUTPUT_DIR/bam/${BASE}.filtered.bam.bai" ]; then
  echo "Indexing BAM..."
  samtools index $OUTPUT_DIR/bam/${BASE}.filtered.bam
else
  echo "BAM index already exists, skipping indexing step..."
fi

# Step 5: Create bigWig file
if [ -f "$OUTPUT_DIR/bigwig/${BASE}.bw" ]; then
  echo "BigWig file already exists, skipping bigWig creation step..."
else
  echo "Creating bigWig for $BASE..."
  bamCoverage \
    --bam $OUTPUT_DIR/bam/${BASE}.filtered.bam \
    --outFileName $OUTPUT_DIR/bigwig/${BASE}.bw \
    --binSize 50 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2652783500 \
    --blackListFileName $BLACKLIST_FILE \
    --extendReads 250 \
    --numberOfProcessors $THREADS
fi

echo "Processing complete for $BASE"