library(DiffBind)
library(tidyverse)
library(BiocParallel)

# Get the experiment name and number of threads from Snakemake
experiment <- snakemake@wildcards[["experiment"]]
threads <- snakemake@threads

# Set up parallel processing
register(MulticoreParam(workers = threads))

# Read sample sheet
# Set DATA_DIR to your data directory containing sample sheets
DATA_DIR <- ""  # e.g., "/path/to/data"
samples <- read.csv(paste0(DATA_DIR, experiment, "_sample_sheet.csv"))

# Create DiffBind object
dba <- dba(sampleSheet=samples)

# Count reads with parallelization
dba <- dba.count(dba, bParallel=TRUE)

# Create a data frame with the number of peaks for each sample
peaks_df <- data.frame(
  Sample = names(dba$peaks),
  Peaks = unlist(dba$peaks)
)

# Write the peaks data frame to a CSV file
write.csv(peaks_df, snakemake@output[["results"]], row.names = FALSE)

# Print number of peaks found for each sample (optional, for console output)
print("Number of peaks found for each sample:")
print(peaks_df)

# Normalize
dba <- dba.normalize(dba)

# Perform differential binding analysis
dba <- dba.contrast(dba, categories=c("Condition", "Tissue"))
dba <- dba.analyze(dba, bParallel=TRUE)

# Extract results
res <- dba.report(dba)

# Write results
write_csv(res, snakemake@output[["results"]])