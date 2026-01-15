# For filtering peaks based on quality metrics
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
peaks_file <- args[1]
frip_file <- args[2]
output_file <- args[3]
min_frip <- as.numeric(args[4])
min_score <- as.numeric(args[5])

# Read FRIP score
frip <- as.numeric(readLines(frip_file)[1])

# Filter peaks if FRIP score is acceptable
if (frip >= min_frip) {
    peaks <- import(peaks_file)
    filtered_peaks <- peaks[peaks$score >= min_score]
    export(filtered_peaks, output_file)
} else {
    warning("FRIP score below threshold")
    file.create(output_file)  # Create empty file
}