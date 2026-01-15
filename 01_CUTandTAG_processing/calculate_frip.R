# For calculating Fraction of Reads in Peaks
library(GenomicAlignments)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
peaks_file <- args[2]
output_file <- args[3]

# Read BAM file
reads <- readGAlignments(bam_file)

# Read peaks
peaks <- import(peaks_file)

# Calculate FRIP
reads_in_peaks <- countOverlaps(reads, peaks)
frip_score <- sum(reads_in_peaks > 0) / length(reads)

# Write output
write(frip_score, output_file)