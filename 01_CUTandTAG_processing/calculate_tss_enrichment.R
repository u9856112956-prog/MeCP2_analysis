# For calculating TSS enrichment
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
gtf_file <- args[2]
score_file <- args[3]
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[5])

# Read GTF and extract TSS
txdb <- makeTxDbFromGFF(gtf_file)
genes <- genes(txdb)
tss <- promoters(genes, upstream=upstream, downstream=downstream)

# Read BAM file
reads <- readGAlignments(bam_file)

# Calculate coverage around TSS
coverage <- coverage(reads)
tss_coverage <- lapply(seqlevels(tss), function(chr) {
    chr_tss <- tss[seqnames(tss) == chr]
    if (length(chr_tss) > 0) {
        chr_coverage <- coverage[[chr]]
        viewSums(Views(chr_coverage, ranges(chr_tss)))
    }
})

# Calculate enrichment score
tss_coverage <- unlist(tss_coverage)
background <- mean(tss_coverage[c(1:100, (length(tss_coverage)-100):length(tss_coverage))])
enrichment_score <- max(tss_coverage) / background

# Write score to file
write(enrichment_score, score_file)