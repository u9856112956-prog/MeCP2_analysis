library(optparse)
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)

# Parse command line arguments
option_list <- list(
    make_option("--peaks", type="character", help="Comma-separated list of peak files"),
    make_option("--output", type="character", help="Output file path"),
    make_option("--genome", type="character", help="Path to genome GTF file")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Read peak files
peak_files <- unlist(strsplit(opt$peaks, " "))
peak_list <- lapply(peak_files, import)
names(peak_list) <- basename(peak_files)

# Create TxDb object from GTF
txdb <- makeTxDbFromGFF(opt$genome)

# Annotate peaks
peak_anno <- lapply(peak_list, function(peaks) {
    annotatePeak(peaks, TxDb=txdb, annoDb="org.Mm.eg.db")
})

# Extract annotations and combine
results <- do.call(rbind, lapply(names(peak_anno), function(sample) {
    anno <- peak_anno[[sample]]@anno
    data.frame(
        Sample = sample,
        PeakID = anno$peak,
        Gene = anno$geneId,
        Distance = anno$distanceToTSS,
        Annotation = anno$annotation
    )
}))

# Write results
write.table(results, opt$output, sep="\t", quote=FALSE, row.names=FALSE) 