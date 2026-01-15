library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Parse command line arguments
library(optparse)

# Create option parser with explicit dest parameter
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory path")
)

# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Working directory path:")
print(opts$work_dir)

# Check if work-dir argument is provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}

# Set working directory
message(sprintf("Setting working directory to: %s", opts$work_dir))
setwd(opts$work_dir)

# Create DE_integrated directory if it doesn't exist
dir.create("DE_integrated", showWarnings = FALSE)

# Function to process a single binding dataset and merge with DE data
process_binding_data <- function(file_path, output_prefix) {
  # Read binding data
  binding_df <- read.csv(file_path)
  
  # Convert to GRanges object
  gr <- GRanges(
    seqnames = binding_df$chr,
    ranges = IRanges(start = binding_df$start, end = binding_df$end),
    enrichment = binding_df$enrichment,
    binding_type = binding_df$binding_type,
    peak_width_exo = binding_df$peak_width_exo,
    peak_width_endo = binding_df$peak_width_endo,
    significant = binding_df$significant
  )
  
  # Get transcript database
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # Annotate peaks with nearest genes
  peakAnno <- annotatePeak(gr, 
                          tssRegion=c(-3000, 3000),
                          TxDb=txdb,
                          annoDb="org.Mm.eg.db")
  
  # Convert to data frame
  anno_df <- as.data.frame(peakAnno@anno)
  
  # Convert ENTREZ IDs to gene symbols
  anno_df$gene_symbol <- mapIds(org.Mm.eg.db,
                               keys=anno_df$geneId,
                               column="SYMBOL",
                               keytype="ENTREZID",
                               multiVals="first")
  
  # Read differential expression data - set path to your DEA results file
  de_data <- read.csv("DEA_NSC.csv")  # Adjust path as needed
  
  # Merge with differential expression data
  merged_data <- merge(anno_df,
                      de_data,
                      by.x="gene_symbol",
                      by.y="gene",
                      all.x=TRUE)
  
  # Add original binding information
  final_data <- merge(merged_data,
                     binding_df,
                     by.x=c("seqnames", "start", "end"),
                     by.y=c("chr", "start", "end"))
  
  # Write output
  output_file <- paste0("DE_integrated/", output_prefix, "_with_DE.csv")
  write.csv(final_data, output_file, row.names=FALSE)
  
  return(final_data)
}

# Process each dataset
both_results <- process_binding_data(
  "mecp2_cpg_enrichment_parallel_both.csv",
  "mecp2_both"
)

endo_results <- process_binding_data(
  "mecp2_cpg_enrichment_parallel_endo_only.csv",
  "mecp2_endo"
)

exo_results <- process_binding_data(
  "mecp2_cpg_enrichment_parallel_exo_only.csv",
  "mecp2_exo"
)