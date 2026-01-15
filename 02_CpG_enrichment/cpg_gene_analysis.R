library(tidyverse)
library(optparse)
library(GenomicRanges)
library(org.Mm.eg.db)
library(clusterProfiler)

# Add package installation checks and installation if needed
required_packages <- c("VennDiagram", "tidyverse", "optparse", "GenomicRanges", 
                      "org.Mm.eg.db", "clusterProfiler")

for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        message(sprintf("Installing package: %s", package))
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install(package, update = FALSE)
    }
    # Load package with error handling
    tryCatch({
        library(package, character.only = TRUE)
    }, error = function(e) {
        message(sprintf("Warning: Failed to load package %s: %s", package, e$message))
    })
}

######################################### Function definitions #########################################

# Function to read and process gene data with error handling
process_gene_data <- function(file_path, category) {
    tryCatch({
        if (!file.exists(file_path)) {
            warning(sprintf("File not found: %s", file_path))
            return(NULL)
        }
        
        data <- read_csv(file_path) %>%
            mutate(category = category) %>%
            rename(chr = seqnames)
        
        # Convert ENTREZ IDs to gene symbols with error handling
        genes_with_symbols <- data %>%
            filter(!is.na(geneId)) %>%
            mutate(
                geneId = as.character(geneId),
                category = category
            )
        
        return(genes_with_symbols)
    }, error = function(e) {
        warning(sprintf("Error processing file %s: %s", file_path, e$message))
        return(NULL)
    })
}

# Function to create Venn diagram
create_venn_diagram <- function(venn_data, output_file) {
    tryCatch({
        venn.plot <- venn.diagram(
            venn_data,
            filename = output_file,
            imagetype = "png",
            height = 3000,
            width = 3000,
            resolution = 300,
            compression = "lzw",
            lwd = 2,
            col = c("#440154FF", "#21908CFF", "#FDE725FF"),
            fill = c(alpha("#440154FF", 0.3), 
                    alpha("#21908CFF", 0.3), 
                    alpha("#FDE725FF", 0.3)),
            main = "Gene Overlaps Between Categories"
        )
        return(TRUE)
    }, error = function(e) {
        warning(sprintf("Failed to create Venn diagram: %s", e$message))
        return(FALSE)
    })
}

# Function to calculate enrichment scores
calculate_enrichment_scores <- function(peaks_file) {
    tryCatch({
        peaks <- read_csv(peaks_file)
        
        enrichment_scores <- peaks %>%
            select(chr, start, end, exo_signal, endo_signal, enrichment) %>%
            mutate(
                log2_enrichment = log2(enrichment),
                normalized_enrichment = (exo_signal - endo_signal) / (exo_signal + endo_signal)
            )
        
        return(enrichment_scores)
    }, error = function(e) {
        warning(sprintf("Error calculating enrichment scores: %s", e$message))
        return(NULL)
    })
}

# Function to create enrichment plots
create_enrichment_plots <- function(data, output_prefix) {
    tryCatch({
        # Distance vs enrichment plot
        p1 <- ggplot(data, aes(x = abs(distanceToTSS), y = normalized_enrichment)) +
            geom_point(alpha = 0.5, color = "#21908CFF") +
            geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#FDE725FF") +
            theme_minimal() +
            scale_x_continuous(trans = "log1p") +
            labs(title = "Enrichment Score vs Distance to TSS",
                 x = "Distance to TSS (bp)",
                 y = "Normalized Enrichment Score")
        
        ggsave(paste0(output_prefix, "_distance_vs_enrichment.pdf"), p1, width = 10, height = 6)
        
        # Enrichment score distribution
        p2 <- ggplot(data, aes(x = normalized_enrichment)) +
            geom_histogram(bins = 30, fill = "#21908CFF", alpha = 0.7, boundary = 0) +
            theme_minimal() +
            labs(title = "Distribution of Normalized Enrichment Scores",
                 x = "Normalized Enrichment Score",
                 y = "Count")
        
        ggsave(paste0(output_prefix, "_enrichment_distribution.pdf"), p2, width = 10, height = 6)
        
        return(TRUE)
    }, error = function(e) {
        warning(sprintf("Error creating enrichment plots: %s", e$message))
        return(FALSE)
    })
}

######################################### Main script #########################################

# Parse command line arguments
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory containing the analysis results")
)

# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
message("Command line arguments received:")
print(opts)
message("Working directory path:")
print(opts$work_dir)

# Check if work-dir argument is provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}

# Set working directory
message(sprintf("Setting working directory to: %s", opts$work_dir))
setwd(opts$work_dir)

# Create output directory for gene analysis
dir.create("gene_analysis", showWarnings = FALSE)

# Read the TSS-proximal peaks data
message("Reading TSS-proximal peaks data...")
tss_peaks_file <- "peaks_annotation/tss_peaks.csv"

if (!file.exists(tss_peaks_file)) {
    stop(sprintf("TSS peaks file not found: %s", tss_peaks_file))
}

genes_data <- process_gene_data(tss_peaks_file, basename(opts$work_dir))

if (is.null(genes_data)) {
    stop("Failed to process gene data")
}

# Calculate basic statistics
message("Calculating gene statistics...")
gene_stats <- genes_data %>%
    group_by(category) %>%
    summarise(
        unique_genes = n_distinct(geneId),
        mean_distance = mean(abs(distanceToTSS)),
        median_distance = median(abs(distanceToTSS))
    )

# Save gene statistics
write_csv(gene_stats, "gene_analysis/gene_statistics.csv")

# Calculate enrichment scores
message("Calculating enrichment scores...")
enrichment_scores <- calculate_enrichment_scores("lists/peaks.csv")

if (!is.null(enrichment_scores)) {
    # Combine enrichment scores with gene information
    genes_with_enrichment <- genes_data %>%
        left_join(
            enrichment_scores,
            by = c("chr", "start", "end")
        ) %>%
        arrange(desc(normalized_enrichment))
    
    # Save detailed enrichment analysis
    write_csv(genes_with_enrichment, "gene_analysis/genes_with_enrichment_scores.csv")
    
    # Create enrichment plots
    message("Creating enrichment plots...")
    create_enrichment_plots(genes_with_enrichment, "gene_analysis/enrichment")
    
    # Generate summary statistics for high-enrichment genes
    high_enrichment_threshold <- quantile(genes_with_enrichment$normalized_enrichment, 0.75)
    high_enrichment_genes <- genes_with_enrichment %>%
        filter(normalized_enrichment >= high_enrichment_threshold) %>%
        arrange(desc(normalized_enrichment))
    
    # Save high-enrichment genes
    write_csv(high_enrichment_genes, "gene_analysis/high_enrichment_genes.csv")
}

# Create summary report
message("Generating summary report...")
report <- c(
    "Gene Analysis Summary",
    "==================",
    "",
    paste("Analysis directory:", basename(opts$work_dir)),
    paste("Total genes analyzed:", nrow(genes_data)),
    "",
    "Gene Statistics:",
    capture.output(print(gene_stats)),
    "",
    if (!is.null(enrichment_scores)) {
        c(
            "Enrichment Score Statistics:",
            paste("Median enrichment score:", median(genes_with_enrichment$normalized_enrichment)),
            paste("Mean enrichment score:", mean(genes_with_enrichment$normalized_enrichment)),
            paste("High-enrichment threshold:", high_enrichment_threshold),
            paste("Number of high-enrichment genes:", nrow(high_enrichment_genes))
        )
    }
)

# Save summary report
writeLines(report, "gene_analysis/analysis_summary.txt")

message("Gene analysis complete. Results saved in gene_analysis directory.")