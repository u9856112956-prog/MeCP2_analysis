library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(gridExtra)

# Initialize txdb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Parse command line arguments
library(optparse)

# Create option parser with explicit dest parameter
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory path"),
    make_option("--data-to-analyze",
                type="character",
                default=NULL,
                dest="data_file",
                help="CSV file to analyze")
)

######################################### Function definitions first #########################################
# Analyze peak distribution relative to TSS
analyze_peaks <- function(peaks, name) {
    # Annotate peaks
    peakAnno <- annotatePeak(peaks, 
                            tssRegion = c(-3000, 3000),  # Define TSS region as ±3kb
                            TxDb = txdb,
                            annoDb = "org.Mm.eg.db")
    
    # Plot peak distribution
    pdf(paste0("peaks_annotation/", name, "_distribution.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    # Return detailed annotation
    return(peakAnno)
}

# Extract peaks near TSS
get_tss_peaks <- function(peaks, anno) {
    anno_df <- as.data.frame(anno@anno)
    tss_peaks <- anno_df %>%
        filter(abs(distanceToTSS) <= 3000) %>%
        select(seqnames, start, end, distanceToTSS, geneId, SYMBOL)
    return(tss_peaks)
}

# Function to read peak file
read_peak_file <- function(file_path) {
    # Read CSV file
    peaks_df <- read.csv(file_path)
    
    # Convert to GRanges object
    gr <- GRanges(
        seqnames = peaks_df$chr,
        ranges = IRanges(
            start = peaks_df$start,
            end = peaks_df$end
        ),
        # Add metadata columns that might be useful for analysis
        exo_signal = peaks_df$exo_signal,
        endo_signal = peaks_df$endo_signal,
        enrichment = peaks_df$enrichment,
        pvalue = peaks_df$pvalue,
        binding_type = peaks_df$binding_type,
        binding_type_by_peaks = peaks_df$binding_type_by_peaks,
        significant = peaks_df$significant,
        cpg_score = peaks_df$cpg_score,
        cpg_name = peaks_df$cpg_name,
        region_length = peaks_df$region_length
    )
    
    return(gr)
}

# Enhanced genomic annotation analysis
analyze_genomic_regions <- function(peaks, name) {
    # Create output directory for plots
    dir.create("peaks_annotation/", showWarnings = FALSE)
    
    # 1. Detailed genomic annotation
    detailed_anno <- annotatePeak(peaks,
                                TxDb = txdb,
                                annoDb = "org.Mm.eg.db",
                                tssRegion = c(-3000, 3000),
                                level = "transcript",
                                verbose = FALSE)
    
    # 2. Generate detailed statistics and map annotations
    anno_stats <- as.data.frame(detailed_anno@anno) %>%
        group_by(annotation) %>%
        summarise(
            count = n(),
            .groups = 'drop' # Drop grouping for further operations
        ) %>%
        mutate(
            percentage = count / sum(count) * 100,
            # Map detailed annotations to broader categories for plotting
            PlotCategory = case_when(
                grepl("Promoter", annotation) ~ "Promoter",
                grepl("UTR|Exon|Intron|Downstream", annotation) ~ "Gene Body",
                grepl("Intergenic", annotation) ~ "Intergenic",
                TRUE ~ "Other" # Fallback category
            ),
            # Ensure PlotCategory is a factor for consistent ordering in plot
            PlotCategory = factor(PlotCategory, levels = c("Promoter", "Gene Body", "Intergenic", "Other"))
        ) %>%
        # Summarize again by the broader PlotCategory
        group_by(PlotCategory) %>%
        summarise(
            TotalCount = sum(count),
            TotalPercentage = sum(percentage),
            .groups = 'drop'
        )

    # 3. Create and save the pie chart using ggplot2
    
    # Define colors similar to the reference script
    plot_colors <- c(
        "Promoter" = "#3182bd",    # Blue for Promoter regions
        "Gene Body" = "#fd8d3c",   # Orange for Gene Body (UTRs, Exons, Introns, Downstream)
        "Intergenic" = "#e7969c",  # Pinkish Red for Intergenic
        "Other" = "#bdbdbd"        # Grey for any other category
    )

    # Create the pie chart
    pie_chart <- ggplot(anno_stats, aes(x = "", y = TotalPercentage, fill = PlotCategory)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = plot_colors) +
        theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5) # Center title
        ) +
        labs(
            title = paste("Genomic Distribution of Peaks:", name),
            fill = "Genomic Region"
        )

    # Save the plot
    ggsave(paste0("peaks_annotation/", name, "_pie_chart.pdf"),
           plot = pie_chart, width = 8, height = 6)
    ggsave(paste0("peaks_annotation/", name, "_pie_chart.png"),
           plot = pie_chart, width = 8, height = 6, dpi = 300)

    # 4. Get promoter vs non-promoter stats (using original detailed stats)
    promoter_stats_detailed <- as.data.frame(detailed_anno@anno) %>%
        group_by(annotation) %>%
        summarise(count = n(), .groups = 'drop')
    
    # 4. Get promoter vs non-promoter stats
    promoter_stats <- data.frame(
        region = c("Promoter", "Non-promoter"),
        count = c(
            sum(promoter_stats_detailed$count[grep("Promoter", promoter_stats_detailed$annotation)]),
            sum(promoter_stats_detailed$count[!grepl("Promoter", promoter_stats_detailed$annotation)])
        )
    )
    promoter_stats$percentage <- promoter_stats$count / sum(promoter_stats$count) * 100
    
    # Return statistics
    return(list(
        detailed_stats = anno_stats,
        promoter_stats = promoter_stats,
        full_annotation = detailed_anno
    ))
}

# Plot distribution of distances to TSS
plot_tss_dist <- function(anno, name) {
  distances <- as.data.frame(anno@anno) %>%
    select(distanceToTSS) %>%
    mutate(Category = name)
  
  ggplot(distances, aes(x = distanceToTSS, fill = Category)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    xlim(-3000, 3000) +
    labs(x = "Distance to TSS (bp)",
         y = "Density",
         title = "Distribution of Peak Distances to TSS")
  
  ggsave(paste0("peaks_annotation/", name, "_tss_distance_distribution.pdf"), 
         width = 10, height = 6)
}

# Modified feature enrichment analysis
calculate_feature_enrichment <- function(peaks, name) {
    # Create output directory if it doesn't exist
    dir.create("peaks_annotation/", showWarnings = FALSE, recursive = TRUE)
    
    # Get peak annotations
    peak_anno <- annotatePeak(peaks,
                            TxDb = txdb,
                            annoDb = "org.Mm.eg.db",
                            tssRegion = c(-3000, 3000))
    
    # Extract gene IDs from peaks
    peak_genes <- unique(peak_anno@anno$geneId)
    
    # Get all genes as background
    all_genes <- keys(org.Mm.eg.db)
    
    # Perform GO enrichment analysis
    ego <- enrichGO(gene = peak_genes,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
    
    # Plot results
    if (!is.null(ego) && nrow(ego) > 0) {
        pdf(paste0("peaks_annotation/", name, "_GO_enrichment.pdf"))
        print(dotplot(ego, showCategory = 15, title = paste(name, "GO Enrichment")))
        dev.off()
        
        # Save results to CSV
        write.csv(as.data.frame(ego), 
                 file = paste0("peaks_annotation/", name, "_GO_enrichment.csv"))
    }
    
    return(ego)
}

#############################################################################################################
# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Working directory path:")
print(opts$work_dir)
print("Data file to analyze:")
print(opts$data_file)

# Check if required arguments are provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}
if (is.null(opts$data_file) || !is.character(opts$data_file)) {
    stop("Invalid or missing data file path")
}

# Set working directory
tryCatch({
    setwd(opts$work_dir)
}, error = function(e) {
    stop(paste("Error setting working directory:", e$message))
})

# Extract the base name without extension for output file naming
file_base <- tools::file_path_sans_ext(opts$data_file)

# Read peak file
file_path <- file.path("lists", opts$data_file)
message("Reading file: ", file_path)
peaks <- read_peak_file(file_path)

# Create peaks_annotation directory if it doesn't exist
dir.create("peaks_annotation", showWarnings = FALSE)

# Analyze peaks
peaks_anno <- analyze_peaks(peaks, file_base)

# Get TSS-proximal peaks
tss_peaks <- get_tss_peaks(peaks, peaks_anno)

# Save TSS-proximal peaks
write_csv(tss_peaks, paste0("peaks_annotation/", file_base, "_tss_peaks.csv"))

# Generate summary statistics
summary_stats <- data.frame(
    Category = file_base,
    Total_Peaks = length(peaks),
    TSS_Proximal = nrow(tss_peaks)
) %>%
    mutate(Percent_TSS = (TSS_Proximal / Total_Peaks) * 100)

write_csv(summary_stats, paste0("peaks_annotation/", file_base, "_peak_location_summary.csv"))

# Enhanced genomic annotation analysis
peaks_analysis <- analyze_genomic_regions(peaks, file_base)

# Save detailed statistics
write_csv(peaks_analysis$detailed_stats, 
         paste0("peaks_annotation/", file_base, "_detailed_stats.csv"))

# Calculate feature enrichment
peaks_enrichment <- calculate_feature_enrichment(peaks, file_base)