# Script to generate volcano plots from DESeq2 results
# with specific thresholds: 50 reads, FC > 0.5, FC > 1
# Uses the same styling as the original generate_volcano_plots.R

library(tidyverse)
library(ggrepel)

# --- Configuration ---
RESULTS_FILE <- "results/DESeq2/DESeq2_results_all.csv"
OUTPUT_DIR <- "results/volcano_plots_filtered"
P_THRESHOLD <- 0.05
BASEMEAN_THRESHOLD <- 50 # Minimum baseMean threshold
GENE_TO_LABEL <- "MECP2" # Gene to highlight (if desired)

# Threshold configurations
THRESHOLDS <- list(
    list(name = "FC05", fc = 0.5, description = "FC > 0.5"),
    list(name = "FC1", fc = 1, description = "FC > 1")
)

# --- Function Definition ---

#' Create a Volcano Plot with Specific Styling
#'
#' Generates a volcano plot with specific styling similar to the original script,
#' with filtering for baseMean and fold change thresholds.
#'
#' @param results_df Data frame containing differential expression results.
#'        Must include columns: 'log2FoldChange', 'padj', 'gene_symbol', 'baseMean'.
#' @param output_file Path where the PDF plot will be saved.
#' @param p_threshold Adjusted p-value significance threshold.
#' @param fc_threshold Absolute log2 fold change significance threshold.
#' @param basemean_threshold Minimum baseMean value for filtering.
#' @param gene_to_label The specific gene symbol to label (optional).
#' @param exclude_target_gene Logical, whether to filter out the gene_to_label before plotting.
#' @param plot_title Title for the plot.
#'
#' @return The ggplot object.
create_styled_volcano <- function(results_df, output_file,
                                   p_threshold, fc_threshold, basemean_threshold,
                                   gene_to_label = NULL,
                                   exclude_target_gene = TRUE,
                                   plot_title = "Differential Expression (Mutant vs Control)") {

    # Ensure required columns exist
    required_cols <- c("gene_symbol", "baseMean", "log2FoldChange", "padj")
    missing_cols <- setdiff(required_cols, colnames(results_df))
    if (length(missing_cols) > 0) {
        stop(paste("Results data frame must contain columns:", paste(missing_cols, collapse = ", ")))
    }

    # --- Filter by baseMean threshold ---
    results_df <- results_df %>% filter(baseMean >= basemean_threshold)
    message(paste("  Filtered to", nrow(results_df), "genes with baseMean >=", basemean_threshold))

    # --- Store original data before potential filtering ---
    original_results_df <- results_df

    # --- Conditionally filter out the specific gene ---
    label_data <- NULL
    if (!is.null(gene_to_label)) {
        if (exclude_target_gene) {
            results_df <- results_df %>% filter(gene_symbol != gene_to_label)
            message(paste("  Filtered out", gene_to_label, "for plot."))
        } else {
            # If including the gene, prepare its data for labeling later
            label_data <- original_results_df %>% filter(gene_symbol == gene_to_label)
            if (nrow(label_data) == 0) {
                warning(paste("Gene to label", gene_to_label, "not found in the results"))
            }
        }
    }

    # Remove rows with NA p-values which cannot be plotted
    results_df <- results_df %>% filter(!is.na(padj))

    # Define significance categories based on original script
    results_df <- results_df %>%
        mutate(
            significance = case_when(
                padj < p_threshold & abs(log2FoldChange) >= fc_threshold ~ "p-value and log2 FC",
                padj < p_threshold & abs(log2FoldChange) < fc_threshold ~ "p-value",
                padj >= p_threshold & abs(log2FoldChange) >= fc_threshold ~ "Log2 FC",
                TRUE ~ "NS" # Not Significant
            ),
            # Ensure significance is a factor with the desired order for the legend
            significance = factor(significance, levels = c("NS", "Log2 FC", "p-value", "p-value and log2 FC"))
        )

    # Count significant genes
    n_sig <- sum(results_df$significance == "p-value and log2 FC")
    n_up <- sum(results_df$significance == "p-value and log2 FC" & results_df$log2FoldChange > 0)
    n_down <- sum(results_df$significance == "p-value and log2 FC" & results_df$log2FoldChange < 0)

    message(sprintf("  Significant genes (p < %.2f & |log2FC| > %.2f): %d (Up: %d, Down: %d)",
                    p_threshold, fc_threshold, n_sig, n_up, n_down))

    # Define colors similar to the original script
    color_map <- c(
        "NS" = "#404040",                  # Very dark grey for non-significant
        "Log2 FC" = "green",               # Intense Green
        "p-value" = "blue",                # Intense Blue
        "p-value and log2 FC" = "red"      # Intense Red
    )

    # --- Calculate y-axis limit ---
    if (is.null(gene_to_label) || exclude_target_gene) {
        # Case 1: Target gene EXCLUDED or no target - Focus on other significant genes
        y_limit_data_df <- results_df # Use the filtered data

        y_limit_data <- y_limit_data_df %>%
            filter(significance != "NS", padj > 0) %>% # Exclude p=0 for quantile calculation
            pull(padj)

        # Handle case where there are few/no other significant genes
        if (length(y_limit_data) < 10) {
            # Fallback: use max -log10(p) of remaining genes, or 10 if none exist
            max_y_base <- max(c(10, -log10(results_df$padj[results_df$padj > 0])), na.rm = TRUE)
        } else {
            # Use 99th percentile of significant genes
            max_y_base <- -log10(quantile(y_limit_data, 0.01, na.rm = TRUE))
        }

        # Add buffer
        max_y <- max_y_base * 1.1

        # Handle Inf cases
        if (is.infinite(max_y)) {
            finite_y_values <- -log10(results_df$padj[results_df$padj > 0])
            max_y <- if (length(finite_y_values) > 0) max(finite_y_values, na.rm = TRUE) * 1.2 else 350
        }
        # Ensure minimum limit and apply cap
        max_y <- max(max_y, 10)
        max_y <- min(max_y, 75) # Apply cap

    } else {
        # Case 2: Target gene INCLUDED - Ensure it's visible
        y_vals_all <- -log10(original_results_df$padj)
        y_vals_all <- y_vals_all[!is.na(y_vals_all) & is.finite(y_vals_all)] # Remove NA/Inf

        if (length(y_vals_all) == 0) {
            max_y_base <- 10 # Default if no valid p-values
        } else {
            max_y_base <- max(y_vals_all, na.rm = TRUE)
        }

        # Add buffer
        max_y <- max_y_base * 1.15

        # Ensure minimum limit
        max_y <- max(max_y, 10)
    }

    # Handle potential remaining Inf case for max_y
    if (is.infinite(max_y)) {
        max_y <- 350 # Set arbitrary high limit if calculation resulted in Inf
    }

    # Create the plot
    p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj)))

    # Add layers sequentially
    p <- p +
        geom_point(aes(color = significance), alpha = 0.6, size = 1.5) # Plot points

    p <- p +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey50") # Significance lines

    p <- p +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "grey50")

    p <- p +
        scale_color_manual(values = color_map) # Set colors

    p <- p +
        theme_bw(base_size = 14) # Set theme

    p <- p +
        theme(
            legend.position = "top",
            legend.title = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    p <- p +
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~italic(P)),
            title = plot_title
        )

    # Use coord_cartesian to zoom without removing data points outside limits
    plot_data_for_lims <- if (is.null(gene_to_label) || exclude_target_gene) results_df else original_results_df
    x_limits <- c(
        min(plot_data_for_lims$log2FoldChange, na.rm = TRUE) - 0.5,
        max(plot_data_for_lims$log2FoldChange, na.rm = TRUE) + 0.5
    )

    p <- p + coord_cartesian(
        ylim = c(0, max_y),
        xlim = x_limits
    )

    # --- Add the specific gene label using ggrepel ONLY if included ---
    if (!is.null(gene_to_label) && !exclude_target_gene && !is.null(label_data) && nrow(label_data) > 0) {
        p <- p + geom_text_repel(
            data = label_data,
            aes(label = gene_symbol),
            size = 4,
            fontface = "bold",
            color = "black",
            nudge_y = (max_y * 0.02),
            nudge_x = (x_limits[2] - x_limits[1]) * 0.05,
            box.padding = 0.5,
            point.padding = 0.5,
            segment.color = 'grey50',
            max.overlaps = Inf
        )
    }

    # Save the plot
    ggsave(output_file, p, width = 8, height = 7, device = "pdf")

    # Also save PNG version
    png_file <- sub("\\.pdf$", ".png", output_file)
    ggsave(png_file, p, width = 8, height = 7, dpi = 300)

    message(paste("  Saved volcano plot to:", output_file))
    message(paste("  Saved volcano plot to:", png_file))

    return(p)
}

# --- Main Script Logic ---

cat("==========================================\n")
cat("Generating Volcano Plots with Thresholds\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Check if results file exists
if (!file.exists(RESULTS_FILE)) {
    stop(paste("Results file not found:", RESULTS_FILE))
}

cat("Reading DESeq2 results from:", RESULTS_FILE, "\n")
res_data <- read.csv(RESULTS_FILE, row.names = 1)
cat("  Total genes:", nrow(res_data), "\n\n")

# Generate plots for each threshold configuration
for (threshold in THRESHOLDS) {
    cat("==========================================\n")
    cat(sprintf("Processing: %s (log2FC threshold: %.2f)\n", threshold$description, threshold$fc))
    cat("==========================================\n")

    # --- Plot 1: Excluding target gene (MECP2) ---
    output_file_excluded <- file.path(
        OUTPUT_DIR,
        sprintf("volcano_plot_%s_50reads_noMECP2.pdf", threshold$name)
    )

    cat("Generating plot EXCLUDING", GENE_TO_LABEL, "\n")
    create_styled_volcano(
        results_df = res_data,
        output_file = output_file_excluded,
        p_threshold = P_THRESHOLD,
        fc_threshold = threshold$fc,
        basemean_threshold = BASEMEAN_THRESHOLD,
        gene_to_label = GENE_TO_LABEL,
        exclude_target_gene = TRUE,
        plot_title = sprintf(
            "Differential Expression (Mutant vs Control)\nbaseMean >= %d, %s, padj < %.2f",
            BASEMEAN_THRESHOLD, threshold$description, P_THRESHOLD
        )
    )
    cat("\n")

    # --- Plot 2: Including target gene (MECP2) ---
    output_file_included <- file.path(
        OUTPUT_DIR,
        sprintf("volcano_plot_%s_50reads_withMECP2.pdf", threshold$name)
    )

    cat("Generating plot INCLUDING", GENE_TO_LABEL, "\n")
    create_styled_volcano(
        results_df = res_data,
        output_file = output_file_included,
        p_threshold = P_THRESHOLD,
        fc_threshold = threshold$fc,
        basemean_threshold = BASEMEAN_THRESHOLD,
        gene_to_label = GENE_TO_LABEL,
        exclude_target_gene = FALSE,
        plot_title = sprintf(
            "Differential Expression (Mutant vs Control)\nbaseMean >= %d, %s, padj < %.2f",
            BASEMEAN_THRESHOLD, threshold$description, P_THRESHOLD
        )
    )
    cat("\n")
}

cat("==========================================\n")
cat("Volcano plot generation complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
