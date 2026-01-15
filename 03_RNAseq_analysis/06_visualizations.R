#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 7: Generate Visualizations
# This script creates publication-quality plots for RNA-seq analysis

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(viridis)
    library(rtracklayer)
})

cat("==========================================\n")
cat("Starting Visualization Generation\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths - configure gtf_file for your annotation
deseq_dir <- "results/DESeq2"
output_dir <- "results/plots"
gtf_file <- ""  # Set to your GTF annotation file (e.g., "gencode.v44.annotation.gtf")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load data
# ============================================================
cat("Loading DESeq2 objects...\n")
dds <- readRDS(file.path(deseq_dir, "dds.rds"))
vsd <- readRDS(file.path(deseq_dir, "vsd.rds"))
res <- readRDS(file.path(deseq_dir, "results.rds"))

# Load metadata
metadata <- read.csv("metadata/sample_info.csv", row.names = 1)

cat("Data loaded successfully\n\n")

# ============================================================
# 2. Create Gene ID to Symbol Mapping
# ============================================================
cat("Loading gene annotations and creating ID to symbol mapping...\n")

# Read GTF file
gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

# Extract gene_id and gene_name mapping (only for genes)
gene_mapping <- gtf_df %>%
    filter(type == "gene") %>%
    select(gene_id, gene_name) %>%
    distinct() %>%
    mutate(gene_id = gsub("\\..*", "", gene_id))  # Remove version numbers

# Create a named vector for easy lookup
gene_id_to_symbol <- setNames(gene_mapping$gene_name, gene_mapping$gene_id)

# Function to convert gene IDs to symbols
convert_to_symbol <- function(gene_ids) {
    gene_ids_clean <- gsub("\\..*", "", gene_ids)  # Remove version numbers
    symbols <- gene_id_to_symbol[gene_ids_clean]
    # If no symbol found, use the cleaned gene ID
    ifelse(is.na(symbols), gene_ids_clean, symbols)
}

cat("  Created mapping for", length(gene_id_to_symbol), "genes\n\n")

# ============================================================
# 3. PCA Plot
# ============================================================
cat("Generating PCA plot...\n")

# Calculate PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_text_repel(size = 4, box.padding = 0.5, point.padding = 0.3) +
    scale_color_manual(values = c("Control" = "#3498db", "Mutant" = "#e74c3c")) +
    labs(
        title = "PCA Plot - Control vs Mutant NPCs",
        x = paste0("PC1: ", percent_var[1], "% variance"),
        y = paste0("PC2: ", percent_var[2], "% variance"),
        color = "Condition"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "PCA_plot.png"),
    pca_plot,
    width = 10,
    height = 8,
    dpi = 300
)

ggsave(
    file.path(output_dir, "PCA_plot.pdf"),
    pca_plot,
    width = 10,
    height = 8
)

cat("  Saved: PCA_plot.png and PCA_plot.pdf\n\n")

# ============================================================
# 4. Sample Distance Heatmap
# ============================================================
cat("Generating sample distance heatmap...\n")

# Calculate sample distances
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

# Clean row and column names (use sample IDs, not condition names)
rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

# Create annotation
annotation_col <- data.frame(
    Condition = vsd$condition,
    row.names = colnames(vsd)
)

# Color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Create heatmap
png(file.path(output_dir, "sample_distance_heatmap.png"), width = 10, height = 9, units = "in", res = 300)
pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    color = colors,
    annotation_col = annotation_col,
    main = "Sample-to-Sample Distance Heatmap",
    fontsize = 12
)
dev.off()

pdf(file.path(output_dir, "sample_distance_heatmap.pdf"), width = 10, height = 9)
pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    color = colors,
    annotation_col = annotation_col,
    main = "Sample-to-Sample Distance Heatmap",
    fontsize = 12
)
dev.off()

cat("  Saved: sample_distance_heatmap.png and sample_distance_heatmap.pdf\n\n")

# ============================================================
# 5. Volcano Plot
# ============================================================
cat("Generating volcano plot...\n")

# Prepare data
volcano_data <- as.data.frame(res) %>%
    mutate(
        gene = rownames(.),
        significant = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
            padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
            TRUE ~ "Not significant"
        )
    ) %>%
    filter(!is.na(padj), !is.na(log2FoldChange))

# Get top genes to label
top_genes <- volcano_data %>%
    filter(significant != "Not significant") %>%
    arrange(padj) %>%
    head(20)

# Create volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    scale_color_manual(
        values = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db", "Not significant" = "gray70"),
        name = "Expression"
    ) +
    labs(
        title = "Volcano Plot - Mutant vs Control",
        x = "log2 Fold Change",
        y = "-log10 (Adjusted P-value)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "volcano_plot.png"),
    volcano_plot,
    width = 12,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(output_dir, "volcano_plot.pdf"),
    volcano_plot,
    width = 12,
    height = 10
)

cat("  Saved: volcano_plot.png and volcano_plot.pdf\n\n")

# ============================================================
# 5b. Additional Volcano Plots with Different Thresholds
# ============================================================
cat("Generating additional volcano plots with different thresholds...\n")

# Volcano plot with 50 reads threshold and FC > 0.5 (log2FC > 0.585)
cat("  Creating volcano plot: baseMean >= 50, FC > 0.5...\n")
volcano_data_50reads_fc05 <- as.data.frame(res) %>%
    filter(baseMean >= 50) %>%
    mutate(
        gene = rownames(.),
        significant = case_when(
            padj < 0.05 & log2FoldChange > log2(1.5) ~ "Upregulated",
            padj < 0.05 & log2FoldChange < -log2(1.5) ~ "Downregulated",
            TRUE ~ "Not significant"
        )
    ) %>%
    filter(!is.na(padj), !is.na(log2FoldChange))

# Count significant genes
n_up_fc05 <- sum(volcano_data_50reads_fc05$significant == "Upregulated")
n_down_fc05 <- sum(volcano_data_50reads_fc05$significant == "Downregulated")

volcano_plot_50reads_fc05 <- ggplot(volcano_data_50reads_fc05,
                                     aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    scale_color_manual(
        values = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db", "Not significant" = "gray70"),
        name = "Expression"
    ) +
    labs(
        title = sprintf("Volcano Plot - Mutant vs Control\n(baseMean >= 50, FC > 0.5, padj < 0.05)\nUp: %d | Down: %d",
                        n_up_fc05, n_down_fc05),
        x = "log2 Fold Change",
        y = "-log10 (Adjusted P-value)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "volcano_plot_50reads_FC05.png"),
    volcano_plot_50reads_fc05,
    width = 12,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(output_dir, "volcano_plot_50reads_FC05.pdf"),
    volcano_plot_50reads_fc05,
    width = 12,
    height = 10
)

cat("    Saved: volcano_plot_50reads_FC05.png and volcano_plot_50reads_FC05.pdf\n")

# Volcano plot with 50 reads threshold and FC > 1 (log2FC > 1)
cat("  Creating volcano plot: baseMean >= 50, FC > 1...\n")
volcano_data_50reads_fc1 <- as.data.frame(res) %>%
    filter(baseMean >= 50) %>%
    mutate(
        gene = rownames(.),
        significant = case_when(
            padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
            padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
            TRUE ~ "Not significant"
        )
    ) %>%
    filter(!is.na(padj), !is.na(log2FoldChange))

# Count significant genes
n_up_fc1 <- sum(volcano_data_50reads_fc1$significant == "Upregulated")
n_down_fc1 <- sum(volcano_data_50reads_fc1$significant == "Downregulated")

volcano_plot_50reads_fc1 <- ggplot(volcano_data_50reads_fc1,
                                    aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    scale_color_manual(
        values = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db", "Not significant" = "gray70"),
        name = "Expression"
    ) +
    labs(
        title = sprintf("Volcano Plot - Mutant vs Control\n(baseMean >= 50, FC > 1, padj < 0.05)\nUp: %d | Down: %d",
                        n_up_fc1, n_down_fc1),
        x = "log2 Fold Change",
        y = "-log10 (Adjusted P-value)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "volcano_plot_50reads_FC1.png"),
    volcano_plot_50reads_fc1,
    width = 12,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(output_dir, "volcano_plot_50reads_FC1.pdf"),
    volcano_plot_50reads_fc1,
    width = 12,
    height = 10
)

cat("    Saved: volcano_plot_50reads_FC1.png and volcano_plot_50reads_FC1.pdf\n\n")

# ============================================================
# 6. MA Plot
# ============================================================
cat("Generating MA plot...\n")

# Prepare data
ma_data <- as.data.frame(res) %>%
    mutate(
        gene = rownames(.),
        significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not significant")
    ) %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(baseMean))

# Create MA plot
ma_plot <- ggplot(ma_data, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black") +
    scale_color_manual(
        values = c("Significant" = "#e74c3c", "Not significant" = "gray70"),
        name = ""
    ) +
    labs(
        title = "MA Plot - Mutant vs Control",
        x = "log10 (Mean Expression)",
        y = "log2 Fold Change"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "MA_plot.png"),
    ma_plot,
    width = 12,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(output_dir, "MA_plot.pdf"),
    ma_plot,
    width = 12,
    height = 10
)

cat("  Saved: MA_plot.png and MA_plot.pdf\n\n")

# ============================================================
# 7. Heatmap of Top Differentially Expressed Genes
# ============================================================
cat("Generating heatmap of top DEGs...\n")

# Get top 50 significant genes by adjusted p-value
top_genes_list <- res %>%
    as.data.frame() %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(50) %>%
    rownames()

# Get VST counts for these genes
top_genes_vst <- assay(vsd)[top_genes_list, ]

# Scale by row (z-score)
top_genes_scaled <- t(scale(t(top_genes_vst)))

# Convert row names to gene symbols
rownames(top_genes_scaled) <- convert_to_symbol(rownames(top_genes_scaled))

# Create annotation (use sample IDs, not condition names)
annotation_col <- data.frame(
    Condition = vsd$condition,
    row.names = colnames(vsd)
)

# Define colors
ann_colors <- list(
    Condition = c("Control" = "#3498db", "Mutant" = "#e74c3c")
)

# Create heatmap
png(file.path(output_dir, "heatmap_top50_DEGs.png"), width = 10, height = 14, units = "in", res = 300)
pheatmap(
    top_genes_scaled,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    fontsize_row = 6,
    main = "Top 50 Differentially Expressed Genes (Z-score)",
    breaks = seq(-3, 3, length.out = 101)
)
dev.off()

pdf(file.path(output_dir, "heatmap_top50_DEGs.pdf"), width = 10, height = 14)
pheatmap(
    top_genes_scaled,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    fontsize_row = 6,
    main = "Top 50 Differentially Expressed Genes (Z-score)",
    breaks = seq(-3, 3, length.out = 101)
)
dev.off()

cat("  Saved: heatmap_top50_DEGs.png and heatmap_top50_DEGs.pdf\n\n")

# ============================================================
# 8. P-value Distribution
# ============================================================
cat("Generating p-value distribution plot...\n")

pval_data <- as.data.frame(res) %>%
    filter(!is.na(pvalue))

pval_plot <- ggplot(pval_data, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "#3498db", color = "black", alpha = 0.7) +
    labs(
        title = "Distribution of P-values",
        x = "P-value",
        y = "Frequency"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
    )

ggsave(
    file.path(output_dir, "pvalue_distribution.png"),
    pval_plot,
    width = 10,
    height = 8,
    dpi = 300
)

ggsave(
    file.path(output_dir, "pvalue_distribution.pdf"),
    pval_plot,
    width = 10,
    height = 8
)

cat("  Saved: pvalue_distribution.png and pvalue_distribution.pdf\n\n")

# ============================================================
# 9. Gene Expression Boxplots for Top Genes
# ============================================================
cat("Generating boxplots for top 9 DEGs...\n")

# Get top 9 genes
top9_genes <- res %>%
    as.data.frame() %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(9) %>%
    rownames()

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Prepare data for plotting
boxplot_data <- norm_counts[top9_genes, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
    left_join(
        rownames_to_column(metadata, "sample"),
        by = "sample"
    ) %>%
    mutate(gene = convert_to_symbol(gene))  # Convert to gene symbols

# Create boxplot
boxplot_top9 <- ggplot(boxplot_data, aes(x = condition, y = log2(count + 1), fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Control" = "#3498db", "Mutant" = "#e74c3c")) +
    labs(
        title = "Expression of Top 9 Differentially Expressed Genes",
        x = "Condition",
        y = "log2(Normalized Count + 1)",
        fill = "Condition"
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
    )

ggsave(
    file.path(output_dir, "boxplot_top9_DEGs.png"),
    boxplot_top9,
    width = 14,
    height = 12,
    dpi = 300
)

ggsave(
    file.path(output_dir, "boxplot_top9_DEGs.pdf"),
    boxplot_top9,
    width = 14,
    height = 12
)

cat("  Saved: boxplot_top9_DEGs.png and boxplot_top9_DEGs.pdf\n\n")

# ============================================================
# 10. Summary Statistics Plot
# ============================================================
cat("Generating summary statistics plot...\n")

# Prepare summary data
summary_data <- data.frame(
    Category = c("Total Genes Tested", "Significant (padj < 0.05)",
                 "Upregulated (FC > 2)", "Downregulated (FC < 0.5)"),
    Count = c(
        sum(!is.na(res$padj)),
        sum(res$padj < 0.05, na.rm = TRUE),
        sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE),
        sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
    )
)

summary_plot <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count)) +
    geom_col(fill = "#3498db", alpha = 0.8) +
    geom_text(aes(label = Count), hjust = -0.2, size = 5) +
    coord_flip() +
    labs(
        title = "Differential Expression Analysis Summary",
        x = "",
        y = "Number of Genes"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave(
    file.path(output_dir, "summary_statistics.png"),
    summary_plot,
    width = 10,
    height = 6,
    dpi = 300
)

ggsave(
    file.path(output_dir, "summary_statistics.pdf"),
    summary_plot,
    width = 10,
    height = 6
)

cat("  Saved: summary_statistics.png and summary_statistics.pdf\n\n")

cat("\nVisualization generation complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nGenerated plots:\n")
cat("  - PCA_plot: Principal Component Analysis\n")
cat("  - sample_distance_heatmap: Sample clustering\n")
cat("  - volcano_plot: Differential expression volcano plot (original)\n")
cat("  - volcano_plot_50reads_FC05: Volcano plot (baseMean >= 50, FC > 0.5)\n")
cat("  - volcano_plot_50reads_FC1: Volcano plot (baseMean >= 50, FC > 1)\n")
cat("  - MA_plot: Mean vs log2FC plot\n")
cat("  - heatmap_top50_DEGs: Heatmap of top 50 genes\n")
cat("  - boxplot_top9_DEGs: Expression boxplots for top 9 genes\n")
cat("  - pvalue_distribution: P-value histogram\n")
cat("  - summary_statistics: Analysis summary bar chart\n")
cat("\nAll plots saved in both PNG (high-res) and PDF formats\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
