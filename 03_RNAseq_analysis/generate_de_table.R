#!/usr/bin/env Rscript

# Generate Complete Differential Expression Table
# This script creates a comprehensive publication-ready table with all detected genes
# and explicit fold change direction annotations

# Using base R - no packages needed

cat("==========================================\n")
cat("Generating Complete Differential Expression Table\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths
input_file <- "results/DESeq2/DESeq2_results_all.csv"
output_dir <- "results/DESeq2"

# ============================================================
# 1. Load DESeq2 results
# ============================================================
cat("Loading DESeq2 results...\n")
results <- read.csv(input_file, stringsAsFactors = FALSE)

cat("Total genes in analysis:", nrow(results), "\n\n")

# ============================================================
# 2. Add comprehensive annotations
# ============================================================
cat("Adding comprehensive annotations...\n")

# Add fold change in non-log scale
results$foldChange <- 2^results$log2FoldChange

# Regulation direction (explicit)
results$regulation <- ifelse(is.na(results$padj), "Not tested",
                      ifelse(results$padj >= 0.05, "Not significant",
                      ifelse(results$log2FoldChange > 0, "Upregulated in Mutant",
                      ifelse(results$log2FoldChange < 0, "Downregulated in Mutant", "No change"))))

# Simplified regulation direction
results$direction <- ifelse(is.na(results$padj), "NA",
                     ifelse(results$padj >= 0.05, "NS",
                     ifelse(results$log2FoldChange > 0, "UP",
                     ifelse(results$log2FoldChange < 0, "DOWN", "NC"))))

# Significance category
results$significance <- ifelse(is.na(results$padj), "Not tested",
                        ifelse(results$padj < 0.001, "p < 0.001",
                        ifelse(results$padj < 0.01, "p < 0.01",
                        ifelse(results$padj < 0.05, "p < 0.05", "Not significant"))))

# Fold change category
results$FC_category <- ifelse(is.na(results$log2FoldChange), "NA",
                       ifelse(abs(results$log2FoldChange) < 0.5, "|FC| < 1.4",
                       ifelse(abs(results$log2FoldChange) < 1, "|FC| < 2",
                       ifelse(abs(results$log2FoldChange) < 2, "|FC| < 4", "|FC| >= 4"))))

# Combined significance and FC category
results$DE_category <- ifelse(is.na(results$padj), "Not tested",
                       ifelse(results$padj >= 0.05, "Not significant",
                       ifelse(results$padj < 0.05 & abs(results$log2FoldChange) < 1, "Significant (FC < 2)",
                       ifelse(results$padj < 0.05 & abs(results$log2FoldChange) >= 1, "Significant (FC >= 2)", "Other"))))

# Reorder columns for better readability
results <- results[, c(
    "gene_id", "gene_symbol", "gene_type",
    "baseMean",
    "log2FoldChange", "foldChange", "direction", "regulation",
    "lfcSE", "stat", "pvalue", "padj",
    "significance", "FC_category", "DE_category"
)]

# ============================================================
# 3. Generate summary statistics
# ============================================================
cat("\nGenerating summary statistics...\n")

# Overall summary
total_genes <- nrow(results)
tested_genes <- sum(!is.na(results$padj))
sig_genes <- sum(results$padj < 0.05, na.rm = TRUE)
sig_up <- sum(results$direction == "UP", na.rm = TRUE)
sig_down <- sum(results$direction == "DOWN", na.rm = TRUE)
sig_fc2 <- sum(results$DE_category == "Significant (FC >= 2)", na.rm = TRUE)

# Summary by gene type
results_tested <- results[!is.na(results$padj), ]
gene_types <- unique(results_tested$gene_type)

gene_type_list <- lapply(gene_types, function(gt) {
    subset_data <- results_tested[results_tested$gene_type == gt, ]
    data.frame(
        gene_type = gt,
        total = nrow(subset_data),
        significant = sum(subset_data$padj < 0.05, na.rm = TRUE),
        upregulated = sum(subset_data$direction == "UP", na.rm = TRUE),
        downregulated = sum(subset_data$direction == "DOWN", na.rm = TRUE),
        stringsAsFactors = FALSE
    )
})

gene_type_summary <- do.call(rbind, gene_type_list)
gene_type_summary$pct_significant <- round(100 * gene_type_summary$significant / gene_type_summary$total, 2)
gene_type_summary <- gene_type_summary[order(-gene_type_summary$total), ]

# ============================================================
# 4. Save results
# ============================================================
cat("Saving complete differential expression table...\n")

# Main comprehensive table
write.csv(
    results,
    file = file.path(output_dir, "Complete_DESeq2_results_annotated.csv"),
    row.names = FALSE
)

# Gene type summary
write.csv(
    gene_type_summary,
    file = file.path(output_dir, "DE_summary_by_gene_type.csv"),
    row.names = FALSE
)

# Create a text summary report
summary_text <- paste0(
    "DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY\n",
    "=========================================\n\n",
    "Comparison: Mutant vs Control NPCs\n\n",
    "OVERALL STATISTICS:\n",
    "-------------------\n",
    "Total genes detected: ", format(total_genes, big.mark = ","), "\n",
    "Genes tested: ", format(tested_genes, big.mark = ","), "\n",
    "Genes with non-zero counts: ", format(tested_genes, big.mark = ","), "\n\n",
    "SIGNIFICANT GENES (padj < 0.05):\n",
    "--------------------------------\n",
    "Total significant: ", format(sig_genes, big.mark = ","), " (",
        round(100 * sig_genes / tested_genes, 2), "%)\n",
    "  - Upregulated in Mutant: ", format(sig_up, big.mark = ","), " (",
        round(100 * sig_up / sig_genes, 2), "%)\n",
    "  - Downregulated in Mutant: ", format(sig_down, big.mark = ","), " (",
        round(100 * sig_down / sig_genes, 2), "%)\n\n",
    "HIGHLY SIGNIFICANT (padj < 0.05, |log2FC| >= 1):\n",
    "-------------------------------------------------\n",
    "Total: ", format(sig_fc2, big.mark = ","), " genes\n\n",
    "TOP GENE TYPES:\n",
    "---------------\n"
)

# Add top gene types to summary
for (i in 1:min(10, nrow(gene_type_summary))) {
    gt <- gene_type_summary[i, ]
    summary_text <- paste0(
        summary_text,
        sprintf("%-30s: %6d total, %6d significant (%5.1f%%)\n",
                gt$gene_type, gt$total, gt$significant, gt$pct_significant)
    )
}

summary_text <- paste0(
    summary_text,
    "\n\nOUTPUT FILES:\n",
    "-------------\n",
    "1. Complete_DESeq2_results_annotated.csv\n",
    "   - All ", format(total_genes, big.mark = ","), " detected genes\n",
    "   - Columns include: gene IDs, symbols, expression levels, fold changes,\n",
    "     statistics, and categorical annotations\n\n",
    "2. DE_summary_by_gene_type.csv\n",
    "   - Summary statistics grouped by gene type\n\n",
    "COLUMN DESCRIPTIONS:\n",
    "--------------------\n",
    "gene_id: Ensembl gene ID\n",
    "gene_symbol: HGNC gene symbol\n",
    "gene_type: Gene biotype (protein_coding, lncRNA, etc.)\n",
    "baseMean: Mean normalized expression across all samples\n",
    "log2FoldChange: Log2 fold change (Mutant vs Control)\n",
    "foldChange: Fold change in linear scale (2^log2FC)\n",
    "direction: UP = upregulated in Mutant, DOWN = downregulated in Mutant, NS = not significant\n",
    "regulation: Detailed description of regulation direction\n",
    "lfcSE: Standard error of log2 fold change\n",
    "stat: Wald test statistic\n",
    "pvalue: Nominal p-value\n",
    "padj: Benjamini-Hochberg adjusted p-value\n",
    "significance: Significance category based on padj\n",
    "FC_category: Fold change magnitude category\n",
    "DE_category: Combined significance and fold change category\n\n",
    "Generated: ", as.character(Sys.time()), "\n"
)

# Save summary text
writeLines(
    summary_text,
    file.path(output_dir, "Complete_DE_analysis_summary.txt")
)

# Print summary to console
cat("\n")
cat(summary_text)
cat("\n")

# ============================================================
# 5. Top genes tables
# ============================================================
cat("Creating top genes tables...\n")

# Top 100 upregulated genes
top_up <- results[results$direction == "UP" & !is.na(results$direction), ]
top_up <- top_up[order(top_up$padj, -abs(top_up$log2FoldChange)), ]
top_up <- head(top_up, 100)

write.csv(
    top_up,
    file = file.path(output_dir, "Top100_upregulated_genes.csv"),
    row.names = FALSE
)

# Top 100 downregulated genes
top_down <- results[results$direction == "DOWN" & !is.na(results$direction), ]
top_down <- top_down[order(top_down$padj, -abs(top_down$log2FoldChange)), ]
top_down <- head(top_down, 100)

write.csv(
    top_down,
    file = file.path(output_dir, "Top100_downregulated_genes.csv"),
    row.names = FALSE
)

# Top genes by absolute fold change
top_fc <- results[!is.na(results$padj) & results$padj < 0.05, ]
top_fc <- top_fc[order(-abs(top_fc$log2FoldChange)), ]
top_fc <- head(top_fc, 100)

write.csv(
    top_fc,
    file = file.path(output_dir, "Top100_by_foldchange.csv"),
    row.names = FALSE
)

cat("\n==========================================\n")
cat("Complete! Output files saved to:", output_dir, "\n")
cat("\nMain files:\n")
cat("  - Complete_DESeq2_results_annotated.csv (all genes with annotations)\n")
cat("  - Complete_DE_analysis_summary.txt (summary statistics)\n")
cat("  - DE_summary_by_gene_type.csv (summary by gene type)\n")
cat("  - Top100_upregulated_genes.csv\n")
cat("  - Top100_downregulated_genes.csv\n")
cat("  - Top100_by_foldchange.csv\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
