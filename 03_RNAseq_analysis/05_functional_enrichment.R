#!/usr/bin/env Rscript

# RNA-seq Pipeline Step 6: Functional Enrichment Analysis (GO and KEGG)
# This script performs Gene Ontology and KEGG pathway enrichment analysis

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(tidyverse)
})

cat("==========================================\n")
cat("Starting Functional Enrichment Analysis\n")
cat("==========================================\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

# Define paths
deseq_results_file <- "results/DESeq2/DESeq2_results_all.csv"
output_dir <- "results/GO_KEGG"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set parameters
padj_cutoff <- 0.05
lfc_cutoff <- 1
pvalue_cutoff <- 0.05
qvalue_cutoff <- 0.2

# ============================================================
# 1. Load DESeq2 results
# ============================================================
cat("Loading DESeq2 results...\n")
deseq_results <- read.csv(deseq_results_file, row.names = 1)

# Remove version numbers from Ensembl IDs (e.g., ENSG00000000003.15 -> ENSG00000000003)
rownames(deseq_results) <- gsub("\\..*", "", rownames(deseq_results))

cat("Loaded results for", nrow(deseq_results), "genes\n\n")

# ============================================================
# 2. Prepare gene lists
# ============================================================
cat("Preparing gene lists...\n")

# Get significant genes
sig_genes <- deseq_results %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
    rownames()

# Upregulated genes
up_genes <- deseq_results %>%
    filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    rownames()

# Downregulated genes
down_genes <- deseq_results %>%
    filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    rownames()

cat("Significant genes (padj <", padj_cutoff, ", |log2FC| >", lfc_cutoff, "):", length(sig_genes), "\n")
cat("  - Upregulated:", length(up_genes), "\n")
cat("  - Downregulated:", length(down_genes), "\n\n")

# Convert Ensembl IDs to Entrez IDs
cat("Converting Ensembl IDs to Entrez IDs...\n")
sig_genes_entrez <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
up_genes_entrez <- bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_genes_entrez <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("Converted", nrow(sig_genes_entrez), "genes to Entrez IDs\n")
cat("  - Upregulated:", nrow(up_genes_entrez), "\n")
cat("  - Downregulated:", nrow(down_genes_entrez), "\n\n")

# ============================================================
# 3. Gene Ontology Enrichment Analysis
# ============================================================
cat("Running Gene Ontology enrichment analysis...\n")

# GO enrichment for all significant genes
cat("  - All significant genes...\n")
go_all <- enrichGO(
    gene = sig_genes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
)

# GO enrichment for upregulated genes
cat("  - Upregulated genes...\n")
go_up <- enrichGO(
    gene = up_genes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
)

# GO enrichment for downregulated genes
cat("  - Downregulated genes...\n")
go_down <- enrichGO(
    gene = down_genes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
)

# ============================================================
# 4. KEGG Pathway Enrichment Analysis
# ============================================================
cat("\nRunning KEGG pathway enrichment analysis...\n")

# KEGG enrichment for all significant genes
cat("  - All significant genes...\n")
kegg_all <- enrichKEGG(
    gene = sig_genes_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
)

# Convert KEGG gene IDs to gene symbols
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
    kegg_all <- setReadable(kegg_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# KEGG enrichment for upregulated genes
cat("  - Upregulated genes...\n")
kegg_up <- enrichKEGG(
    gene = up_genes_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
)

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# KEGG enrichment for downregulated genes
cat("  - Downregulated genes...\n")
kegg_down <- enrichKEGG(
    gene = down_genes_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
)

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# ============================================================
# 5. Save results
# ============================================================
cat("\nSaving enrichment results...\n")

# Save GO results
if (!is.null(go_all) && nrow(go_all) > 0) {
    write.csv(as.data.frame(go_all), file.path(output_dir, "GO_enrichment_all_genes.csv"), row.names = FALSE)
    cat("  - GO all genes:", nrow(go_all), "terms\n")
}

if (!is.null(go_up) && nrow(go_up) > 0) {
    write.csv(as.data.frame(go_up), file.path(output_dir, "GO_enrichment_upregulated.csv"), row.names = FALSE)
    cat("  - GO upregulated:", nrow(go_up), "terms\n")
}

if (!is.null(go_down) && nrow(go_down) > 0) {
    write.csv(as.data.frame(go_down), file.path(output_dir, "GO_enrichment_downregulated.csv"), row.names = FALSE)
    cat("  - GO downregulated:", nrow(go_down), "terms\n")
}

# Save KEGG results
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
    write.csv(as.data.frame(kegg_all), file.path(output_dir, "KEGG_enrichment_all_genes.csv"), row.names = FALSE)
    cat("  - KEGG all genes:", nrow(kegg_all), "pathways\n")
}

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    write.csv(as.data.frame(kegg_up), file.path(output_dir, "KEGG_enrichment_upregulated.csv"), row.names = FALSE)
    cat("  - KEGG upregulated:", nrow(kegg_up), "pathways\n")
}

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    write.csv(as.data.frame(kegg_down), file.path(output_dir, "KEGG_enrichment_downregulated.csv"), row.names = FALSE)
    cat("  - KEGG downregulated:", nrow(kegg_down), "pathways\n")
}

# Save R objects
saveRDS(go_all, file.path(output_dir, "go_all.rds"))
saveRDS(go_up, file.path(output_dir, "go_up.rds"))
saveRDS(go_down, file.path(output_dir, "go_down.rds"))
saveRDS(kegg_all, file.path(output_dir, "kegg_all.rds"))
saveRDS(kegg_up, file.path(output_dir, "kegg_up.rds"))
saveRDS(kegg_down, file.path(output_dir, "kegg_down.rds"))

# ============================================================
# 6. Generate enrichment plots
# ============================================================
cat("\nGenerating enrichment plots...\n")

# GO dotplots
if (!is.null(go_all) && nrow(go_all) > 0) {
    p <- dotplot(go_all, showCategory = 20, title = "GO Enrichment - All Significant Genes")
    ggsave(file.path(output_dir, "GO_dotplot_all.png"), p, width = 12, height = 10, dpi = 300)
}

if (!is.null(go_up) && nrow(go_up) > 0) {
    p <- dotplot(go_up, showCategory = 20, title = "GO Enrichment - Upregulated Genes")
    ggsave(file.path(output_dir, "GO_dotplot_upregulated.png"), p, width = 12, height = 10, dpi = 300)
}

if (!is.null(go_down) && nrow(go_down) > 0) {
    p <- dotplot(go_down, showCategory = 20, title = "GO Enrichment - Downregulated Genes")
    ggsave(file.path(output_dir, "GO_dotplot_downregulated.png"), p, width = 12, height = 10, dpi = 300)
}

# KEGG dotplots
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
    p <- dotplot(kegg_all, showCategory = 20, title = "KEGG Pathway Enrichment - All Significant Genes")
    ggsave(file.path(output_dir, "KEGG_dotplot_all.png"), p, width = 12, height = 10, dpi = 300)
}

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    p <- dotplot(kegg_up, showCategory = 20, title = "KEGG Pathway Enrichment - Upregulated Genes")
    ggsave(file.path(output_dir, "KEGG_dotplot_upregulated.png"), p, width = 12, height = 10, dpi = 300)
}

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    p <- dotplot(kegg_down, showCategory = 20, title = "KEGG Pathway Enrichment - Downregulated Genes")
    ggsave(file.path(output_dir, "KEGG_dotplot_downregulated.png"), p, width = 12, height = 10, dpi = 300)
}

cat("\nFunctional enrichment analysis complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("\nKey output files:\n")
cat("  - GO_enrichment_*.csv: Gene Ontology enrichment results\n")
cat("  - KEGG_enrichment_*.csv: KEGG pathway enrichment results\n")
cat("  - *_dotplot_*.png: Visualization of enrichment results\n")
cat("  - *.rds: R objects for downstream analysis\n")
cat("\nEnd time:", as.character(Sys.time()), "\n")
cat("==========================================\n")
