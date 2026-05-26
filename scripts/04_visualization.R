rm(list = ls())
options(stringsAsFactors = FALSE)

# ==============================================================================
# LIBRARIES & SETUP
# ==============================================================================
library(pheatmap)
library(ggplot2)

dir_processed <- "data/processed/"
dir_raw <- "data/raw/"
dir_plots <- "results/plots/"
dir_tables <- "results/tables/"

thr_FC <- 2
thr_pval <- 0.05

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("Loading DEG matrix and results...\n")
results <- read.table(paste0(dir_tables, "DEG.txt"), header = TRUE, sep = "\t")
data_DEG <- read.table(paste0(dir_processed, "matrix_DEG.txt"), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

normal_samples <- read.table(paste0(dir_raw, "normal.txt"), header = FALSE)$V1
tumor_samples <- read.table(paste0(dir_raw, "tumor.txt"), header = FALSE)$V1

# ==============================================================================
# STEP 2: POST-FILTER VOLCANO PLOT
# ==============================================================================
cat("Generating post-filter Volcano Plot...\n")
pdf(paste0(dir_plots, "Volcano_PostFilter.pdf"))
plot(results$logFC, -log10(results$pval_adj), main = "Volcano Plot (After Filtering)",
     xlab = "log2 Fold Change", ylab = "-log10(Adjusted P-value)", pch = 20, 
     col = ifelse(results$direction == "UP", "magenta", "cyan"))
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = 'red')
abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = 'blue')
dev.off()

# ==============================================================================
# STEP 3: BOXPLOTS OF TOP UP AND DOWN DEGs
# ==============================================================================
cat("Generating boxplots for the most UP and DOWN regulated genes...\n")

# Subset results by direction
results_up <- results[results$direction == "UP", ]
results_down <- results[results$direction == "DOWN", ]

# Ensure we actually have genes in both directions before trying to plot
genes_to_plot <- list()

if (nrow(results_up) > 0) {
    genes_to_plot$UP <- list(gene = results_up$geneSymbol[1], pval = results_up$pval_adj[1])
}
if (nrow(results_down) > 0) {
    genes_to_plot$DOWN <- list(gene = results_down$geneSymbol[1], pval = results_down$pval_adj[1])
}

for (dir in names(genes_to_plot)) {
    gene <- genes_to_plot[[dir]]$gene
    pval <- genes_to_plot[[dir]]$pval
    
    df_boxplot <- data.frame(
        expression = as.numeric(data_DEG[gene, ]),
        condition = factor(ifelse(colnames(data_DEG) %in% normal_samples, "Normal", "Cancer"), levels = c("Normal", "Cancer"))
    )
    
    pdf(paste0(dir_plots, "Boxplot_", dir, "_", gene, ".pdf"))
    boxplot(expression ~ condition, data = df_boxplot,
            main = paste0(gene, " Expression (Top ", dir, ")\nAdj p-val = ", format(pval, scientific = TRUE)),
            ylab = "Log2 Expression", col = c("green", "orange"), notch = TRUE)
    dev.off()
}

# ==============================================================================
# STEP 4: HEATMAP
# ==============================================================================
cat("Generating DEG heatmap...\n")
# Sample Annotations (Columns)
annotation_col <- data.frame(condition = ifelse(colnames(data_DEG) %in% normal_samples, "Normal", "Cancer"))
rownames(annotation_col) <- colnames(data_DEG)

# Gene Annotations (Rows)
annotation_row <- data.frame(direction = results$direction)
rownames(annotation_row) <- results$geneSymbol

# Legend Colors
annotation_colors <- list(
    condition = c(Normal = "green", Cancer = "orange"),
    direction = c(UP = "magenta", DOWN = "cyan")
)

pheatmap(data_DEG, scale = "row", border_color = NA,
         cluster_cols = TRUE, cluster_rows = TRUE,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         clustering_method = "average",
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "blue3", "black", "yellow3", "yellow"))(100),
         show_rownames = FALSE, show_colnames = FALSE,
         cutree_cols = 2, cutree_rows = 2,
         filename = paste0(dir_plots, "Heatmap_DEGs.pdf"), width = 10, height = 10) # 

# ==============================================================================
# STEP 5: OPTIONAL PCA 
# ==============================================================================
cat("Computing PCA of DEGs...\n")
# PCA expects samples as rows and genes as columns
pca_res <- prcomp(t(data_DEG), scale. = TRUE)
pca_data <- data.frame(Sample = rownames(pca_res$x),
                       PC1 = pca_res$x[, 1],
                       PC2 = pca_res$x[, 2],
                       Condition = annotation_col$condition)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Normal" = "green", "Cancer" = "orange")) +
    theme_minimal() +
    labs(title = "PCA of Differentially Expressed Genes",
         x = paste0("PC1: ", round(summary(pca_res)$importance[2, 1] * 100, 1), "% variance"),
         y = paste0("PC2: ", round(summary(pca_res)$importance[2, 2] * 100, 1), "% variance"))

pdf(paste0(dir_plots, "PCA_DEGs.pdf"))
print(pca_plot)
dev.off()

cat("Visualization complete!\n")
