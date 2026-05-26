rm(list = ls())
options(stringsAsFactors = FALSE)

# ==============================================================================
# LIBRARIES & SETUP
# ==============================================================================
dir_processed <- "data/processed/"
dir_raw <- "data/raw/"
dir_plots <- "results/plots/"
dir_tables <- "results/tables/"

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

# Parameters
thr_FC <- 2
thr_pval <- 0.05
paired_test <- TRUE

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("Loading filtered matrix and sample lists...\n")
data <- read.table(paste0(dir_processed, "filtered_matrix.txt"), sep = "\t", 
                   header = TRUE, row.names = 1, check.names = FALSE)

normal_samples <- read.table(paste0(dir_raw, "normal.txt"), header = FALSE)$V1
tumor_samples <- read.table(paste0(dir_raw, "tumor.txt"), header = FALSE)$V1

dataN <- data[, normal_samples]
dataC <- data[, tumor_samples]

# ==============================================================================
# STEP 2: LOG FOLD CHANGE & T-TEST
# ==============================================================================
cat("Calculating LogFC and performing paired t-test...\n")
# Data is log-scale: logFC = mean(case) - mean(control) 
logFC <- rowMeans(dataC, na.rm = TRUE) - rowMeans(dataN, na.rm = TRUE)

# Paired t-test assuming normal distribution due to n > 30 
pval <- apply(data, 1, function(x) {
    t.test(x[tumor_samples], x[normal_samples], paired = paired_test)$p.value
})

# Adjust p-values for multiple testing 
pval_adj <- p.adjust(pval, method = "fdr")

# ==============================================================================
# STEP 3: PRE-FILTERING VOLCANO PLOT
# ==============================================================================
cat("Generating pre-filter Volcano Plot...\n")
pdf(paste0(dir_plots, "Volcano_PreFilter.pdf"))
plot(logFC, -log10(pval_adj), main = "Volcano Plot (Before Filtering)",
     xlab = "log2 Fold Change", ylab = "-log10(Adjusted P-value)", pch = 20, col = "darkgrey")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = 'red')
abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = 'blue')
dev.off()

# ==============================================================================
# STEP 4: FILTERING DEGs
# ==============================================================================
cat("Filtering DEGs based on thresholds...\n")
keep_deg <- abs(logFC) >= log2(thr_FC) & pval_adj < thr_pval

data_DEG <- data[keep_deg, ]
genes_DEG <- rownames(data_DEG)
logFC_DEG <- logFC[keep_deg]
pval_DEG <- pval[keep_deg]
pval_adj_DEG <- pval_adj[keep_deg]

cat("Differentially Expressed Genes found:", nrow(data_DEG), "\n")

# ==============================================================================
# STEP 5: EXPORT RESULTS
# ==============================================================================
cat("Exporting DEG results...\n")
direction <- ifelse(logFC_DEG > 0, "UP", "DOWN")
results <- data.frame(geneSymbol = genes_DEG,
                      pvalue = pval[keep_deg], 
                      pval_adj = pval_adj[keep_deg],
                      logFC = logFC_DEG, 
                      direction = direction)

# Sort by absolute logFC descending
results <- results[order(abs(results$logFC), decreasing = TRUE), ]

write.table(results, file = paste0(dir_tables, "DEG.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(data_DEG, file = paste0(dir_processed, "matrix_DEG.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

cat("Differential expression analysis complete!\n")
