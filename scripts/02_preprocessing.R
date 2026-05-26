rm(list = ls())
options(stringsAsFactors = FALSE)

# ==============================================================================
# LIBRARIES & SETUP
# ==============================================================================
# Define directories
dir_raw <- "data/raw/"
dir_processed <- "data/processed/"
dir_plots <- "results/plots/"

# Create directories safely
dir.create(dir_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

# Parameters
prc_IQR <- 0.25 # 25th percentile threshold

# ==============================================================================
# STEP 1: LOAD RAW DATA
# ==============================================================================
cat("Loading raw matrix...\n")
raw_matrix <- read.table(paste0(dir_raw, "matrix.txt"), sep = "\t", 
                         header = TRUE, row.names = 1, check.names = FALSE)

# ==============================================================================
# STEP 2: LOG2 TRANSFORMATION
# ==============================================================================
cat("Applying log2(data + 1) transformation...\n")
log_matrix <- log2(raw_matrix + 1)

# ==============================================================================
# STEP 3: IQR FILTERING
# ==============================================================================
cat("Calculating IQR and filtering genes...\n")
gene_iqr <- apply(log_matrix, 1, IQR)

# Plot IQR distribution to justify the threshold 
pdf(paste0(dir_plots, "IQR_distribution.pdf"))
hist(gene_iqr, breaks = 100, col = "steelblue", 
     main = "Frequency of IQR Distribution", 
     xlab = "Interquartile Range (IQR)", ylab = "Frequency")
abline(v = quantile(gene_iqr, prc_IQR), col = "red", lwd = 2, lty = 2)
dev.off()

# Determine threshold and filter 
thr_prc <- quantile(gene_iqr, prc_IQR)
keep_iqr <- gene_iqr > thr_prc
filtered_matrix <- log_matrix[keep_iqr, ]

cat("Initial genes:", nrow(raw_matrix), "\n")
cat("Genes after IQR filter (>25th percentile):", nrow(filtered_matrix), "\n")

# ==============================================================================
# STEP 4: EXPORT PROCESSED DATA
# ==============================================================================
cat("Exporting filtered matrix...\n")
write.table(filtered_matrix, paste0(dir_processed, "filtered_matrix.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Free up memory
rm(raw_matrix, log_matrix, gene_iqr, keep_iqr)
cat("Preprocessing complete!\n")
