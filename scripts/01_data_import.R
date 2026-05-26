rm(list = ls())
options(stringsAsFactors = FALSE)

# ==============================================================================
# LIBRARIES & SETUP
# ==============================================================================
library("GEOquery")

# Define dataset parameters and relative paths
series <- "GSE22058"
platform <- "GSE22058-GPL6793_series_matrix.txt.gz"
dir_raw <- "data/raw/"

# Create data directory safely
dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STEP 1: DOWNLOADING DATA
# ==============================================================================
cat("Downloading and extracting GEO series:", series, "...\n")
readr::local_edition(1)

tmp <- getGEO(GEO = series)
set <- tmp[[platform]]

pData <- phenoData(set)
metadata <- pData@data

aData <- assayData(set)
exprs_matrix <- data.frame(aData$exprs)

# Free up memory
rm(tmp, pData, aData)

# ==============================================================================
# STEP 2: PREPARING & ANNOTATING DATA
# ==============================================================================
cat("Annotating probes to Gene Symbols...\n")
annotation <- fData(set)
geneSymbol <- annotation$GeneSymbol

# Map matrix rows to annotation IDs
exprs_matrix <- exprs_matrix[annotation$ID, ]

# Aggregate multiple probes for the same gene using the mean
exprs_matrix <- aggregate(exprs_matrix, list(geneSymbol), FUN = mean, na.rm = TRUE)

# Remove rows where Gene Symbol is missing
ind_empty <- which(exprs_matrix$Group.1 == "")
if (length(ind_empty) > 0) {
    exprs_matrix <- exprs_matrix[-ind_empty, ]
}

# Set row names to Gene Symbols and remove the grouping column
rownames(exprs_matrix) <- exprs_matrix$Group.1
exprs_matrix <- exprs_matrix[, -1]

# Free up memory
rm(set, annotation, ind_empty, geneSymbol)

# ==============================================================================
# STEP 3: EXTRACTING TUMOR AND NORMAL SAMPLES (PAIRED)
# ==============================================================================
cat("Splitting into paired normal and tumor samples...\n")
metadata <- metadata[, c("geo_accession", "individual:ch1", "tissue:ch1")]

# Split metadata based on tissue type
list_tissue <- split(metadata, metadata$`tissue:ch1`)
normal <- list_tissue$`adjacent liver non-tumor`
tumor <- list_tissue$`liver tumor`

# Deduplicate to keep only one sample per patient (individual)
normal <- normal[!duplicated(normal$`individual:ch1`), ]
tumor <- tumor[!duplicated(tumor$`individual:ch1`), ]

# Find common patients
common_patients <- intersect(normal$`individual:ch1`, tumor$`individual:ch1`)

# Filter both sets to only include patients present in both conditions
normal <- normal[normal$`individual:ch1` %in% common_patients, ]
tumor <- tumor[tumor$`individual:ch1` %in% common_patients, ]

# Order by individual to ensure the paired t-test works correctly downstream
normal <- normal[order(normal$`individual:ch1`), "geo_accession"]
tumor <- tumor[order(tumor$`individual:ch1`), "geo_accession"]

# Subset matrix and metadata to only include these paired samples
final_data <- exprs_matrix[, c(normal, tumor)]
final_metadata <- metadata[c(normal, tumor), ]

# Free up memory
rm(list_tissue, exprs_matrix)

# ==============================================================================
# STEP 4: EXPORTING DATA
# ==============================================================================
cat("Exporting parsed data to", dir_raw, "...\n")

write.table(final_data, paste0(dir_raw, "matrix.txt"), sep= "\t",
            col.names = NA, row.names = TRUE, quote = FALSE)
write.table(normal, paste0(dir_raw, "normal.txt"), sep= "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(tumor, paste0(dir_raw, "tumor.txt"), sep= "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(final_metadata, paste0(dir_raw, "metadata.txt"), sep= "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Clean up downloaded GEO files to save local space
unlink(series, recursive = TRUE)

cat("Data import complete!\n")
