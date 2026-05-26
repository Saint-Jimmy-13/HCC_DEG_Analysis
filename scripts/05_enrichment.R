rm(list = ls())
options(stringsAsFactors = FALSE)

# ==============================================================================
# LIBRARIES & SETUP
# ==============================================================================
library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)

dir_tables <- "results/tables/"
dir_enrich <- "results/plots/Functional_Enrichment/"
dir.create(dir_enrich, recursive = TRUE, showWarnings = FALSE)

# Parameters
top_term <- 10
thr_pval <- 0.05

# ==============================================================================
# HELPER FUNCTION: PLOT ENRICHMENT
# ==============================================================================
getEnrichmentPlot <- function(annotation, type, suffix, top_term, thr_pval, dirEnrich) {
    
    if (is.null(annotation) || nrow(annotation) == 0) return(NULL)
    
    annotation <- annotation[annotation$Adjusted.P.value < thr_pval, ]
    if (nrow(annotation) == 0) {
        cat("No significant terms found for", type, "(", suffix, ")\n")
        return(NULL)
    }
    
    annotation <- annotation[order(annotation$Adjusted.P.value, decreasing = FALSE), ]
    annotation$Gene_count <- sapply(annotation$Genes, function(x) length(unlist(strsplit(x, split = ";"))))
    annotation$Gene_ratio <- unlist(lapply(annotation$Overlap, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        parts[1] / parts[2]
    }))
    
    annotation_top <- if(top_term > 0 && top_term <= nrow(annotation)) annotation[1:top_term, ] else annotation
    file_prefix <- paste0(dirEnrich, type, "_", suffix)
    
    # Bar Plot
    g1 <- ggplot(annotation_top, aes(x = Gene_count, y = fct_reorder(Term, Gene_count), fill = Adjusted.P.value)) +
        geom_bar(stat = "identity") +
        scale_fill_continuous(low = "red", high = "blue", name = "Adjusted P-value", guide = guide_colorbar(reverse = TRUE)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme_bw(base_size = 10) + ylab(NULL) + ggtitle(paste(type, "-", suffix, "Terms"))
    
    pdf(paste0(file_prefix, "_barplot.pdf"))
    print(g1)
    dev.off()
    
    # Dot Plot
    g2 <- ggplot(annotation_top, aes(x = Gene_count, y = fct_reorder(Term, Gene_count))) +
        geom_point(aes(size = Gene_ratio, color = Adjusted.P.value)) +
        scale_colour_gradient(limits = c(min(annotation_top$Adjusted.P.value), max(annotation_top$Adjusted.P.value)), low = "red", high = "blue") +
        theme_bw(base_size = 10) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        ylab(NULL) + ggtitle(paste(type, "-", suffix, "Terms"))
    
    pdf(paste0(file_prefix, "_dotplot.pdf"))
    print(g2)
    dev.off()
    
    # Export table
    write.table(annotation[, c("Term", "Overlap", "P.value", "Adjusted.P.value", "Gene_count", "Gene_ratio", "Genes")],
                paste0(file_prefix, "_adj_pval_", thr_pval, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

# ==============================================================================
# STEP 1: LOAD & SPLIT GENE LISTS
# ==============================================================================
cat("Loading DEG list and splitting by direction...\n")
input_list <- read.table(paste0(dir_tables, "DEG.txt"), sep = "\t", header = TRUE, check.names = FALSE, quote = "")

genes_UP <- input_list$geneSymbol[input_list$direction == "UP"]
genes_DOWN <- input_list$geneSymbol[input_list$direction == "DOWN"]

cat("Found", length(genes_UP), "UP-regulated genes and", length(genes_DOWN), "DOWN-regulated genes.\n")

# Export separate lists for MIENTURNET
write.table(genes_UP, file = paste0(dir_tables, "gene_symbols_UP_mienturnet.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_DOWN, file = paste0(dir_tables, "gene_symbols_DOWN_mienturnet.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# ==============================================================================
# STEP 2: ENRICHR QUERY
# ==============================================================================
cat("Querying EnrichR databases for UP and DOWN genes separately...\n")
dbs <- c("Reactome_Pathways_2024",
         "GO_Molecular_Function_2025",
         "GO_Biological_Process_2025",
         "KEGG_2026",
         "DisGeNET")

df_UP <- enrichr(genes_UP, dbs)
df_DOWN <- enrichr(genes_DOWN, dbs)

# ==============================================================================
# STEP 3: GENERATE PLOTS
# ==============================================================================
cat("Generating enrichment plots...\n")

# UP Regulated Enrichment
getEnrichmentPlot(df_UP$GO_Biological_Process_2025, "GO_BP", "UP", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_UP$GO_Molecular_Function_2025, "GO_MF", "UP", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_UP$KEGG_2026, "KEGG", "UP", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_UP$Reactome_Pathways_2024, "Reactome", "UP", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_UP$DisGeNET, "DisGeNET", "UP", top_term, thr_pval, dir_enrich)

# DOWN Regulated Enrichment
getEnrichmentPlot(df_DOWN$GO_Biological_Process_2025, "GO_BP", "DOWN", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_DOWN$GO_Molecular_Function_2025, "GO_MF", "DOWN", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_DOWN$KEGG_2026, "KEGG", "DOWN", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_DOWN$Reactome_Pathways_2024, "Reactome", "DOWN", top_term, thr_pval, dir_enrich)
getEnrichmentPlot(df_DOWN$DisGeNET, "DisGeNET", "DOWN", top_term, thr_pval, dir_enrich)

cat("Directional functional enrichment complete!\n")

# ==============================================================================
# miRNA ENRICHMENT (MIENTURNET)
# ==============================================================================
# Perform a microRNA-target enrichment analysis using MIENTURNET.
# MIENTURNET is a web-based tool. To complete this requirement:
# 1. Navigate to: http://userver.bio.uniroma1.it/apps/mienturnet/
# 2. Upload your 'results/tables/gene_symbols_[DOWN|UP]mienturnet.txt' file.
# 3. Export the resulting network diagrams and tables to your 'results/' folder.
