# HCC Gene Expression Data Analysis

This repository contains a complete gene expression analysis pipeline for hepatocellular carcinoma (HCC), based on the GEO dataset **GSE22058**. The aim is to identify differentially expressed genes between liver tumor samples and adjacent non-tumor liver tissues, and to functionally characterize the results using pathway and miRNA-target enrichment analysis.

## Project Aim

The main objectives are:

- import and preprocess gene expression data;
- retain paired tumor and adjacent non-tumor samples;
- identify differentially expressed genes (DEGs);
- visualize the DEG expression profile;
- perform functional enrichment analysis;
- perform miRNA-target enrichment analysis using MIENTURNET;
- interpret the biological meaning of the results using pathway analysis and literature validation.

## Analysis Pipeline

| Script | Description |
|---|---|
| `01_data_import.R` | Downloads and parses GEO data, extracts metadata, selects paired tumor/non-tumor samples, and exports raw matrices. |
| `02_preprocessing.R` | Applies `log2(data + 1)` transformation and removes low-variability genes using IQR filtering. |
| `03_diff_expression.R` | Performs paired differential expression analysis using logFC, paired Student's t-test, and FDR correction. |
| `04_visualization.R` | Generates volcano plots, boxplots, heatmap, and PCA plots. |
| `05_enrichment.R` | Performs functional enrichment using EnrichR and exports gene lists for MIENTURNET. |

## Dataset

- GEO series: **GSE22058**
- Study design: paired tumor vs adjacent non-tumor liver tissue
- Final samples:
  - 96 adjacent non-tumor samples
  - 96 liver tumor samples
  - 192 total samples

The paired design allows tumor and adjacent non-tumor samples from the same patient to be compared directly.

## Differential Expression Criteria

Genes were considered differentially expressed using the following thresholds:

- adjusted p-value `< 0.05`
- `|log2FC| >= 1`
- equivalent to fold-change `>= 2`

Final DEG results:

- **1667 differentially expressed genes**
- **602 UP-regulated genes**
- **1065 DOWN-regulated genes**

## Main Biological Findings

UP-regulated genes were mainly enriched for:

- cell cycle;
- DNA replication;
- mitosis;
- chromosome segregation;
- mitotic checkpoints;
- p53 signaling;
- microRNAs in cancer.

DOWN-regulated genes were mainly enriched for:

- cytochrome P450 metabolism;
- xenobiotic and drug metabolism;
- retinol metabolism;
- bile secretion;
- biological oxidations;
- cytokine and chemokine signaling;
- immune/inflammatory pathways.

Overall, the results suggest that HCC tumor tissue activates proliferative cell-cycle programs while losing part of the normal liver metabolic and detoxification phenotype.

## Additional Analysis and Literature-Based Validation

Additional investigation was performed by integrating the DEG visualizations, representative gene boxplots, functional enrichment results, and MIENTURNET miRNA-target enrichment.

Selected findings were compared with previous HCC/cancer literature:

- **ZIC2**, the top UP-regulated gene, has been reported as a prognostic and immune-response marker in liver cancer/HCC.
- **CLEC1B**, the top DOWN-regulated gene, has been reported as an HCC prognostic biomarker related to immune infiltration.
- The down-regulation of **cytochrome P450** and drug-metabolism genes is consistent with published evidence of impaired or deregulated CYP expression/activity in HCC.
- The MIENTURNET findings are biologically plausible: **miR-146a** is linked to inflammatory Toll-like receptor/NF-kB signaling, while **miR-34a** is linked to p53, cell-cycle arrest, DNA-damage response, and apoptosis.

These observations support the biological plausibility of the computational results, although they remain hypothesis-generating and would require independent validation.

## miRNA-Target Enrichment

MIENTURNET was used manually with the exported DEG lists:

```text
results/tables/gene_symbols_UP_mienturnet.txt
results/tables/gene_symbols_DOWN_mienturnet.txt
```

Main candidate miRNAs:

- UP-regulated gene modules:
  - `hsa-let-7b-5p`
  - `hsa-miR-192-5p`
  - `hsa-miR-193b-3p`
  - `hsa-miR-34a-5p`
  - `hsa-miR-524-5p`

- DOWN-regulated gene modules:
  - `hsa-miR-146a-5p`
  - `hsa-miR-335-5p`

These miRNAs should be interpreted as candidate regulators inferred from target enrichment, not as directly measured differentially expressed miRNAs.

## Repository Structure

```text
.
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── plots/
│   └── tables/
├── 01_data_import.R
├── 02_preprocessing.R
├── 03_diff_expression.R
├── 04_visualization.R
├── 05_enrichment.R
└── README.md
```

## Requirements

Main R packages used:

```r
GEOquery
ggplot2
pheatmap
enrichR
forcats
stringr
```

Install missing packages before running the scripts.

## How to Run

Run the scripts in order from the project root directory:

```bash
Rscript 01_data_import.R
Rscript 02_preprocessing.R
Rscript 03_diff_expression.R
Rscript 04_visualization.R
Rscript 05_enrichment.R
```

After running `05_enrichment.R`, upload the exported UP and DOWN gene lists to MIENTURNET to complete the miRNA-target enrichment analysis.

## Output

The pipeline generates:

- processed expression matrices;
- DEG tables;
- volcano plots;
- representative gene boxplots;
- DEG heatmap;
- PCA plot;
- functional enrichment tables and barplots;
- gene lists for MIENTURNET.

## Conclusion

This project identifies a coherent HCC gene expression signature characterized by increased proliferation-related gene expression in tumor tissue and reduced expression of normal liver metabolic and detoxification genes. Functional enrichment, miRNA-target enrichment, and literature validation provide additional biological interpretation of the DEG signature.
