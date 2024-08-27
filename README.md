# DE-Analysis

Script Overview

Data Loading
The RNA-seq count data is loaded into R using read.table() for both WT and mutant samples.
Data Preprocessing

Counts are filtered to remove rows with very low counts (sum across all samples < 20).
Specific sample columns are selected depending on the comparison of interest (WT vs. mutants, WT vs. rescue clones).
Differential Expression Analysis

DESeq2 is used to perform differential expression analysis, comparing conditions specified in the condition vector.
Results are saved into an Excel file (WTvsRescue_DE_Results1.xlsx).
Significant Genes Filtering

Significant genes (padj < 0.05) are filtered and saved into a separate Excel file.
Separate lists of upregulated and downregulated genes are also generated and saved.
Visualization

Heatmaps of significant genes and top 50 DE genes are generated using ComplexHeatmap and pheatmap respectively.
Enhanced Volcano plots are created to visualize the DE genes, focusing on the most significant ones.
Specific genes of interest are highlighted in the volcano plot.
Gene Ontology (GO) Analysis

GO enrichment analysis is performed separately for upregulated and downregulated genes.
Results are visualized using barplots and dotplots.


How to Run the Script
Install Required Packages:

Install any missing R packages using install.packages() or BiocManager::install() for Bioconductor packages.
Execute the Script:

Run the R script line by line or source it as a whole. Ensure that the data files are in the working directory or specify the correct file paths.
Review Outputs:

Check the Excel files for DE results and significant gene lists.
Review the generated heatmaps, volcano plots, and GO analysis plots.

Notes

Ensure that your RNA-seq count data files are properly formatted and contain only the necessary columns before running the script.
Adjust file paths in the script to match your working directory structure.
