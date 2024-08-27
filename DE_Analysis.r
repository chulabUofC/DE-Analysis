#MDSC 528 12Sep2023
#HES7 mutation (R25W aa or C73T nt)  
#Sub789: 23 samples, (WT x 4 , homo mutants, rescue) 
#Sub789_23samples_EC rename_rescueClones

#TASK 1:  DE on RNA-seq data Compare 1: WT vs. Homozygous mutated clones 
#Upregulated gene list (mutant vs. WT) 
#downregulated gene list (mutant vs. WT) 
#DE on RNA-seq data Compare 2: WT vs. rescued clones 
#Q1: are there signaling pathways or families of genes are enriched in DE gene lists (up and down) 
#Q2: specifically look at NOTCH, WNT, FGF signaling pathways  
#Q3: anything elses? (Transcriptions factors families? circadian?) 
#Q4: we will want to compare to time course osc data  



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("EnhancedVolcano")
library(dada2)
library(dplyr)
library(DESeq2)
library(ggplot2)





#load bul RNA SEQ data

counts <- read.table("Sub789_23samples_EC rename.txt", header=TRUE, row.names=1, sep="\t")

counts <- read.table("Sub793_noP_vs_Resync_EC_7samples.txt", header=TRUE, row.names=1, sep="\t")
counts <- counts %>% select(starts_with("Control_rep"), starts_with("Resync_rep"))

#skip the steps below if your file has all the columns you will use in your analysis. 
#if you have extra columns like different mutation/sample, then specify the columns as below

# Select ONLY the WT and SPECIFIED sample columns
counts <- counts %>% select(starts_with("WT"), starts_with("r11rM1c5B"), starts_with("r12rM1c1H"), starts_with("r12rM1c8C"))
# Select only the WT and r7M1c7F, r7M1c1G columns 
counts <- counts %>% select(starts_with("WT_"), starts_with("r7M1c7F"), starts_with("r7M1c1G"))

counts <- counts %>% select(starts_with("r7M1c7F"),starts_with("r7M1c1G"),starts_with("r11rM1c5B"), starts_with("r12rM1c1H"),starts_with("r12rM1c8C"))





# Preprocessing - filter out rows with very low counts
counts <- counts[which(rowSums(counts) > 20),]
counts

# Define the condition factor for DESeq2  
#(if you have 4 control then 4 Cs if you have 10 samples then 10 S, adjust accordingly)
condition <- factor(c("C","C","C","C","S","S","S","S","S","S","S","S"))
coldata <- data.frame(row.names = colnames(counts), condition)
coldata


# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)

dds <-DESeq(dds)
vsdata <- vst(dds, blind=FALSE)


#PCA PLOT QC
plotPCA(vsdata, intgroup = "condition")

#PLOT DISPERTION
plotDispEsts(dds)

library(openxlsx)
res <- results(dds,contrast = c("condition", "S", "C"))
res
result_df2 <- data.frame(Gene = rownames(res), res)
coldata

# Write the data frame to an Excel file
write.xlsx(result_df2, "WTvsRescue_DE_Results1.xlsx", row.names = FALSE)

#Correct for padj value usually <0.05 to get significant genes

sigs <- na.omit(res)
sigs <- sigs[sigs$padj<0.05,]
sigs

result_df3 <- data.frame(Gene = rownames(sigs), sigs)

# Write the data frame to an Excel file
write.xlsx(result_df3, "WTvsRescue_SigResults.xlsx", row.names = FALSE)



# Create separate dataframes for upregulated and downregulated genes for your reference and save those lists
upregulated <- sigs[sigs$log2FoldChange > 0, ]
downregulated <- sigs[sigs$log2FoldChange < 0, ]

uplregulated_df <- data.frame(Gene = rownames(upregulated), upregulated)
downregulated_df <- data.frame(Gene = rownames(downregulated), downregulated)
# Specify the file names for Excel files
# Write the data frame to an Excel file
write.xlsx(uplregulated_df, "WTVSRESCUEUpregulated Genes1.xlsx", row.names = FALSE)
write.xlsx(downregulated_df, "WTVSRESCUEDownregulated Genes1.xlsx", row.names = FALSE)


#     SIMPLE HEATMAP
library(ComplexHeatmap)
library("org.Hs.eg.db")

## WT vs SCDO4 
sigs.df <- as.data.frame(sigs)

sigs.df <-sigs.df[(sigs.df$baseMean>70)& (abs(sigs.df$log2FoldChange) >1),]

sigs.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df), keytype = "SYMBOL", column = "SYMBOL")
sigs.df
#normalized dds
mat <-counts(dds, normalized = T)[rownames(sigs.df),]
mat.z <- t(apply(mat,1,scale))   
colnames(mat.z) <- rownames(coldata)

ht<- Heatmap(mat.z, cluster_rows = T, cluster_columns =T, column_labels = colnames(mat.z), name = "Z-score", row_labels =rownames(sigs.df), 
            )

############  TOP 50 DE GENES HEATMAP

top_genes <- rownames(sigs[order(sigs$padj),])[1:50]

# Subset the counts matrix to include only the top 50 genes
top_counts <- counts[top_genes,]

# Normalize the counts
top_counts_normalized <- t(apply(top_counts, 1, scale))
pheatmap(top_counts_normalized, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         fontsize_row = 7, fontsize_col = 9,
         border_color = NA,
         main = "WT VS RESCUE - Top 50 DE Genes Heatmap",
         labels_col = colnames(top_counts))



# ENHANCED VOLCANO for significant genes
library(org.Hs.eg.db)
sigs.df <- as.data.frame(sigs)
sigs.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df), keytype = "SYMBOL", column = "SYMBOL")
sigs.df


EnhancedVolcano(sigs.df,x ="log2FoldChange", y = "padj", lab = rownames(sigs.df), pCutoff =  1e-4, FCcutoff = 1)

# EHNANCED VOLCANO FOR SELECTED GENES
selected= c("HES7", "LFNG,","NOTCH1", "DLL1","HES1","JAG1", "HEY1", "HEYL","SNAI2", "BMP7", "TMEME100", "HES5")
EnhancedVolcano(sigs.df,x ="log2FoldChange", y = "padj", lab = rownames(sigs.df), pCutoff =  1e-4, FCcutoff = 1, selectLab = selected)


## GO ANALYSIS

##### UPREGULATED GENES  > 0.05

### we saved significant genes list in sigs 

library(enrichGO)
library(clusterProfiler)
library(org.Hs.eg.db)
#filter for genes that have log2FoldChnage > 0.5 which will give us positive values
# Positive values mean upregulated genes and Upregulated GO terms
genes_to_test<- rownames(sigs[sigs$log2FoldChange>0.5,])
genes_to_test

####  Run GO analysis 
# ketype is symbol- our genes to test as gene names, this can be GeneID or others depending on your dataset
# BP IS biological process, can use other arguments for molecular function MF or CC
GO_results1 <- enrichGO(genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont ="BP" )
as.data.frame(GO_results1)

# See top 20 GO terms, adjust number of terms
plot(barplot(GO_results1, showCategory = 20))


# Perform GO enrichment analysis
GO_results1 <- enrichGO(gene         = genes_to_test, 
                        OrgDb        = org.Hs.eg.db, 
                        keyType      = "SYMBOL", 
                        ont          = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05,
                        readable     = TRUE)

# Convert results to a data frame
GO_results_df <- as.data.frame(GO_results1)

# Plot the results
# If you are using clusterProfiler's barplot function, ensure that it is loaded or use the full namespace
dotplot <- dotplot(GO_results1, showCategory = 20)
print(dotplot)


###########DOWNREGULATED


#filtering genes from sigs that have a Log2FoldChnage < 0.5 for negative values
#therefore downregulated GOterms
genes_to_test2<- rownames(sigs[sigs$log2FoldChange<0.5,])
genes_to_test2
GO_results2 <- enrichGO(genes_to_test2, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont ="BP" )
as.data.frame(GO_results2)
plot(barplot(GO_results2, showCategory = 20))

GO_results2 <- enrichGO(gene         = genes_to_test2, 
                        OrgDb        = org.Hs.eg.db, 
                        keyType      = "SYMBOL", 
                        ont          = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05,
                        readable     = TRUE)

# Convert results to a data frame
GO_results_df <- as.data.frame(GO_results2)

# Plot the results
# If you are using clusterProfiler's barplot function, ensure that it is loaded or use the full namespace
dotplot <- dotplot(GO_results2, showCategory = 20)
print(dotplot)





