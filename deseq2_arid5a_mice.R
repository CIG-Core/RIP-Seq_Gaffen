# Differential Expression Analysis (all the sample)

## Libraries
suppressMessages({
  library(openxlsx)
  library(biomaRt)
  library(stringr)
  library(matrixStats)
  library(pheatmap)
  library(dplyr)
  library(ggfortify)
  library(org.Hs.eg.db)
  library(WebGestaltR)
  library(DESeq2)
  library(apeglm)
  library(EnhancedVolcano)
  library(ReactomeContentService4R)
  library('janitor') #use the function: row_to_names
})

## Sourcing
source("/ix/cigcore/utils/code/de.heatmap.R")

## Set paths
proj_path <- "/ix/cigcore/proj/gaffen"
dataPath <- file.path(proj_path, "data", "Arid5a_CLIP-Seq_2023")
resPath <- gsub(pattern="data", replacement="result", dataPath)
experiment <- "11-16-2023-DE_Analysis_deseq2_RSEM_remove_IgG_1_IL17"
experimentPath <- file.path(resPath, experiment)
dir.create(experimentPath, showWarnings = TRUE, recursive = TRUE)
# webgestalt_path <- file.path(experimentPath, "10-30-2023-pathway_analysis_GSEA")
#dir.create(webgestalt_path, recursive=TRUE)
sig_DEG_YL <- read.xlsx(file.path(experimentPath, "sig_DEG_YL.xlsx"))
colnames(sig_DEG_YL)[1] <- "geneid"

## Read data file
coldata <- read.csv(file.path(dataPath, "Sissi_MetaData_DEseq2.csv"), row.names = 1) #Metadata
coldata$sample <- rownames(coldata)  # new variable:sample = rownames in 'coldata'
# Remove the outlier: "IgG_1_IL17"
coldata <- coldata[!(coldata$sample %in% c("IgG_1_IL17")), ]
coldata$groups <- gsub("\\.", "_", coldata$groups) # set groups as: IgG_Med; Arid5a_Med; IgG_IL17; Arid5a_IL17
coldata$groups <- as.factor(coldata$groups) # convert groups as factor variables
coldata <- coldata[ , c("sample", "groups")]

count <- read.table(file.path(resPath, "mapped_star", "rsem_counts",
                              "expected_counts.txt"), row.names = 1, header = TRUE)
count[, 1:ncol(count)] <- sapply(count[, 1:ncol(count)], as.numeric) #convert it to numeric factors
rawData <- count[ ,rownames(coldata)]
head(rawData, n=10)



all(colnames(rawData) %in% rownames(coldata)) # checking to see if sample names match between the two files
all(colnames(rawData) == rownames(coldata)) # checking to see if the sample names are in the same order in both files

## Create DESeqDataSet object
# DESeqDataSetFromMatrix: create a data metrix to store the input values, intermediate calculations and results of an analysis of differential expression. 
#The DESeqDataSet class enforces non-negative integer values in the "counts" matrix stored as the first element in the assay list.
dds <- DESeqDataSetFromMatrix(countData = round(rawData), colData = coldata, design = ~groups) # design is our conditions for comparisons
dds$groups 

## Make list with conditions and 4 comparisons
# List sequence: front --> treatment; back--> control
comparisons <- list(c("Arid5a_IL17", "IgG_IL17"), #whether Arid5a enhances or inhibits gene expression under stimulated (IL-17) conditions 
                    c("Arid5a_Med", "IgG_Med"), #whether Arid5a enhances or inhibits gene expression under unstimulated (Med) conditions  baseline effects of Arid5a
                    c("IgG_IL17", "IgG_Med"), #how IL-17 affects gene expression in IgG (control)
                    c("Arid5a_IL17", "Arid5a_Med")) #how IL-17 affects gene expression in Arid5a (treatment) 

de_results <- list() 

#data frames to store the summary of DE genes
DEgene_summary_all_biotype <- data.frame(Comparison = character(0), Sum_DE_Genes = numeric(0), Sum_Upregulated_DE_Genes = numeric(0), Sum_Downregulated_DE_Genes = numeric(0))
DEgene_summary_mRNA <- data.frame(Comparison = character(0), Sum_mRNA_DE_Genes = numeric(0), Sum_mRNA_Upregulated = numeric(0), Sum_mRNA_Downregulated = numeric(0))

for(c1 in 1:length(comparisons)){
  names <- paste0(comparisons[[c1]][1], "_vs_", comparisons[[c1]][2], sep ="") #create new variable 'names' for 4 comparisons: concatenates the first and second elements of the c1-th comparison using the paste0 function, separated by "vs"
  metadata <- subset(colData(dds), groups %in% comparisons[[c1]])  
  #Metadata: subset the colData with the condition that groups match the comparisons
  rawcount <- rawData[, metadata$sample] 
  #Rawcount: subsets the count data (gene expression counts) from 'dds': including all rows (genes) and only the columns that correspond to the samples specified in metadata$sample.
  
  # Create DESeqDataSet object with colData and countData subset
  dds_subset <- DESeqDataSetFromMatrix(countData = round(rawcount), #replace as.matrix to round since got the error: Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are not integers
                                       colData = metadata, design = ~groups)
  dds_subset <- estimateSizeFactors(dds_subset)
  idx <- rowSums( counts(dds_subset, normalized=TRUE) >= 10 ) >= 6 
  
  dds_subset <- dds_subset[idx,]
  
  new_level <- print(comparisons[[c1]][2]) #set up a reference group or level for the DESeq2
  dds_subset$groups <- relevel(dds_subset$groups, ref = new_level)#set the reference group in dds_subset
  
  dds1 <- DESeq(dds_subset) #DE analysis
  res <- results(dds1) # See the result from DE analysis
  print("DESeq done")
  
  ## Volcano Plot
  # Adjust 'padj' values and 'log2FoldChange' values
  res_volcano <- res  # Create a new data frame 'res_volcano' to store the adjusted values
  res_volcano$padj <- ifelse(res$padj < 10^(-10), 10^(-10), ifelse(res$padj > 10^0.4, 10^0.4, res$padj))
  res_volcano$log2FoldChange <- ifelse(res$log2FoldChange < -2.5, -2.5, ifelse(res$log2FoldChange > 3, 3, res$log2FoldChange))
  
  filtered_labels <- rownames(res_volcano)[rownames(res_volcano) %in% sig_DEG_YL$geneid]
  
  # Create a vector of NA values with the same length as res_volcano
  na_labels <- rep(NA, nrow(res_volcano))
  
  # Replace the NA values with the filtered_labels for matching genes
  na_labels[match(filtered_labels, rownames(res_volcano))] <- filtered_labels
  
  
  volcano_plot <- EnhancedVolcano(
    res_volcano,
    lab = na_labels,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    labSize = 3.0,
    title = paste0(names),
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    xlim = c(-2.5, 3),
    ylim = c(-0.4, 10.0)
  )
  
  
  
  ggsave(file = file.path(experimentPath, paste0(names, "_volcano_plot_label.png")), width = 10, height = 7) #delete #
  
  volcano_plot_no_label <- EnhancedVolcano(
    res_volcano, lab = NA, x = 'log2FoldChange', y = 'padj',
    pCutoff = 0.05, labCol = 'black', labFace = 'bold',
    boxedLabels = TRUE, labSize = 3.0, title = paste0(names),
    drawConnectors = TRUE, widthConnectors = 1.0,
    colConnectors = 'black',
    xlim = c(-2.5, 3),  # Adjust x-axis limits
    ylim = c(-0.4, 10.0)  # Adjust y-axis limits
  )
  
  
  ggsave(file = file.path(experimentPath, paste0(names, "_volcano_plot_no_label.png")), width = 10, height = 7) #delete #
  
  ## DE p.adj <0.05
  res2 <- res[complete.cases(res), ] #removing rows with missing data (NAs) from res
  res2 <- res2[order(res2$padj), ] #reorder the based on the "padj" column in ascending order.
  
  resDE <- as.data.frame(res2[which(res2$padj < 0.05), ]) #resDE:include only rows where the "padj" value is less than 0.05
  #padj < 0.05: identifying statistically significant DE genes
  ensemblIDs <- rownames(resDE) # ensemblIDs: statistically significant DE genes
  de_results[[paste0(names, "_", "DeG")]] <- resDE # new dataframe with data from resDE, index is $comparison(names)_DeG
  resDE$ENSEMBL <- rownames(resDE) #create new column in resDE called ENSEMBL (same value as rownames in resDE)
  

  ## Get biotype information
  ensemblMart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")#access and retrieve data from Ensembl for mouse
  g2g <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "gene_biotype",
                            "entrezgene_id", "ensembl_gene_id_version"),
               values = ensemblIDs, mart = ensemblMart) 
  #retrieve specific attributes (external_gene_name, ensembl_gene_id, gene_biotype, entrezgene_id, ensembl_gene_id_version) for a list of gene identifiers (ensemblIDs) 
  write.csv(g2g, file.path(experimentPath, paste0(names, "_g2g.csv")), row.names = F) #delete #
  # g2g <- read.csv(file.path(experimentPath, paste0(names, "_g2g.csv")) 
  
  duplicates <- g2g$ensembl_gene_id_version[duplicated(g2g$ensembl_gene_id_version)]
  g2g2 <- g2g[!duplicated(g2g$ensembl_gene_id_version), ]#remove duplicate rows based on the "ensembl_gene_id_version" 
  annotateData <- merge(resDE, g2g2, by.x = "ENSEMBL", by.y = "external_gene_name") #merge two data by matching "ENSEMBL" in resDE and "external_gene_name" in g2g2. 
  colnames(annotateData)[colnames(annotateData) == "Row.names"] <- "ENSEMBL"
  de_results[[paste0(names, "AllBT")]] <- annotateData #store the annotateData under de_results
  
  
  # Count the number of significant DE genes and store in de_counts_df
  annotateData_no_na <- annotateData[!is.na(annotateData$ENSEMBL), ]
  annotateData_nodup <- annotateData_no_na[!duplicated(annotateData_no_na$ENSEMBL), ]
  DEg_Up <- resDE[annotateData_nodup$log2FoldChange > 0, ] # upregulated DE genes in the condition of interest.
  DEg_Dn <- resDE[annotateData_nodup$log2FoldChange < 0, ] # downregulated DE genes in the condition of interest.
  Sum_significant_genes <- nrow(annotateData_nodup)
  Sum_upregulated_genes <- nrow(DEg_Up)
  Sum_downregulated_genes <- nrow(DEg_Dn)
  DEgene_summary_all_biotype <- rbind(DEgene_summary_all_biotype, data.frame(Comparison = names, Sum_DE_Genes = Sum_significant_genes, Sum_Upregulated_DE_Genes = Sum_upregulated_genes, Sum_Downregulated_DE_Genes = Sum_downregulated_genes))
  

  write.table(annotateData[ , c(1, 3)], #extract 1st and 3rd column from annotateDate: gene id(ENSEMBL) and log2FoldChange
              file = file.path(experimentPath, paste0(names, "_gsea.rnk")),#delete #
              quote=FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)#delete #
  
  biotypes <- unique(g2g2$gene_biotype) #unique gene biotypes
  pseudoBiotypes <- biotypes[biotypes %in% (str_subset(biotypes, pattern = "_pseudogene"))]
  #pseudoBiotypes: contain the substring "_pseudogene." 
  otherBiotypes <- biotypes[!biotypes %in% c(pseudoBiotypes, "protein_coding", "lncRNA", "rRNA", "miRNA")]
  # otherBiotypes variable: contain biotypes of "pseudoBiotypes," "protein_coding," "lncRNA," "rRNA," and "miRNA." 
  
  #create new columns "grouped_biotypes" with the assigned categorie in biotypes:
  annotateData[which(annotateData$gene_biotype %in% pseudoBiotypes), "grouped_biotypes"] <- "pseudogenes"
  annotateData[which(annotateData$gene_biotype %in% otherBiotypes), "grouped_biotypes"] <- "miscellinous"
  annotateData[which(annotateData$gene_biotype %in% c("protein_coding")), "grouped_biotypes"] <- "mRNA"
  annotateData[which(!annotateData$gene_biotype %in% c("protein_coding", pseudoBiotypes, otherBiotypes)),
               "grouped_biotypes"] <- annotateData$gene_biotype[!annotateData$gene_biotype %in% c("protein_coding", pseudoBiotypes, otherBiotypes)]
  
  # Filter DE genes within the biotype "mRNA"
  #question
  mRNA_resDE_raw <- as.data.frame(annotateData[which(annotateData$gene_biotype == "protein_coding"), ])

#  de_results[[paste0(names, "_mDeU")]] <- mRNA_resDE_Up
#  de_results[[paste0(names, "_mDeD")]] <- mRNA_resDE_Dn
  
  # Count the number of mRNA DE genes and their upregulated and downregulated subsets
  mRNA_resDE_no_na <- mRNA_resDE_raw[!is.na(mRNA_resDE_raw$ENSEMBL), ]
  mRNA_resDE <- mRNA_resDE_no_na[!duplicated(mRNA_resDE_no_na$ENSEMBL), ]
  mRNA_resDE_Up <- resDE[mRNA_resDE$log2FoldChange > 0, ] # upregulated DE genes in the condition of interest.
  mRNA_resDE_Dn <- resDE[mRNA_resDE$log2FoldChange < 0, ]
  Sum_mRNA_DE_genes <- nrow(mRNA_resDE)
  Sum_mRNA_Upregulated <- nrow(mRNA_resDE_Up)
  Sum_mRNA_Downregulated <- nrow(mRNA_resDE_Dn)
  
  de_results[[paste0(names, "_mDeG")]] <- mRNA_resDE
  
  # Append the summary to the mRNA_summary data frame
  DEgene_summary_mRNA <- rbind(DEgene_summary_mRNA, data.frame(Comparison = names, Sum_mRNA_DE_Genes = Sum_mRNA_DE_genes, Sum_mRNA_Upregulated = Sum_mRNA_Upregulated, Sum_mRNA_Downregulated = Sum_mRNA_Downregulated))
  
  ## Plot pie chart
  freqBiotypes <- as.data.frame(table(annotateData$grouped_biotypes))
  freqBiotypes$percent <- round(freqBiotypes$Freq/sum(freqBiotypes$Freq) * 100, digits=2)
  freqBiotypes$label = paste(freqBiotypes$Var1," (", freqBiotypes$percent,"%)", sep = "")
  rownames(freqBiotypes) <- freqBiotypes$Var1
  orderBiotypes <- c("mRNA", "pseudogenes", "miscellinous", "lncRNA", "rRNA", "miRNA")
  freqBiotypes = freqBiotypes[orderBiotypes, ]
  cols <- rainbow(length(orderBiotypes))
  freqBiotypes$Col <- cols
  freqBiotypes <- na.omit(freqBiotypes)
  
   png(file = file.path(experimentPath, paste0(names, "_geneBiotypes.png")),#delete #
       width=1500, height=1200, res = 300)#delete #
   pie(freqBiotypes$Freq, labels = freqBiotypes$label, col = freqBiotypes$Col, cex=0.5)#delete #
   dev.off()#delete #
  
  ## Keeping only mRNAs for heatmaps
  print("Keeping only mRNAs for heatmaps")
  keepBiotypes <- biotypes[biotypes %in% c("protein_coding")]
  subAnnoData <- annotateData[annotateData$gene_biotype %in% keepBiotypes, ]
  orderedAnnoData2 <- subAnnoData[order(abs(subAnnoData$log2FoldChange), decreasing = TRUE), ]
  orderedAnnoData3 <- subAnnoData[order(subAnnoData$log2FoldChange, decreasing = TRUE), ]
#  de_results[[paste0(names, "onlymrna")]] <- orderedAnnoData2
  
  ## Normalized Count matrix
  vsd <- vst(dds1, blind=FALSE)
  vsDE <- vsd[rownames(vsd) %in% rownames(resDE), ]
  
  inPCA <- as.data.frame(t(assay(vsDE)))
  colnames(inPCA) <- rownames(assay(vsDE))
  inPCA$groups <- dds_subset$groups
  
  samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)
  a1 <- autoplot(samplePCA, data = inPCA, colour = "groups", size = 3)#delete #
  a1 + theme(legend.position = "right", axis.text=element_text(size=9))
  a1 + geom_text(aes(label = rownames(inPCA)), nudge_x = 0.05, nudge_y = 0.1) + ggtitle(paste0(names))
   ggsave(file = file.path(experimentPath, paste0(names, "_samplePCA.png")), width = 10, height = 7) #delete #
  
   samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)#delete #
   a1 <- autoplot(samplePCA, data = inPCA, colour = "groups", size = 3)#delete #
   a1 + theme(legend.position = "right", axis.text=element_text(size=9))#delete #
   a1 + ggtitle(paste0(names))
   ggsave(file = file.path(experimentPath, paste0(names, "_samplePCA_no_labels.png")), width = 10, height = 7)#delete #
  
  ### Heatmaps 
  ## top 40 genes
   # Filter out rows with missing ENSEMBL values
   orderedAnnoData4 <- orderedAnnoData3[!is.na(orderedAnnoData3$ENSEMBL), ]
   # Order the data by padj
   orderedAnnoData4 <- orderedAnnoData4[order(orderedAnnoData4$padj), ]
   #only keep upregualted genes
   orderedAnnoData5 <- orderedAnnoData4[orderedAnnoData4$log2FoldChange > 0, ] 
   
   # Extract normalized counts for genes in the top 40 based on Z-scores
   mrnaMat <- assay(vsDE)[rownames(vsDE) %in% orderedAnnoData5$ENSEMBL, ]
   
   # Calculate standard deviations for each gene across samples
   row_std = rowSds(mrnaMat)
   
   # Calculate mean expression values for each gene across samples
   row_means = rowMeans(mrnaMat)
   
   # Standardize the gene expression matrix by subtracting the mean and dividing by the standard deviation
   zMatrix <- (mrnaMat - row_means) / row_std
   
   # Subset the standardized matrix based on the specified gene order
   zMatrix <- zMatrix[orderedAnnoData5$ENSEMBL, ]
   
   # Assign external gene names as row names in the standardized matrix
   rownames(zMatrix) <- orderedAnnoData5$external_gene_name
   
   # Create a data frame for sample grouping information
   colAnnot <- data.frame(groups = dds_subset$groups)
   
   # Assign column names in the annotation data frame to match the column names in the standardized matrix
   row.names(colAnnot) <- colnames(zMatrix)
   
   # Print the number of rows in zMatrix before heatmap generation
   print(paste("Number of rows in zMatrix for", names, ":", nrow(zMatrix)))
  
   # Generate a heatmap visualization for the top 40 genes (if more than 40 genes exist)
   if(nrow(zMatrix) > 40){
     # Define the file path for saving the heatmap image
     hmapFile <- file.path(experimentPath, paste0(names, "_top40.pdf"))
     
     # Create a heatmap using the top 40 genes, with additional clustering and annotation options
     p1 <- pheatmap(zMatrix[1:40,], show_rownames=TRUE, cluster_cols = TRUE, 
                    cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                    cutree_cols = 2, fontsize_row = 7)
     
     # Save the generated heatmap as an image file
     
     save.heatmap.pdf(x = p1, filename = hmapFile)
     
   }
   

}

write.csv(DEgene_summary_all_biotype, file.path(experimentPath, "DEgene_summary_all_biotype.csv"), row.names = FALSE)
write.csv(DEgene_summary_mRNA, file.path(experimentPath, "DEgene_summary_mRNA.csv"), row.names = FALSE)

write.xlsx(x=de_results, file = file.path(experimentPath, paste0("sig_DEG.xlsx")),rowNames = TRUE)

