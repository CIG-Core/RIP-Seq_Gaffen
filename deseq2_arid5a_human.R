# Differential Expression Analysis 

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
})

## Sourcing
source("/ix/cigcore/utils/code/save.heatmap.R")

## Set paths
dataPath <- "/ix/cigcore/proj/gaffen/data"
resPath <- gsub(pattern="data", replacement="result", dataPath)
experiment <- "03-21-2023-DE_analysis_deseq2_0.05"
experimentPath <- file.path(resPath, "RIP_seq_Arid5a", "RIPseq_arid5aVSigg", experiment)
dir.create(experimentPath, showWarnings = TRUE, recursive = TRUE)
webgestalt_path <- file.path(experimentPath, "03-21-2023-pathway_analysis_GSEA")
dir.create(webgestalt_path, recursive=TRUE)

## Read data file 
coldata <- read.csv(file.path(dataPath, "MetaData_Arid5aRIP.csv"), row.names = 1)
count <- read.table(file.path(resPath, "RIP_seq_Arid5a", "RIPseq_arid5aVSigg", "mapped_star", "rsem", "expected_counts.txt"))

coldata$sample <- rownames(coldata) 
coldata$groups <- gsub("\\.", "_", coldata$groups)
coldata$groups <- as.factor(coldata$groups)
rawData <- count[, rownames(coldata)]

all(colnames(rawData) %in% rownames(coldata)) # checking to see if sample names match between the two files
all(colnames(rawData) == rownames(coldata)) # checking to see if the sample names are in the same order in both files

## Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(rawData), colData = coldata, design = ~groups)
dds$groups 

## Make list with conditions and 4 comparisons
comparisons <- list(c("IgG_IL17","Arid5a_IL17"), 
                    c("IgG_Med","Arid5a_Med"), 
                    c("IgG_Med", "IgG_IL17"), 
                    c("Arid5a_Med", "Arid5a_IL17"))

de_results <- list() 

for(c1 in 1:length(comparisons)){
  names <- paste0(comparisons[[c1]][1], "_vs_", comparisons[[c1]][2], sep ="")
  metadata <- subset(colData(dds), groups %in% comparisons[[c1]])  #subset the colData
  rawcount <- counts(dds)[, metadata$sample] #subset the countData
  
  # Create DESeqDataSet object with colData and countData subset
  dds_subset <- DESeqDataSetFromMatrix(countData = rawcount, colData = metadata, design = ~groups)
  
  new_level <- print(comparisons[[c1]][1])
  dds_subset$groups <- relevel(dds_subset$groups, ref = new_level)
  
  dds1 <- DESeq(dds_subset)
  res <- results(dds1)
  print("DESeq done")
  
  ## Volcano Plot
  EnhancedVolcano(res, lab = rownames(res),  x = 'log2FoldChange', y = 'pvalue', 
                  labCol = 'black', labFace = 'bold', boxedLabels = TRUE,
                  labSize = 3.0, title = paste0(names),
                  drawConnectors = TRUE, widthConnectors = 1.0,
                  colConnectors = 'black')
  
  # ggsave(file = file.path(experimentPath, paste0(names, "_volcano_plot.png")))
  
  ## DE p.adj <0.05
  res2 <- res[complete.cases(res), ]
  res2 <- res2[order(res2$padj), ]
  
  resDE <- as.data.frame(res2[which(res2$padj < 0.05), ])
  ensemblIDs <- rownames(resDE)
  de_results[[paste0(names, "_", "DeG")]] <- resDE
  
  resDE_Up <- resDE[resDE$log2FoldChange > 0, ]
  resDE_Dn <- resDE[resDE$log2FoldChange < 0, ]
  # de_results[[paste0(names,"_", "DeU")]] <- resDE_Up
  # de_results[[paste0(names,"_", "DeD")]] <- resDE_Dn

  ## Get biotype information
  ensemblMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  g2g <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "gene_biotype",
                            "entrezgene_id", "ensembl_gene_id_version"), filters = "ensembl_gene_id_version",
               values = ensemblIDs, mart = ensemblMart)

  # write.csv(g2g, file.path(experimentPath, paste0(names, "_g2g.csv")), row.names = F)
  # g2g <- read.csv(file.path(experimentPath, "Arid5a_Med_vs_Arid5a_IL17_g2g.csv"))
  
  g2g2 <- g2g[!duplicated(g2g$ensembl_gene_id_version), ]
  annotateData <- merge(resDE, g2g2, by.x=0, by.y="ensembl_gene_id_version")
  colnames(annotateData)[colnames(annotateData) == "Row.names"] <- "ENSEMBL"
  de_results[[paste0(names, "wthBtype")]] <- annotateData
  
  # write.table(annotateData[ , c(8, 3)],
  #             file = file.path(experimentPath, paste0(names, "_gsea.rnk")),
  #             quote=FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)
  
  biotypes <- unique(g2g2$gene_biotype)
  pseudoBiotypes <- biotypes[biotypes %in% (str_subset(biotypes, pattern = "_pseudogene"))]
  otherBiotypes <- biotypes[!biotypes %in% c(pseudoBiotypes, "protein_coding", "lncRNA", "rRNA", "miRNA")]

  annotateData[which(annotateData$gene_biotype %in% pseudoBiotypes), "grouped_biotypes"] <- "pseudogenes"
  annotateData[which(annotateData$gene_biotype %in% otherBiotypes), "grouped_biotypes"] <- "miscellinous"
  annotateData[which(annotateData$gene_biotype %in% c("protein_coding")), "grouped_biotypes"] <- "mRNA"
  annotateData[which(!annotateData$gene_biotype %in% c("protein_coding", pseudoBiotypes, otherBiotypes)),
               "grouped_biotypes"] <- annotateData$gene_biotype[!annotateData$gene_biotype %in% c("protein_coding", pseudoBiotypes, otherBiotypes)]

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
  
  # png(file = file.path(experimentPath, paste0(names, "_geneBiotypes.png")),
  #     width=1500, height=1200, res = 300)
  # pie(freqBiotypes$Freq, labels = freqBiotypes$label, col = freqBiotypes$Col, cex=0.5)
  # dev.off()

  ## Keeping only mRNAs for heatmaps
  print("Keeping only mRNAs for heatmaps")
  keepBiotypes <- biotypes[biotypes %in% c("protein_coding")]
  subAnnoData <- annotateData[annotateData$gene_biotype %in% keepBiotypes, ]
  orderedAnnoData2 <- subAnnoData[order(abs(subAnnoData$log2FoldChange), decreasing = TRUE), ]
  orderedAnnoData3 <- subAnnoData[order(subAnnoData$log2FoldChange, decreasing = TRUE), ]
  de_results[[paste0(names, "onlymrna")]] <- orderedAnnoData2
  
  ## Normalized Count matrix
  vsd <- vst(dds1, blind=FALSE)
  vsDE <- vsd[rownames(vsd) %in% rownames(resDE), ]
  
  inPCA <- as.data.frame(t(assay(vsDE)))
  colnames(inPCA) <- rownames(assay(vsDE))
  inPCA$groups <- dds_subset$groups
   
  samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)
  a1 <- autoplot(samplePCA, data = inPCA, colour = "groups", size = 3)
  a1 + theme(legend.position = c(0.85, 0.1), axis.text=element_text(size=12))
  a1 + geom_text(aes(label = rownames(inPCA)), nudge_y = 0.1) + ggtitle(paste0(names))
  # ggsave(file = file.path(experimentPath, paste0(names, "_samplePCA.png")))

  # samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)
  # a1 <- autoplot(samplePCA, data = inPCA, colour = "groups", size = 3)
  # a1 + theme(legend.position = "right", axis.text=element_text(size=12))
  # ggsave(file = file.path(experimentPath, paste0(names, "_samplePCA_no_labels.png")))
  
  ## Heatmaps 
  ## top 50 genes
  orderedAnnoData4 <- orderedAnnoData3[orderedAnnoData3$external_gene_name !="" , ]
  #orderedAnnoData4 <- data.frame(orderedAnnoData4[order(orderedAnnoData4[,"padj"]) , ])
  
  mrnaMat <- assay(vsDE)[rownames(vsDE) %in% orderedAnnoData4$ENSEMBL, ]
  row_std = rowSds(mrnaMat)
  row_means = rowMeans(mrnaMat)
  zMatrix <- (mrnaMat - row_means) / row_std
  zMatrix <- zMatrix[orderedAnnoData4$ENSEMBL, ]
  rownames(zMatrix) <- orderedAnnoData4$external_gene_name
  colAnnot <- data.frame(groups = dds_subset$groups)
  row.names(colAnnot) <- colnames(zMatrix)
  if(nrow(zMatrix) > 50){
  hmapFile <- file.path(experimentPath, paste0(names, "_top50.png"))
  pOrder <- pheatmap(zMatrix[1:50,], show_rownames=TRUE, cluster_cols = TRUE, 
                     cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                     cutree_cols = 2, fontsize_row = 7)
  o1 <- pOrder$tree_col[["order"]]
  o1rev <- c(o1[4:6], o1[1:3])
  p1 <- pheatmap(zMatrix[1:50, o1rev], show_rownames=TRUE, cluster_cols = FALSE,
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2,
                 fontsize_row = 8)
  
  #save.heatmap(x = p1, filename = hmapFile)
  
  gene_list <- annotateData[ , c(8, 3)]
  # wgr_out <- WebGestaltR(enrichMethod = "GSEA", organism = "hsapiens", enrichDatabase = "pathway_Reactome", 
  #                        interestGene = gene_list, interestGeneType = "genesymbol", 
  #                        referenceSet = "illumina_humanref_8_v3", minNum = 5, maxNum = 2000,
  #                        reportNum = 40, isOutput = TRUE, outputDirectory = webgestalt_path,
  #                        projectName = names)
  
  if(names == "Arid5a_Med_vs_Arid5a_IL17"){
  # Arid5a_Med versus Arid5a_IL17
  ## Identifying Translation "GO:0006412" and inflammatory response genes "GO:0002437"
  print("Identifying translation related genes")
  
  ## Translation "GO:0006412"
  transGeneList <- mget(c("GO:0006412"),org.Hs.egGO2ALLEGS)
  transGenes <- unique(transGeneList[[1]])
  inBoth <- intersect(transGenes, orderedAnnoData2$entrezgene_id)
  subTransData <- orderedAnnoData2[orderedAnnoData2$entrezgene_id %in% inBoth, ]
  
  
  ## Reactome "R-HSA-72766"
  transGeneList2 <- event2Ids(event.id = "R-HSA-72766") # can only input id here
  transGenes2 <- unique(transGeneList2$geneSymbol)
  transInBoth <- intersect(transGenes2, orderedAnnoData2$external_gene_name)
  subTransData2 <- orderedAnnoData2[orderedAnnoData2$external_gene_name %in% transInBoth, ]
  
  print("Identifying inflammatory response related genes")
  inflamGeneList <- mget(c("GO:0002437"), org.Hs.egGO2ALLEGS) 
  inflamGenes <- unique(inflamGeneList[[1]])
  inBothInflam <- intersect(inflamGenes, orderedAnnoData2$entrezgene_id)
  subInflamData <- orderedAnnoData2[orderedAnnoData2$entrezgene_id %in% inBothInflam, ]
  
  print("Identifying IL-17 inflammatory response related genes")
  il17Genes <- c("CEBPB", "CEBPD", "CXCL2", "IGF2BP1", "IL13RA1", "IL6", "IL6ST", "NFKBIZ",
                 "ZC3H12A", "ZC3H12C", "ZC3H3", "ZC3H7A")
  subIL17Data <- orderedAnnoData2[orderedAnnoData2$external_gene_name %in% il17Genes, ]
 
  ## Heatmaps 
  ## Translation related genes, GO:0006412
  subTransData$row_std = rowSds(assay(vsDE)[rownames(vsDE) %in% subTransData$ENSEMBL, ])
  subTransData$row_means = rowMeans(assay(vsDE)[rownames(vsDE) %in% subTransData$ENSEMBL, ])
  zMatrixTrans <- (assay(vsDE)[rownames(vsDE) %in% subTransData$ENSEMBL, ] - subTransData$row_means) / subTransData$row_std
  all(subTransData$ENSEMBL == rownames(zMatrixTrans)) 
  zMatrixTrans <- zMatrixTrans[subTransData$ENSEMBL, ]
  all(subTransData$ENSEMBL == rownames(zMatrixTrans)) 
  rownames(zMatrixTrans) <- subTransData$external_gene_name
  removeGenes <- c("ZNF540", "ZC3H12D", "DDX25", "EIF5B")
  zMatrixTrans <- zMatrixTrans[!(rownames(zMatrixTrans)%in% removeGenes), ]
  
  hmapFile2 <- file.path(experimentPath, "Arid5a_Med_vs_Arid5a_IL17_translation.png")
  pOrder <- pheatmap(zMatrixTrans, show_rownames=TRUE, cluster_cols = TRUE, 
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 cutree_cols = 2, fontsize_row = 5)
  o1 <- pOrder$tree_col[["order"]]
  o1rev <- c(o1[4:6], o1[1:3])
  p2 <- pheatmap(zMatrixTrans[ , o1rev], show_rownames=TRUE, cluster_cols = FALSE, 
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 fontsize_row = 7)
  
  #save.heatmap(x = p2, filename = hmapFile2)
  
  ## Inflammatory Response Genes, GO:0002437
  subInflamData$row_std = rowSds(assay(vsDE)[rownames(vsDE) %in% subInflamData$ENSEMBL, ])
  subInflamData$row_means = rowMeans(assay(vsDE)[rownames(vsDE) %in% subInflamData$ENSEMBL, ])
  zMatrixInflam <- (assay(vsDE)[rownames(vsDE) %in% subInflamData$ENSEMBL, ] - subInflamData$row_means) / subInflamData$row_std
  all(subInflamData$ENSEMBL == rownames(zMatrixInflam)) 
  zMatrixInflam <- zMatrixInflam[subInflamData$ENSEMBL, ]
  all(subInflamData$ENSEMBL == rownames(zMatrixInflam)) 
  rownames(zMatrixInflam) <- subInflamData$external_gene_name
  
  hmapFile3 <- file.path(experimentPath, "Arid5a_Med_vs_Arid5a_IL17_inflammation.png")
  pOrder <- pheatmap(zMatrixInflam, show_rownames=TRUE, cluster_cols = TRUE, 
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 cutree_cols = 2, fontsize_row = 12)
  o1 <- pOrder$tree_col[["order"]]
  o1rev <- c(o1[4:6], o1[1:3])
  p3 <- pheatmap(zMatrixInflam[ , o1rev], show_rownames=TRUE, cluster_cols = FALSE, 
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 fontsize_row = 12)
  
  
  #save.heatmap(x = p3, filename = hmapFile3)
  
  ## IL17 related Inflammatory Response Genes
  subIL17Data$row_std = rowSds(assay(vsDE)[rownames(vsDE) %in% subIL17Data$ENSEMBL, ])
  subIL17Data$row_means = rowMeans(assay(vsDE)[rownames(vsDE) %in% subIL17Data$ENSEMBL, ])
  zMatrixIL17 <- (assay(vsDE)[rownames(vsDE) %in% subIL17Data$ENSEMBL, ] - subIL17Data$row_means) / subIL17Data$row_std
  all(subIL17Data$ENSEMBL == rownames(zMatrixIL17))
  zMatrixIL17 <- zMatrixIL17[subIL17Data$ENSEMBL, ]
  all(subIL17Data$ENSEMBL == rownames(zMatrixIL17))
  rownames(zMatrixIL17) <- subIL17Data$external_gene_name
  
  hmapFile4 <- file.path(experimentPath, "Arid5a_Med_vs_Arid5a_IL17_il17.png")
  pOrder <- pheatmap(zMatrixIL17, show_rownames=TRUE, cluster_cols = TRUE, cellwidth = 24,
                 cellheight = 8, cluster_rows = TRUE, annotation_col = colAnnot, 
                 cutree_rows = 2, cutree_cols = 2, fontsize_row = 8)
  o1 <- pOrder$tree_col[["order"]]
  o1rev <- c(o1[4:6], o1[1:3])
  p4 <- pheatmap(zMatrixIL17[ , o1rev], show_rownames=TRUE, cluster_cols = FALSE, cellwidth = 24,
                 cellheight = 8, cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 fontsize_row = 7)
  
  #save.heatmap(x = p4, filename = hmapFile4)
  
  ## Translation related genes, Reactome:R-HSA-72766
  subTransData2$row_std = rowSds(assay(vsDE)[rownames(vsDE) %in% subTransData2$ENSEMBL, ])
  subTransData2$row_means = rowMeans(assay(vsDE)[rownames(vsDE) %in% subTransData2$ENSEMBL, ])
  zMatrixTrans2 <- (assay(vsDE)[rownames(vsDE) %in% subTransData2$ENSEMBL, ] - subTransData2$row_means) / subTransData2$row_std
  all(subTransData2$ENSEMBL == rownames(zMatrixTrans2)) 
  zMatrixTrans2 <- zMatrixTrans2[subTransData2$ENSEMBL, ]
  all(subTransData2$ENSEMBL == rownames(zMatrixTrans2)) 
  rownames(zMatrixTrans2) <- subTransData2$external_gene_name
  
  hmapFile5 <- file.path(experimentPath, "Arid5a_Med_vs_Arid5a_IL17_translation_reactome.png")
  pOrder <- pheatmap(zMatrixTrans2, show_rownames=TRUE, cluster_cols = TRUE, 
                     cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                     cutree_cols = 2, fontsize_row = 7)
  o1 <- pOrder$tree_col[["order"]]
  o1rev <- c(o1[4:6], o1[1:3])
  
  p5 <- pheatmap(zMatrixTrans2[ , o1rev], show_rownames=TRUE, cluster_cols = FALSE, 
                 cluster_rows = TRUE, annotation_col = colAnnot, cutree_rows = 2, 
                 fontsize_row = 7)
  
  #save.heatmap(x = p5, filename = hmapFile5)
  
    }
  }
}

# names(de_results)[c(11, 12)] <- c("Arid5a_Med_vs_Arid5a_IL17wiBtyp", "Arid5a_Med_vs_Arid5a_IL17_mrna")
# write.xlsx(x=de_results, file = file.path(experimentPath, paste0("sig_DEG.xlsx")),rowNames = TRUE)



