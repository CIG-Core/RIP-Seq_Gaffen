# Load the required libraries
library(readxl)
library(writexl)
library(dplyr)
library('orthogene')
# Venn plot
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(grid)
library(openxlsx)
library(GeneOverlap) #overlap significance analysis

# Define the file path to your Excel file
proj_path <- "/ix/cigcore/proj/gaffen/result/Arid5a_CLIP-Seq_2023/11-10-2023-DE_Analysis_deseq2_RSEM_remove_IgG_1_IL17"
sissi_path <- file.path(proj_path, "sig_DEG.xlsx")
oge_path <- "/ix/cigcore/proj/gaffen/result/RIP_seq_Arid5a/RIPseq_arid5aVSigg/03-21-2023-DE_analysis_deseq2_0.05/sig_DEG.xlsx"
all_biotype_path <- file.path(proj_path, "intersection", "Intersect_DEgenes_All_Biotypes.xlsx")
mRNA_path <- file.path(proj_path, "intersection", "Intersect_DEgenes_Only_mRNA.xlsx")
save_path <- file.path(proj_path, "intersection")
intersect_mRNA_path <- file.path(proj_path, "intersection", "Intersect_DEgene_List_mRNA.xlsx")
intersect_all_biotype_path <- file.path(proj_path, "intersection", "Intersect_DEgene_List_All_Biotypes.xlsx")
mouse_mRNA_path <- file.path(proj_path, "intersection", "mouse_mRNA.xlsx")
human_mRNA_path <- file.path(proj_path, "intersection", "human_mRNA.xlsx")
mouse_all_biotype_path <- file.path(proj_path, "intersection", "mouse_All_Biotypes.xlsx")
significant_path <- file.path(save_path, "GeneOverlapResults.rds")
# dir.create(significant_path, recursive=TRUE)
intersect_up_all_biotype_path <- file.path(proj_path, "intersection", "Intersect_UP_DEgene_List_All_Biotypes.xlsx")
###Extract data (only mRNA)
# Read the data from human's the specific sheet and column
human_ar17_ig17_mRNA <- read_excel(oge_path, sheet = "IgG_IL17_vs_Arid5a_IL17onlymrna", col_names = TRUE)
human_armed_igmed_mRNA <- read_excel(oge_path, sheet = "IgG_Med_vs_Arid5a_Medonlymrna", col_names = TRUE)
human_ig17_igmed_mRNA <- read_excel(oge_path, sheet = "IgG_Med_vs_IgG_IL17onlymrna", col_names = TRUE)
human_ar17_armed_mRNA <- read_excel(oge_path, sheet = "Arid5a_Med_vs_Arid5a_IL17_mrna", col_names = TRUE)

# Read the data from mouse's the specific sheet and column
mouse_ar17_ig17_mRNA <- read_excel(sissi_path, sheet = "Arid5a_IL17_vs_IgG_IL17_mDeG", col_names = TRUE)
mouse_armed_igmed_mRNA <- read_excel(sissi_path, sheet = "Arid5a_Med_vs_IgG_Med_mDeG", col_names = TRUE)
mouse_ig17_igmed_mRNA <- read_excel(sissi_path, sheet = "IgG_IL17_vs_IgG_Med_mDeG", col_names = TRUE)
mouse_ar17_armed_mRNA <- read_excel(sissi_path, sheet = "Arid5a_IL17_vs_Arid5a_Med_mDeG", col_names = TRUE)

###Intersect for All Biotypes
# Read the data from Oge's the specific sheet and column
human_ar17_ig17 <- read_excel(oge_path, sheet = "IgG_IL17_vs_Arid5a_IL17wthBtype", col_names = TRUE)
human_armed_igmed <- read_excel(oge_path, sheet = "IgG_Med_vs_Arid5a_MedwthBtype", col_names = TRUE)
human_ig17_igmed <- read_excel(oge_path, sheet = "IgG_Med_vs_IgG_IL17wthBtype", col_names = TRUE)
human_ar17_armed <- read_excel(oge_path, sheet = "Arid5a_Med_vs_Arid5a_IL17wiBtyp", col_names = TRUE)


# Read the data from Sissi's the specific sheet and column
mouse_ar17_ig17 <- read_excel(sissi_path, sheet = "Arid5a_IL17_vs_IgG_IL17AllBT", col_names = TRUE)
mouse_armed_igmed <- read_excel(sissi_path, sheet = "Arid5a_Med_vs_IgG_MedAllBT", col_names = TRUE)
mouse_ig17_igmed <- read_excel(sissi_path, sheet = "IgG_IL17_vs_IgG_MedAllBT", col_names = TRUE)
mouse_ar17_armed <- read_excel(sissi_path, sheet = "Arid5a_IL17_vs_Arid5a_MedAllBT", col_names = TRUE)


### Only mRNA ####

# Create a list of data frames for mRNA
data_frames_mRNA <- list(
  "Arid5a_IL17_vs_IgG_IL17" = list(mouse_ar17_ig17_mRNA, human_ar17_ig17_mRNA),
  "Arid5a_Med_vs_IgG_Med" = list(mouse_armed_igmed_mRNA, human_armed_igmed_mRNA),
  "IgG_IL17_vs_IgG_Med" = list(mouse_ig17_igmed_mRNA, human_ig17_igmed_mRNA),
  "Arid5a_IL17_vs_Arid5a_Med" = list(mouse_ar17_armed_mRNA, human_ar17_armed_mRNA)
)

# Initialize a list to store common gene counts for mRNA
common_DEg_counts_mRNA <- list()

# Create a list to store summary tables
summary_tables_mRNA <- list()
venn_diagrams_mRNA <- list()
intersection_genes_mRNA <- list()
mouse_mRNA <- list()
human_mRNA <- list()

# Perform the significance analysis using GeneOverlap
go_objs <- list()

# Loop through the data frames for mRNA and perform the necessary operations
for (df_name_mRNA in names(data_frames_mRNA)) {
  mouse_mRNA <- data_frames_mRNA[[df_name_mRNA]][[1]]
  human_mRNA <- data_frames_mRNA[[df_name_mRNA]][[2]]
  #  names <- paste0(data_frames_mRNA[[c1]])
  
  # Convert mm for mRNA
  mouse_mRNA_orthologs <- orthogene::convert_orthologs(mouse_mRNA, gene_input = "ENSEMBL", 
                                                   gene_output = "columns", input_species = "mmusculus",
                                                   output_species = "hsapiens", non121_strategy = "keep_both_species")
  
  # Perform an inner join on the common column "ensembl_gene_id" for mRNA
  ## Notes: remove NA: check the data, whther they have length-space
  # check the dataset before/after removing NAs
  mouse_mRNA_no_na <- mouse_mRNA[!is.na(mouse_mRNA$ENSEMBL), ]
  mouse_mRNA2 <- mouse_mRNA_no_na[!duplicated(mouse_mRNA$ENSEMBL), ]

  
  # Filter out rows where ortholog_gene is not missing (not NA)
  mouse_mRNA_orthologs_no_na <- mouse_mRNA_orthologs[!is.na(mouse_mRNA_orthologs$ortholog_gene), ]
  # Subsetting based on unique rows
  mouse_mRNA_orthologs2 <- mouse_mRNA_orthologs_no_na[!duplicated(mouse_mRNA_orthologs_no_na$ortholog_gene), ]
  
  # Filter out rows where external_gene_name is not missing (not NA)
  human_mRNA_no_na <- human_mRNA[!is.na(human_mRNA$external_gene_name), ]
  # Subsetting based on unique rows
  human_mRNA2 <- human_mRNA_no_na[!duplicated(human_mRNA_no_na$external_gene_name), ]
  
  
  common_genes_df_mRNA <- mouse_mRNA_orthologs2[mouse_mRNA_orthologs2$ortholog_gene %in% human_mRNA2$external_gene_name, ]
  
  # Count the number of common genes for mRNA
  common_DEg_count_mRNA <- nrow(common_genes_df_mRNA)
  common_DEg_counts_mRNA[[df_name_mRNA]] <- common_DEg_count_mRNA
  
  # Create a summary table for this data frame for mRNA
  summary_table_mRNA <- data.frame(
    Data_Frame = c("RIPseq_arid5aVSigg", "Arid5a_CLIP-Seq_2023", "Common Genes"),
    Number_of_Genes = c(nrow(human_mRNA2), nrow(mouse_mRNA2), common_DEg_count_mRNA)
  )
  
  summary_tables_mRNA[[df_name_mRNA]] <- summary_table_mRNA
  
  
  # List the intersection genes for mRNA
  intersection_genes_mRNA[[df_name_mRNA]] <- data.frame(Intersection_Genes = common_genes_df_mRNA$ortholog_gene)
  
  mouse_mRNA[[df_name_mRNA]] <- mouse_mRNA_orthologs2
  human_mRNA[[df_name_mRNA]] <- human_mRNA2
  
  #Perform the significance analysis using GeneOverlap
  # Notes: check the number later
  total_universe_genes <- 21694 # The Ensembl gene set consists of 21,694 protein coding genes
  go_obj <- newGeneOverlap(mouse_mRNA_orthologs2$ortholog_gene, human_mRNA$external_gene_name, genome.size = total_universe_genes)
  go_obj <- testGeneOverlap(go_obj)
  
  # Extract the relevant information for each comparison
  listA_df <- data.frame(ListA = go_obj@listA)
  listB_df <- data.frame(ListB = go_obj@listB)
  intersection_df <- data.frame(Intersection = go_obj@intersection)
  union_df <- data.frame(Union = go_obj@union)
  genome_size_df <- data.frame(GenomeSize = go_obj@genome.size)
  pvalue_df <- data.frame(PValue = go_obj@pval)
  odds_ratio_df <- data.frame(OddsRatio = go_obj@odds.ratio)
  jaccard_index_df <- data.frame(JaccardIndex = go_obj@Jaccard)
  
  # Store the data frames in the list
  go_objs[[df_name_mRNA]] <- list(
    ListA = listA_df,
    ListB = listB_df,
    Intersection = intersection_df,
    Union = union_df,
    GenomeSize = genome_size_df,
    PValue = pvalue_df,
    OddsRatio = odds_ratio_df,
    JaccardIndex = jaccard_index_df
  )

  
  # Venn plot
  de_mrna_mouse <- mouse_mRNA_orthologs2$ortholog_gene
  de_mrna_human <- human_mRNA2$external_gene_name
  
  # Create Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(de_mrna_human),
    area2 = length(de_mrna_mouse),
    cross.area = common_DEg_count_mRNA,
    category = c("Human Data (mRNA)", "Mouse Data (mRNA)"),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 550,
    width = 550,
    resolution = 300,
    col = c("lightpink", 'lightblue'),  # Color for the category titles
    fontface = "bold",  # Font style for category titles
    cex = 1.5,  # Font size for category titles
    fill = c("lightpink", 'lightblue'),  # Fill color for the category titles
    label.fill = "white",  # Label color inside the category titles
    cross.category = "purple",  # Fill color for the set intersections
    cat.pos = 2,
    cat.cex = 1.5,
    cat.dist = 0.03
  )
  
  # Store the Venn diagram in the list
  venn_diagrams_mRNA[[df_name_mRNA]] <- venn.plot
  
  # Save the Venn diagram
  file_name <- paste0(df_name_mRNA, "_venn_mRNA.png")
  ggsave(file.path(save_path, file_name), plot = venn.plot, width = 10, height = 7)
  
  # show individual plot in the end of codes
  dev.off()
}

# Create an Excel workbook and write the summary tables as separate sheets
write_xlsx(summary_tables_mRNA, path = mRNA_path)
write_xlsx(intersection_genes_mRNA, path = intersect_mRNA_path)
write_xlsx(mouse_mRNA, path = mouse_mRNA_path)
write_xlsx(human_mRNA, path = human_mRNA_path)


# save as rds
saveRDS(go_objs, file = file.path(significant_path))



### ALL Biotypes###
# Create a list of data frames and comparisons
data_frames_comparisons <- list(
  "Arid5a_IL17_vs_IgG_IL17" = list(mouse_ar17_ig17, human_ar17_ig17),
  "Arid5a_Med_vs_IgG_Med" = list(mouse_armed_igmed, human_armed_igmed),
  "IgG_IL17_vs_IgG_Med" = list(mouse_ig17_igmed, human_ig17_igmed),
  "Arid5a_IL17_vs_Arid5a_Med" = list(mouse_ar17_armed, human_ar17_armed)
)

# Initialize a list to store common gene counts for mRNA
common_DEg_counts_comparison <- list()

# Create a list to store summary tables
summary_tables_comparison <- list()
venn_diagrams_comparison <- list()
intersection_genes_comparison <- list()
mouse_comparison <- list()
intersection_genes_upregulated <- list()
# Loop through the specified comparisons and perform the necessary operations
for (df_name_comparison in names(data_frames_comparisons)) {
  mouse_all <- data_frames_comparisons[[df_name_comparison]][[1]]
  human_all <- data_frames_comparisons[[df_name_comparison]][[2]]
  
  
  # Convert mm for the specified comparison
  #EnsemblBioMart --> check the mapping 
  mouse_all_orthologs <- orthogene::convert_orthologs(mouse_all, gene_input = "ENSEMBL", 
                                                         gene_output = "columns", input_species = "mmusculus",
                                                         output_species = "hsapiens", non121_strategy = "keep_both_species")
  
  # Perform an inner join on the common column "ensembl_gene_id" for the specified comparison
  mouse_all_no_na <- mouse_all[!is.na(mouse_all$ENSEMBL), ]
  mouse_all2 <- mouse_all_no_na[!duplicated(mouse_all$ENSEMBL), ]
  
  mouse_all_orthologs_no_na <- mouse_all_orthologs[!is.na(mouse_all_orthologs$ortholog_gene), ]
  mouse_all_orthologs2 <- mouse_all_orthologs_no_na[!duplicated(mouse_all_orthologs_no_na$ortholog_gene), ]
  
  # Filter out rows where external_gene_name is not missing (not NA)
  human_all_no_na <- human_all[!is.na(human_all$external_gene_name), ]
  human_all2 <- human_all_no_na[!duplicated(human_all_no_na$external_gene_name), ]
  
  common_genes_df_comparison <- mouse_all_orthologs2[mouse_all_orthologs2$ortholog_gene %in% human_all2$external_gene_name, ]
  
  
  # Count the number of common genes for all biotypes
  common_DEg_count_comparison <- nrow(common_genes_df_comparison)
  common_DEg_counts_comparison [[df_name_comparison]] <- common_DEg_count_comparison
  
  # Create a summary table for this comparison
  summary_table_comparison <- data.frame(
    Data_Frame = c("RIPseq_arid5aVSigg", "Arid5a_CLIP-Seq_2023", "Common Genes"),
    Number_of_Genes = c(nrow(human_all2), nrow(mouse_all2), common_DEg_count_comparison)
  )
  
  summary_tables_comparison [[df_name_comparison]] <- summary_table_comparison
  
  # List the intersection genes for mRNA
  intersection_genes_comparison[[df_name_comparison]] <- data.frame(Intersection_Genes = common_genes_df_comparison$ortholog_gene)
  
  # Mouse data with original geneid and orthlog geneid
  mouse_comparison[[df_name_comparison]] <- mouse_all_orthologs2

  # common gene for upregulated
  mouse_upregulated <- as.data.frame(mouse_all2[which(mouse_all2$log2FoldChange > 0), ])
  human_upregulated <- as.data.frame(human_all2[which(human_all2$log2FoldChange > 0), ])
  
  mouse_up_comparison <- orthogene::convert_orthologs(mouse_upregulated, gene_input = "ENSEMBL", 
                                                         gene_output = "columns", input_species = "mmusculus",
                                                         output_species = "hsapiens", non121_strategy = "keep_both_species")
  
  common_genes_df_upregulated <- mouse_up_comparison[mouse_up_comparison$ortholog_gene %in% human_upregulated$external_gene_name, ]
  
  intersection_genes_upregulated[[df_name_comparison]] <- data.frame(Intersection_Genes = common_genes_df_upregulated$ortholog_gene)
  
  
  # Venn plot
  de_comparison_mouse <- mouse_all_orthologs2$ortholog_gene
  de_comparison_human <- human_all2$external_gene_name
  
  # Create Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(de_comparison_human),
    area2 = length(de_comparison_mouse),
    cross.area = common_DEg_count_comparison,
    category = c("Human Data (All Biotypes)", "Mouse Data (All Biotypes)"),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 550,
    width = 550,
    resolution = 300,
    col = c("lightpink", 'lightblue'),  # Color for the category titles
    fontface = "bold",  # Font style for category titles
    cex = 1.5,  # Font size for category titles
    fill = c("lightpink", 'lightblue'),  # Fill color for the category titles
    label.fill = "white",  # Label color inside the category titles
    cross.category = "purple",  # Fill color for the set intersections
    cat.pos = 2,
    cat.cex = 1.5,
    cat.dist = 0.03
  )
  

  # Store the Venn diagram in the list
  venn_diagrams_comparison[[df_name_comparison]] <- venn.plot
  
  # Save the Venn diagram
  file_name <- paste0(df_name_comparison, "_venn_all_biotypes.png")
  ggsave(file.path(save_path, file_name), plot = venn.plot, width = 10, height = 7)
  
  # show individual plot in the end of codes
  dev.off()
}

write_xlsx(summary_tables_comparison, path = all_biotype_path)
write_xlsx(intersection_genes_comparison, path = intersect_all_biotype_path)
write_xlsx(mouse_comparison, path = mouse_all_biotype_path)
write_xlsx(intersection_genes_upregulated, path = intersect_up_all_biotype_path)

