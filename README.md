# RIP-Seq_Gaffen
Mouse data: analyzed by Sissi from 2023/10/12 to 2023/11/16:
1) Preprocessing analysis including 1st QC Raw data (MultiQC), trimming, and 2nd QC trimmed data (MultiQC). Need to mention that Sissi modified trimmomatic adapter to TruSeq3-PE-2.fa in trimmomatic.job, in order to remove the "Illumina Universal Adapter" and "SOLID small RNA Adapter" from the raw sequence read.
2) Kallisto the trimmed data, the index/annotation used mouse database (mm10), since it didn't have good alignment, switch to the STAR-DEseq2.
3) Counted reads using both RSEM and FeatureCount, using RSEM for DEseq2 since it had better alignment scores.
4) DEseq2 using code 'deseq2_arid5a_mice.R', the final result is under cluster "11-16-2023-DE_Analysis_deseq2_RSEM_remove_IgG_1_IL17" after removing one outlier (IgG_1_IL17).
5) Intersect for significant DE genes between Human data (data: "RIPseq_arid5aVSigg) and Mouse data  (data: "Arid5a_CLIP-Seq_2023) after orthologs conversions, and code can be found in 'intersect.R'.
6) GSEA analysis was done through WebGestalt-Web version (Minimum number of IDs in the category: 15 and Significance Level: FDR < 0.25).

All results and PowerPoints for this analysis were shared with Sarah and Yang. 

Result pathways:
1) DE analysis result:
/ix/cigcore/proj/gaffen/result/Arid5a_CLIP-Seq_2023/11-16-2023-DE_Analysis_deseq2_RSEM_remove_IgG_1_IL17
2) Intersect gene summary, gene lists, and venn diagram: /ix/cigcore/proj/gaffen/result/Arid5a_CLIP-Seq_2023/11-10-2023-DE_Analysis_deseq2_RSEM_remove_IgG_1_IL17/intersection
3) rnk.file for GSEA analysis for genes (all upregulated genes):
/ix/cigcore/proj/gaffen/result/Arid5a_CLIP-Seq_2023/11-13-2023-GSEA-ORA
4) txt.file for ORA analysis for intersect genes (significant upregulated intersect genes):/ix/cigcore/proj/gaffen/result/Arid5a_CLIP-Seq_2023/11-13-2023-ORA-intersect