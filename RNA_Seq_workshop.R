####################
#RNA_seq_workshop  #
#Windows 11        #
#R version 4.3.3   #
#script by Omar K. #
####################

#Install library
#Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#Install DESeq2 using BiocManager
BiocManager::install("DESeq2", force = TRUE)


#Install tidyverse and ggplot2 from CRAN
install.packages("tidyverse")
install.packages("ggplot2")


#Install pheatmap from CRAN
install.packages("pheatmap")


#import_library
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
#Load_SRA_RNA_matrix
RNA = read.csv("/R/transcriptome_SRA_name.csv", header = TRUE)
#check_the_dimension
dim(RNA)
#Remove_duplicates
RNA = RNA[!duplicated(RNA$gene_name), ]
#Assign_each_row_to_gene_name
row.names(RNA) = RNA$gene_name
#remove_gene_name_column
RNA$gene_name = NULL
#obtain_the_data_into_matrix
RNA = as.matrix(RNA)
#convert_each_count_value_to_integer_number
genes = row.names(RNA)
RNA = apply(RNA,2,as.integer)
row.names(RNA) = genes

#Load_meta_data
meta = read.csv("/R/first_meta_data.csv", row.names = "X")

#Filter_other_virus_Cases
meta = meta %>% filter(infection != "other_virus", Age_cat == "child")
meta = meta[order(meta$Age_cat), ]

#sort_RNA_data_with_meta_data_sample_id
RNA_data = subset(as.matrix(RNA), select = as.character(meta$Run))



#Assign_condition_one_to_SARS_COV2"SC2"patient_and_condtion2_is_control
cond1 = "SC2"
cond2 = "no_virus"

#Create_deseq2_object
dds = DESeqDataSetFromMatrix(countData = RNA_data, colData = meta, design = ~ infection)

#Run_Deseq
dds = DESeq(dds)



#########
#Differential_analysis_statistics
res = results(dds, contrast = c("infection", cond1 ,cond2))
#Remove_NaN
res=as.data.frame(res[complete.cases(res), ])
#Add_threeshould
deseq_deg = res[res$padj < 0.05 & abs(res$log2FoldChange) > 2, ]
#Create_norm_transformation_df 
ntd = normTransform(dds)


#extrace_sigDif_genes
significant_gene_names = rownames(deseq_deg)
test_matrix = as.matrix(assay(ntd)[significant_gene_names[1:25],])


#Create_df_object_contain_samples_annotation
df <- as.data.frame(colData(dds)[,c("infection")])
colnames(df) = "Infection"
rownames(df) <- colnames(test_matrix)


#Creat_heatmap_using_ntd_object
pheatmap(test_matrix,cluster_rows = TRUE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         cluster_cols = TRUE, annotation_col= df)

#Plot_PCA_depending_on_Age_infection_and_both
plotPCA(ntd, intgroup=c("infection"))


#Plot_volcano_plot
par(mfrow=c(1,1))
res.k = res
with(res.k, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,4)))
with(subset(res.k, pvalue< 0.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.k, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
