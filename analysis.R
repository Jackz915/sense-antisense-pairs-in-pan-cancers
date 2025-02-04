
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BRCA",             
  data.category = "Transcriptome Profiling",  
  data.type = "Gene Expression Quantification",  
)

GDCdownload(query)
data <- GDCprepare(query)
sample_info_list=colData(data)@listData
sample_info=as.data.frame(sample_info_list[1:10])
sample_info=sample_info[sample_info$definition %in% c('Primary solid Tumor', 'Solid Tissue Normal'),]
sample_info=sample_info[order(sample_info$definition),]
table(sample_info$definition)

gene_info=as.data.frame(rowRanges(data))
data_expr_tpm <- assay(data,i = "tpm_unstrand")
data_expr_tpm <- data_expr_tpm[,sample_info$barcode]
rownames(data_expr_tpm) <- substring(rownames(data_expr_tpm),1,15)
data_expr_tpm <- log2(data_expr_tpm + 1)


putative_bidirectional %>%
  group_by(PairID) %>%
  summarise(
    sense_gene = Gene_ID[Strand == "1"],
    antisense_gene = Gene_ID[Strand == "-1"]
  ) -> pair_info


result_matrix <- matrix(NA, nrow = nrow(pair_info), ncol = ncol(data_expr_tpm))

for (i in 1:nrow(pair_info)) {
  sense_gene <- pair_info$sense_gene[[i]]
  antisense_gene <- pair_info$antisense_gene[[i]]
  
  if (sense_gene %in% rownames(data_expr_tpm) &&
      antisense_gene %in% rownames(data_expr_tpm)) {
    for (j in 1:ncol(data_expr_tpm)) {
      sense_gene_exp <- data_expr_tpm[sense_gene, j]
      antisense_gene_exp <- data_expr_tpm[antisense_gene, j]
      sum_ = sense_gene_exp + antisense_gene_exp
      if (sum_ == 0) {
        result_matrix[i, j] <- 1
      }
      else {
        result_matrix[i, j] <- sense_gene_exp / sum_
      }
    }
  }
}
rownames(result_matrix) <- pair_info$PairID
colnames(result_matrix) <- colnames(data_expr_tpm)
result_matrix_filter <- na.omit(result_matrix)



library(limma)
library(dplyr)
colnames(result_matrix_filter) <- c(rep("tumor", table(sample_info$sample_type)[1]),
                                    rep("normal", table(sample_info$sample_type)[2]))
list <- colnames(result_matrix_filter) %>% factor(., levels = c("tumor", "normal"), ordered = F)
head(list)
list <- model.matrix(~factor(list)+0)  
colnames(list) <- c("tumor", "normal")
df.fit <- lmFit(result_matrix_filter, list)  
df.matrix <- makeContrasts("tumor-normal", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,n = Inf, adjust = "fdr")
head(DEG)


DEG$neg_log_pval <- -log10(DEG$P.Value)
DEG$neg_log_adj_pval <- -log10(DEG$adj.P.Val)

ggplot(DEG, aes(x = logFC, y = neg_log_pval)) +
  geom_point(aes(color = (adj.P.Val < 0.05)), size = 3) +  
  scale_color_manual(values = c("gray", "red")) +  
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log Fold Change (logFC)",
    y = "-log10(P-value)",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )




pair_info <- pair_info[pair_info$PairID %in% rownames(result_matrix_filter), ]
correlation_change_matrix <- matrix(NA, nrow = nrow(pair_info), ncol = 4)
normal_samples <- sample_info$barcode[sample_info$definition == "Solid Tissue Normal"]
tumor_samples <- sample_info$barcode[sample_info$definition == "Primary solid Tumor"]

for (i in 1:nrow(pair_info)) {
  sense_gene <- pair_info$sense_gene[[i]]
  antisense_gene <- pair_info$antisense_gene[[i]]
  
  protein_coding_exp_normal <- data_expr_tpm[sense_gene, normal_samples]
  lncRNA_exp_normal <- data_expr_tpm[antisense_gene, normal_samples]
  
  protein_coding_exp_tumor <- data_expr_tpm[sense_gene, tumor_samples]
  lncRNA_exp_tumor <- data_expr_tpm[antisense_gene, tumor_samples]
  
  cor_normal <- cor.test(protein_coding_exp_normal, lncRNA_exp_normal, method = "pearson")
  cor_normal_p <- as.numeric(cor_normal[3])
  cor_normal_correlation <- as.numeric(cor_normal[4])
  
  cor_tumor <- cor.test(protein_coding_exp_tumor, lncRNA_exp_tumor, method = "pearson")
  cor_tumor_p <- as.numeric(cor_tumor[3])
  cor_tumor_correlation <- as.numeric(cor_tumor[4])
  
  correlation_change_matrix[i, 1] <- cor_normal_correlation
  correlation_change_matrix[i, 2] <- cor_normal_p
  correlation_change_matrix[i, 3] <- cor_tumor_correlation
  correlation_change_matrix[i, 4] <- cor_tumor_p

}

colnames(correlation_change_matrix) <- c('normal_cor', 'normal_p', 'tumor_cor', 'tumor_p')
result <- cbind(pair_info, correlation_change_matrix, DEG[pair_info$PairID,])
result <- na.omit(result)

result$PairType <- NA 
for (i in 1:nrow(result)) {
  normal_p <- result$normal_p[i]
  tumor_p <- result$tumor_p[i]
  normal_cor <- result$normal_cor[i]
  tumor_cor <- result$tumor_cor[i]

  if (normal_p < 0.05 && tumor_p < 0.05) {
    if (sign(normal_cor) == sign(tumor_cor)) {
      if (abs(tumor_cor) > abs(normal_cor)) {
        result$PairType[i] <- "T_strong"
      } else {
        result$PairType[i] <- "N_strong"
      }
    } else  {
      result$PairType[i] <- "Discordant"
    } 
  } else if (normal_p < 0.05 && tumor_p > 0.05) {
    result$PairType[i] <- "N_specific"
  } else if (tumor_p < 0.05 && normal_p > 0.05) {
    result$PairType[i] <- "T_specific"
  } else  {
    result$PairType[i] <- "Discordant"
  }
}

# result$DEType <- ifelse(result$adj.P.Val < 0.05 & result$logFC > 0, 'Up', 
#                                                       ifelse(result$adj.P.Val < 0.05 & result$logFC < 0, 'Down', 'Not'))

result$DEType <- ifelse(result$adj.P.Val < 0.05,  'Sign', 'Not')

library(ggplot2)
library(dplyr)
library(tidyr)

df_summary <- result %>%
  group_by(PairType, DEType) %>%
  summarise(Count = n()) %>%
  ungroup()

ggplot(df_summary, aes(x = PairType, y = Count, fill = DEType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Not" = "gray", "Sign" = "red")) +
  labs(x = "Pair Type", y = "Count", fill = "DEType") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  ggtitle("Stacked Bar Chart of Pair Types and DE Types")


write.table(result,'GDCdata/TCGA-BRCA/sense_antisense_result.csv', quote = F,sep = ",", row.names = F,col.names = T)
write.table(sample_info,'GDCdata/TCGA-BRCA/sample_info.csv',quote = F,sep = ",", row.names = F,col.names = T)
colnames(result_matrix_filter)  <- sample_info$barcode
write.table(result_matrix_filter,'GDCdata/TCGA-BRCA/sense_antisense_ratio.csv',quote = F,sep = ",",row.names = T,col.names = T)
