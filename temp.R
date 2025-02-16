
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BRCA",             
  data.category = "Transcriptome Profiling",  
  data.type = "Gene Expression Quantification",  
)



GDCdownload(query, files.per.chunk = 50)
data <- GDCprepare(query)
sample_info_list=colData(data)@listData
sample_info=as.data.frame(sample_info_list[1:10])
sample_info=sample_info[sample_info$definition %in% c('Primary solid Tumor', 'Solid Tissue Normal'),]
sample_info=sample_info[order(sample_info$definition),]
table(sample_info$definition)

data_expr_tpm <- assay(data,i = "tpm_unstrand")
data_expr_tpm <- data_expr_tpm[,sample_info$barcode]
rownames(data_expr_tpm) <- substring(rownames(data_expr_tpm),1,15)
data_expr_tpm <- log2(data_expr_tpm + 1)







pair_info <- putative_bidirectional[,c(1:3)]
result_matrix <- matrix(NA, nrow = nrow(pair_info), ncol = ncol(data_expr_tpm))

for (i in 1:nrow(pair_info)) {
  sense_gene <- pair_info$Gene_ID_sense[[i]]
  antisense_gene <- pair_info$Gene_ID_antisense[[i]]
  
  if (sense_gene %in% rownames(data_expr_tpm) &&
      antisense_gene %in% rownames(data_expr_tpm)) {
    for (j in 1:ncol(data_expr_tpm)) {
      sense_gene_exp <- data_expr_tpm[sense_gene, j]
      antisense_gene_exp <- data_expr_tpm[antisense_gene, j]
      sum_ = sense_gene_exp + antisense_gene_exp
      if (sum_ == 0) {
        result_matrix[i, j] <- 0.5
      } else if (sense_gene_exp == 0) {
        result_matrix[i, j] <- 0
      } else if (antisense_gene_exp == 0) {
        result_matrix[i, j] <- 1
      } else {
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

















pair_info_filter <- pair_info[pair_info$PairID %in% rownames(result_matrix_filter), ]



correlation_change_matrix <- matrix(NA, nrow = nrow(pair_info_filter), ncol = 4)
normal_samples <- sample_info$barcode[sample_info$definition == "Solid Tissue Normal"]
tumor_samples <- sample_info$barcode[sample_info$definition == "Primary solid Tumor"]

# 对每个PairID计算表达相关性变化
for (i in 1:nrow(pair_info_filter)) {
  sense_gene <- pair_info_filter$Gene_ID_sense[[i]]
  antisense_gene <- pair_info_filter$Gene_ID_antisense[[i]]
  
  # 提取肿瘤和正常样本的表达数据
  protein_coding_exp_normal <- data_expr_tpm[sense_gene, normal_samples]
  lncRNA_exp_normal <- data_expr_tpm[antisense_gene, normal_samples]
  
  protein_coding_exp_tumor <- data_expr_tpm[sense_gene, tumor_samples]
  lncRNA_exp_tumor <- data_expr_tpm[antisense_gene, tumor_samples]
  
  # 计算正常样本的表达相关性
  cor_normal <- cor.test(protein_coding_exp_normal, lncRNA_exp_normal, method = "pearson")
  cor_normal_p <- as.numeric(cor_normal[3])
  cor_normal_correlation <- as.numeric(cor_normal[4])
  
  # 计算肿瘤样本的表达相关性
  cor_tumor <- cor.test(protein_coding_exp_tumor, lncRNA_exp_tumor, method = "pearson")
  cor_tumor_p <- as.numeric(cor_tumor[3])
  cor_tumor_correlation <- as.numeric(cor_tumor[4])
  
  correlation_change_matrix[i, 1] <- cor_normal_correlation
  correlation_change_matrix[i, 2] <- cor_normal_p
  correlation_change_matrix[i, 3] <- cor_tumor_correlation
  correlation_change_matrix[i, 4] <- cor_tumor_p

}

colnames(correlation_change_matrix) <- c('normal_cor', 'normal_p', 'tumor_cor', 'tumor_p')
result <- cbind(pair_info_filter, correlation_change_matrix, DEG[pair_info_filter$PairID,])
result <- na.omit(result)

result$PairType <- NA 
for (i in 1:nrow(result)) {
  normal_p <- result$normal_p[i]
  tumor_p <- result$tumor_p[i]
  normal_cor <- result$normal_cor[i]
  tumor_cor <- result$tumor_cor[i]

  if (normal_p < 0.05 && tumor_p < 0.05 && abs(normal_cor) >= 0.1  && abs(tumor_cor) >= 0.1) {
    if (sign(normal_cor) == sign(tumor_cor)) {
      # 如果相关性同号，比较绝对值
      if (abs(tumor_cor) > abs(normal_cor)) {
        result$PairType[i] <- "T_strong"
      } else {
        result$PairType[i] <- "N_strong"
      }
    } else  {
      result$PairType[i] <- "Discordant"
    } 
  } else if (normal_p < 0.05 && abs(normal_cor) >= 0.1 && tumor_p > 0.05) {
    result$PairType[i] <- "N_specific"
  } else if (tumor_p < 0.05 && abs(tumor_cor) >= 0.1 && normal_p > 0.05) {
    result$PairType[i] <- "T_specific"
  } else  {
    result$PairType[i] <- "Discordant"
  }
}

# 依据 adj.P.Val 是否小于 0.05 来分类
# result$DEType <- ifelse(result$adj.P.Val < 0.05 & result$logFC > 0, 'Up', 
#                                                       ifelse(result$adj.P.Val < 0.05 & result$logFC < 0, 'Down', 'Not'))

result$DEType <- ifelse(result$adj.P.Val < 0.05,  'Sign', 'Not')

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 整理数据
df_summary <- result %>%
  group_by(PairType, DEType) %>%
  summarise(Count = n()) %>%
  ungroup()

# 2. 绘制堆叠柱状图
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


write.table(data_expr_tpm,'GDCdata/TCGA-BRCA/expr_tpm.csv', quote = F,sep = ",", row.names = T,col.names = T)
write.table(result,'GDCdata/TCGA-BRCA/sense_antisense_result.csv', quote = F,sep = ",", row.names = F,col.names = T)
write.table(sample_info,'GDCdata/TCGA-BRCA/sample_info.csv',quote = F,sep = ",", row.names = F,col.names = T)
colnames(result_matrix_filter)  <- sample_info$barcode
write.table(result_matrix_filter,'GDCdata/TCGA-BRCA/sense_antisense_ratio.csv',quote = F,sep = ",",row.names = T,col.names = T)


library(dplyr)

# 设置 GDCdata 目录路径
root_dir <- "GDCdata"

# 获取所有 TCGA- 开头的文件夹
folders <- list.files(root_dir, pattern = "^TCGA-", full.names = TRUE)

# 初始化空数据框
merged_df <- data.frame()

# 遍历文件夹
for (folder in folders) {
  csv_file <- file.path(folder, "sense_antisense_result.csv")
  
  # 确保 CSV 文件存在
  if (file.exists(csv_file)) {
    # 读取 CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE,  header = T)
    
    # 提取肿瘤类别（文件夹名）
    tumor_type <- basename(folder)
    
    # 添加肿瘤类别列
    df$Tumor_Type <- tumor_type
    
    # 合并数据（按行）
    merged_df <- bind_rows(merged_df, df)  # 也可用 rbind()
  }
}

# 保存合并后的 CSV
write.table(merged_df, "GDCdata/merged_sense_antisense_result.csv", quote = F,sep = ",", row.names = F,col.names = T)






# 仅保留 PairType 为 T_specific 或 T_strong 的数据
merged_df <- read.csv("GDCdata/merged_sense_antisense_result.csv",header = T,sep = ",")
filtered_df <- merged_df %>%
  filter(PairType %in% c("T_specific", "T_strong"))

# 计算各个 Tumor_Type 下 DEType == Not 和 非Not 的 hallmark_genes 占比
library(msigdbr)
hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")

tumor_stats <- filtered_df %>%
  group_by(Tumor_Type, DEType) %>%
  summarise(
    hallmark_count = length(intersect(hallmark_genes$ensembl_gene,
                                      c(Gene_ID_sense, Gene_ID_antisense))),
    total_genes = n()
  ) %>%
  mutate(percentage = hallmark_count / total_genes * 100) %>%
  ungroup()


ggplot(tumor_stats, aes(x = Tumor_Type, y = percentage, fill = DEType)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL,
       x = "Tumor Type", y = "Hallmark Genes Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Not" = "blue", "Sign" = "red"))  # 颜色可调整




library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)

# 计算每个 Tumor_Type 的 sense_gene + antisense_gene 集合
tumor_gene_sets <- filtered_df %>%
  group_by(Tumor_Type) %>%
  summarise(gene_set = list(unique(c(Gene_ID_sense, Gene_ID_antisense)))) %>%
  ungroup()

# 计算 Jaccard 距离矩阵
tumor_types <- tumor_gene_sets$Tumor_Type
num_tumors <- length(tumor_types)
jaccard_matrix <- matrix(0, nrow = num_tumors, ncol = num_tumors, dimnames = list(tumor_types, tumor_types))

for (i in 1:num_tumors) {
  for (j in 1:num_tumors) {
    if (i != j) {
      genes_i <- tumor_gene_sets$gene_set[[i]]
      genes_j <- tumor_gene_sets$gene_set[[j]]
      
      intersection_size <- length(intersect(genes_i, genes_j))
      union_size <- length(unique(c(genes_i, genes_j)))
      
      jaccard_matrix[i, j] <- 1 - (intersection_size / union_size)  # 计算 Jaccard 距离
    }
  }
}

# 绘制热图
pheatmap(jaccard_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Jaccard Distance Heatmap of Tumor Types",
         display_numbers = TRUE)




putative_bidirectional <- read.table('GDCdata/putative_bidirectional.txt',header = T,sep = "\t",check.names = F)
merged_df <- read.csv("GDCdata/merged_sense_antisense_result.csv",header = T,sep = ",")
filtered_df <- merged_df %>%
  filter(PairType %in% c("T_specific", "T_strong"))

filtered_df <- merge(filtered_df, putative_bidirectional[,c('PairID','putative_TF')], by.x.y = 'PairID')
# filtered_df <- filtered_df[filtered_df$putative_TF != '' & filtered_df$DEType == 'Sign',]
filtered_df$putative_TF <- toupper(filtered_df$putative_TF)
# filtered_df <- filtered_df[filtered_df$Tumor_Type == 'TCGA-BRCA',]

# 初始化一个空列表来保存数据框
ratio_dict <- list()
sample_dict <- list()
expr_dict <- list()
clinic_dict <- list()

# 设置 GDCdata 目录路径
root_dir <- "GDCdata"
# 获取所有 TCGA- 开头的文件夹
folders <- list.files(root_dir, pattern = "^TCGA-", full.names = TRUE)

# 遍历文件夹

for (folder in folders) {
  csv_file <- file.path(folder, "sense_antisense_ratio.csv")
  
  # 确保 CSV 文件存在
  if (file.exists(csv_file)) {
    # 读取 CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = F, header = T)
    
    # 获取文件夹名作为 key (如 TCGA-COAD)
    folder_name <- basename(folder)
    
    # 将数据框保存到字典中
    ratio_dict[[folder_name]] <- df
  }
}

for (folder in folders) {
  csv_file <- file.path(folder, "sample_info.csv")
  
  # 确保 CSV 文件存在
  if (file.exists(csv_file)) {
    # 读取 CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = F, header = T)
    
    # 获取文件夹名作为 key (如 TCGA-COAD)
    folder_name <- basename(folder)
    
    # 将数据框保存到字典中
    sample_dict[[folder_name]] <- df
  }
}

for (folder in folders) {
  # 获取文件夹名作为 key (如 TCGA-COAD)
  folder_name <- basename(folder)
  csv_file <- file.path(folder, paste(folder_name, ".survival.tsv",sep = ""))
  
  # 确保 CSV 文件存在
  if (file.exists(csv_file)) {
    # 读取 CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE, sep = '\t', check.names = F, header = T)

    
    # 将数据框保存到字典中
    clinic_dict[[folder_name]] <- df
  }
}


for (folder in folders) {
  csv_file <- file.path(folder, "expr_tpm.csv")
  
  # 确保 CSV 文件存在
  if (file.exists(csv_file)) {
    # 读取 CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = F, header = T, row.names = NULL)
    colnames(df)[1] <- 'Gene_id'
    df <- df[!duplicated(df$Gene_id), ] 
    rownames(df) <- df$Gene_id
    df <- df[,-1]
    
    # 获取文件夹名作为 key (如 TCGA-COAD)
    folder_name <- basename(folder)
    
    # 将数据框保存到字典中
    expr_dict[[folder_name]] <- df
  }
}



# 获取Primary Tumor条形码的函数
get_primary_tumor_barcodes <- function(tumor_type, sample_dict) {
  sample_data <- sample_dict[[tumor_type]]
  primary_tumor_data <- sample_data %>% filter(definition == "Primary solid Tumor")
  primary_tumor_data <- primary_tumor_data[!duplicated(primary_tumor_data$sample_submitter_id),]
  primary_tumor_barcodes <- primary_tumor_data$barcode
  primary_tumor_short_barcodes <- primary_tumor_data$sample_submitter_id
  
  # 合并两个向量并返回
  return(list(primary_tumor_barcodes, primary_tumor_short_barcodes))
}

# 获取TF表达量的函数
get_tf_expression <- function(tf, expr_dict, gene_info, tumor_type, barcodes) {
  tf_ensembl <- unique(gene_info %>% filter(hgnc_symbol == tf) %>% pull(ensembl_gene_id))
  if (length(tf_ensembl) == 0) {
    return(NULL)  # 如果没有找到对应的TF，返回NULL
  }
  
  if (length(tf_ensembl) > 1) {
    tf_ensembl <- tf_ensembl[1]
  }
  
  # 获取对应表达值并只保留Primary Tumor的条形码
  expr_values <- expr_dict[[tumor_type]][tf_ensembl, barcodes]
  
  if (is.list(expr_values)) {  # 如果是list类型，需要转换为数值型向量
    expr_values <- unlist(expr_values)
  }
  return(expr_values)
}

# 获取ratio值的函数
get_ratio_values <- function(pair_id, ratio_dict, tumor_type, barcodes) {
  ratio_values <- ratio_dict[[tumor_type]][pair_id, barcodes]
  if (is.list(ratio_values)) {  # 如果是list类型，转换为数值型向量
    ratio_values <- unlist(ratio_values)
  }
  return(ratio_values)
}

# 计算相关性的函数
calculate_correlation <- function(expr_values, ratio_values) {
  expr_values <- as.numeric(expr_values)
  ratio_values <- as.numeric(ratio_values)
  
  # 删除任何NA值
  valid_indices <- !is.na(expr_values) & !is.na(ratio_values)
  expr_values <- expr_values[valid_indices]
  ratio_values <- ratio_values[valid_indices]
  
  # 如果长度不一致，或者都没有有效值，返回NA
  if (length(expr_values) == 0 || length(ratio_values) == 0 || length(expr_values) != length(ratio_values)) {
    return(c(NA, NA))
  }
  
  # 计算Pearson相关性和p值
  cor_result <- cor.test(expr_values, ratio_values, method = "pearson")
  return(c(round(cor_result$estimate, 4), round(cor_result$p.value, 4)))
}

calculate_survival_pvalue_from_data <- function(ratio_values, clinic) {
  # 可以使用生存分析的统计方法，如 Cox 比例风险模型或 Log-rank 检验
  # 这里用 Cox 模型作为示例
  library(survival)
  surv_data <- data.frame(OS = clinic$OS, time = clinic$OS.time, ratio = ratio_values)
  cox_model <- coxph(Surv(time, OS) ~ ratio, data = surv_data)
  
  # 返回 p 值
  return(summary(cox_model)$sctest[3])  # 获取 cox 模型的 p 值
}

# 计算putative_TF_cor和putative_TF_pvalue的函数
calculate_putative_TF <- function(filtered_df, expr_dict, ratio_dict, gene_info, sample_dict) {
  filtered_df$putative_TF_cor <- ""
  filtered_df$putative_TF_pvalue <- ""
  
  for (i in 1:nrow(filtered_df)) {
    if (filtered_df$putative_TF[i] == "") {
      next 
    }
    else {
    putative_TFs <- strsplit(filtered_df$putative_TF[i], ",")[[1]]
    pair_id <- filtered_df$PairID[i]
    tumor_type <- filtered_df$Tumor_Type[i]
    
    primary_barcodes <- get_primary_tumor_barcodes(tumor_type, sample_dict)[1]
    
    # 获取ratio值
    ratio_values <- get_ratio_values(pair_id, ratio_dict, tumor_type, primary_barcodes)
    
    # 获取每个TF的表达量
    expr_values_list <- lapply(putative_TFs, function(tf) {
      get_tf_expression(tf, expr_dict, gene_info, tumor_type, primary_barcodes)
    })
    
    # 计算相关性和p值
    results <- mapply(function(expr_values) {
      calculate_correlation(expr_values, ratio_values)
    }, expr_values_list)
    
    # 每个TF的相关性和p值
    cor_vals <- results[1, , drop = FALSE]
    p_vals <- results[2, , drop = FALSE]
    
    # 将每个TF的cor和p值按逗号连接起来
    filtered_df$putative_TF_cor[i] <- paste(cor_vals, collapse = ",")
    filtered_df$putative_TF_pvalue[i] <- paste(p_vals, collapse = ",")
    }
  }
  
  return(filtered_df)
}

calculate_survival_pvalue <- function(filtered_df, ratio_dict, clinic_dict, sample_dict) {
  # 添加新列用于存储生存p值和背景p值
  filtered_df$survival_pvalue <- NA
  filtered_df$background_pvalue <- NA
  
  for (i in 1:nrow(filtered_df)) {
    pair_id <- filtered_df$PairID[i]
    tumor_type <- filtered_df$Tumor_Type[i]
    
    # 获取对应的clinic数据和barcodes
    clinic_data <- clinic_dict[[tumor_type]]
    barcodes_list <- get_primary_tumor_barcodes(tumor_type, sample_dict)
    barcodes <- data.frame('long' = unlist(barcodes_list[1]), 'short' = unlist(barcodes_list[2]))
    
    # 获取匹配的clinic数据
    clinic <- clinic_data %>%
      filter(sample %in% barcodes$short)
    
    # order barcodes 
    rownames(barcodes) <- barcodes$short
    barcodes <- barcodes[clinic$sample, ]
    
    # 获取符合条件的ratio值
    ratio_values <- get_ratio_values(pair_id, ratio_dict, tumor_type, barcodes$long)
    
    # 计算生存p值
    survival_pvalue <- calculate_survival_pvalue_from_data(ratio_values, clinic)
    
    # 计算背景p值：随机生成的ratio_values
    random_ratio_values <- runif(length(ratio_values), 0, 1)
    background_pvalue <- calculate_survival_pvalue_from_data(random_ratio_values, clinic)
    
    # 将结果存入新列
    filtered_df$survival_pvalue[i] <- survival_pvalue
    filtered_df$background_pvalue[i] <- background_pvalue
  }
  
  return(filtered_df)
}


# 运行计算
filtered_df <- calculate_putative_TF(filtered_df, expr_dict, ratio_dict, gene_info, sample_dict)
filtered_df <- calculate_survival_pvalue(filtered_df, ratio_dict, clinic_dict, sample_dict)

# 使用 ggplot 绘制 DEType 分类的 survival_pvalue 密度图
ggplot(filtered_df, aes(x = survival_pvalue, color = DEType, fill = DEType)) +
  geom_density(alpha = 0.4) +  # alpha 调节透明度
  labs(title = "Density Plot of Survival p-values by DEType",
       x = "Survival p-value", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top")  # 可根据需要调整图例位置
