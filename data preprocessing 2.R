# 1. 设置工作目录
setwd("G:\\187geneMR\\03.download\\GSE96804")

# 2. 检查并加载HTA-2.0平台的注释包
if (!require("hta20transcriptcluster.db")) {
  BiocManager::install("hta20transcriptcluster.db")  # HTA-2.0专用注释包
}
library(hta20transcriptcluster.db)

# 3. 读取表达矩阵文件（假设已下载的矩阵是TXT格式）
exp_file <- "GSE96804_series_matrix.txt"  # 替换为实际文件名
exp_data <- read.table(exp_file, header=TRUE, sep="\t", row.names=1)

# 检查数据格式（示例）
head(exp_data[, 1:4])  # 查看前4列

# 4. 探针ID转换为基因符号（使用hta20transcriptcluster.db）
probe_ids <- rownames(exp_data)  # 提取探针ID（如"2824546_st"）

# 转换为基因符号（处理多对一关系）
gene_symbols <- mapIds(
  hta20transcriptcluster.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"  # 多探针对应一基因时取第一个
)

# 5. 处理缺失值和重复基因名
# 去除无基因符号的探针
keep <- !is.na(gene_symbols)
exp_data <- exp_data[keep, ]
gene_symbols <- gene_symbols[keep]

# 对重复基因取表达均值
gene_exp_matrix <- aggregate(
  exp_data,
  by = list(Gene = gene_symbols),
  FUN = mean
)
rownames(gene_exp_matrix) <- gene_exp_matrix$Gene
gene_exp_matrix$Gene <- NULL

# 6. 查看结果
dim(gene_exp_matrix)  # 基因数×样本数
head(gene_exp_matrix[, 1:4])

# 7. 保存结果
#write.csv(gene_exp_matrix, file = "GSE96804_gene_expression_matrix.csv")

# 导出为TXT（兼容后续分析）
geneMatrix <- cbind(Gene = rownames(gene_exp_matrix), gene_exp_matrix)
write.table(
  geneMatrix,
  file = "GSE96804_Matrix.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
