

# 1. 加载RDS文件（不要赋值，直接查看结构）
#str(readRDS("GSM3823939_control.s1.dgecounts.rds"), max.level = 2)
# 2. 或者查看对象名称
#names(readRDS("GSM3823939_control.s1.dgecounts.rds"))
# 3. 如果是Seurat或SingleCellExperiment对象
#class(readRDS("GSM3823939_control.s1.dgecounts.rds"))
# 加载测试安装是否成功
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)
library(Matrix)
rm(list=ls())
options(stringsAsFactors = F)

# 创建目录
setwd("G:\\187geneMR\\25.GSE131882\\0.download&process\\GSE131882_RAW")
## 导入数据

data <- readRDS("G:/187geneMR/25.GSE131882/0.download&process/GSE131882_RAW/GSM3823939_control.s1.dgecounts.rds")

#1. 正确提取表达矩阵
# 检查数据结构
names(data$umicount)  # 应该显示"exon", "inex", "intron"等
names(data$umicount$exon)  # 检查是否有"all"或"downsampled"

# 推荐先使用外显子矩阵（标准分析）
counts <- data$umicount$exon$all
if(is.null(counts)) counts <- data$umicount$exon$downsampled

#基因ID转换（优化版）
library(org.Hs.eg.db)

# 处理ENSEMBL ID（去除版本号）
rownames(counts) <- sub("\\..*", "", rownames(counts))

# 转换ID
gene_info <- select(
  org.Hs.eg.db,
  keys = rownames(counts),
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "ENSEMBL"
)
# 过滤蛋白编码基因
protein_coding <- gene_info[
  gene_info$GENETYPE == "protein-coding" &  # 条件1：基因类型为蛋白编码
    !is.na(gene_info$SYMBOL) &             # 条件2：存在有效基因符号（非NA）
    !duplicated(gene_info$ENSEMBL),        # 条件3：去除ENSEMBL ID重复项
]

# 重命名矩阵
counts <- counts[protein_coding$ENSEMBL, ]
rownames(counts) <- protein_coding$SYMBOL
# 1. 检查原始数据中的重复基因
dup_genes <- protein_coding$SYMBOL[duplicated(protein_coding$SYMBOL)]
dup_genes  # 查看前几个重复基因


#创建Seurat对象
library(Seurat)

#直接使用make.unique()处理重复基因名（最快且内存友好）
rownames(counts) <- make.unique(rownames(counts))  # 例如：MALAT1 → MALAT1.1, MALAT1.2

# 创建Seurat对象
controls1 <- CreateSeuratObject(
  counts = counts,
  project = "Control01",
  min.cells = 3,      # 宽松的基因过滤
  min.features = 200  # 核通常检出基因较少
)

# 3. 检查关键标记基因是否被拆分（如CD3D、VEGFA）
markers <- c("CD3D", "VEGFA", "PECAM1","COL18A1","VIM")
lapply(markers, function(x) grep(paste0("^", x), rownames(Control1_qc), value = TRUE))


# 线粒体基因（注意核可能含更多mtRNA）
controls1[["percent.mt"]] <- PercentageFeatureSet(controls1, pattern = "^MT-")
controls1 <- subset(controls1, 
               subset = nFeature_RNA > 300 & 
                 nFeature_RNA < 5000 & 
                 percent.mt < 15)
# 重命名细胞
colnames(controls1) <- paste0("Control01", "_", seq_len(ncol(controls1)))
# 关键步骤：同步基因名到counts层

gene_names <- rownames(controls1)
if("counts" %in% names(controls1@assays[["RNA"]]@layers)) {
  rownames(controls1@assays[["RNA"]]@layers[["counts"]]) <- gene_names
}
# 验证同步结果
stopifnot(identical(
  rownames(controls1),
  rownames(LayerData(controls1, layer = "counts"))
))
# 保存质控后的独立样本
saveRDS(controls1, file = "controls1_qc.rds")


########处理6各样本##########
########处理6各样本##########
library(Seurat)
library(org.Hs.eg.db)

# 设置工作目录和文件路径
setwd("G:\\187geneMR\\25.GSE131882\\0.download&process\\GSE131882_RAW")

# 定义样本信息
sample_files <- c(
  "GSM3823939_control.s1.dgecounts.rds",
  "GSM3823940_control.s2.dgecounts.rds",
  "GSM3823941_control.s3.dgecounts.rds",
  "GSM3823942_diabetes.s1.dgecounts.rds",
  "GSM3823943_diabetes.s2.dgecounts.rds",
  "GSM3823944_diabetes.s3.dgecounts.rds"
)

sample_names <- c(
  "Control1", "Control2", "Control3",
  "Early_DKD1", "Early_DKD2", "Early_DKD3"
)

### 第一阶段：独立处理每个样本 ###
processed_samples <- list()

for(i in seq_along(sample_files)){
  cat("\n>>> 处理样本", i, ":", sample_names[i], "\n")
  
  # 1. 读取RDS文件
  data <- readRDS(sample_files[i])
  
  # 2. 提取外显子矩阵
  counts <- data$umicount$exon$all %||% data$umicount$exon$downsampled
  stopifnot(!is.null(counts))
  
  # 3. 基因ID转换
  rownames(counts) <- sub("\\..*", "", rownames(counts))
  gene_info <- select(
    org.Hs.eg.db,
    keys = rownames(counts),
    columns = c("SYMBOL", "GENETYPE"),
    keytype = "ENSEMBL"
  )
  
  # 过滤蛋白编码基因
  protein_coding <- gene_info[
    gene_info$GENETYPE == "protein-coding" & 
      !is.na(gene_info$SYMBOL),
  ]
  
  # 将基因信息与表达矩阵合并
  counts <- counts[protein_coding$ENSEMBL, ]
  protein_coding$total_counts <- rowSums(counts)
  
  # 对于重复基因名，保留总表达量最高的那个
  protein_coding <- protein_coding[order(protein_coding$SYMBOL, -protein_coding$total_counts), ]
  protein_coding <- protein_coding[!duplicated(protein_coding$SYMBOL), ]
  
  # 筛选最终基因
  counts <- counts[protein_coding$ENSEMBL, ]
  rownames(counts) <- protein_coding$SYMBOL
  
  # 4. 创建Seurat对象
  sobj <- CreateSeuratObject(
    counts = counts,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  
  # 5. 质控过滤
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj <- subset(sobj, 
                 subset = nFeature_RNA > 300 & 
                   nFeature_RNA < 5000 & 
                   percent.mt < 15)
  
  # 6. 重命名细胞
  colnames(sobj) <- paste0(sample_names[i], "_", seq_len(ncol(sobj)))
  
  # 7. 同步基因名到counts层
  gene_names <- rownames(sobj)
  if("counts" %in% names(sobj@assays[["RNA"]]@layers)) {
    rownames(sobj@assays[["RNA"]]@layers[["counts"]]) <- gene_names
  }
  
  # 验证同步结果
  stopifnot(identical(
    rownames(sobj),
    rownames(LayerData(sobj, layer = "counts"))
  ))
  
  # 8. 保存并存储
  saveRDS(sobj, file = paste0(sample_names[i], "_qc.rds"))
  processed_samples[[i]] <- sobj
  cat(">>> 已保存", paste0(sample_names[i], "_qc.rds"), "\n")
}

### 第二阶段：获取共同基因 ###
cat("\n>>> 获取质控后的共同基因...\n")
gene_list <- lapply(processed_samples, rownames)
common_genes <- Reduce(intersect, gene_list)
cat(">>> 共同基因数量:", length(common_genes), "\n")

### 第三阶段：裁剪并合并 ###
final_list <- list()
group_names <- c(rep("Control", 3), rep("EarlyDKD", 3))

for(i in seq_along(processed_samples)){
  cat("\n>>> 裁剪样本", i, "到共同基因...\n")
  sobj <- processed_samples[[i]]
  
  # 裁剪到共同基因
  sobj <- subset(sobj, features = common_genes)
  
  # 添加分组信息
  sobj$group <- group_names[i]
  
  # 再次确保基因名同步
  if("counts" %in% names(sobj@assays[["RNA"]]@layers)) {
    rownames(sobj@assays[["RNA"]]@layers[["counts"]]) <- rownames(sobj)
  }
  
  # 保存裁剪后的版本
  saveRDS(sobj, file = paste0(sample_names[i], "_final.rds"))
  final_list[[i]] <- sobj
}

# 合并对象
cat("\n>>> 合并所有样本...\n")
merged_seurat <- merge(
  final_list[[1]], 
  y = final_list[2:6],
  add.cell.ids = sample_names
)

# 同步所有层的基因名
for(layer_name in Layers(merged_seurat)){
  rownames(merged_seurat@assays$RNA@layers[[layer_name]]) <- rownames(merged_seurat)
}

# 验证合并结果
cat("\n>>> 合并后对象信息：\n")
print(merged_seurat)
cat(">>> 细胞数:", ncol(merged_seurat), "\n")
cat(">>> 基因数:", nrow(merged_seurat), "\n")

# 保存最终结果
saveRDS(merged_seurat, "merged_seurat_final.rds")
cat("\n>>> 分析流程完成！结果已保存为 merged_seurat_final.rds\n")