setwd("G:\\CKDurine\\04.urine_ekidney")    
###################################04.04.数据前期处理和矫正###################################
#GITR (TNFRSF18) IL2RA仅存在/ILRα/Tac/IL2R/p55/ (CD25); IL7R (CD127); SELL (CD62L);
#"ITGAE"="CD103" ;TIM-3(CD366/HAVCR2); 
#读取文件，并对重复基因取均值
gc()
####################从这一部开始分析#####################
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(harmony)
library(cowplot)
library(patchwork)
library(devtools)
library(tidyverse)
#install.packages('devtools')
#devtools::install_github('immunogenomics/harmony')
#devtools::install_github('junjunlab/scRNAtoolVis')
#devtools::install_github('sajuukLyu/ggunchull')
#devtools::install_github('lydiaMyr/ImmuCellAI')
#devtools::install_github('mojaveazure/seurat-disk')
library(scRNAtoolVis) #为了clusterCornerAxes函数
library(ggunchull)
#本地安装包,小红书教程
#install.packages("F:\\R-4.5.1\\library\\scRNAtoolVis-master.zip", repos = NULL, type = "source")
#install.packages("F:\\R-4.5.1\\library\\scRNAtoolVis-master", repos = NULL, type = "source") 

#c("PSCA", "PLAT", "KRT13", "KRT4", "KRT17", "TSPAN1", "FXYD3")


library(Seurat)
library(ggplot2)
library(patchwork)

Early <- subset(merged_all, subset = disease_stage == "Early_DKD")


# 1. 数据合并与预处理 ---------------------------------------------------
urine_merged <- JoinLayers(merged_all)

rm(merged_all)
# 添加样本信息（使用Advanced替代Late）
# 1. 提取纯净样本名（去除序号）

# 3. 验证分组结果
table(urine_merged$disease_stage)


# 2. 数据整合 --------------------------------------------------------
# 基因筛选
# 1. 确保使用完全相同的基因集
shared_genes <- intersect(rownames(urine_merged), rownames(Early))
bladder_markers <- c("PSCA", "PLAT", "KRT13", "FXYD3", "KRT4", "KRT19")
use_genes <- union(shared_genes, bladder_markers[bladder_markers %in% shared_genes])

# 2. 统一预处理流程（关键步骤）
preprocess <- function(obj) {
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData(features = use_genes) %>%  # 限定基因集
    RunPCA(features = use_genes, npcs = 30, verbose = FALSE)
  return(obj)
}

urine_processed <- preprocess(subset(urine_merged, features = use_genes))
kidney_processed <- preprocess(subset(Early, features = use_genes))

# 3. 检查PCA维度一致性
stopifnot(ncol(urine_processed[["pca"]]) == ncol(kidney_processed[["pca"]]))
gc()


# 重命名细胞群
Idents(urine_processed) <- paste0("U_", Idents(urine_processed))
Idents(kidney_processed) <- paste0("K_", Idents(kidney_processed))#这样是否可行
# 4. 重新运行整合
anchors <- FindIntegrationAnchors(
  object.list = list(urine = urine_processed, kidney = kidney_processed),
  anchor.features = 2000,
  reduction = "rpca",
  dims = 1:30,
  k.anchor = 20  # 增加锚点数量提高稳定性
)

gc()
rm(Early)
rm(urine_merged)
# 整合数据
combined <- IntegrateData(anchorset = anchors)

# 标准分析流程
combined <- ScaleData(combined) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)%>% 
  RunTSNE(dims = 1:30)

table(combined@active.ident)
table(urine_processed@active.ident)
table(kidney_processed@active.ident)
# ===================== 绘图部分（保持您的原始风格）=====================


table(combined@active.ident)
# 绘图1：UMAP展示样本来源
pdf("05.pcaHeatmap.pdf", width=6, height=5)
DimPlot(combined, group.by = "orig.ident", label = TRUE,repel = TRUE)
dev.off()

# 绘图2：UMAP基础可视化
pdf("06.UMAP.pdf", width=6, height=5)
DimPlot(combined, reduction = "umap", label = TRUE,repel = TRUE)
dev.off()


table(combined@active.ident)

# 绘图3：带分组的UMAP
# 复制 active.ident 到 simple_label
combined$simple_label <- combined@active.ident


table(combined$simple_label)
pdf("06.umap_simple_label.pdf", width=6, height=5)
DimPlot(combined,
        reduction = "umap",
        group.by = "simple_label",
        cols = c(rep("#FAC858", 9), rep("#91CC75", 13)), # 您的原始配色
        pt.size = 0.5,
        label = TRUE,
        label.size = 5,
        repel = TRUE) +
  ggtitle("")  # 仅去除标题，保留图例和其他设置
dev.off()


# 4. 后续分析 --------------------------------------------------------
# 假设肾脏数据的注释列是 "celltype"
anchors <- FindTransferAnchors(
  reference = kidney_processed,
  query = urine_processed,
  reference.reduction = "pca",
  dims = 1:30
)

# 转移注释
predictions <- TransferData(
  anchorset = anchors,
  refdata = kidney_processed$celltype,  # 替换为您的肾脏注释列名
  dims = 1:30
)

# 添加到尿液数据
urine_sub <- AddMetaData(urine_processed, metadata = predictions)

# 检查预测结果
table(urine_sub$predicted.id)



# 8. 专业级可视化（保持您的风格）------------------------------------------
# UMAP标注
pdf("Final_UMAP.pdf", width = 8, height = 6)
DimPlot(urine_sub, 
        group.by = "predicted.id",
        cols = c("#FAC858", "#91CC75", "#EE6666","#FAC858", "#91CC75", "#EE6666","#FAC858", "#91CC75", "#EE6666"),  # 您原有的配色
        pt.size = 0.5,
        label = TRUE,
        label.size = 5,
        repel = TRUE) + NoLegend()
dev.off()

# 标记基因表达
pdf("Bladder_Markers.pdf", width = 10, height = 4)
FeaturePlot(urine_sub, 
            features = bladder_markers[bladder_markers %in% rownames(urine_sub)],
            pt.size = 0.2,
            ncol = 3,
            order = TRUE)
dev.off()


#4. 识别尿液中的膀胱/尿道上皮细胞######

# 方法1：低置信度预测可能是非肾脏细胞
low_conf_cells <- urine_sub$prediction.score.max < 0.5
table(low_conf_cells)  # 查看数量

# 方法2：膀胱标记基因高表达
urine_sub <- AddModuleScore(
  urine_sub,
  features = list(bladder_markers),
  name = "Bladder_Score"
)

# 标记膀胱细胞
bladder_cells <- WhichCells(urine_sub, expression = Bladder_Score1 > 0.5)
urine_sub$final_annotation <- ifelse(
  colnames(urine_sub) %in% bladder_cells,
  "Bladder/Urothelial",
  urine_sub$predicted.id
)

# 可视化
pdf("07.umap_final_annotation.pdf", width=10, height=6)
DimPlot(urine_sub, 
        reduction= "umap",
        group.by = "final_annotation", 
        label = TRUE,
        repel = TRUE) +
  FeaturePlot(urine_sub, features = bladder_markers)
dev.off()

# 9. 结果保存 ----------------------------------------------------------
saveRDS(combined, file = "Integrated_Urine_aKidney.rds")



cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542") #,,
