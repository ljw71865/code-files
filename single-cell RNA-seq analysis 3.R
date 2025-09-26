
setwd("G:\\CKDurine\\03.early_DKD\\6.monocle")                 #设置工作目录

#导入注释好的seurat对象（已注释）\\advanced


responder <- subset(Early, idents = c("PTC","Inj-PTC"))
table(responder@active.ident)  # 检查细胞数量
#"Inj-PTC","Prolif-PTC"
library(monocle)
library(dplyr)
library(Matrix)
library(Seurat)
library(ggplot2)  

library(libcoin)
library(partykit)

gc()
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
#expr_matrix <- as(as.matrix(responder@assays[["RNA"]]@counts), 'sparseMatrix')
expr_matrix <- as.sparse(responder@assays[["RNA"]]@layers[["counts"]])
##提取表型信息到p_data(phenotype_data)里面 
p_data <- responder@meta.data 
p_data$celltype <- responder@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(responder),row.names = row.names(responder))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
# 在运行任何Monocle函数前执行

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit=.5,
                      expressionFamily = negbinomial.size())
#save(cds, file = "before_estimateSizeFactors.rda")

#伪时间分析流程分支1和分支2
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#,cores=4,relative_expr=T
cds=detectGenes(cds,min_expr = 0.1) #计算每个基因在多少细胞中表达
print(head(fData(cds)))  

#save(cds, file = "select_gene_for_monocle.rda")

#load("select_gene_for_monocle.rda")
#下一步选择要构建轨迹所用的基因,对轨迹影响很大
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 3)) 
#过滤掉在小于10个细胞中表达的基因
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(responder) 

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype") 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)
##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
##差异基因的结果文件保存
write.table(deg,file="gene_for_trajectory.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
save(cds, file = "cds-for_plot.rda")
#load("G:/RCC/sc_immune_therapy/6.monocle/cds-for_plot.rda")


pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()

cds=reduceDimension(cds,reduction_method = "DDRTree",max_components = 2) 

#首次运行orderCells()生成State信息（不指定root_state）
cds <- orderCells(cds)  # 此时会自动分配State


# 手动指定root_state（假设PCT是State1）
#cds <- orderCells(cds, root_state = 3)  # 替换为实际的State编号

# 加载ggplot2包（如果尚未加载）
library(ggplot2)

# 1. 绘制伪时间轨迹图 -----------------------------------------------------------------
pdf("train.monocle.pseudotime.pdf", width = 3, height = 3)  # 创建PDF文件，设置宽高为7英寸
plot_cell_trajectory(
  cds, 
  color_by = "Pseudotime",  # 按伪时间值着色
  size = 1,                 # 点的大小
  show_backbone = TRUE      # 显示主干轨迹线
) + 
  theme(
    text = element_text(size = 12),        # 全局文本大小
    axis.title = element_text(size = 14),   # 坐标轴标题大小
    axis.text = element_text(size = 12),    # 坐标轴刻度标签大小
    legend.title = element_text(size = 14), # 图例标题大小
    legend.text = element_text(size = 12)   # 图例项目文本大小
  )
dev.off()  # 关闭图形设备


# 2. 绘制细胞类型轨迹图 --------------------------------------------------------------
pdf("train.monocle.celltype.pdf", width = 3, height = 3)
plot_cell_trajectory(
  cds,
  color_by = "celltype",  # 按细胞类型着色
  size = 1,
  show_backbone = TRUE
) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(nrow = 2))  # 关键修改：图例分两列显示
dev.off()


# 3. 绘制状态轨迹图 ------------------------------------------------------------------
pdf("train.monocle.state.pdf", width = 4, height = 4)
plot_cell_trajectory(
  cds, 
  color_by = "State",  # 按Monocle状态着色
  size = 1,
  show_backbone = TRUE
) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()


# 5. 分面绘制状态轨迹图（横向排列）--------------------------------------------------
pdf("train.monocle.state.faceted.pdf", width = 6, height = 3)
plot_cell_trajectory(
  cds, 
  color_by = "celltype"
) + 
  facet_wrap("~celltype", nrow = 1) +  # 按State分面，单行排列
  theme(
    text = element_text(size = 14),         # 增大基础字号
    axis.title = element_text(size = 16),    # 增大坐标轴标题
    axis.text = element_text(size = 14),     # 增大刻度标签
    legend.title = element_text(size = 16),  # 增大图例标题
    legend.text = element_text(size = 14),   # 增大图例文本
    strip.text = element_text(size = 14)     # 分面标签大小
  )
dev.off()

# 加载包
library(monocle)
library(ggplot2)
library(tidyr)


# 拟时序基因分析 ==============================================================


# 3. 关键基因拟时序表达曲线 ---------------------------------------------------

# 定义关键标记基因（需替换为实际关注的基因）
#key_markers <- c("VCAM1", "HAVCR1","CENPF",  "TOP2A")

key_markers <- c("PDK4", "RHCG", "FBP1")
# 筛选数据集中存在的基因
valid_markers <- key_markers[key_markers %in% fData(cds)$gene_short_name]

pdf("魔性基因_expression.pdf", width = 5, height = 3)
plot_genes_in_pseudotime(
  cds[valid_markers,], 
  color_by = "celltype",
  ncol = 2,                 # 分两列排列
  panel_order = valid_markers  # 按指定顺序排列
) + 
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"))  # 自定义颜色
dev.off()


# 目标基因"FN1", "ALDH2"拟时序分析####
target_genes <- c("FN1", "ALDH2")  # 确保基因名与数据匹配

# 检查基因是否存在
missing_genes <- setdiff(target_genes, rownames(cds))
if (length(missing_genes) > 0) {
  stop(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 1. 拟时序热图 (PDF输出)
pdf("FN1_ALDH2_pseudotime_heatmap.pdf", width = 8, height = 5)
plot_genes_in_pseudotime(cds[target_genes, ], 
                         color_by = "celltype",
                         ncol = 2) +
  ggtitle("FN1 and ALDH2 Expression along Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

library(tidyr)
# 2. 拟时序趋势曲线 (PDF输出)
pdf("FN1_ALDH2_trend_curves.pdf", width = 3, height = 3)
exprs <- exprs(cds)[target_genes, ]
pdata <- pData(cds)
df <- data.frame(
  Pseudotime = pdata$Pseudotime,
  State = pdata$State,
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df, aes(x = Pseudotime, y = Expression, color = Gene)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "Pseudotime", y = "Expression", 
       title = "Expression Trends") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#####3. 拟时序趋势曲线 (PDF输出)
pdf("FN1_ALDH2_trend_curves_by_celltype.pdf", width = 6, height = 4)
exprs <- exprs(cds)[target_genes, ]
pdata <- pData(cds)
df <- data.frame(
  Pseudotime = pdata$Pseudotime,
  CellType = pdata$celltype,  # 确保列名与您的metadata一致
  State = pdata$State,
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df, aes(x = Pseudotime, y = Expression, color = Gene)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "Pseudotime", y = "Expression", 
       title = "FN1 and ALDH2 Expression Trends by Cell Type") +
  facet_wrap(~ CellType, ncol = 3) +  # 按细胞类型分面
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
dev.off()




# 5. 小提琴图 (PDF输出)
pdf("FN1_ALDH2_violin_by_celltype.pdf", width = 8, height = 5)
df_violin <- data.frame(
  celltype = as.factor(pdata$celltype),
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df_violin, aes(x = celltype, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("FN1" = "#d9352a", "ALDH2" = "#0072B2")) +
  labs(x = "celltype", y = "Expression", 
       title = "FN1 and ALDH2 Expression by celltype") +
  theme_minimal(base_size = 12) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# 5. 小提琴图 (PDF输出)
# 5. 小提琴图 (PDF输出)
pdf("FN1_ALDH2_violin_by_state.pdf", width = 8, height = 5)
df_violin <- data.frame(
  State = as.factor(pdata$State),
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df_violin, aes(x = State, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "State", y = "Expression", 
       title = "FN1 and ALDH2 Expression by State") +
  theme_minimal(base_size = 12) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


## 目标基因ALDOB_LRP2_SLC34A1拟时序小提琴图####

target_genes <- c("ALDOB", "LRP2", "SLC34A1")
stopifnot(all(target_genes %in% rownames(exprs(cds))))

# 创建小提琴图数据
df_violin <- data.frame(
  celltype = as.factor(pData(cds)$celltype),
  ALDOB = exprs(cds)["ALDOB", ],
  LRP2 = exprs(cds)["LRP2", ],
  SLC34A1 = exprs(cds)["SLC34A1", ]
) %>% 
  pivot_longer(cols = all_of(target_genes),
               names_to = "Gene", 
               values_to = "Expression")

# 绘制小提琴图
pdf("ALDOB_LRP2_SLC34A1_violin_by_celltype.pdf", width = 8, height = 6)
ggplot(df_violin, aes(x = celltype, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("ALDOB" = "#4979b6", "LRP2" = "#d98d2a", "SLC34A1" = "#d9352a")) +
  labs(x = "Cell Type", y = "Expression") +
  theme_minimal(base_size = 14) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
dev.off()


######可视化基因集分数沿伪时间的变化,以celltype#######

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 提取关键数据
library(stringr)
plot_data <- pData(cds) %>%
  select(Pseudotime, celltype, 
         EPITHELIAL_MESENCHYMAL_TRANSITION_1:INFLAMMATORY_RESPONSE_1) %>%
  pivot_longer(cols = -c(Pseudotime, celltype),
               names_to = "Geneset",
               values_to = "Score")


# 移除Geneset名称中的"_1"

# 查看数据结构
head(plot_data)
plot_data <- plot_data %>%
  mutate(Geneset = stringr::str_remove(Geneset, "_1$"))


pdf("Geneset_Scores_along_Pseudotime_by_Celltype.pdf", width = 8, height = 6)

ggplot(plot_data, aes(x = Pseudotime, y = Score, color = celltype)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  facet_wrap(~ Geneset, scales = "free_y", ncol = 3) +
  labs(x = "Pseudotime", y = "Geneset Score", 
       title = "Geneset Scores along Pseudotime by Cell Type") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("PTC−LOH" = "#1F77B4", 
                                "Inj−PTC" = "#FF7F0E",
                                "Prolif−PTC" = "#91CC75",
                                "Inj−PTC−LOH" = "red"))

dev.off()

pdf("Geneset_Scores_vs_Pseudotime_Scatter.pdf", width = 12, height = 8)

ggplot(plot_data, aes(x = Pseudotime, y = Score, color = celltype)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_grid(Geneset ~ celltype, scales = "free_y") +
  labs(x = "Pseudotime", y = "Score", 
       title = "Geneset Scores Distribution by Cell Type") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8)) +
  scale_color_manual(values = c("PTC−LOH" = "#1F77B4", 
                                         "Inj−PTC" = "#FF7F0E",
                                         "Prolif−PTC" = "#91CC75",
                                         "Inj−PTC−LOH" = "red"))

dev.off()

#####旧代码########
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
# 加载必要包
library(monocle)
library(dplyr)
library(tidyr)

# 1. 拟时序差异分析（你的原有代码）
stares_de <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
stares_de <- stares_de[,c(5,2,3,4,1,6,7)] # 把gene放前面
stares_de <- stares_de[order(stares_de$qval), ]
write.csv(stares_de, "states_diff_gene.csv", row.names = F)  

# 2. 绘制热图并获取对象（你的原有代码）
Time_genes <- stares_de %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, 
                             show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

# 只选择q值最小的前100个基因（可根据需要调整数量）
top_genes <- head(stares_de$gene_short_name, 50) %>% as.character()
p <- plot_pseudotime_heatmap(cds[top_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapTop50.pdf", p, width = 5, height = 5)

# 只选择q值最小的前100个基因（可根据需要调整数量）
top_genes <- head(stares_de$gene_short_name, 100) %>% as.character()
p <- plot_pseudotime_heatmap(cds[top_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapTop100.pdf", p, width = 5, height = 10)
# 2. 获取差异基因列表
Time_genes <- stares_de %>% pull(gene_short_name) %>% as.character()

# 3. 提取表达矩阵和伪时间信息（关键修正部分）
expr_matrix <- exprs(cds[Time_genes,])
pseudotime <- pData(cds)$Pseudotime

# 4. 按伪时间排序细胞
cell_order <- order(pseudotime)
ordered_expr <- expr_matrix[, cell_order]
ordered_pseudotime <- pseudotime[cell_order]

# 5. 定义极值区域（前10%和后10%伪时间）
left_cutoff <- quantile(ordered_pseudotime, 0.05)
right_cutoff <- quantile(ordered_pseudotime, 0.95)

# 6. 识别极值高表达基因
left_cells <- which(ordered_pseudotime <= left_cutoff)
right_cells <- which(ordered_pseudotime >= right_cutoff)

# 标准化表达矩阵（Z-score）
zscore_matrix <- t(scale(t(ordered_expr)))

left_high_genes <- rownames(zscore_matrix)[
  rowMeans(zscore_matrix[, left_cells]) > 0.8  # Z-score > 1
]
right_high_genes <- rownames(zscore_matrix)[
  rowMeans(zscore_matrix[, right_cells]) > 0.8  # Z-score > 1
]

# 7. 输出基因列表
write.table(left_high_genes, "m2-mac_high_genes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(right_high_genes, "inf-mac_high_genes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 8. 保存带标注的完整结果
result_df <- data.frame(
  gene = Time_genes,
  Prolif_PTC_high = Time_genes %in% left_high_genes,
  Inj_PTC_high = Time_genes %in% right_high_genes
) %>% left_join(stares_de, by = c("gene" = "gene_short_name"))

write.csv(result_df, "annotated_time_dependent_genes.csv", row.names = FALSE)


# 10. 打印结果摘要
cat("=== 分析结果摘要 ===\n")
cat("Prolif-PTC高表达基因数量:", length(left_high_genes), "\n")
cat("示例基因:", paste(head(left_high_genes, 3), collapse = ", "), "\n\n")
cat("Inj-PTC高表达基因数量:", length(right_high_genes), "\n")
cat("示例基因:", paste(head(right_high_genes, 3), collapse = ", "), "\n")






library(ggpubr)
df <- pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
#View(df)

pdf("geom_density.pdf",width = 4,height = 3)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
dev.off()


######报错##########

#寻找以依赖于分支的方式调控的基因
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res_branch_point_one.csv", row.names = F)
#branch_point = 1产生的是node1处的差异基因
pdf("genes_branched_heatmap",width = 10,height = 7)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了
dev.off()

BEAM_res <- BEAM(cds[ordergene,], branch_point = 2, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res_branch_point_two.csv", row.names = F)



#手动设置颜色
ClusterName_color_panel <- c(
  "Epithelial Luminal" = "#DC143C", "Epithelial Other" = "#0000FF", "Endothelial" = "#20B2AA",
  "Fibroblast" = "#FFA500", "Monocyte" = "#9370DB", "B cell" = "#98FB98",
  "Epithelial Basal" = "#F08080", "Epithelial NE" = "#0000FF")  #, "Platelet" = "#20B2AA"
  pdf("geom_density2.pdf",width = 10,height = 7)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ scale_fill_manual(name = "", values = ClusterName_color_panel)+scale_color_manual(name = "", values = ClusterName_color_panel)
dev.off()
