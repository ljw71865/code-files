setwd("G:\\CKDurine\\06.urine_annotation\\8.gsva")

library(Seurat)
 library(msigdbr) 
 library(GSVA) 
 library(tidyverse)
  library(clusterProfiler) 
  library(patchwork) 
  library(limma) 
  #rm(list=ls())

# 1. 提取目标细胞数据（严格使用celltype列）
prolif_cells <- colnames(merged_all)[merged_all$celltype == "Prolif-PTC"]
inj_cells <- colnames(merged_all)[merged_all$celltype == "Inj-PTC"]

# 2. 获取表达矩阵（Seurat v5.3正确方式）
expr_prolif <- LayerData(merged_all, assay="RNA", layer="data")[, prolif_cells]
expr_inj <- LayerData(merged_all, assay="RNA", layer="data")[, inj_cells]
combined_expr <- cbind(expr_prolif, expr_inj)

# 3. GSVA分析（KEGG基因集）
# 读取 .gmt 文件（最简写法）
gmt <- read.gmt("c2.cp.kegg.Hs.symbols.gmt")
# 2. 准备基因集列表（确保格式正确）
gene_sets <- split(gmt$gene, gmt$term)



gmt2 <- read.gmt("h.all.v2025.1.Hs.symbols.gmt")
# 2. 准备基因集列表（确保格式正确）
gene_sets <- split(gmt2$gene, gmt2$term)

# 1. 准备参数对象（适用于GSVA 1.40.0+）
param <- ssgseaParam(
  exprData = as.matrix(combined_expr),  # 确保是基因×细胞的矩阵
  geneSets = gene_sets,            # 基因集列表
  minSize = 5,                     # 单细胞建议降低至5-10（因基因检出率低）
  maxSize = 1000,                  # 可适当提高以保留更多通路
  alpha = 0.25,                    # 保持默认权重（0.25-0.75间调试）
  normalize = TRUE)               # 推荐TRUE（尤其当两亚群细胞数差异大时）


# 2. 运行分析
gsva_res <- gsva(param, verbose = TRUE)

# 4. 差异分析（精确匹配您的参数）
group <- factor(c(rep("Prolif-PTC", length(prolif_cells)), 
                  rep("Inj-PTC", length(inj_cells))),
                levels = c("Inj-PTC", "Prolif-PTC"))
design <- model.matrix(~ group)
fit <- lmFit(gsva_res, design)
fit <- eBayes(fit)
diff <- topTable(fit, coef=2, number=Inf, adjust.method="BH")
write.csv(diff, "diff_raw_results.csv", row.names = TRUE)

diff <- read.csv("diff_modified.csv", row.names = 1)
# 1. kegg差异分析结果处理（严格匹配您的要求）
# 使用FDR校正的p值
diff$group <- ifelse(diff$logFC > 0 & diff$adj.P.Val < 0.05, "up",
                     ifelse(diff$logFC < 0 & diff$adj.P.Val < 0.05, "down", "noSig"))
#diff$group <- ifelse(diff$logFC > 0.04 & diff$P.Value < 0.05, "up",
 #                    ifelse(diff$logFC < -0.02 & diff$P.Value < 0.05, "down", "noSig"))

# 2. 准备绘图数据（完全匹配您的处理方式）
plot_data <- diff %>%
  filter(group != "noSig") %>%
  mutate(
    hjust = ifelse(t > 0, 1, 0),
    nudge_y = ifelse(t > 0, -0.1, 0.1)
  ) %>%
  arrange(t) %>%
  rownames_to_column("ID") %>%
  mutate(ID = factor(ID, levels = ID))

# 3. 计算坐标轴范围
limt <- ceiling(max(abs(plot_data$t))) * 1.1  # 增加10%的边距
                
                # 4. 双向柱状图绘制（完全匹配您的可视化参数）
                pdf("Prolif_vs_Inj_PTC_GSVA_Result.pdf", width=7.5, height=6)
                ggplot(plot_data, aes(x=ID, y=t, fill=group)) +
                  geom_col(width=0.7, alpha=0.7) +
                  scale_fill_manual(values=c("up"="#08519C", "down"="#008020")) +
                  geom_text(
                    aes(label=ID, y=nudge_y),
                    hjust=plot_data$hjust,
                    size=3,
                    nudge_x=0,
                    nudge_y=0
                  ) +
      
                  labs(
                    x = "KEGG pathways",
                    y = 't value of GSVA score\nProlif-PTC | Inj-PTC',  # 正确字符串写法
                    title = "GSVA: Prolif vs Inj-PTC"
                  ) +
                  scale_y_continuous(limits=c(-limt, limt)) +
                  coord_flip() +
                  theme_bw() +
                  theme(
                    panel.grid=element_blank(),
                    panel.border=element_rect(size=0.6),
                    plot.title=element_text(hjust=0.5, size=18, face="bold"),
                    axis.text.y=element_blank(),
                    axis.title=element_text(hjust=0.5, size=14),
                    axis.line=element_blank(),
                    axis.ticks.y=element_blank(),
                    legend.position="none"
                  ) +
                  geom_vline(xintercept=0, linetype="dashed", color="gray40")
                dev.off()
  
#logFC > 0, 表示该通路在 Prolif-PTC 细胞中显著高表达（上调） 
 #logFC < 0,表示该通路在 Inj-PTC 细胞中显著高表达（上调）
                
     
                
  ##########hallmark##########
                setwd("G:\\CKDurine\\06.urine_annotation\\8.gsva\\hallmark_PTC")
                
                library(Seurat)
                library(msigdbr) 
                library(GSVA) 
                library(tidyverse)
                library(clusterProfiler) 
                library(patchwork) 
                library(limma) 
                #rm(list=ls())
                
                # 1. 提取目标细胞数据（严格使用celltype列）
                prolif_cells <- colnames(merged_all)[merged_all$celltype == "Prolif-PTC"]
                inj_cells <- colnames(merged_all)[merged_all$celltype == "Inj-PTC"]
                
                # 2. 获取表达矩阵（Seurat v5.3正确方式）
                expr_prolif <- LayerData(merged_all, assay="RNA", layer="data")[, prolif_cells]
                expr_inj <- LayerData(merged_all, assay="RNA", layer="data")[, inj_cells]
                combined_expr <- cbind(expr_prolif, expr_inj)
                
                # 3. GSVA分析（KEGG基因集）
              
                
                gmt2 <- read.gmt("h.all.v2025.1.Hs.symbols.gmt")
                # 2. 准备基因集列表（确保格式正确）
                gene_sets <- split(gmt2$gene, gmt2$term)
                
                # 1. 准备参数对象（适用于GSVA 1.40.0+）
                param <- ssgseaParam(
                  exprData = as.matrix(combined_expr),  # 确保是基因×细胞的矩阵
                  geneSets = gene_sets,            # 基因集列表
                  minSize = 5,                     # 单细胞建议降低至5-10（因基因检出率低）
                  maxSize = 1000,                  # 可适当提高以保留更多通路
                  alpha = 0.25,                    # 保持默认权重（0.25-0.75间调试）
                  normalize = TRUE)               # 推荐TRUE（尤其当两亚群细胞数差异大时）
                
                
                # 2. 运行分析
                gsva_res <- gsva(param, verbose = TRUE)
                
                # 4. 差异分析（精确匹配您的参数）
                group <- factor(c(rep("Prolif-PTC", length(prolif_cells)), 
                                  rep("Inj-PTC", length(inj_cells))),
                                levels = c("Inj-PTC", "Prolif-PTC"))
                design <- model.matrix(~ group)
                fit <- lmFit(gsva_res, design)
                fit <- eBayes(fit)
                diff <- topTable(fit, coef=2, number=Inf, adjust.method="BH")
                write.csv(diff, "hallamrk_diff_raw_results.csv", row.names = TRUE)
                
                diff <- read.csv("hallamrk_diff_modified.csv", row.names = 1)
                # 1. kegg差异分析结果处理（严格匹配您的要求）
                diff$group <- ifelse(diff$logFC > 0 & diff$P.Value < 0.05, "up",
                                     ifelse(diff$logFC < 0 & diff$P.Value < 0.05, "down", "noSig"))
                
                # 2. 准备绘图数据（完全匹配您的处理方式）
                plot_data <- diff %>%
                  filter(group != "noSig") %>%
                  mutate(
                    hjust = ifelse(t > 0, 1, 0),
                    nudge_y = ifelse(t > 0, -0.1, 0.1)
                  ) %>%
                  arrange(t) %>%
                  rownames_to_column("ID") %>%
                  mutate(ID = gsub("HALLMARK_", "", ID),  # 新增这行删除KEGG_前缀
                         ID = factor(ID, levels = ID))
                
                # 3. 计算坐标轴范围
                limt <- ceiling(max(abs(plot_data$t))) * 1.1  # 增加10%的边距
                
                # 4. 双向柱状图绘制（完全匹配您的可视化参数）
                pdf("hallmark_Prolif_vs_Inj_PTC_GSVA_Result.pdf", width=5, height=5)
                ggplot(plot_data, aes(x=ID, y=t, fill=group)) +
                  geom_col(width=0.7, alpha=0.7) +
                  scale_fill_manual(values=c("up"="#08519C", "down"="#008020")) +
                  geom_text(
                    aes(label=ID, y=nudge_y),
                    hjust=plot_data$hjust,
                    size=3,
                    nudge_x=0,
                    nudge_y=0
                  ) +
                  
                  labs(
                    x = "Hallmark pathways",
                    y = 't value of GSVA score\nProlif-PTC | Inj-PTC',  # 正确字符串写法
                    title = "GSVA: Prolif vs Inj-PTC"
                  ) +
                  scale_y_continuous(limits=c(-limt, limt)) +
                  coord_flip() +
                  theme_bw() +
                  theme(
                    panel.grid=element_blank(),
                    panel.border=element_rect(size=0.6),
                    plot.title=element_text(hjust=0.5, size=18, face="bold"),
                    axis.text.y=element_blank(),
                    axis.title=element_text(hjust=0.5, size=14),
                    axis.line=element_blank(),
                    axis.ticks.y=element_blank(),
                    legend.position="none"
                  ) +
                  geom_vline(xintercept=0, linetype="dashed", color="gray40")
                dev.off()
                
                #logFC > 0, 表示该通路在 Prolif-PTC 细胞中显著高表达（上调） 
                #logFC < 0,表示该通路在 Inj-PTC 细胞中显著高表达（上调）
                
                