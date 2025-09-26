setwd("G:\\CKDurine\\03.early_DKD\\1.Seurat")    
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

Early <- readRDS("G:/187geneMR/25.GSE131882/1.Seurat/merged_seurat_final.rds")


# 合并前提取第一个样本的基因顺序（作为基准）
#pre_merge_genes <- rownames(Early@assays$RNA@layers[["counts.Control1"]])


# 2. 合并数据层（关键步骤）
Early <- JoinLayers(Early)
# 确保assay的行名与counts层一致
rownames(Early@assays$RNA@layers$counts) <- rownames(Early)



Early[["percent.mt"]] <- PercentageFeatureSet(object = Early, pattern = "^MT-")
pdf(file="featureViolin.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = Early, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


##美化,pt.size = 0
pdf(file="featureViolin2.pdf",width=10,height=6)
Early[["percent.mt"]] <- PercentageFeatureSet(object = Early, pattern = "^MT-")
VlnPlot(Early,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0,group.by = "orig.ident")&NoLegend()&labs(x = '') #&geom_hline(yintercept = 2500)
dev.off()



pdf(file="featureCor.pdf",width=10,height=6)             #保存基因特征相关性图
plot1 <- FeatureScatter(object = Early, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = Early, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

Early <- NormalizeData(object = Early, normalization.method = "LogNormalize", scale.factor = 10000)
Early <- FindVariableFeatures(object = Early, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = Early), 10)  #在图片上标记出前十个
pdf(file="featureVar.pdf",width=10,height=6)                 #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = Early)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



table(Early$orig.ident)  #查看各类细胞数目
#  Control1   Control2   Control3 Early_DKD1 Early_DKD2 Early_DKD3 
#3732       2342       5219       3744       2682       2501 


###################################05.PCA主成分分析###################################
Early=ScaleData(Early)                   #PCA降维之前的标准预处理步骤
#Early=RunPCA(object= Early,npcs = 30,pc.genes=VariableFeatures(object = Early))     #PCA分析
Early <- RunPCA(Early,npcs = 30,verbose = FALSE)
#ElbowPlot(Early,ndims=50)
Early@meta.data$patient <- ifelse(
  grepl("^Control", Early$orig.ident),
  "Control",
  "Early_DKD"
)

# 验证分组
table(Early@meta.data[["patient"]])
#  Control Early_DKD 
#11293      8927

save(Early, file = "pbmc_after_pca.RDA")

options(repr.plot.height = 5, repr.plot.width = 12)

pdf(file="05.before harmony.pdf",width=6.5,height=6)
DimPlot(object = Early, reduction = "pca", pt.size = .1, group.by = "patient")
dev.off()

pdf(file="05.before harmony2.pdf",width=6.5,height=6)
VlnPlot(object = Early, features = "PC_1", group.by = "patient", pt.size = .1)
dev.off()

#
options(repr.plot.height = 2.5, repr.plot.width = 6)
Early <- Early %>% 
  RunHarmony("patient", plot_convergence = TRUE)

#获取Harmony 矫正之后的信息，使用Embeddings()函数
harmony_embeddings <- Embeddings(Early, 'harmony')
harmony_embeddings[1:5, 1:5]

##查看数据Harmony整合之后的前两个维度上数据是不是很好的整合，最好是很好的整合结果。
pdf(file="05.after harmony.pdf",width=12,height=5)
#options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = Early, reduction = "harmony", pt.size = .1, group.by = "patient")
p2 <- VlnPlot(object = Early, features = "harmony_1", group.by = "patient", pt.size = .1)
plot_grid(p1,p2)
dev.off()




#绘制每个PCA成分的相关基因,分为20个pc?#
pdf(file="051.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = Early, dims = 1:4, reduction = "pca",nfeatures = 20) 
#只显示出四个pc图#每个图展示20个基因
dev.off()

#主成分分析图形
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = Early, reduction = "pca")
dev.off()

#主成分分析热图#只展示前四个热图#每张图显示30个基因#一行两个图#
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = Early, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#每个PC的p值分布和均匀分布
Early <- JackStraw(object = Early, num.replicate = 100)
Early <- ScoreJackStraw(object = Early, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = Early, dims = 1:20)  #解螺旋是15,生信自学网20
dev.off()
save(Early, file = "pbmc_ScoreJackStraw_20220.rds")

###################################06.TSNE聚类分析和marker基因###################################
setwd("G:\\187geneMR\\25.GSE131882\\1.Seurat")
####UMAP/tSNE聚类分析#综合解螺旋的代码#比较好用#
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(Seurat)
library(patchwork) 
pcSelect=30
Early <- FindNeighbors(Early,reduction = "harmony", dims = 1:pcSelect)       #计算邻接距离 #解螺旋取得值为10,生信自学网为20          
Early <- FindClusters(Early, reduction = "harmony",resolution = 0.4)        #对细胞分组,优化标准模块化 
####选择不同的resolution值可以获得不同的cluster数目，值越大cluster数目越多，默认值是0.5.         
Early <- RunTSNE(Early, reduction = "harmony",dims = 1:pcSelect)      #TSNE聚类 #解螺旋取得值为10,生信自学网为20               

TSNEPlot(Early, label = TRUE)
ggsave(filename = "04.TSNE2.pdf",width=6.5,height=6)
#write.table(Early$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)
Early <- RunUMAP(Early,reduction = "harmony", dims = 1:pcSelect)    #解螺旋取得值为10,生信自学网为20     


DimPlot(Early,reduction = "umap",label = TRUE)
 #UAMP可视化#标注序号# pt.size = 2为点的大小#pt.size = 2
ggsave(filename = "06.UMAP.pdf",width=3,height=4)
dev.off()

DimPlot(Early,pt.size = 1,reduction = "umap",group.by="orig.ident",label = TRUE,label.size = 3)
 #UAMP可视化#标注序号# pt.size = 2为点的大小#
ggsave(filename = "06.UMAP_sample.pdf",width=6.5,height=6)
dev.off()

DimPlot(Early,pt.size = 1,reduction = "umap",group.by="group",label = TRUE,label.size = 3)
#UAMP可视化#标注序号# pt.size = 2为点的大小#
ggsave(filename = "06.UMAP_group.pdf",width=6.5,height=6)
dev.off()


# 寻找所有以聚类的差异基因`
cluster.markers <- FindAllMarkers(
  object = Early,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  only.pos = TRUE,  # 只保留上调基因
  verbose = TRUE
)
significant.markers <- subset(cluster.markers, 
                              p_val_adj < 0.05 & abs(avg_log2FC) > 2)

# 3. 按cluster和log2FC排序输出
significant.markers <- significant.markers[order(
  significant.markers$cluster, 
  -abs(significant.markers$avg_log2FC)
), ]

sig50_markers <- subset(cluster.markers, p_val_adj < 0.05 & avg_log2FC > 1)

# 2. 提取每个cluster前50基因
top50 <- sig50_markers %>% 
  group_by(cluster) %>% 
  top_n(50, avg_log2FC) %>% 
  arrange(cluster, -avg_log2FC)

# 3. 输出到CSV
write.csv(top50, "top50_markers_per_cluster.csv", row.names = FALSE)
# 4. 输出完整结果和筛选结果
write.table(cluster.markers, 
            file = "all_cluster_markers.xls",
            sep = "\t", 
            row.names = TRUE, 
            quote = FALSE)

write.table(significant.markers,
            file = "significant_cluster_markers.xls",
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)

###################################07.注释细胞类型###############################
save(Early, file = "Early_annotation_20000.rda")

setwd("G:\\CKDurine\\15.early_DKD\\3.annotation\\新建文件夹")
pdf(file="02CD4_CD8.pdf",width=17,height=6)
#CD4_CD8
cluster4Marker=c("PDGFRB","FN1","COL4A2","NPHS1","PODXL","NPHS2","SLC4A11","AQP9","SLC34A1", "CUBN", "LRP2", "SLC22A6", "GPX3", "ALDOB", "SLC22A8","ALDH1A2","ALDH2") #uster4中logFC最正和最负值的五个基因
DotPlot(object = Early, features = cluster4Marker)
dev.off()


pdf(file="1损伤PCT.pdf",width=17,height=6)
cluster4Marker=c("FN1", "VCAM1", "CD74", "TYROBP", "CTSB", 
                 "CD44", "HMOX1", "SERPINE1", "ANXA1", 
                 "S100A11", "S100A10", "VIM", "LGALS1",
                 "CHI3L1", "LCN2", "TIMP1", "B2M","HAVCR1")
DotPlot(object = Early, features = cluster4Marker)
dev.off()

pdf(file="2损伤PCT.pdf",width=17,height=6)
cluster4Marker=c("SLC4A11", "AQP9", "ALDH1A2", "LRP2", "CCN2", "VIM", 
                 "THBS1", "COL4A1", "COL4A2", "FNBP1", "PLOD2", "ITGB3", 
                 "ITGA3", "HAVCR1", "PDGFB", "ZEB2", "NFKB1", "C3", 
                 "CXCL1", "IL34", "PCNA", "PROM1", "EPCAM")
DotPlot(object = Early, features = cluster4Marker)
dev.off()



pdf(file="代谢活跃_PCT.pdf",width=18,height=6)
cluster4Marker=c("CYP4A11", "PCK1", "FABP1", "SLC5A12", "HAO2", 
                 "ACSM2A", "SLC13A3", "UGT1A8", "DPEP1", "SLC22A8",
                 "SLC4A11","ALDOB","LRP2", "CUBN","SLC34A1", "SLC22A6", "GPX3")
DotPlot(object = Early, features = cluster4Marker)
dev.off()


# table(Early@active.ident)
pdf(file="损伤pct.pdf",width=19,height=6)  #width=10.5,height=8
cluster4Marker=c("SLC4A11","AQP9","ALDH1A2","LRP2","CCN2","VIM","THBS1","COL4A1","COL4A2","FNBP1","PLOD2","ITGB3","ITGA3","HAVCR1","PDGFB","ZEB2","NFKB1","C3","CXCL1","IL34","PCNA","PROM1","EPCAM")
DotPlot(object = Early, features = cluster4Marker)
dev.off()

#######分组展示
DimPlot(Early,
        pt.size = 1,
        reduction = "tsne",
        group.by = "seurat_clusters",
        split.by = "group",  # 按group列拆分
        label = TRUE,
        label.size = 3) +
  plot_annotation(title = "Normal vs DKD (Clusters 0-14)")

ggsave("06.UMAP_Normal_vs_DKD_split.pdf", width = 12, height = 6)
dev.off()


#重新对注释后的细胞可视化
#Early <- subset(Early, idents = c("21"), invert = TRUE)#去掉低质量细胞群
new.cluster.ids <- c("0"="CD-PC", 
                     "1"="PTC", 
                     "2"="DCT", 
                     "3"="LOH", 
                     "4"="CD-ICA", 
                     "5"="LOH", 
                     "6"="CD-PC", 
                     "7"="Inj-PTC", 
                     "8"="PTC", 
                     "9"="EC", 
                     "10"="Inj-PTC", 
                     "11"="Podo",
                     "12"="CD-ICB",
                     "13"="Mes",
                     "14"="Imm")

#save(Early, file = "pbmc_for_cell_markers.RDA") #按照上面注释,生成此文件


Early <- RenameIdents(Early, new.cluster.ids)                        
Early$celltype <- Early@active.ident #将注释结果添加到metadata
#也可以从这一步开始进行monocle分析,参见简书代码

my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8542')

cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
                "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542") #,,
DimPlot(Early, reduction = "umap",group.by = "celltype",label = T,pt.size = 0.2,cols=cbPalette,repel=T)#
#DimPlot(Early,pt.size = 0.2,reduction = "umap",label = TRUE)
#UAMP可视化#标注序号# pt.size = 2为点的大小#
ggsave(filename = "注释后.sample_cellUMAP.pdf",width=3.5,height=3)
dev.off()
TSNEPlot(Early, pt.size = 0.2, label = TRUE,cols=cbPalette)
ggsave(filename = "注释后.ample_cellTSNE.pdf",width=7.5,height=6)
dev.off()
save(Early, file = "2Early_annotated.RDA")
#按照sample对cluster进行注释

#biomamama
## umap/tsne
library(RColorBrewer)
color_ct=c(brewer.pal(12, "Set3"),"#b3b3b3",
           brewer.pal(5, "Set1"),
           brewer.pal(3, "Dark2"),
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
clusterCornerAxes(Early,reduction = 'umap',clusterCol = 'celltype',pSize = 0.05,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(cbPalette,0.65)) + NoLegend() +
  scale_fill_manual(values = alpha(cbPalette,0.65))
#pSize = 0.5为细胞点的大小;cellLabelSize = 5文字大小
ggsave("bioma_umap.pdf",width = 3,height = 3,dpi = 600)

clusterCornerAxes(object = Early,reduction = 'tsne',clusterCol = 'celltype',cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(cbPalette,0.65)) +  NoLegend() +
  scale_fill_manual(values = alpha(cbPalette,0.65))
ggsave("bioma_tsne.pdf",width = 3,height = 3,dpi = 600)

clusterCornerAxes(object = Early,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("bioma_umap_indi.pdf",width = 5.5,height = 5,dpi = 600)

clusterCornerAxes(object = Early,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("bioma_tsne_indi.pdf",width = 5.5,height = 5,dpi = 600)

cluster4Marker=c(
  #集合管主细胞
  "SCNN1G", "SCNN1B", "AQP2",  
  # PCT         
  "ALDOB","LRP2" ,"CUBN",   
  #DCT           
  "SLC12A3","CNNM2", "TRPM6", 
      #LOH       
   "SLC12A1", "UMOD", "CLDN16", 
     
      #集合管细胞-A型闰细胞          
 "SLC26A7", "AQP6",  
  
  #内皮         
   "PECAM1", "EMCN", 
        #足细胞         
  "NPHS1", "NPHS2",  
  # 集合管细胞-B型闰细胞      
   "FOXI1",  "ATP6V1G3","SLC26A4",
          #系膜细胞     
   "PDGFRB","ACTA2",   
           #免疫细胞      
   "PTPRC")

#,"FBP1","PDK4","RHCG"
# 使用更美观的渐变色方案（深蓝-浅蓝-白-橙-红）
pdf("2细胞注释气泡图.pdf", width=9, height=5)

DotPlot(
  object = Early,
  features = cluster4Marker,
  cols = c("#ffffbf", "#d73027"),  # 浅黄→深红渐变
  col.min = 0,                     # 匹配图例最小值
  col.max = 10,                    # 匹配图例最大值
  dot.scale = 10                    # 适当调大气泡
) + 
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=12),
    panel.grid.major = element_line(color="grey90"),)  # 添加浅灰色网格线

dev.off()
# 寻找聚类4和聚类7的差异基因
clusterL_NE.markers <- FindMarkers(Early, ident.1 = "Luminal", ident.2 = "Luminal/NE", min.pct = 0.25,logfc.threshold = 0.5)
head(clusterL_NE.markers, n = 5)
write.table(clusterL_NE.markers,file="clusterL_NE.markers.xls",sep="\t",row.names=T,quote=F)
# 寻找聚类1和所以聚类的差异基因
cluster3.markers <- FindMarkers(Early, ident.1 = "EX CD8+T", min.pct = 0.25,logfc.threshold = 0.5)
head(cluster3.markers, n = 5)
write.table(cluster3.markers,file="EX CD8+T_markers.xls",sep="\t",row.names=T,quote=F)


#NEPC_CRPC比较
marker <- FindMarkers(Early, ident.1 = c(12,21), ident.2= c(7,10),min.pct = 0.25,only.pos = T,logfc.threshold = 0.5)
head(marker, n = 5)
marker<- subset(marker,p_val_adj<0.05)
marker=rbind(id=colnames(marker),marker)
write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",quote=F,col.names=F)
#write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",row.names=T,quote=F)
###min.pct表示基因在多少细胞中表达的阈值，only.pos = TRUE表示只求高表达的基因，

#NEPC_NS比较
marker <- FindMarkers(Early, ident.1 = c(12,21), ident.2= 3,min.pct = 0.25,only.pos = T,logfc.threshold = 0.5)
head(marker, n = 5)
marker<- subset(marker,p_val_adj<0.05)`
marker=rbind(id=colnames(marker),marker)
write.table(marker,file="NEPC_NS_Luminal_0.25.markers.xls",sep="\t",quote=F,col.names=F)
#write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",row.names=T,quote=F)
###min.pct表示基因在多少细胞中表达的阈值，only.pos = TRUE表示只求高表达的基因，

#ͬ同上绘制marker的小提琴图

VlnPlot(Early,features = c("ALDH2", "FN1"),pt.size = 0,group.by = "celltype")&NoLegend()&labs(x = '') #&geom_hline(yintercept = 2500)
ggsave(filename = "ALDH2_FN1_violin.pdf",width=10,height=6)  

VlnPlot(Early, features = c("ALDH2", "FN1"), slot = "counts", log = TRUE)
ggsave(filename = "ALDH2_FN1_violin_Violin log(count).pdf",width=10,height=6) #取log值


#同上,选择性绘制不同cluster中某些基因的散点图#和上面没啥区别?
pdf(file="ALDH2_FN1_scatter.pdf",width=10,height=6)
cluster4Marker=c("FN1", "ALDH2") 
FeaturePlot(Early, reduction="umap",features =  cluster4Marker)#cols = c("green", "red")
dev.off()

###组内差异分析##############
setwd("G:\\RCC\\GSE181061_ATAC_TCR_RNA\\scRNAseq_Tcell\\3.annotation")
library(dplyr)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
Early.markers <- FindAllMarkers(Early, min.pct = 0.25,logfc.threshold = 1) 
#因基因太多logfc.threshold = 1,基因少改为0.5#only.pos = TRUE,
Early.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)   ##n=2与n=10运行结果没多大区别?##
write.table(Early.markers,file="all_1.xls",sep="\t",row.names=F,quote=F)
Early = Early[, Idents(Early) %in% c( "EX CD8+T" , "C1-RM CD8+T","C2-RM CD8+T","Cytotoxic CD8+T")]#

###此代码出的结果过有部分PadjValue值不满足小于等于0.05,需要自己手工筛选#
Early.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#输出tsneTop10差异基因热图
DoHeatmap(Early, features = top10$gene) + NoLegend()
ggsave(filename = "02Top5gene_heatmap.pdf",width=15,height=15)
#save(Early, file = "pbmc_for_markers.rda")



#结合解螺旋的代码进行14个cluster注释#可成功运行
setwd("G:\\RCC\\GSE181061_ATAC_TCR_RNA\\scRNAseq_Tcell\\3.annotation") 
test <- GetAssayData(Early, slot="data") 
clusters<-Early@meta.data$seurat_clusters   
library(celldex)
library(SingleR)
HumanCell <- celldex::HumanPrimaryCellAtlasData()  #加载人类的
singler <- SingleR(test, ref = HumanCell,
                   labels = HumanCell$label.main,clusters = clusters)  
#test为需要注释的矩阵,ref = HumanCell参考矩阵为人类细胞文件#label.main主标签#  clusters为之前得到的聚类结果14个              
write.table(singler,file="03.17个cluster聚类注释2.txt",quote=F,sep="\t",col.names=F)
#包含微调前（first.labels）、微调后（labels）以及修剪后（pruned.labels),选择labels即可###

#对4504个细胞进行细胞注释
singler2 <- SingleR(test, ref = HumanCell, labels = HumanCell$label.main) #对每个细胞进行细胞注释
write.table(singler2,file="07.各个细胞类型注释-解.txt",quote=F,sep="\t",col.names=F)  
#需比对生信自学网产生的数据类型,看是否一致#

# 比较Seurat聚类和SignleR聚类的结果
table(singler2$labels,clusters)
write.table(table(singler2$labels,clusters),file="07.SingleR-Seurat聚类和SignleR聚类的结果对比.txt",quote=F,sep="\t",col.names=T)


#SingleR-绘制细胞类型结果热图Heatmap
pdf(file="07.SingleR-绘制细胞类型结果Heatmap.pdf",width=10,height=8)
plotScoreHeatmap(singler)
dev.off()



#组间差异分析########
logFCfilter=0.5
adjPvalFilter=0.05

Early.markers=FindMarkers(Early, ident.1 = "Responder", ident.2 = "Non.Responder", group.by = 'patient')
Early.markers=cbind(Gene=row.names(Early.markers), Early.markers)
write.table(Early.markers,file="0.5allGene.txt",sep="\t",row.names=F,quote=F)
###vg_log2FC：记录两组之间平均表达的倍数变化
#。正值表示该基因在第一组中表达更高。


sig.markers=Early.markers[(abs(as.numeric(as.vector(Early.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(Early.markers$p_val_adj))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="0.5diffGene.txt",sep="\t",row.names=F,quote=F)

