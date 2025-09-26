#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05     #p值过滤条件
adjPvalFilter=1       #矫正后的p值过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

setwd("G:\\CKDurine\\19.go\\125 genes\\2.kegg")      #设置工作目录
rt=read.csv("aging_genes_diff_results.csv", header=T, sep=",", check.names=F)     #读取输入文件

#提取疾病相关基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,"exposure"]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=entrezIDs[!is.na(entrezIDs)]
gene=as.character(entrezIDs)
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
#kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)

#基因id转换成基因名字
kk = setReadable(kk, #前面分析的结果
                 OrgDb = "org.Hs.eg.db",    #人类数据库
                 keyType = "ENTREZID")      #要转换的基因类型
KEGG=as.data.frame(kk)
KEGG=KEGG[KEGG$category != "Human Diseases",]

#输出显著富集的结果
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#######修改
KEGG <- read.delim(
  "KEGG.txt",
  header = TRUE,
  row.names = 3,       # 用第3列（ID）作为行名
  stringsAsFactors = FALSE
)


#设置需要展示的通路
showTerm=as.vector(KEGG[,"Description"])
if(length(showTerm)>30){
	showTerm=showTerm[1:30]
}

#柱状图
pdf(file="barplot.pdf", width=5, height=5)
barplot(kk, drop=TRUE, showCategory=showTerm, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=5, height=5)
dotplot(kk, showCategory=showTerm, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#绘制基因和通路的关系图
pdf(file="cnetplot.pdf", width=10, height=7)
cnet=cnetplot(kk, circular=TRUE, showCategory=showTerm, colorEdge=TRUE)
print(cnet)
dev.off()

#绘制基因和通路的关系图
pdf(file="cnetplot.pdf", width=9, height=6.5)
cnet=cnetplot(kk, circular=TRUE, showCategory=showTerm, colorEdge=TRUE)
print(cnet)
dev.off()

