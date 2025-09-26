
#install.packages("rms")
#install.packages("rmda")


#引用包
library(rms)
library(rmda)

inputFile="data.test.txt"             #表达数据文件
geneFile="importanceGene.XGB.txt"     #基因列表文件
setwd("G:\\CKDurine\\11.group1_model\\23.Nomo")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#读取基因列表文件,提取疾病特征基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息
data=t(data)
# 建议改为：
group=ifelse(grepl("_DKD$", row.names(data)), "DKD", "Control")
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type~ PDK4+FBP1+RHCG, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.1, 0.9),
	lp=F, funlabel="DKD Risk Probability")
#输出列线图
pdf("Nomo.pdf", width=5, height=4)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.1, height=4)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

#绘制决策曲线
rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ PDK4+FBP1+RHCG, data=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)
#输出DCA图形
pdf(file="DCA.pdf", width=5, height=5)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()

pdf(file = "DCA.pdf", width = 5, height = 5)
plot_decision_curve(dc,
                    curve.names = "Model",
                    xlab = "Threshold probability",
                    cost.benefit.axis = TRUE,
                    col = "red",
                    confidence.intervals = FALSE,
                    standardize = FALSE,
                    legend.x = 0.7,  # 水平位置 (0-1)
                    legend.y = 0.2)  # 垂直位置 (0-1)
dev.off()
