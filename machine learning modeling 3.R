# 加载包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)
setwd("G:\\CKDurine\\11.group1_model\\22.model")
set.seed(123)

# 读取数据
train_data <- read.table("data.train.txt", header=T, sep="\t", row.names=1)
test_data <- read.table("data.test.txt", header=T, sep="\t", row.names=1)

# 读取LASSO筛选的8个基因
lasso_genes <- read.table("interGenes.txt")[,1]

# 转置数据（基因在列，样本在行）
train_data <- t(train_data)
test_data <- t(test_data)

# 提取分组信息函数
get_group <- function(ids) {
  sapply(strsplit(ids, "_"), function(x) x[length(x)])
}

# 获取分组信息
train_groups <- get_group(rownames(train_data))
test_groups <- get_group(rownames(test_data))

# 创建数据框（保留LASSO基因）
train_df <- data.frame(train_data[, lasso_genes], Type = factor(train_groups))
test_df <- data.frame(test_data[, lasso_genes], Type = factor(test_groups))

# 检查数据分布
table(train_df$Type)  # 训练集类别分布
table(test_df$Type)   # 测试集类别分布


# 数据标准化（如果尚未标准化）
train_df[, -ncol(train_df)] <- scale(train_df[, -ncol(train_df)])
test_df[, -ncol(test_df)] <- scale(test_df[, -ncol(test_df)])

# 定义交叉验证（针对小样本优化）
# 针对不平衡数据的优化控制参数
ctrl <- trainControl(
  method = "repeatedcv", 
  number = 5,          # 5折交叉验证
  repeats = 3,         # 重复3次提高稳定性
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  sampling = "up"      # 上采样处理类别不平衡
)

# 训练模型（简化版，减少计算量）
models <- list()
# 1. 随机森林（限制树的数量）
models$rf <- train(Type ~ ., data = train_df, 
                   method = "rf",
                   ntree = 100,
                   trControl = ctrl,
                   metric = "ROC")

# 2. SVM（简化参数调优）
# 更精细的核参数搜索
# 2. SVM精细调参
models$svm <- train(
  Type ~ .,
  data = train_df,
  method = "svmRadial",
  tuneGrid = expand.grid(
    sigma = seq(0.05, 0.2, length.out = 5),
    C = c(0.5, 1, 2, 5)
  ),
  class.weights = c(Control = 1, DKD = 0.8),
  trControl = ctrl,
  metric = "ROC"
)

# 3. XGBoost调参（增强泛化能力）
models$xgb <- train(
  Type ~ .,
  data = train_df,
  method = "xgbTree",
  tuneGrid = expand.grid(
    nrounds = 100,
    max_depth = 3,              # 限制树深度
    eta = 0.1,                  # 适中学习率
    gamma = 0.1,                # 增加正则化
    colsample_bytree = 0.8,
    min_child_weight = 2,
    subsample = 0.8
  ),
  trControl = trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE
  ),
  metric = "ROC"
)

# 4. 逻辑回归（验证数据质量）
models$glm <- train(
  Type ~ .,
  data = train_df,
  method = "glm",
  family = "binomial",
  trControl = trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE
  ),
  metric = "ROC"
)
# 更新模型评估函数以计算训练集和测试集的ROC及置信区间
eval_model <- function(model, train_data, test_data) {
  # 训练集评估
  train_pred_prob <- predict(model, train_data, type = "prob")[, "DKD"]
  train_roc <- roc(response = train_data$Type, predictor = train_pred_prob)
  train_ci <- ci.auc(train_roc)
  
  # 测试集评估
  test_pred_prob <- predict(model, test_data, type = "prob")[, "DKD"]
  test_roc <- roc(response = test_data$Type, predictor = test_pred_prob)
  test_ci <- ci.auc(test_roc)
  
  list(
    train_roc = train_roc,
    train_auc = auc(train_roc),
    train_ci_lower = train_ci[1],
    train_ci_upper = train_ci[3],
    test_roc = test_roc,
    test_auc = auc(test_roc),
    test_ci_lower = test_ci[1],
    test_ci_upper = test_ci[3]
  )
}

# [前面的数据加载和模型训练部分保持不变...]
# 评估模型后，立即锁定结果（不变部分省略...）
results <- lapply(models, eval_model, train_df, test_df)

# 生成表格时直接从results提取数据，确保与绘图完全一致
auc_results <- data.frame(
  Model = rep(names(results), each=2),
  Dataset = rep(c("Train", "Test"), times=length(results)),
  AUC = c(
    results$rf$train_auc, results$rf$test_auc,
    results$svm$train_auc, results$svm$test_auc,
    results$xgb$train_auc, results$xgb$test_auc,
    results$glm$train_auc, results$glm$test_auc
  ),
  CI = sprintf(
    "[%.3f-%.3f]",
    c(
      results$rf$train_ci_lower, results$rf$test_ci_lower,
      results$svm$train_ci_lower, results$svm$test_ci_lower,
      results$xgb$train_ci_lower, results$xgb$test_ci_lower,
      results$glm$train_ci_lower, results$glm$test_ci_lower
    ),
    c(
      results$rf$train_ci_upper, results$rf$test_ci_upper,
      results$svm$train_ci_upper, results$svm$test_ci_upper,
      results$xgb$train_ci_upper, results$xgb$test_ci_upper,
      results$glm$train_ci_upper, results$glm$test_ci_upper
    )
  )
)

# 绘图时使用表格中的数据进行标注（关键修改！）
pdf("4.5ROC_train_CI.pdf", width=4.5, height=4.5)
plot(results$rf$train_roc, col="#1F77B4", legacy.axes=T, main="Training Set ROC Curves")
plot(results$svm$train_roc, col="#FF7F0E", add=TRUE)
plot(results$xgb$train_roc, col="#2CA02C", add=TRUE)
plot(results$glm$train_roc, col="#9467BD", add=TRUE)


# 从表格提取数据生成图例
train_auc_data <- subset(auc_results, Dataset == "Train")
legend_text <- sprintf(
  "%s (AUC=%.3f %s)",
  c("RF", "SVM", "XGB", "GLM"),
  train_auc_data$AUC,
  train_auc_data$CI
)

legend("bottomright", legend=legend_text, 
       col=c("#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD"), lwd=2, cex=0.8)
dev.off()

# 测试集绘图同理修改
pdf("4ROC_test_CI.pdf", width=4.5, height=4.5)
plot(results$rf$test_roc, col="#1F77B4", legacy.axes=T, main="Validation ROC Curves")
plot(results$svm$test_roc, col="#FF7F0E", add=TRUE)
plot(results$xgb$test_roc, col="#2CA02C", add=TRUE)
plot(results$glm$test_roc, col="#9467BD", add=TRUE)

test_auc_data <- subset(auc_results, Dataset == "Test")
legend_text <- sprintf(
  "%s (AUC=%.3f %s)",
  c("RF", "SVM", "XGB", "GLM"),
  test_auc_data$AUC,
  test_auc_data$CI
)

legend("bottomright", legend=legend_text,
       col=c("#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD"), lwd=2, cex=0.8)
dev.off()

# 打印最终表格
print(auc_results)
write.csv(auc_results, "AUC_Results_with_CI.csv", row.names=FALSE)

# [后面的测试集ROC绘图和基因重要性分析保持不变...]

# 基因重要性分析（保持不变）
pdf("Feature_Importance.pdf", width=4, height=3)
par(mfrow=c(1,2))
topnum=3
plot(varImp(models$rf), top=topnum, main="RF - Top 3 Genes")
plot(varImp(models$xgb), top=topnum, main="XGB - Top 3 Genes")
plot(varImp(models$glm), top=topnum, main="GLM - Top 3 Genes")
plot(varImp(models$svm), top=topnum, main="SVM - Top 3 Genes")
dev.off()

# 保存重要基因（保持不变）
important_genes <- list(
  RF = rownames(varImp(models$rf)$importance)[1:5],
  XGB = rownames(varImp(models$xgb)$importance)[1:5],
  svm = rownames(varImp(models$svm)$importance)[1:5],
  glm = rownames(varImp(models$glm)$importance)[1:5]
)
write.table(important_genes$RF, "Important_Genes_RF.txt", quote=F, row.names=F)
write.table(important_genes$XGB, "Important_Genes_XGB.txt", quote=F, row.names=F)
write.table(important_genes$svm,"Important_Genes_svm.txt", quote=F, row.names=F)
write.table(important_genes$glm, "Important_Genes_glm.txt", quote=F, row.names=F)

