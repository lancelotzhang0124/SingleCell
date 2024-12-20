# Load required libraries
library(survival)
library(caret)
library(survminer)
library(pROC)
library(cowplot)

# Define result directory
res_home <- "/home/lfZhang/sc/Figure/"

# Load dataset and select specific proteins
df <- read.csv(file = 'Cox.csv', sep = ",")
# 33 features with non-zero coefficients were selected.
specific_proteins <- c("ATP1B1", "CLU", "CNTN5", "CRHBP", "CXCL14", "DCC", "DPP10", "ERBB4", "GAPDH", "GPC5", "IDS", "KIT", "NCAM2", "NEFL", "NPTXR", "NRGN", "NXPH1", "PCDH7", "PDE5A", "PEBP1", "PILRB", "PKD1", "PRKG1", "PVALB", "ROBO2", "SAT2", "SLIT2", "SNAP25", "SPARCL1", "THY1", "TMSB10", "TSPYL1", "VSNL1")

trainData[sapply(trainData, is.character)] <- lapply(trainData[sapply(trainData, is.character)], as.factor)
testData[sapply(testData, is.character)] <- lapply(testData[sapply(testData, is.character)], as.factor)

formula_str_basic <- "Surv(follow_Duration, MDD==1) ~ Age + Sex_cate + family_history + Townsend_index_cate + Qualifications_cate + Smoking_status_cate + Alcohol_status_cate + BMI_cate + IPAQ_cate + chronic_num + Drug +"
formula_str <- paste(formula_str_basic, paste(specific_proteins, collapse = " + "))
cox_formula <- as.formula(formula_str)

x_train <- model.matrix(cox_formula, data = trainData)[, -1]  
x_test <- model.matrix(cox_formula, data = testData)[, -1]

coxboost_fit <- CoxBoost(
  time = trainData$follow_Duration,
  status = trainData$MDD == 1,
  x = x_train,
  stepno = 500
)
summary(coxboost_fit)

# Five features with non-zero coefficients were selected.
specific_proteins <- c("NEFL", "PKD1", "TMSB10", "CNTN5", "NCAM2")
df_subset <- df[, c("MDD", specific_proteins)]

# Split dataset into training (70%) and testing (30%) sets
set.seed(123)
trainIndex <- createDataPartition(df_subset$MDD, p = 0.7, list = FALSE)
trainData <- df[trainIndex, ]
testData <- df[-trainIndex, ]

# Define and fit the Cox model
formula_str <- paste(
  "Surv(follow_Duration, MDD==1) ~ Sex_cate + Townsend_index_cate + Qualifications_cate +",
  "Smoking_status_cate + Alcohol_status_cate + BMI_cate + Age + IPAQ_cate + family_history +",
  "chronic_num + Drug +", paste(specific_proteins, collapse = " + ")
)
cox_model <- coxph(as.formula(formula_str), data = trainData)

# Predict risk scores and calculate C-index
testData$risk_score <- predict(cox_model, newdata = testData, type = "risk")
c_index <- survConcordance(Surv(follow_Duration, MDD == 1) ~ risk_score, data = testData)$concordance
print(paste("C-index:", c_index))

# Kaplan-Meier survival analysis
testData$risk_group <- ifelse(testData$risk_score > median(testData$risk_score), "High", "Low")
surv_object <- Surv(testData$follow_Duration, testData$MDD == 1)
fit <- survfit(surv_object ~ risk_group, data = testData)
surv_plot <- ggsurvplot(
  fit, data = testData, pval = TRUE, title = "Survival Analysis",
  xlab = "Time (Months)", ylab = "Survival Probability",
  legend.labs = c("High Risk", "Low Risk"), palette = c("#b33d47", "#2a7306"),
  ggtheme = theme_bw(), conf.int = FALSE, risk.table = TRUE
)

# Save Kaplan-Meier plot
ggsave(filename = paste0(res_home, "curve_cox.png"), plot = surv_plot$plot, height = 10, width = 8, dpi = 300)

# ROC analysis
time_point <- median(testData$follow_Duration)
testData$event_binary <- ifelse(testData$follow_Duration <= time_point & testData$MDD == 1, 1, 0)
roc_result <- roc(testData$event_binary, testData$risk_score)
print(paste("AUC at time point", time_point, ":", auc(roc_result)))

# Save ROC curve plot
png(paste0(res_home, "roc_curve_cox.png"), width = 6, height = 6, units = "in", res = 300)
plot(roc_result, main = "ROC Curve", col = "#cc3c22", lwd = 3)
dev.off()
