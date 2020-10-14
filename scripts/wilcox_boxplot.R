setwd("../summary_data/*/")
#read the datasets
celltypes<-read.table('../sepsis/summary_data/cell_types/CombinedLM22.tsv', 
                      sep = '\t', header = TRUE)
head(celltypes)
gene<-read.table('../sepsis/summary_data/expression_data/symb_all.tsv', 
                 sep = '\t', header = TRUE)
head(gene)
dim(celltypes)
dim(gene)

#statistical analysis

#Wilcoxon test
survivor = celltypes[celltypes$survivor==0,]
death = celltypes[celltypes$survivor==1,]
names(celltypes)
#pull out columns that are statistically significant with survivor
for (i in 9:35){
  print(names(celltypes)[i])
  result = wilcox.test(survivor[,i],death[,i])
  print(result)}
for (i in 5:6){
  print(names(celltypes)[i])
  result = wilcox.test(survivor[,i],death[,i])
  print(result)}

#logistic regression
install.packages('devtools')
devtools::install_github("hadley/ggplot2")
devtools::install_github("sachsmc/plotROC")
library(MASS) # for Pima data sets
library(ggplot2)
library(plotROC)

# train model on training data
model_glm <- glm(survivor~T_cells+T_cells_CD4+Macrophages+Mast_cells_resting+Neutrophils+age, 
                 data = celltypes, family = 'binomial',
                 na.action = na.exclude)


# ROC curves
library(pROC)
test_prob = predict(model_glm, newdata = celltypes, type = "response")
test_roc = roc(celltypes$survivor ~ test_prob, plot = TRUE, print.auc = TRUE)

# combine linear predictor and known truth for training and test datasets into one data frame
df <- rbind(data.frame(predictor = predict(model_glm, celltypes),
                       known.truth = celltypes$survivor,
                       model = "train"))
# ,
# data.frame(predictor = predict(model_glm, Pima.te),
#            known.truth = Pima.te$type,
#            model = "test"))

ggplot(df, aes(d = known.truth, m = predictor, color = model)) + 
  geom_roc(n.cuts = 0)

coef(model_glm)
head(predict(model_glm, type = 'response'))
model_glm_pred = ifelse(predict(model_glm, type = "response") > 0.5, "Yes", "No")
calc_class_err = function(actual, predicted) {
  mean(actual != predicted)
}

calc_class_err(actual = celltypes$survivor, predicted = model_glm_pred)
train_tab = table(predicted = model_glm_pred, actual = celltypes$survivor)


get_logistic_error = function(mod, data, res = "y", pos = 0, neg = 1, cut = 0.5) {
  probs = predict(mod, newdata = celltypes, type = "response", na.action = na.exclude)
  preds = ifelse(probs > cut, pos, neg)
  calc_class_err(actual = data[, res], predicted = preds)
}
get_logistic_error(model_glm, data = celltypes, 
                   res = "survivor", pos = "Yes", neg = "No", cut = 0.5)

plot(survivor ~ T_cells+age, data = celltypes, 
     col = "darkorange", pch = "|", ylim = c(-0.2, 1),
     main = "Using Logistic Regression for Classification")
abline(h = 0, lty = 3)
abline(h = 1, lty = 3)
abline(h = 0.5, lty = 2)
curve(predict(model_glm, data.frame(T_cells+age = x),  type = "response"), 
      add = TRUE, lwd = 3, col = "dodgerblue")
abline(v = -coef(model_glm)[1] / coef(model_glm)[2], lwd = 2)



#plot
par(mfcol = c(2,3))
boxplot(celltypes$age~celltypes$survivor, xlab = 'Survivor', ylab = 'Age')
boxplot(celltypes$T_cells~celltypes$survivor, xlab = 'Survivor', ylab = 'T cells')
boxplot(celltypes$T_cells_CD4~celltypes$survivor, xlab = 'Survivor', ylab = 'CD4+ T cells')
boxplot(celltypes$Macrophages~celltypes$survivor, xlab = 'Survivor', ylab = 'Macrophages')
boxplot(celltypes$Mast_cells_resting~celltypes$survivor, xlab = 'Survivor', ylab = 'Mast cells')
boxplot(celltypes$Neutrophils~celltypes$survivor, xlab = 'Survivor', ylab = 'Neutrophils')


install.packages("ggpubr")
library(ggpubr)

p <- ggboxplot(lm22, x = "Survivor", y = "T_cells_CD4",
                               order = c("YES", "NO"),
                               color = "#2E86C1", add = "jitter",
                               bxp.errorbar = TRUE, xlab = FALSE,
                               ylab = "T cells (CD4+)") + 
       font("ylab", size = 16) +
       font("xy.text", size = 16)+
       scale_x_discrete(labels = c("YES" = "Survived", "NO" = "Died")) + 
       stat_compare_means(label = "p.format", label.x = 1.5, size = 6)
ggsave("lm22CD4Tcells.png", path = "../figures/", dpi = 300)