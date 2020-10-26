setwd("../summary_data/*/")
#read the datasets
celltypes<-read.table('../summary_data/cell_types/CombinedLM22.tsv', 
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

# split training and test data
table(celltypes$survivor)
trainingData <- celltypes[celltypes$use == 'discovery',]

# train model on training data
model_glm <- glm(survivor~Mast_cells_resting+Macrophages+T_cells_CD4, 
                 data = trainingData, family = 'binomial',
                 na.action = na.exclude)


# add an age group column
celltypes$AgeGroup <- cut(celltypes$age, breaks = c(0, 18, 65, Inf), 
                          labels = c('Children', 'Adults', 'Elderly'))

EMEXP3567.test <- celltypes[celltypes$accession == 'E-MEXP-3567'|celltypes$AgeGroup == 'Adults',]
EMTAB4421.test <- celltypes[celltypes$accession == 'E-MTAB-4421'|celltypes$AgeGroup == 'Adults',]
gse10474.test <- celltypes[celltypes$accession == 'GSE10474'|celltypes$AgeGroup == 'Adults',]
gse21802.test <- celltypes[celltypes$accession == 'GSE21802'|celltypes$AgeGroup == 'Adults',]
gse27131.test <- celltypes[celltypes$accession == 'GSE27131'|celltypes$AgeGroup == 'Adults',]
gse32707.test <- celltypes[celltypes$accession == 'GSE32707'|celltypes$AgeGroup == 'Adults',]
gse33341.test <- celltypes[celltypes$accession == 'GSE33341'|celltypes$AgeGroup == 'Adults',]
gse40586.test <- celltypes[celltypes$accession == 'GSE40586'|celltypes$AgeGroup == 'Adults',]
gse54514.test <- celltypes[celltypes$accession == 'GSE54514'|celltypes$AgeGroup == 'Adults',]
gse63042.test <- celltypes[celltypes$accession == 'GSE63042'|celltypes$AgeGroup == 'Adults',]
gse63990.test <- celltypes[celltypes$accession == 'GSE63990'|celltypes$AgeGroup == 'Adults',]
gse66890.test <- celltypes[celltypes$accession == 'GSE66890'|celltypes$AgeGroup == 'Adults',]

EMEXP3567.elder <- celltypes[celltypes$accession == 'E-MEXP-3567'|celltypes$AgeGroup == 'Elderly',]
EMTAB4421.elder <- celltypes[celltypes$accession == 'E-MTAB-4421'|celltypes$AgeGroup == 'Elderly',]
gse10474.elder <- celltypes[celltypes$accession == 'GSE10474'|celltypes$AgeGroup == 'Elderly',]
gse21802.elder <- celltypes[celltypes$accession == 'GSE21802'|celltypes$AgeGroup == 'Elderly',]
gse27131.elder <- celltypes[celltypes$accession == 'GSE27131'|celltypes$AgeGroup == 'Elderly',]
gse32707.elder <- celltypes[celltypes$accession == 'GSE32707'|celltypes$AgeGroup == 'Elderly',]
gse33341.elder <- celltypes[celltypes$accession == 'GSE33341'|celltypes$AgeGroup == 'Elderly',]
gse40586.elder <- celltypes[celltypes$accession == 'GSE40586'|celltypes$AgeGroup == 'Elderly',]
gse54514.elder <- celltypes[celltypes$accession == 'GSE54514'|celltypes$AgeGroup == 'Elderly',]
gse63042.elder <- celltypes[celltypes$accession == 'GSE63042'|celltypes$AgeGroup == 'Elderly',]
gse63990.elder <- celltypes[celltypes$accession == 'GSE63990'|celltypes$AgeGroup == 'Elderly',]
gse66890.elder <- celltypes[celltypes$accession == 'GSE66890'|celltypes$AgeGroup == 'Elderly',]

EMEXP3567.chil <- celltypes[celltypes$accession == 'E-MEXP-3567'|celltypes$AgeGroup == 'Children',]
EMTAB4421.chil <- celltypes[celltypes$accession == 'E-MTAB-4421'|celltypes$AgeGroup == 'Children',]
gse10474.chil <- celltypes[celltypes$accession == 'GSE10474'|celltypes$AgeGroup == 'Children',]
gse21802.chil <- celltypes[celltypes$accession == 'GSE21802'|celltypes$AgeGroup == 'Children',]
gse27131.chil <- celltypes[celltypes$accession == 'GSE27131'|celltypes$AgeGroup == 'Children',]
gse32707.chil <- celltypes[celltypes$accession == 'GSE32707'|celltypes$AgeGroup == 'Children',]
gse33341.chil <- celltypes[celltypes$accession == 'GSE33341'|celltypes$AgeGroup == 'Children',]
gse40586.chil <- celltypes[celltypes$accession == 'GSE40586'|celltypes$AgeGroup == 'Children',]
gse54514.chil <- celltypes[celltypes$accession == 'GSE54514'|celltypes$AgeGroup == 'Children',]
gse63042.chil <- celltypes[celltypes$accession == 'GSE63042'|celltypes$AgeGroup == 'Children',]
gse63990.chil <- celltypes[celltypes$accession == 'GSE63990'|celltypes$AgeGroup == 'Children',]
gse66890.chil <- celltypes[celltypes$accession == 'GSE66890'|celltypes$AgeGroup == 'Children',]

# ROC curves
library(pROC)
EMEXP3567.pred = predict(model_glm, newdata = EMEXP3567.test, type = "response")
EMTAB4421.pred = predict(model_glm, newdata = EMTAB4421.test, type = "response")
gse10474.pred = predict(model_glm, newdata = gse10474.test, type = "response")
gse21802.pred = predict(model_glm, newdata = gse21802.test, type = "response")
gse27131.pred = predict(model_glm, newdata = gse27131.test, type = "response")
gse32707.pred = predict(model_glm, newdata = gse32707.test, type = "response")
gse33341.pred = predict(model_glm, newdata = gse33341.test, type = "response")
gse40586.pred = predict(model_glm, newdata = gse40586.test, type = "response")
gse54514.pred = predict(model_glm, newdata = gse54514.test, type = "response")
gse63042.pred = predict(model_glm, newdata = gse63042.test, type = "response")
gse63990.pred = predict(model_glm, newdata = gse63990.test, type = "response")
gse66890.pred = predict(model_glm, newdata = gse66890.test, type = "response")

EMEXP3567.pred1 = predict(model_glm, newdata = EMEXP3567.elder, type = "response")
EMTAB4421.pred1 = predict(model_glm, newdata = EMTAB4421.elder, type = "response")
gse10474.pred1 = predict(model_glm, newdata = gse10474.elder, type = "response")
gse21802.pred1 = predict(model_glm, newdata = gse21802.elder, type = "response")
gse27131.pred1 = predict(model_glm, newdata = gse27131.elder, type = "response")
gse32707.pred1 = predict(model_glm, newdata = gse32707.elder, type = "response")
gse33341.pred1 = predict(model_glm, newdata = gse33341.elder, type = "response")
gse40586.pred1 = predict(model_glm, newdata = gse40586.elder, type = "response")
gse54514.pred1 = predict(model_glm, newdata = gse54514.elder, type = "response")
gse63042.pred1 = predict(model_glm, newdata = gse63042.elder, type = "response")
gse63990.pred1 = predict(model_glm, newdata = gse63990.elder, type = "response")
gse66890.pred1 = predict(model_glm, newdata = gse66890.elder, type = "response")

EMEXP3567.pred2 = predict(model_glm, newdata = EMEXP3567.chil, type = "response")
EMTAB4421.pred2 = predict(model_glm, newdata = EMTAB4421.chil, type = "response")
gse10474.pred2 = predict(model_glm, newdata = gse10474.chil, type = "response")
gse21802.pred2 = predict(model_glm, newdata = gse21802.chil, type = "response")
gse27131.pred2 = predict(model_glm, newdata = gse27131.chil, type = "response")
gse32707.pred2 = predict(model_glm, newdata = gse32707.chil, type = "response")
gse33341.pred2 = predict(model_glm, newdata = gse33341.chil, type = "response")
gse40586.pred2 = predict(model_glm, newdata = gse40586.chil, type = "response")
gse54514.pred2 = predict(model_glm, newdata = gse54514.chil, type = "response")
gse63042.pred2 = predict(model_glm, newdata = gse63042.chil, type = "response")
gse63990.pred2 = predict(model_glm, newdata = gse63990.chil, type = "response")
gse66890.pred2 = predict(model_glm, newdata = gse66890.chil, type = "response")


EMEXP3567.roc = roc(EMEXP3567.test$survivor ~ EMEXP3567.pred, plot = TRUE, print.auc = TRUE)
EMTAB4421.roc = roc(EMTAB4421.test$survivor ~ EMTAB4421.pred, plot = TRUE, print.auc = TRUE)
gse10474.roc = roc(gse10474.test$survivor ~ gse10474.pred, plot = TRUE, print.auc = TRUE)
gse21802.roc = roc(gse21802.test$survivor ~ gse21802.pred, plot = TRUE, print.auc = TRUE)
gse27131.roc = roc(gse27131.test$survivor ~ gse27131.pred, plot = TRUE, print.auc = TRUE)
gse32707.roc = roc(gse32707.test$survivor ~ gse32707.pred, plot = TRUE, print.auc = TRUE)
gse33341.roc = roc(gse33341.test$survivor ~ gse33341.pred, plot = TRUE, print.auc = TRUE)
gse40586.roc = roc(gse40586.test$survivor ~ gse40586.pred, plot = TRUE, print.auc = TRUE)
gse54514.roc = roc(gse54514.test$survivor ~ gse54514.pred, plot = TRUE, print.auc = TRUE)
gse63042.roc = roc(gse63042.test$survivor ~ gse63042.pred, plot = TRUE, print.auc = TRUE)
gse63990.roc = roc(gse63990.test$survivor ~ gse63990.pred, plot = TRUE, print.auc = TRUE)
gse66890.roc = roc(gse66890.test$survivor ~ gse66890.pred, plot = TRUE, print.auc = TRUE)

EMEXP3567.roc1 = roc(EMEXP3567.elder$survivor ~ EMEXP3567.pred1)
EMTAB4421.roc1 = roc(EMTAB4421.elder$survivor ~ EMTAB4421.pred1)
gse10474.roc1 = roc(gse10474.elder$survivor ~ gse10474.pred1)
gse21802.roc1 = roc(gse21802.elder$survivor ~ gse21802.pred1)
gse27131.roc1 = roc(gse27131.elder$survivor ~ gse27131.pred1)
gse32707.roc1 = roc(gse32707.elder$survivor ~ gse32707.pred1)
gse33341.roc1 = roc(gse33341.elder$survivor ~ gse33341.pred1)
gse40586.roc1 = roc(gse40586.elder$survivor ~ gse40586.pred1)
gse54514.roc1 = roc(gse54514.elder$survivor ~ gse54514.pred1)
gse63042.roc1 = roc(gse63042.elder$survivor ~ gse63042.pred1)
gse63990.roc1 = roc(gse63990.elder$survivor ~ gse63990.pred1)
gse66890.roc1 = roc(gse66890.elder$survivor ~ gse66890.pred1)

EMEXP3567.roc2 = roc(EMEXP3567.chil$survivor ~ EMEXP3567.pred2)
EMTAB4421.roc2 = roc(EMTAB4421.chil$survivor ~ EMTAB4421.pred2)
gse10474.roc2 = roc(gse10474.chil$survivor ~ gse10474.pred2)
gse21802.roc2 = roc(gse21802.chil$survivor ~ gse21802.pred2)
gse27131.roc2 = roc(gse27131.chil$survivor ~ gse27131.pred2)
gse32707.roc2 = roc(gse32707.chil$survivor ~ gse32707.pred2)
gse33341.roc2 = roc(gse33341.chil$survivor ~ gse33341.pred2)
gse40586.roc2 = roc(gse40586.chil$survivor ~ gse40586.pred2)
gse54514.roc2 = roc(gse54514.chil$survivor ~ gse54514.pred2)
gse63042.roc2 = roc(gse63042.chil$survivor ~ gse63042.pred2)
gse63990.roc2 = roc(gse63990.chil$survivor ~ gse63990.pred2)
gse66890.roc2 = roc(gse66890.chil$survivor ~ gse66890.pred2)

str(EMEXP3567.roc)

# plot ROC curves

png('../figures/roc_AgeGroups.png', height = 900, width = 2700, res = 300)

par(mfcol = c(1,3))

plot(EMEXP3567.roc, lwd = 1, col = '#808080', main = "ROC Curves - Adults")
plot(EMTAB4421.roc, add = TRUE, lwd = 1, col = '#000000')
plot(gse10474.roc, col = '#FF0000', add = TRUE, lwd = 1)
plot(gse21802.roc, col = '#800000', add = TRUE, lwd = 1)
plot(gse27131.roc, col = '#FFFF00', add = TRUE, lwd = 1)
plot(gse32707.roc, col = '#008000', add = TRUE, lwd = 1)
plot(gse33341.roc, col = '#00FFFF', add = TRUE, lwd = 1)
plot(gse40586.roc, col = '#008080', add = TRUE, lwd = 1)
plot(gse54514.roc, col = '#0000FF', add = TRUE, lwd = 1)
plot(gse63042.roc, col = '#000080', add = TRUE, lwd = 1)
plot(gse63990.roc, col = '#FF00FF', add = TRUE, lwd = 1)
plot(gse66890.roc, col = '#800080', add = TRUE, lwd = 1)

legend("bottomright", 
       paste(c(sprintf("EMEXP3567 AUC = %0.3f", EMEXP3567.roc$auc), 
               sprintf("EMTAB4421 AUC = %0.3f", EMTAB4421.roc$auc), 
               sprintf("GSE10474 AUC = %0.3f", gse10474.roc$auc), 
               sprintf("GSE21802 AUC = %0.3f", gse21802.roc$auc),
               sprintf("GSE27131 AUC = %0.3f", gse27131.roc$auc),
               sprintf("GSE32707 AUC = %0.3f", gse32707.roc$auc),
               sprintf("GSE33341 AUC = %0.3f", gse33341.roc$auc),
               sprintf("GSE40586 AUC = %0.3f", gse40586.roc$auc),
               sprintf("GSE54514 AUC = %0.3f", gse54514.roc$auc),
               sprintf("GSE63042 AUC = %0.3f", gse63042.roc$auc),
               sprintf("GSE63990 AUC = %0.3f", gse63990.roc$auc),
               sprintf("GSE66890 AUC = %0.3f", gse66890.roc$auc))), 
       col = c("#808080", "#000000", "#FF0000", "#800000", "#FFFF00",
               "#008000", "#00FFFF", "#008080", "#0000FF", "#000080", 
               "#FF00FF", "#800080"),
       lty = 1, lwd = 2, cex = .5)

plot(EMEXP3567.roc1, lwd = 1, col = '#808080', main = "ROC Curves - Elderly")
plot(EMTAB4421.roc1, add = TRUE, lwd = 1, col = '#000000')
plot(gse10474.roc1, col = '#FF0000', add = TRUE, lwd = 1)
plot(gse21802.roc1, col = '#800000', add = TRUE, lwd = 1)
plot(gse27131.roc1, col = '#FFFF00', add = TRUE, lwd = 1)
plot(gse32707.roc1, col = '#008000', add = TRUE, lwd = 1)
plot(gse33341.roc1, col = '#00FFFF', add = TRUE, lwd = 1)
plot(gse40586.roc1, col = '#008080', add = TRUE, lwd = 1)
plot(gse54514.roc1, col = '#0000FF', add = TRUE, lwd = 1)
plot(gse63042.roc1, col = '#000080', add = TRUE, lwd = 1)
plot(gse63990.roc1, col = '#FF00FF', add = TRUE, lwd = 1)
plot(gse66890.roc1, col = '#800080', add = TRUE, lwd = 1)

legend("bottomright", 
       paste(c(sprintf("EMEXP3567 AUC = %0.3f", EMEXP3567.roc1$auc), 
               sprintf("EMTAB4421 AUC = %0.3f", EMTAB4421.roc1$auc), 
               sprintf("GSE10474 AUC = %0.3f", gse10474.roc1$auc), 
               sprintf("GSE21802 AUC = %0.3f", gse21802.roc1$auc),
               sprintf("GSE27131 AUC = %0.3f", gse27131.roc1$auc),
               sprintf("GSE32707 AUC = %0.3f", gse32707.roc1$auc),
               sprintf("GSE33341 AUC = %0.3f", gse33341.roc1$auc),
               sprintf("GSE40586 AUC = %0.3f", gse40586.roc1$auc),
               sprintf("GSE54514 AUC = %0.3f", gse54514.roc1$auc),
               sprintf("GSE63042 AUC = %0.3f", gse63042.roc1$auc),
               sprintf("GSE63990 AUC = %0.3f", gse63990.roc1$auc),
               sprintf("GSE66890 AUC = %0.3f", gse66890.roc1$auc))), 
       col = c("#808080", "#000000", "#FF0000", "#800000", "#FFFF00",
               "#008000", "#00FFFF", "#008080", "#0000FF", "#000080", 
               "#FF00FF", "#800080"),
       lty = 1, lwd = 2, cex = .5)

plot(EMEXP3567.roc2, lwd = 1, col = '#808080', main = "ROC Curves - Children")
plot(EMTAB4421.roc2, add = TRUE, lwd = 1, col = '#000000')
plot(gse10474.roc2, col = '#FF0000', add = TRUE, lwd = 1)
plot(gse21802.roc2, col = '#800000', add = TRUE, lwd = 1)
plot(gse27131.roc2, col = '#FFFF00', add = TRUE, lwd = 1)
plot(gse32707.roc2, col = '#008000', add = TRUE, lwd = 1)
plot(gse33341.roc2, col = '#00FFFF', add = TRUE, lwd = 1)
plot(gse40586.roc2, col = '#008080', add = TRUE, lwd = 1)
plot(gse54514.roc2, col = '#0000FF', add = TRUE, lwd = 1)
plot(gse63042.roc2, col = '#000080', add = TRUE, lwd = 1)
plot(gse63990.roc2, col = '#FF00FF', add = TRUE, lwd = 1)
plot(gse66890.roc2, col = '#800080', add = TRUE, lwd = 1)

legend("bottomright", 
       paste(c(sprintf("EMEXP3567 AUC = %0.3f", EMEXP3567.roc2$auc), 
               sprintf("EMTAB4421 AUC = %0.3f", EMTAB4421.roc2$auc), 
               sprintf("GSE10474 AUC = %0.3f", gse10474.roc2$auc), 
               sprintf("GSE21802 AUC = %0.3f", gse21802.roc2$auc),
               sprintf("GSE27131 AUC = %0.3f", gse27131.roc2$auc),
               sprintf("GSE32707 AUC = %0.3f", gse32707.roc2$auc),
               sprintf("GSE33341 AUC = %0.3f", gse33341.roc2$auc),
               sprintf("GSE40586 AUC = %0.3f", gse40586.roc2$auc),
               sprintf("GSE54514 AUC = %0.3f", gse54514.roc2$auc),
               sprintf("GSE63042 AUC = %0.3f", gse63042.roc2$auc),
               sprintf("GSE63990 AUC = %0.3f", gse63990.roc2$auc),
               sprintf("GSE66890 AUC = %0.3f", gse66890.roc2$auc))), 
       col = c("#808080", "#000000", "#FF0000", "#800000", "#FFFF00",
               "#008000", "#00FFFF", "#008080", "#0000FF", "#000080", 
               "#FF00FF", "#800080"),
       lty = 1, lwd = 2, cex = .5)




dev.off()


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
library(ggplot2)
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