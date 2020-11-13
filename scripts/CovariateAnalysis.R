setwd("../summary_data/*/")


#read the datasets
celltypes<-read.table('../summary_data/cell_types/CombinedLM22.tsv', 
                      sep = '\t', header = TRUE)
head(celltypes)
gene<-read.table('../summary_data/expression_data/symb_all.tsv', 
                 sep = '\t', header = TRUE)
head(gene)
dim(celltypes)
dim(gene)


# clinical score
table(celltypes$score_type)
summary(celltypes$severity_score)

# add a column of APACHE scores
library(dplyr)
# celltypes
celltypes=celltypes%>%
  mutate(APACHE.scores=ifelse(score_type %in% c("APACHE","APACHE_II"),
                              severity_score,NA))

table(celltypes$APACHE.scores)
summary(celltypes$APACHE.scores)
# gene expression
gene=gene%>%
  mutate(APACHE.scores=ifelse(score_type %in% c("APACHE","APACHE_II"),severity_score,NA))
summary(gene$APACHE.scores)

# add an age group column
celltypes$AgeGroup <- cut(celltypes$age, breaks = c(0, 18, 65, Inf), 
                          labels = c('Children', 'Adults', 'Elderly'))

# data distribution
install.packages("ggpubr")
library(ggplot2)
library(ggpubr)

score <- ggboxplot(celltypes, x = "survivor", y = "APACHE.scores",
                  order = c('1', '0'),
                  color = "#2E86C1", add = "jitter",
                  # add.params = list(size=.2),
                  bxp.errorbar = TRUE, xlab = FALSE,
                  ylab = "Clinical Severity Scores") + 
  font("ylab", size = 16) +
  font("xy.text", size = 16)+
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .9,
                     size = 5)

ggsave("scores_boxplot.png", plot = score, path = "../figures/", dpi = 300,)


# there are many outliers, the data doesn't follow normal distribution
# so we'll use wilcoxon test to see if there is statistically significant
# difference between mortality and clinical score

# statistical analysis
survivor = celltypes[celltypes$survivor==1,]
death = celltypes[celltypes$survivor==0,]
names(celltypes)

for (i in 8){
  result = wilcox.test(survivor[,i],death[,i])
  print(result)
}

# split training and test data
table(celltypes$survivor)

trainingData <- celltypes[celltypes$use == 'discovery',]
trainingAdult <- celltypes[celltypes$AgeGroup == 'Adults',]
trainingElder <- celltypes[celltypes$AgeGroup == 'Elderly',]


# logistic regression model for clinical score
model.apache <- glm(survivor~APACHE.scores,trainingData,family = 'binomial')
summary(model.apache)


### plot roc curves for clinical score
# create a empty list, full accessions and corresponding colors
list<-list()
full_name<-c('E-MEXP-3567','E-MTAB-4421','GSE10474','GSE21802','GSE27131',
             'GSE32707','GSE33341','GSE40586','GSE54514','GSE63042',
             'GSE63990','GSE66890')
col<-c("#808080","#000000", "#FF0000", "#800000", "#FFFF00",
       "#008000", "#00FFFF", "#008080", "#0000FF", "#000080",
       "#FF00FF", "#800080")

# combine accessions that have APACHE scores
name<-c('E-MTAB-4421','GSE10474','GSE32707','GSE54514','GSE66890')

### plot
library(pROC)

png('../figures/roc_ClinicalScore.png', height = 1800, width = 1800, res = 300)

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i],]
  pred<-predict(model.apache, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves: Predicting Mortality Using Clinical Scores',
            lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .9)
dev.off()






# pull out columns that are statistically significant with survivor

# cell types

for (i in c(5:6,9:35)){
  print(names(celltypes)[i])
  result = wilcox.test(survivor[,i],death[,i])
  print(result)}


# gene expression
for (i in 9:67){
  print(names(gene)[i])
  result = wilcox.test(survivor[,i],death[,i])
  print(result)}
library(limma)


#boxplot for statistically significant columns
Ages <- ggboxplot(celltypes, x = "survivor", y = "age",
                  order = c('1', '0'),
                  color = "#2E86C1", add = "jitter",
                  add.params = list(size=.2),
                  bxp.errorbar = TRUE, xlab = FALSE,
                  ylab = "Age") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .1)

Tcell <- ggboxplot(celltypes, x = "survivor", y = "T_cells",
                   order = c('1', '0'),
                   color = "#2E86C1", add = "jitter",
                   add.params = list(size=.2),
                   bxp.errorbar = TRUE, xlab = FALSE,
                   ylab = "T cells") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .9)

CD4Tcell <- ggboxplot(celltypes, x = "survivor", y = "T_cells_CD4",
               order = c('1', '0'),
               color = "#2E86C1", add = "jitter",
               add.params = list(size=.2),
               bxp.errorbar = TRUE, xlab = FALSE,
               ylab = "T cells (CD4+)") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .9)

macrophages <- ggboxplot(celltypes, x = "survivor", y = "Macrophages",
                      order = c('1', '0'),
                      color = "#2E86C1", add = "jitter",
                      add.params = list(size=.2),
                      bxp.errorbar = TRUE, xlab = FALSE,
                      ylab = "Macrophages") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .9)

mastcells <- ggboxplot(celltypes, x = "survivor", y = "Mast_cells_resting",
                      order = c('1', '0'),
                      color = "#2E86C1", add = "jitter",
                      add.params = list(size=.2),
                      bxp.errorbar = TRUE, xlab = FALSE,
                      ylab = "Mast cells") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = 0.5,label.y.npc = 0.9)

neutrophils <- ggboxplot(celltypes, x = "survivor", y = "Neutrophils",
                         order = c('1', '0'),
                         color = "#2E86C1", add = "jitter", 
                         add.params = list(size=.2),
                         bxp.errorbar = TRUE, xlab = FALSE,
                         ylab = "Neutrophils") + 
  scale_x_discrete(labels = c('1' = "Survived", '0' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .9)
# multipanel figure
install.packages('patchwork')
library(patchwork)

figure1<-Ages+Tcell+CD4Tcell+neutrophils+macrophages+mastcells

ggsave("significant_cells.png", plot = figure1, path = "../figures/", dpi = 300)

################
################
# logistic regression model on cell types that are statistically significant 
# associated with mortality
model_glm <- glm(survivor~Mast_cells_resting+Macrophages+T_cells_CD4, 
                 data = trainingData, family = 'binomial')

# overall prediction
par(mfcol = c(3,2))
# standard predictors
test <- gene[gene$score_type %in% c("APACHE","APACHE_II"),]

pred <- predict(model.apache, newdata = test, type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc =TRUE,
           main="ROC Curves: Clinical Severity Scores")
# severity+celltypes
pred <- predict(model_glm, newdata = test, type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves: Clinical Scores and Cell Types Composition")
# Adult
test <- celltypes[celltypes$AgeGroup == "Adults",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Adults")
# Elderly
test <- celltypes[celltypes$AgeGroup == "Elderly",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Elderly")
# Children
test <- celltypes[celltypes$AgeGroup == "Children",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Children")
# All
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(gene$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for All Ages")

summary(model_glm)
# This will make predictions on the training data that you use to fit the model 
# and give me a vector of fitted probabilities.
glm.prob <- predict(model_glm,trainingData,type = 'response')
glm.prob[1:5]

# turn the probabilities into classifications by thresholding at 0.5
glm.pred <- ifelse(glm.prob > 0.5, 1, 0)
attach(trainingData)

table(glm.pred,trainingData$survivor)
mean(glm.pred == trainingData$survivor)

# Concordance
# higher the concordance, the better is the quality of model
Concordance(trainingData$survivor, glm.prob)

# sensitivity and specificity
sensitivity(trainingData$survivor, glm.prob)
specificity(trainingData$survivor, glm.prob)

# Misclassification Error
# percentage mismatch of predcited vs actuals, irrespective of 1’s or 0’s. 
# The lower the misclassification error, the better is your model.
misClassError(trainingData$survivor,glm.prob)

# plot roc curves for training and test data
png('../figures/roc_train.test.png', height = 2300, width = 3300, res = 300)

par(mfrow = c(2,3))
# training
# cell types
list <- list()
train <- c("E-MEXP-3567","GSE10474","GSE27131","GSE32707","GSE40586","GSE63042",
           "GSE66890")
legend_list<-c()
col_plot<-c()
for (i in 1:length(train)){
  test<-celltypes[celltypes$accession == train[i],]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(train[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Cell types',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(train[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)

# severity alone
model.apache <- glm(survivor~APACHE.scores,trainingData,family = 'binomial')
list<-list()
train<-c('GSE10474','GSE32707','GSE66890')
legend_list<-c()
col_plot<-c()
for (i in 1:length(train)){
  test<-celltypes[celltypes$accession == train[i],]
  pred<-predict(model.apache, newdata = test, type = "response")
  k=match(train[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Clinical scores alone',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(train[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)

# severity+cell types
model_joint <- glm(survivor~Mast_cells_resting+Macrophages+T_cells_CD4+APACHE.scores, 
                 data = trainingData, family = 'binomial')
summary(model_joint)
list<-list()
legend_list<-c()
col_plot<-c()
for (i in 1:length(train)){
  test<-celltypes[celltypes$accession == train[i],]
  pred<-predict(model_joint, newdata = test, type = "response")
  k=match(train[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Cell types + clinical scores',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(train[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)

# test
# cell types
list <- list()
testing <- c("E-MTAB-4421","GSE21802","GSE33341","GSE54514","GSE63990")
legend_list<-c()
col_plot<-c()
for (i in 1:length(testing)){
  test<-celltypes[celltypes$accession == testing[i],]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(testing[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Cell types',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(testing[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)

# severity alone
list <- list()
testing <- c("E-MTAB-4421", "GSE54514")
legend_list<-c()
col_plot<-c()
for (i in 1:length(testing)){
  test<-celltypes[celltypes$accession == testing[i],]
  pred<-predict(model.apache, newdata = test, type = "response")
  k=match(testing[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Clinical scores alone',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(testing[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)

# severity+cell types
list <- list()
testing <- c("E-MTAB-4421", "GSE54514")
legend_list<-c()
col_plot<-c()
for (i in 1:length(testing)){
  test<-celltypes[celltypes$accession == testing[i],]
  pred<-predict(model_joint, newdata = test, type = "response")
  k=match(testing[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'Cell types + clinical scores',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(testing[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .8)
dev.off()

# open a png file
png('../figures/roc_AgeGroups.png', height = 2300, width = 2300, res = 300)
# set ncols and nrows
par(mfcol = c(2,2))

### plot roc curves for adults
list<-list()

name<-c('E-MTAB-4421','GSE10474','GSE27131','GSE32707','GSE33341','GSE40586',
        'GSE54514','GSE66890')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i]&celltypes$AgeGroup == 'Adults',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Adults',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for elderly
list<-list()

name<-c('E-MTAB-4421','GSE10474','GSE32707','GSE33341','GSE54514','GSE66890')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i]&celltypes$AgeGroup == 'Elderly',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Elderly',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for children
list<-list()

name<-c('E-MEXP-3567')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i]&celltypes$AgeGroup == 'Children',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Children',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for all ages
list<-list()

name<-full_name

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i],]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for All Ages',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

dev.off()


str(roc)

###############################
###############################
# logistic regression model on cell types, which are statistically significant 
# associated with mortality, and clinical scores
model_glm <- glm(survivor~Mast_cells_resting+Macrophages+T_cells_CD4+APACHE.scores, 
                 data = trainingData, family = 'binomial')
summary(model_glm)
# This will make predictions on the training data that you use to fit the model 
# and give me a vector of fitted probabilities.
glm.prob <- predict(model_glm,trainingData,type = 'response')
glm.prob[1:5]

# turn the probabilities into classifications by thresholding at 0.5
glm.pred <- ifelse(glm.prob > 0.5, 1, 0)
attach(trainingData)

table(glm.pred,trainingData$survivor)
mean(glm.pred == trainingData$survivor)

# Concordance
# higher the concordance, the better is the quality of model
Concordance(trainingData$survivor, glm.prob)

# sensitivity and specificity
sensitivity(trainingData$survivor, glm.prob)
specificity(trainingData$survivor, glm.prob)

# Misclassification Error
# percentage mismatch of predcited vs actuals, irrespective of 1’s or 0’s. 
# The lower the misclassification error, the better is your model.
misClassError(trainingData$survivor,glm.prob)

# open a png file
png('../figures/roc_cell+score.png', height = 1800, width = 1800, res = 300)

### plot roc curves
list<-list()

name<-c('E-MTAB-4421','GSE10474','GSE32707','GSE54514','GSE66890')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-celltypes[celltypes$accession == name[i],]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves: Clinical Scores and Cell Types Data',
            lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .9)

dev.off()

#####################################################
gene$Mast_cells_resting <- celltypes$Mast_cells_resting
gene$Macrophages <- celltypes$Macrophages
gene$T_cells_CD4 <- celltypes$T_cells_CD4
# add an age group column
gene$AgeGroup <- cut(gene$age, breaks = c(0, 18, 65, Inf), 
                     labels = c('Children', 'Adults', 'Elderly'))

# logistic regression model on gene expression that are statistically significant 
# associated with mortality
trainingData.gene <- gene[gene$use == 'discovery',]
model.apache <- glm(survivor~APACHE.scores, data=trainingData.gene,family="binomial")
model_glm <- glm(survivor~CFD+DDIT4+DEFA4+IFI27+IL1R2+MAFF+AIM2+APH1A+CCR2+EIF5A+RAB40B+VNN3+Mast_cells_resting+Macrophages+T_cells_CD4+APACHE.scores, 
                 data = trainingData.gene, family = 'binomial')
summary(model_glm)

# standard predictors
test <- gene[gene$score_type %in% c("APACHE","APACHE_II"),]

pred <- predict(model.apache, newdata = test, type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc =TRUE,
           main="ROC Curves: Clinical Severity Scores")
# severity+celltypes+gene
pred <- predict(model_glm, newdata = test, type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves: Clinical Scores + Cell Types + Gene")
# Adult
test <- gene[gene$AgeGroup == "Adults",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Adults")
# Elderly
test <- gene[gene$AgeGroup == "Elderly",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Elderly")
# Children
test <- gene[gene$AgeGroup == "Children",]
pred <- predict(model_glm,newdata = test,type = "response")
roc <- roc(test$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for Children")
# All
pred <- predict(model_glm,newdata = gene,type = "response")
roc <- roc(gene$survivor~pred,plot=TRUE,print.auc=TRUE,
           main="ROC Curves for All Ages")

# open a png file
png('../figures/roc_AgeGroups.png', height = 2300, width = 2300, res = 300)
# set ncols and nrows
par(mfcol = c(2,2))
full_name<-c('EMEXP3567','E-MTAB-4421','GSE10474','GSE21802','GSE27131',
             'GSE32707','GSE33341','GSE40586','GSE54514','GSE63042',
             'GSE63990','GSE66890')
### plot roc curves for adults
list<-list()

name<-c('E-MTAB-4421','GSE10474','GSE27131','GSE32707','GSE33341','GSE40586',
        'GSE54514','GSE66890')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-gene[gene$accession == name[i]&gene$AgeGroup == 'Adults',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Adults',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for elderly
list<-list()

name<-c('E-MTAB-4421','GSE10474','GSE32707','GSE33341','GSE54514','GSE66890')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-gene[gene$accession == name[i]&gene$AgeGroup == 'Elderly',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Elderly',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for children
list<-list()

name<-c('EMEXP3567')

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-gene[gene$accession == name[i]&gene$AgeGroup == 'Children',]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for Children',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

### plot roc curves for all ages
list<-list()

name<-full_name

legend_list<-c()
col_plot<-c()
for (i in 1:length(name)){
  test<-gene[gene$accession == name[i],]
  pred<-predict(model_glm, newdata = test, type = "response")
  k=match(name[i],full_name)
  col_plot<-c(col_plot,col[k])
  if (i==1)
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,
            main = 'ROC Curves for All Ages',lwd = 1, col = col[k])
  }
  else
  {
    roc=roc(test$survivor ~ pred, plot = TRUE,col = col[k], add = TRUE,
            lwd = 1)
  }
  legend_list<-paste(c(legend_list,paste(name[i],sprintf("AUC = %0.3f", roc$auc))))
  
}
legend("bottomright", 
       legend_list, 
       col = col_plot,
       lty = 1, lwd = 2, cex = .5)

dev.off()

