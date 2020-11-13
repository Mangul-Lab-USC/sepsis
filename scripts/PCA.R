# packages used
library(survival)
library(survminer)
# install.packages("factoextra")
library(factoextra)
install.packages("survivalROC")
library(survivalROC)
library(ggfortify)

# read processed data
celltypes <- read.table("../summary_data/cell_types/CombinedLM22.tsv",
                        sep = "\t", header = TRUE)
# preview
head(celltypes)
dim(celltypes)

# PCA
result.PCA <- prcomp(celltypes[,9:19], center=TRUE, scale=FALSE)
# result of PCA
summary(result.PCA)
str(result.PCA)

# Visualize eigenvalues (scree plot). Show the percentage of 
# variances explained by each principal component.
fviz_eig(result.PCA)

# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(result.PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)




autoplot(result.PCA, data = celltypes, colour = "survivor", loadings = TRUE,
         loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3.5)



# Cox proportional hazards regression model
dim(celltypes)
celltypes<-celltypes[!is.na(celltypes$age),]
dim(celltypes)
celltypes<-celltypes[celltypes$sex %in% c(0,1),]
dim(celltypes)


test1 <- list(age=celltypes$age, 
              status=celltypes$survivor, 
              x=result.PCA$x[as.integer(rownames(celltypes)),'PC1'], 
              sex=celltypes$sex,
              cohort=celltypes$accession
)

table(celltypes$accession)
junk<-data.frame(test1)
junk[junk$cohort=='GSE66890',]

t.test(junk$x~junk$status)

model_cox <- coxph(Surv(time=age, event=status) ~ x + cohort + strata(sex), test1)
summary(model_cox)
plot(cox.zph(model_cox))
ggsurvplot(survfit(model_cox), col = "#2E9FDF", data = test1, xlab="Age",
           ggtheme = theme_minimal())
# cox_fit <- survfit(model_cox)
# autoplot(cox_fit)
# aa_fit <-aareg(Surv(age, status) ~ x + cohort, test1)
# summary(aa_fit)
# autoplot(aa_fit)

fun.survivalROC <- function(lp, t) {
  res <- with(celltypes,
              survivalROC(Stime        = age,
                          status       = survivor,
                          marker       = get(lp),
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  
  ## Plot ROCs
  with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC)))
  abline(a = 0, b = 1, lty = 2)
  
  res
}
fun.survivalROC(lp = "model_cox", t)
