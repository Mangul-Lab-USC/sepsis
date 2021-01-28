# load packages 
library(survival)
library(survminer)
install.packages("factoextra")
library(factoextra)
install.packages("survivalROC")
library(survivalROC)
library(ggfortify)
library("corrplot")

# read processed data
celltypes <- read.table("../summary_data/cell_types/CombinedLM22.tsv",
                        sep = "\t", header = TRUE)
# preview table 
head(celltypes)

# drop rows without ages
dim(celltypes)
celltypes<-celltypes[!is.na(celltypes$age),]
dim(celltypes)
# preview sex column
table(celltypes$sex)
# drop rows without sex and with undefined sex
celltypes<-celltypes[celltypes$sex %in% c(0,1),] 
dim(celltypes)

# compute PCA (for cell types)
result.PCA <- prcomp(celltypes[,9:19], center=TRUE, scale=FALSE)
# result of PCA
summary(result.PCA)
str(result.PCA)

# Eigenvalues
eig.val <- get_eigenvalue(result.PCA)
eig.val

# Visualize eigenvalues (scree plot). Show the percentage of 
# variances explained by each principal component.
fviz_eig(result.PCA, addlabels = TRUE)


# Results for Variables
res.var <- get_pca_var(result.PCA)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# visualize the contribution of variables on all the dimensions
corrplot(res.var$contrib, is.corr = FALSE)


# Graph of variables. Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(result.PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# contribution of variables to PC1-3
fviz_contrib(result.PCA, choice = "var", axes = 1:3)

# correlations between variables and PC1-3
result.PCA$rotation[,1:3]

# Biplot of individuals and variables
autoplot(result.PCA, data = celltypes, colour = "survivor", loadings = TRUE,
         loadings.colour = 'blue', 
         loadings.label = TRUE, loadings.label.size = 3.5)




# create test data list
test1 <- list(age=celltypes$age, 
              status=celltypes$survivor, 
              PC1=result.PCA$x[,1],
              PC2=result.PCA$x[,2],
              PC3=result.PCA$x[,3],
              sex=celltypes$sex,
              cohort=celltypes$accession
)

# check cohort name
table(celltypes$accession)

# save test data as a data frame
test.df<-data.frame(test1)

# statistical test between PCs and status
t.test(test.df$PC1~test.df$status)
t.test(test.df$PC2~test.df$status)
t.test(test.df$PC3~test.df$status)

# build Cox regression model
model_cox <- coxph(Surv(age, status) ~ PC1 + PC2 + PC3 + cohort + sex, 
                   test1)

summary(model_cox)

# logit model
mylogit <- glm(status ~ PC1 + PC2 + PC3 + cohort + age, data = test.df, 
               family = "binomial")
summary(mylogit)

## odds ratios and 95% CI
exp(cbind(OR = coef(mylogit), confint(mylogit)))

test.df$prob <- predict(mylogit, newdata = test.df, type = "response")
summary(test.df$prob)


# roc(test.df$status~test.df$prob,plot=TRUE,print.auc=TRUE,)


