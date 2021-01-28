# load packages
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(lubridate) # extract date (year,month ....)
library(tidyr) # unite function

#read the datasets
celltypes<-read.table('../summary_data/cell_types/CombinedHPCA.tsv', 
                      sep = '\t', header = TRUE)
head(celltypes)
dim(celltypes)

summary(celltypes$age)

# statistical analysis
t.test(celltypes$Platelets~celltypes$Survivor)

survivor = celltypes[celltypes$Survivor==0,]
death = celltypes[celltypes$Survivor==1,]
names(celltypes)

for (i in 9:26){
  print(names(celltypes)[i])
  result = wilcox.test(survivor[,i],death[,i])
  print(result)
}

# data distribution
platelets <- ggboxplot(celltypes, x = "Survivor", y = "Platelets",
                  order = c('0', '1'),
                  color = "#2E86C1", add = "jitter",
                  add.params = list(size=.2),
                  bxp.errorbar = TRUE, xlab = FALSE,
                  ylab = "Platelet abundances") + 
  scale_x_discrete(labels = c('0' = "Survived", '1' = "Died")) +
  stat_compare_means(label = "p.format", label.x.npc = .5, label.y.npc = .8,
                     size = 5)

# grouped by age&sex
attach(celltypes)
celltypes$age_sex[age <= 35 & sex == 0] <- "Young men"
celltypes$age_sex[age > 35 & sex == 0] <- "Old men"
celltypes$age_sex[age <= 35 & sex == 1] <- "Young women"
celltypes$age_sex[age > 35 & sex == 1] <- "Old women"
# drop NA
celltypes<-celltypes[!is.na(celltypes$age_sex),]
dim(celltypes)
table(celltypes$age_sex)
# statistical analysis

# create boxplot
ggplot(celltypes, aes(x=factor(Survivor), y=Platelets, fill = age_sex)) +
  geom_boxplot() +
  stat_compare_means() +
  ylab(label = "Platelet abundances") +
  xlab(label = "") +
  scale_x_discrete(labels = c('0' = "Survivor", '1' = "Non-survivor"))

