setwd("../summary_data/*/")
df = read.csv("CellType*.tsv", sep = '\t', header = TRUE)
yes = df[df$Survivor == "YES",]
no = df[df$Survivor == "NO",]
names(df)
for (i in 3:20){
  print(names(df)[i])
  result = wilcox.test(yes[,i],no[,i])
  print(result)}
boxplot(#column names
  ~Survivor,df, xlab = "", ylab = "")

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