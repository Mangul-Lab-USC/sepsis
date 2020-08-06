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