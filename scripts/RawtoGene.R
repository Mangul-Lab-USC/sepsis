# First, we want to convert CEL files which are downloaded from GEO
# to either csv or tsv files, we do this using Rstudio

setwd('') # the path of working directory where the CEL files are

library(affy)
data = ReadAffy()
normalized = rma(data)
write.exprs(normalized, file = "OutputFileName.tsv", sep =  "\t")

# On GEO accession website, we can download the conversion table
# through "Platforms" to convert probe id to gene symbols
# I use a python script to do this for individual cohort 
# ('IDtoSymb.py'), you can find the script under "scripts" folder in the 
# GitHub repository I sent to you days ago, then you will get
# gene expression data for this single cohort