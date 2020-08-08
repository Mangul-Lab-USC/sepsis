mydict = {}

df = open("../raw_data/gse54514gpl6947PHENO.csv","r")
for line in df:
  fields = line.strip().split(",")
  sample = fields[15]
  survivor = fields[7]
  mydict[sample] = survivor

with open("../summary_data/GSE54514/LM22prediction54514.tsv","r") as celltype:
     for line in celltype:
          fields = line.strip().split("\t")
          sample = fields[0]
          if sample in mydict:
              row = "\t".join([sample] + fields[1:])
              print(row)     
