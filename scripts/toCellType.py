mydict = {}

df = open("../raw_data/E-MTAB-4421/EMTAB4421.51PHENO.csv","r")
for line in df:
  fields = line.strip().split(",")
  sample = fields[0]
  survivor = fields[4]
  mydict[sample] = survivor

with open("../summary_data/E-MTAB-4421/GEDITPredictions4421.tsv","r") as celltype:
     for line in celltype:
          fields = line.strip().split("\t")
          sample = fields[0]
          if sample in mydict:
              row = "\t".join([sample] + fields[1:])
              print(row)     
