mydict = {}

df = open("../raw_data/gse63042gpl9115PHENO.csv","r", encoding="utf8", errors='ignore')
for line in df:
  fields = line.strip().split(",")
  sample = fields[0]
  survivor = fields[28]
  mydict[sample] = survivor

with open("../summary_data/symb.tsv","r") as gene:
     for line in gene:
          fields = line.strip().split("\t")
          sample = fields[0]
          if sample in mydict:
              survivor = mydict[sample]
              row = "\t".join([sample] + [survivor] + fields[1:])
              print(row)     
