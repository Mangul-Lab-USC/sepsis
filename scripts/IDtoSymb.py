cDict = {}
#print(df.columns)

fstream = open("../summary_data/ConversionTable/affyhg1.0ST.tsv","r")
for line in fstream:
    fields = line.strip().split("\t")
    probeID = fields[1]
    geneSymb = fields[4]
    if fields[3] == "protein_coding" and probeID not in cDict:
              cDict[probeID] = geneSymb
#print(cDict)
with open("../summary_data/GSE27131/GSE27131.tsv","r") as data:
     for line in data:
          fields = line.strip().split("\t")
          probeID = fields[0]
          if probeID in cDict:
              geneSymb = cDict[probeID]
              #print(geneSymb)
              new = "\t".join([geneSymb] + fields[1:])
              print(new)
