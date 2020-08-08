cDict = {}
#print(df.columns)

fstream = open("../summary_data/ConversionTable/IlluminaHT12v3.0.tsv","r")
for line in fstream:
    fields = line.strip().split("\t")
    probeID = fields[0]
    geneSymb = fields[13]
    if probeID not in cDict:
    	cDict[probeID] = geneSymb
#print(cDict)
with open("../summary_data/GSE54514/GSE54514.tsv","r") as data:
     for line in data:
          fields = line.strip().split("\t")
          probeID = fields[0]
          if probeID in cDict:
              geneSymb = cDict[probeID]
              #print(geneSymb)
              new = "\t".join([geneSymb] + fields[1:])
              print(new)
