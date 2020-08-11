def getGeneLine(fileList, geneName):
  outDict = {}
  #outDict['placeholder']=['test','test2']

  for fileName in fileList:
        geneLine = False
        Fstream = open(fileName, 'r')
        First = True
        for line in Fstream:
           if First:
             firstline = line
             First = False
             continue
           else:
             curGene = line.strip().split('\t')[0]
             #print(curGene)
             #print(geneName)
             #print(curGene == geneName)
             if curGene == geneName:
                #print(line)
                geneLine = line
                outDict[fileName] = [firstline, geneLine]
                break
           

  return outDict

expressionFileList = [\
        '../summary_data/E-MEXP-3567/GeneSymbEMEXP3567.tsv',\
        '../summary_data/E-MTAB-4421/GeneSymb4421.tsv',\
        '../summary_data/GSE10474/GeneSymb10474.tsv',\
        '../summary_data/GSE21802/GeneSymb21802.tsv',\
        '../summary_data/GSE27131/GeneSymb27131.tsv',\
        '../summary_data/GSE32707/GeneSymb32707.tsv',\
        '../summary_data/GSE33341/GeneSymb33341.tsv',\
        '../summary_data/GSE40586/GeneSymb40586.tsv',\
        '../summary_data/GSE54514/GeneSymb54514.tsv',\
        '../summary_data/GSE63042/GeneSymb63042.tsv',\
        '../summary_data/GSE63990/GeneSymb63990.tsv',\
        '../summary_data/GSE66890/GeneSymb66890.tsv']

outDict = getGeneLine(expressionFileList, "UBA7")

for file in outDict:
     print(outDict[file][0].strip())
     print(outDict[file][1])
