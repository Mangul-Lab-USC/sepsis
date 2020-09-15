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
        '../summary_data/expression_data/GeneSymbEMEXP3567.tsv',\
        '../summary_data/expression_data/GeneSymb4421.tsv',\
        '../summary_data/expression_data/GeneSymb10474.tsv',\
        '../summary_data/expression_data/GeneSymb21802.tsv',\
        '../summary_data/expression_data/GeneSymb27131.tsv',\
        '../summary_data/expression_data/GeneSymb32707.tsv',\
        '../summary_data/expression_data/GeneSymb33341.tsv',\
        '../summary_data/expression_data/GeneSymb40586.tsv',\
        '../summary_data/expression_data/GeneSymb54514.tsv',\
        '../summary_data/expression_data/GeneSymb63042.tsv',\
        '../summary_data/expression_data/GeneSymb63990.tsv',\
        '../summary_data/expression_data/GeneSymb66890.tsv']

outDict = getGeneLine(expressionFileList, "CEACAM8")

for file in outDict:
     print(outDict[file][0].strip())
     print(outDict[file][1])
