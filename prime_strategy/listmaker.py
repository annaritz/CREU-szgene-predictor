import itertools
import time
import math
import statistics
import matplotlib.pyplot as plt
import astropy.stats as ast
import numpy as np
badGenes=[]
allfoldchanges=[]
allps=[]
allEffectConfidences=[]

def main():
    print('Making GeneMap')
    GeneMapfile=open('Homo_sapiens.txt','r')
    GeneMap, allgenenames=createGeneMap(GeneMapfile)

    print('Reading 108 Loci')
    LociFile=open('108Loci.csv', 'r')
    LociList=LociGeneLister(LociFile,GeneMap, allgenenames)

    print('Reading DeadExpression')
    DeadFile=open('DeadExpression.csv','r')
    DeadExpressionList=DeadGeneLister(DeadFile,GeneMap, allgenenames)

    # print('Reading hiPSC')
    # hiPSCFile=open('hiPSC.csv', 'r')
    # hiPSCList=hiPSCLister(hiPSCFile, GeneMap, allgenenames)

    print('Reading szgene')
    szgeneFile=open('szgene.csv', 'r')
    szgeneList=szgeneLister(szgeneFile, GeneMap, allgenenames)

    # print('Reading chromoconfo')
    # chromoconfoFile=open('chromoconfo.csv', 'r')
    # chromoconfoList=chromoLister(chromoconfoFile, GeneMap, allgenenames)




    BigList=LociList+DeadExpressionList+szgeneList
    collectionList=[LociList,DeadExpressionList,szgeneList]

    sortExpressionList=sorted(BigList, key=lambda x: x[2])
    # for i in sortExpressionList:
    #     print(i)

    foldDistribution=[]
    pDistribution=[]
    effectConfidenceDistribution=[]

    

    for collection in collectionList:
        newfoldList=[]
        newpList=[]
        newECList=[]
        for gene in collection:
            if gene[3]<1:

                print('nrlanfnalfnldsnfanfljnsdlnf')
                print(collectionList.index(collection))
                print(gene)
                print('nrlanfnalfnldsnfanfljnsdlnf')


            newfoldList.append(gene[3])
            newpList.append(gene[4])
            newECList.append(gene[2])
        foldDistribution.append(newfoldList)
        pDistribution.append(newpList)
        effectConfidenceDistribution.append(newECList)

    
    deleteList=[]
    print(len(BigList))
    overlaps=0
    already_done=[]
    for i in range(len(BigList)):
        for j in range(len(BigList)):
            if BigList[i][1]==BigList[j][1]:
                if i!=j and (i not in already_done or j not in already_done):
                    already_done.append(i)
                    already_done.append(j)
                    if BigList[i][4]>BigList[j][4]:
                        deleteList.append(i)
                        overlaps=overlaps+1
                        BigList[j][4]=BigList[i][4]*BigList[j][4]
                    if BigList[i][4]<BigList[j][4]:
                        overlaps=overlaps+1
                        deleteList.append(j)
                        BigList[i][4]=BigList[i][4]*BigList[j][4]
                    BigList[i][2]=BigList[i][2]+BigList[j][2]
                    BigList[j][2]=BigList[i][2]
    for i in range(len(deleteList)):
        BigList.pop(deleteList[i])
        for j in range(len(deleteList)):
            if deleteList[j]>deleteList[i]:
                deleteList[j]=deleteList[j]-1

    BigList=sorted(BigList, key=lambda x: x[4])
    outputFile1=open('SZ_positives.txt', 'w')
    # outputFile2=open('E2_positives.txt', 'w')
    # outputFile3=open('E3_positives.txt', 'w')
    overlapControlledECdistribution=[]

    E1=[]
    E2=[]
    E3=[]
    Actual_positives=[]

    for gene in BigList:
        if gene[4]<0.05:
            Actual_positives.append(gene)

    for gene in Actual_positives:
        outputFile1.write(str(gene[1])+'\n')
    #     if Actual_positives.index(gene)<len(Actual_positives)//3:
    #         E1.append(gene)
    #     elif Actual_positives.index(gene)<(len(Actual_positives)//3)*2:
    #         E2.append(gene)
    #     else:
    #         E3.append(gene)

    # for gene in E1:
    #     outputFile1.write(str(gene[1])+'\t'+str(gene[4])+'\n')
    # for gene in E2:
    #     outputFile1.write(str(gene[1])+'\t'+str(gene[4])+'\n')
    # for gene in E3:
    #     outputFile3.write(str(gene[1])+'\t'+str(gene[4])+'\n')









        overlapControlledECdistribution.append(gene[2])

    collectionNames=['108 Loci','Common Mind Consortium','SZGene Metaanalysis']
    print('foldDistribution optimal bins:', max(allfoldchanges)/ast.knuth_bin_width(allfoldchanges))
    print('allps optimal bins:', max(allps)/ast.knuth_bin_width(allps))
    print('EC optimal bins:', max(allEffectConfidences)/ast.knuth_bin_width(allEffectConfidences))

    plt.hist(foldDistribution, bins=143,histtype='barstacked', label=collectionNames)
    plt.legend(prop={'size': 10})
    plt.title('all fold-changes')
    plt.show()

    plt.hist(pDistribution, bins=np.logspace(np.log10(0.00000000000001),np.log10(0.05)), stacked=True, label=collectionNames)
    plt.gca().set_xscale("log")
    plt.legend(prop={'size': 10})
    plt.title('all ps')
    plt.show()

    plt.hist(effectConfidenceDistribution, bins=35, stacked=True, label=collectionNames)
    plt.legend(prop={'size': 10})
    plt.title('all Effect Confidences')
    plt.show()

    plt.hist(overlapControlledECdistribution, bins=35, stacked=True, label=collectionNames)
    plt.title('overlap controlled Effect Confidences')
    plt.show()


    return

class Gene:
    def __init__(self, name, aliases=None, entrez=None, ens=None,locustype=None):
        self.name=name
        self.aliases=aliases
        self.entrez=entrez
        self.ens=ens
        self.locustype=locustype
        self.effectconfidence=0

class ProtoGene:
    def __init__(self, name=None, entrez=None, ens=None, effectConfidence=None):
        self.name=name
        self.entrez=entrez
        self.ens=ens
        self.effectConfidence=effectConfidence     

class Counter:
    def __init__(self, count,badGenes):
        self.count=0
        self.genesChanged=[]
        self.badGenes=[]

def calculateEffectConfidence(foldchange, p):
    allfoldchanges.append(foldchange)
    allps.append(math.log(p)*-1)
    # x=(foldchange**2)*math.log(p)*-1
    # x=foldchange**(math.log(p)*-1)
    # x=foldchange*math.log(p)*-1
    x=math.log(foldchange, 2)*math.log(p,2)*-1
    # x=(foldchange**2)*(math.log(p)*-1)**2
    allEffectConfidences.append(x)
    return x



def createGeneMap(GeneMapfile):
    GeneMap=[]
    allgenenames=set()
    start=time.time()

    j=0

    for i in GeneMapfile:
        i=i.split('\t')
        if i[0]=='9606' and i[9]!='biological-region' and i[9]!='tRNA' and i[9]!='rRNA':
            ientrez=i[1]
            iname=i[2]
            ialiases=set()
            if iname!=i[10] and i[10]!= '-':
                ialiases.add(i[10])
            if i[4]!='-':
                if '|' in i[4]:
                    synonyms=set(i[4].split('|'))
                    for item in synonyms:
                        ialiases.add(item)
                else:
                    ialiases.add(i[4])
            database=i[5].split('|')

            iens=None
            for identifier in database:
                if 'Ensembl:' in identifier:
                    iens= identifier[8:]
            locustype=i[9]
            igene=Gene(iname,ialiases,ientrez,iens, locustype)
            GeneMap.append(igene)
            allgenenames.add(iname)
    return GeneMap, allgenenames






def LociAliasDestroy(GeneList, GeneMap, allgenenames):
    for h in range(len(GeneList)):
        found=False

        for i in GeneMap:

            if GeneList[h][0] in i.aliases:

                GeneList[h][0]=i.name
                found=True

            if GeneList[h][0] == i.name:
                GeneList[h]=[GeneList[h][0], i.entrez, GeneList[h][1], GeneList[h][2], GeneList[h][3]]
                found=True


        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])
            print(badGenes)

    return GeneList



    # Loci=open('108Loci.csv','r')
    # chromoconfo=open('chromoconfo.csv','r')
    # DeadExpression=open('DeadExpression.csv','r')
    # hiPSC=open('hiPSC.csv','r')
    # szgene=open('szgene.csv','r')

def csvLister(csv):
    csvList=[]
    for line in csv:
        line=line.split(',')
        csvList.append(line)
    return csvList



def LociGeneLister(LociFile,GeneMap, allgenenames):
    rawList= csvLister(LociFile)
    GeneList=[]
    for i in range(len(rawList)):
        #We need a measure of confidence as to how much effect it has on SZ. 
        #It's as simple as the proportional frequency change between SZ
        # and control squared times the p value.
        pList1=rawList[i][11].split('E')
        p1=float(pList1[0])**float(pList1[1])
        pList2=rawList[i][12].split('E')
        p2=float(pList2[0])**float(pList2[1])
        p=(p1+p2)/2
        foldchange=max(float(rawList[i][8])/float(rawList[i][9]),float(rawList[i][9])/float(rawList[i][8]))
        effectConfidence=calculateEffectConfidence(foldchange,p)
        GeneList.append([rawList[i][3],effectConfidence, foldchange, p])
    GeneList=LociAliasDestroy(GeneList, GeneMap, allgenenames)

    return GeneList


    # Loci=csvLister(Loci)
    # Loci=LociGeneLister(Loci)

def DeadAliasDestroy(GeneList, GeneMap, allgenenames):
    for h in range(len(GeneList)):
        found=False
        effectConfidence=GeneList[h][1]
        foldchange=GeneList[h][2]
        p=GeneList[h][3]
        

        for i in GeneMap:

            if GeneList[h][0] in i.aliases:

                GeneList[h][0]=i.name
                found=True

            if GeneList[h][0] == i.name:
                GeneList[h]=[GeneList[h][0], i.entrez, effectConfidence, foldchange, p]
        
                found=True

        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])



    return GeneList

def DeadGeneLister(deadFile,GeneMap, allgenenames):
    rawList=csvLister(deadFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i!=0:

            FClist=rawList[i][2].split('E')
            logFC=float(FClist[0])**float(FClist[1])
            pList=rawList[i][7].split('E')
            p=float(pList[0])**float(pList[1])
            foldchange= 2.0**abs(logFC)
            if foldchange<=0:
                foldchange=foldchange*-1
            if foldchange<1:
                foldchange=1/foldchange
            effectConfidence=calculateEffectConfidence(foldchange, p)


            if rawList[i][1] not in allgenenames or rawList[i][1]=='.':
                
                for gene in GeneMap:
                    if gene.ens==rawList[i][0]:
                        # print
                        rawList[i][1]=gene.name

            if rawList[i][1]!='.':
                GeneList.append([rawList[i][1], effectConfidence, foldchange, p])


    GeneList=DeadAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList


def hiPSCAliasDestroy(GeneList, GeneMap, allgenenames):
    newGeneList=[]
    for h in range(len(GeneList)):
        found=False
        effectConfidence=GeneList[h][1]
        foldchange=GeneList[h][2]
        p=GeneList[h][3]
        

        for i in GeneMap:


            if GeneList[h][0] in i.aliases:

                GeneList[h][0]=i.name
                found=True

            if GeneList[h][0] == i.name:

                newGeneList.append([GeneList[h][0], i.entrez, effectConfidence, foldchange, p])
        
                found=True

        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])



    return newGeneList


def hiPSCLister(hiPSCFile,GeneMap, allgenenames):
    rawList=csvLister(hiPSCFile)
    GeneList=[]
    for i in range(len(rawList)):



        if rawList[i][0]!='' and i>1:
            foldchange=abs(float(rawList[i][2]))
            p=float(rawList[i][3])
            if p==0:
                p=p+0.001
            effectconfidence=calculateEffectConfidence(foldchange,p)
            GeneList.append([rawList[i][0],effectconfidence, foldchange, p])

    GeneList.pop(0)
    for i in range(len(GeneList)):
        try:
            


            if GeneList[i][0][:3]=='LOC':
                GeneList.remove(GeneList[i])
        except:

            pass
    GeneList=hiPSCAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList


    # hiPSC=hiPSCLister(hiPSC)

    # szgene=csvLister(szgene)

def szgeneAliasDestroy(GeneList, GeneMap, allgenenames):
    newGeneList=[]
    for h in range(len(GeneList)):
        found=False
        effectConfidence=GeneList[h][1]
        foldchange=GeneList[h][2]
        p=GeneList[h][3]
        for i in GeneMap:
            if GeneList[h][0] in i.aliases:
                GeneList[h][0]=i.name
                found=True
            if GeneList[h][0] == i.name:
                newGeneList.append([GeneList[h][0], i.entrez, effectConfidence, foldchange, p])
                found=True
        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])
    return newGeneList


def szgeneLister(szgeneFile,GeneMap,allgenenames):
    rawList=csvLister(szgeneFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i>0:
            I2=1.2+((100.0-float(rawList[i][3]))*0.002)
            pList= rawList[i][2].split('(')
            p=float(pList[1][:-1])+0.000001
            rawList[i][2]=p
            if rawList[i][2]<0.05:
                effectConfidence=calculateEffectConfidence(I2,p)
                GeneList.append([rawList[i][0], effectConfidence, I2, p])

    GeneList=szgeneAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList

def chromoAliasDestroy(GeneList, GeneMap, allgenenames):
    newGeneList=[]
    for h in range(len(GeneList)):
        found=False
        effectConfidence=GeneList[h][1]
        foldchange=GeneList[h][2]
        p=GeneList[h][3]
        for i in GeneMap:
            if GeneList[h][0] in i.aliases:
                GeneList[h][0]=i.name
                found=True
            if GeneList[h][0] == i.name:
                newGeneList.append([GeneList[h][0], i.entrez, effectConfidence, foldchange, p])
                found=True
        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])

    return newGeneList




def chromoLister(chromoconfoFile,GeneMap,allgenenames):
    rawList=csvLister(chromoconfoFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i>0:
            foldchange=2.0**(abs(float(rawList[i][8])))
            if foldchange<1:
                foldchange=1.0/foldchange

            FDR=abs(float(rawList[i][9]))
            effectConfidence=calculateEffectConfidence(foldchange,FDR)
            GeneList.append([rawList[i][3],effectConfidence,foldchange,FDR])
            
    GeneList=chromoAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList








    # combinedList=Loci+szgene+hiPSC+DeadExpression+chromoconfo
    # print len(combinedList)
    # TotalMarks=combinedList


    # combinedList=list(set(combinedList))
    # print len(combinedList)




    # class GeneList:
    #     def __init__(self, name, genes):
    #         self.name=name
    #         self.genes=genes

    # class Combo:
    #     def __init__(self, GeneLists):
    #         self.GeneLists=GeneLists
    #         self.name=''
    #         self.Genes=[]
    #         self.GeneCount=0
    #         for collection in GeneLists:
    #             if self.name!='':
    #                 self.name=self.name+', '+ collection.name
    #             else:
    #                 self.name=collection.name

    #     def isin(self,gene):
    #         yes=True
    #         for GeneList in self.GeneLists:
    #             # try:


    #             if gene in GeneList.genes:
    #                 pass
    #             else:
    #                 return False
    #             # except:
    #             #     if gene in GeneList.genes:
    #             #         return True
    #         return True




    # Loci=GeneList('108 Loci', Loci)
    # szgene=GeneList('SZGene', szgene)
    # hiPSC=GeneList('hiPSC', hiPSC)
    # DeadExpression=GeneList('DeadExpression', DeadExpression)
    # chromoconfo=GeneList('Chromosome Conformation', chromoconfo)

    # x=[Loci,szgene,hiPSC,DeadExpression,chromoconfo]

    # ListOfCombos=[]

    # for i in range(len(x)):
    #     i=i+1
    #     for comb in itertools.combinations(x,i):
    #         ListOfCombos.append(comb)

    # for i in range(len(ListOfCombos)):
    #     ListOfCombos[i]=Combo(ListOfCombos[i])

    # combinedList.pop(0)


    # for i in range(len(combinedList)):
    #     theChosenCombo=''


    #     for j in range(len(ListOfCombos)):

    #         if ListOfCombos[j].isin(combinedList[i]):
    #             theChosenCombo=ListOfCombos[j]

    #     if theChosenCombo=='':
    #         faniwonfiabn
    #     theChosenCombo.GeneCount=theChosenCombo.GeneCount+1
    #     theChosenCombo.Genes.append(combinedList[i])
            
    # outputFile=open('GeneOverlap.txt','w')

    # for combo in ListOfCombos:
    #     outputfileline=combo.name+';'+str(combo.Genes)+'\n'
    #     outputFile.write(outputfileline)
    #     print
    #     print combo.name
    #     print combo.GeneCount
    #     print combo.Genes


    # for thing in x:
    #     print thing.name, len(thing.genes)
    # print 'gene names changed:', AliasDestructionCount.count
    # print AliasDestructionCount.genesChanged
    # print 'genes not found in the human gene map:', len(AliasDestructionCount.badGenes)
    # print AliasDestructionCount.badGenes

main()

