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

    print('Reading CMC Expression')
    DeadFile=open('common_mind/Genes-Table.csv','r')
    DeadExpressionList=DeadGeneLister(DeadFile,GeneMap, allgenenames)

    # print('Reading CMC Isoform Expression')
    # IsoFile=open('common_mind/Isoforms-Table.csv','r')
    # DeadExpressionList2=DeadGeneLister2(IsoFile,GeneMap, allgenenames)

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

    # sortExpressionList=sorted(BigList, key=lambda x: x[2])
    # for i in sortExpressionList:
    #     print(i)



    
    deleteList=[]
    print(len(BigList))
    overlaps=0
    already_done=[]
    BigSet=set()
    for i in range(len(BigList)):
        BigSet.add(BigList[i][1])
        for j in range(len(BigList)):
            if BigList[i][1]==BigList[j][1]:
                if i!=j and (i not in already_done or j not in already_done):
                    already_done.append(i)
                    already_done.append(j)

                    deleteList.append(i)
                    overlaps=overlaps+1
                    if BigList[j][2]!=BigList[i][2]:
                        BigList[j][2]=BigList[i][2]+' and '+BigList[j][2]
      


    for i in range(len(deleteList)):
        BigList.pop(deleteList[i])
        for j in range(len(deleteList)):
            if deleteList[j]>deleteList[i]:
                deleteList[j]=deleteList[j]-1
    print(len(BigList), 'BigList')
    print(len(BigSet), 'BigSet')
    BigList=sorted(BigList, key=lambda x: x[2])
    outputFile1=open('SZ_positives.txt', 'w')
    outputFile2=open('SZ_positives_source.txt','w')
    # outputFile2=open('E2_positives.txt', 'w')
    # outputFile3=open('E3_positives.txt', 'w')
    overlapControlledECdistribution=[]

    E1=[]
    E2=[]
    E3=[]
    Actual_positives=[]



    for gene in BigList:
        outputFile1.write(str(gene[1])+'\n')
        outputFile2.write(str(gene[1])+'\t'+str(gene[2]+'\n'))


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
                GeneList[h]=[GeneList[h][0], i.entrez, 'PGC']
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
        GeneList.append([rawList[i][3]])
    GeneList=LociAliasDestroy(GeneList, GeneMap, allgenenames)

    return GeneList


    # Loci=csvLister(Loci)
    # Loci=LociGeneLister(Loci)

def DeadAliasDestroy(GeneList, GeneMap, allgenenames):
    for h in range(len(GeneList)):
        found=False


        for i in GeneMap:

            if GeneList[h][0] in i.aliases:

                GeneList[h][0]=i.name
                found=True

            if GeneList[h][0] == i.name:
                GeneList[h]=[GeneList[h][0], i.entrez, 'CMC Expression']
        
                found=True

        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])



    return GeneList

def DeadGeneLister(deadFile,GeneMap, allgenenames):
    rawList=csvLister(deadFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i>1:


            if rawList[i][1] not in allgenenames or rawList[i][1]=='.':
                
                for gene in GeneMap:
                    if gene.ens==rawList[i][0]:
                        # print
                        rawList[i][1]=gene.name

            if rawList[i][1]!='.':
                GeneList.append([rawList[i][1]])


    GeneList=DeadAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList

def DeadAliasDestroy2(GeneList, GeneMap, allgenenames):
    for h in range(len(GeneList)):
        found=False


        for i in GeneMap:

            if GeneList[h][0] in i.aliases:

                GeneList[h][0]=i.name
                found=True

            if GeneList[h][0] == i.name:
                GeneList[h]=[GeneList[h][0], i.entrez, 'CMC Isoform Expression']
        
                found=True

        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])

    return GeneList

def DeadGeneLister2(deadFile,GeneMap, allgenenames):
    rawList=csvLister(deadFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i>1 and rawList[i]!=['', '', '', '', '', '', '', '', '', '', '\n']:

            if rawList[i][2] not in allgenenames or rawList[i][2]=='.':
                
                for gene in GeneMap:
                    if gene.ens==rawList[i][0]:
                        
                        rawList[i][2]=gene.name

            if rawList[i][2]!='.' and rawList[i][2]!='':
                GeneList.append([rawList[i][2]])


    GeneList=DeadAliasDestroy2(GeneList, GeneMap, allgenenames)
    return GeneList


def szgeneAliasDestroy(GeneList, GeneMap, allgenenames):
    newGeneList=[]
    for h in range(len(GeneList)):
        found=False

        for i in GeneMap:
            if GeneList[h][0] in i.aliases:
                GeneList[h][0]=i.name
                found=True
            if GeneList[h][0] == i.name:
                newGeneList.append([GeneList[h][0], i.entrez, 'SZGene'])
                found=True
        if GeneList[h][0] not in allgenenames and not found:
            badGenes.append(GeneList[h])
    return newGeneList


def szgeneLister(szgeneFile,GeneMap,allgenenames):
    rawList=csvLister(szgeneFile)
    GeneList=[]
    for i in range(len(rawList)):
        if i>0:
            pList= rawList[i][2].split('(')
            p=float(pList[1][:-1])+0.000001
            rawList[i][2]=p
            if rawList[i][2]<0.05:
                GeneList.append([rawList[i][0]])

    GeneList=szgeneAliasDestroy(GeneList, GeneMap, allgenenames)
    return GeneList

main()

