import astropy.stats as ast
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
cellgeneFile=open('cellgene_rankings.txt','r')
szgeneFile=open('gene_rankings.txt','r')

class Gene:
    def __init__(self, entrez, name, SZscore, SZlabel):
        self.entrez=entrez
        self.name=name
        self.SZscore=SZscore
        self.SZlabel=SZlabel
        self.Cell_score=0
        self.Cell_label=None
        self.totalScore=0



SZlist=[]
SZnegativeentrez=set()
totalList=[]


for line in szgeneFile:
    line=line[1:len(line)-2]
    line=line.split(',')
    for i in range(len(line)):
        line[i]=line[i].lstrip(' ')
        line[i]=line[i].strip('\'')
    SZlist.append(Gene(line[0],line[3],line[1],line[2]))
    if 'Negative' in line[2]:
        SZnegativeentrez.add(line[0])



for line in cellgeneFile:
    line=line[1:len(line)-2]
    line=line.split(',')
    for i in range(len(line)):
        line[i]=line[i].lstrip(' ')
        line[i]=line[i].strip('\'')
    for gene in SZlist:
        if line[0]==gene.entrez:
            gene.Cell_score=line[1]
            gene.Cell_label=line[2]


for gene in SZlist:
    if gene.Cell_label==None:
        gene.Cell_label=' \'Unlabeled\' '
        gene.Cell_score=0.0


for gene in SZlist:
    score=(float(gene.SZscore)) * (float(gene.Cell_score))
    gene.totalScore=score
    totalList.append(gene)


x=sorted(totalList, key=lambda x: x.totalScore, reverse=True)

overlapScore=open('overlapScore.txt','w')
overlapScore.write('HGNC Symbol'+'\t'+'Total Score'+'\t'+'entrez'+'\t'+'SZ label'+'\t'+'Cell Motility Label'+'\n\n')

scoreList=[]

for gene in x:
    print(gene.name)
    print(gene.totalScore)
    print(gene.entrez)
    print(gene.SZlabel)
    print(gene.Cell_label+'\n')
    overlapScore.write(gene.name+'\t'+str(gene.totalScore)+'\t'+ gene.entrez+'\t'+gene.SZlabel+'\t'+gene.Cell_label+'\n')
    scoreList.append(gene.totalScore)

print(np.logspace(np.log10(0.0001),np.log10(1.0)))
x=sorted(scoreList, reverse=True)
# plt.hist(x, bins=np.logspace(np.log10(0.0001),np.log10(1.0)))
# pl.gca().set_xscale("log")
plt.hist(x, bins=50)
plt.title('Scores')
plt.show()












