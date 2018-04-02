import astropy.stats as ast
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
cellgeneFile=open('cellgene_rankings.txt','r')
szgeneFile=open('gene_rankings.txt','r')

class Gene:
    def __init__(self, entrez, SZscore, SZlabel):
        self.entrez=entrez
        self.name='None'
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
    if line[1]!= '0.5':
        SZlist.append(Gene(line[0],line[1],line[2]))
    else:
        SZlist.append(Gene(line[0],'0.0',line[2]))

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
            if line[1]!= '0.5':
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


GeneMapfile=open('Homo_sapiens.txt', 'r')
GeneMap=[]
for i in GeneMapfile:
    i=i.split('\t')
    iname=i[2]
    
    if iname!='Approved Symbol':
        
        ientrez=i[1]
        igene=[ientrez,iname]
        GeneMap.append(igene)




for gene in x:
    counter=0
    found=False
    while found==False:
        if gene.entrez==GeneMap[counter][0]:
            print(GeneMap[counter])
            gene.name=GeneMap[counter][1]
            found=True
        counter=counter+1
        if counter==60566:
            found=True










overlapScore=open('overlapScore.txt','w')
overlapScore.write('HGNC Symbol'+'\t'+'SZ Score'+'\t'+'Cell Score'+'\t'+'Total Score'+'\t'+'entrez'+'\t'+'SZ label'+'\t'+'Cell Motility Label'+'\n\n')

scoreList=[]
for gene in x:
    # print(str(gene.name)+'\t'+str(gene.SZscore)+'\t'+str(gene.Cell_score)+'\t'+str(gene.totalScore)+'\t'+ gene.entrez+'\t'+gene.SZlabel+'\t'+gene.Cell_label+'\n')

    overlapScore.write(str(gene.name)+'\t'+str(gene.SZscore)+'\t'+str(gene.Cell_score)+'\t'+str(gene.totalScore)+'\t'+ gene.entrez+'\t'+gene.SZlabel+'\t'+gene.Cell_label+'\n')
    scoreList.append(gene.totalScore)


print(np.logspace(np.log10(0.0001),np.log10(1.0)))
x=sorted(scoreList, reverse=True)
# plt.hist(x, bins=np.logspace(np.log10(0.0001),np.log10(1.0)))
# pl.gca().set_xscale("log")
plt.hist(x, bins=50)
plt.title('Scores')
plt.show()












