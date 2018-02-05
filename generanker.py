import time
import networkx as nx
from operator import itemgetter

iterations=1000

GWASoverlap=open('GeneOverlap.txt','r')

GeneMapfile=open('hugogenelist.txt','r')

GeneMap=[]
allgenenames=[]
DataCombos=[]

posEntrez={}
posEntrez=set()

G = nx.Graph()

# class Gene:
#     def __init__(self, name, aliases=None, entrez=None, UID=None, ens=None):
#         self.name=name
#         self.aliases=aliases
#         self.entrez=entrez
#         self.UID=UID
#         self.ens=ens
#         if self.ens==None:
#             self.ens=''

class DataCombo:
    def __init__(self, name, genes=None):
        self.name=name
        self.genes=genes

for i in GeneMapfile:

    i=i.split('\t')
    iname=i[0]

    if iname!='Approved Symbol':
        iUID=i[4]
        ientrez=i[3]

        ialiases=i[1]+','+i[2]

        ialiases=ialiases.split(',')
        for j in range(len(ialiases)):
            if ialiases[j]!='':
                if ialiases[j][0]==' ':
                    ialiases[j]=ialiases[j][1:]
        if ialiases[0]=='' and ialiases[1]=='':
            ialiases=['']
        if len(ialiases)>1:
            if ialiases[0]=='' and ialiases[1]!='':
                ialiases=ialiases[1:]
        iens=i[5]


        print(iname,ialiases,ientrez,iUID, iens)
        igene=[iname,ientrez]
        GeneMap.append(igene)
        G.add_node(ientrez, SZconfidence=0, positive=False, name=iname)


print(len(GeneMap))
def geneLister(superstring):
    length=len(superstring)
    superstring=superstring[1:length-2]
    superstring=superstring.replace('\'','')
    superstring=superstring.split(',')
    GeneList=superstring
    return GeneList

def geneobjectifier(GeneList):
    SuperGeneList=[]
    for i in range(len(GeneList)):
        found=False
        if GeneList[i]=='':
            continue

        if GeneList[i][0]==' ':
            GeneList[i]=GeneList[i][1:]
        for supergene in GeneMap:
            if supergene[0]==GeneList[i]:
                found=True
                GeneList[i]=supergene
                posEntrez.add(supergene[1])
                G.nodes[supergene[1]]['SZconfidence']=1
                G.nodes[supergene[1]]['positive']=True

    return GeneList




for combo in GWASoverlap:
    combo=combo.split(';')
    comboname=combo[0]
    print('----------------')
    print(comboname)
    genes=combo[1]
    # print genes
    genes=geneLister(genes)
    genes=geneobjectifier(genes)
    DataCombos.append(DataCombo(comboname,genes))

for DataCombo in DataCombos:
    print(DataCombo.name)
    print(DataCombo.genes)

GIANTbrain=open('brain_top.txt','r')

start=time.time()
# GIANTentrez={}
# GIANTentrez=set()
# GIANTentrez5={}
# GIANTentrez5=set()
# GIANTentrez75={}
# GIANTentrez75=set()
# GIANTentrez9={}
# GIANTentrez9=set()
# GIANTentrez25={}
# GIANTentrez25=set()
# GIANTentrez10={}
# GIANTentrez10=set()
x=0
for line in GIANTbrain:

    timepassedSinceStart=time.time()-start
    x=x+1
    done=x*100.0/41668864.0
    if x%100000==0:
        print() 
        print(done, 'percent completed')
        print('time remaining:', (100-done)*timepassedSinceStart/done, 'seconds')

    line=line.split('\t')
    line[2]=float(line[2][:len(line[2])-2])

    if line[2]>0.2:

        G.add_edge(line[0],line[1], weight=line[2])


GIANTbrain.close()

print(G.number_of_edges())
print(G.number_of_nodes())
print()

for node in G.nodes():
    if 'positive' not in G.nodes[node]:
        G.nodes[node]['positive']=False
    if 'name' not in G.nodes[node]:
        G.nodes[node]['name']='no name'



for iteration in range(iterations):
    print('iteration', iteration+1)
    dicv={}
    nodes=G.nodes()
    
    for node in nodes:
        

        try:
            dicv[node]=G.nodes[node]['SZconfidence']
        except:
            dicv[node]=0.5
            G.nodes[node]['SZconfidence']=0.5


    for node in nodes:
        newConfidence=0
        sumofweights=0
        for neighbor, datadict  in G.adj[node].items():
            newConfidence=newConfidence+datadict['weight']*G.nodes[neighbor]['SZconfidence']
            sumofweights=sumofweights+datadict['weight']
        if sumofweights>0:
            if G.nodes[node]['positive']==False:
                G.nodes[node]['SZconfidence']=float(newConfidence)/float(sumofweights)

print('initializing list')
nodeValues=set()
for node in G.nodes:
    
    nodeValues.add((node, G.nodes[node]['SZconfidence'], G.nodes[node]['positive'],G.nodes[node]['name']))
    

nodeValues=sorted(nodeValues, key=itemgetter(1))


x=open('gene_rankings.txt', 'w')
for node in nodeValues:
    # if G.nodes[node]['positive']==False:
    #     print(node)
    entrez,value,positive,name=node
    if not positive:
        print(node)
        x.write(str(node)+'\n')


















# start=time.time()
# x=0
# for node in G:
#     x=x+1
#     done=x*100.0/G.number_of_nodes()
#     if x%1==0:
#         print() 
#         print(done, 'percent completed')
#         print('time remaining:', (100-done)*timepassedSinceStart/done, 'seconds')
#     if len(G.nodes[node])>0:
#         if G.nodes[node]['positive']==False:

#             positiveNeighborCount=0
#             weightedCount=0.0
#             for adjacent,datadict in G.adj[node].items():
#                 if len(G.nodes[adjacent])>0:
#                     if G.nodes[adjacent]['positive']==True:
#                         positiveNeighborCount=positiveNeighborCount+1
#                         weightedCount=weightedCount+datadict['weight']
#             if positiveNeighborCount>0:
#                 HasPositiveNeighbors.append([node,positiveNeighborCount])
#                 weightedHasPositiveNeighbors.append([node,weightedCount])



#     # for gene in GeneMap:

# SortedHasPositiveNeighbors=sorted(HasPositiveNeighbors, key=lambda student: student[1], reverse=True)
# SortedWeighted=sorted(weightedHasPositiveNeighbors, key=lambda student: student[1], reverse=True)

# for gene in SortedHasPositiveNeighbors:
#     for supergene in GeneMap:
#         if supergene[1]==gene[0]:
#             gene[0]=supergene[0]


# for gene in SortedWeighted:
#     for supergene in GeneMap:
#         if supergene[1]==gene[0]:
#             gene[0]=supergene[0]
















# output=open('genesWithAdjacentPositivesSansPositives0.5.txt','w')

# for i in SortedHasPositiveNeighbors:    
#     output.write(str(i)+'\n')

# output.write('\n\n\n')

# for i in SortedWeighted:    
#     output.write(str(i)+'\n')



