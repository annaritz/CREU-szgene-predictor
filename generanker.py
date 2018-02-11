print('Initializing Modules')
import time
import networkx as nx
from operator import itemgetter


def main():
    timesteps = 150
    print('Opening Files')
    SZnegativeFile=open('SZnegatives.csv','r')
    SZpositiveFile=open('SZpositives.txt','r')
    
    GeneMapfile=open('hugogenelist.txt','r')
    G = nx.Graph()
    GeneMap, nodeset=read_gene_map_file(GeneMapfile,G)
    positiveReader(SZpositiveFile, G, nodeset)
    negativeReader(SZnegativeFile, G, nodeset)
    edgeFile=open('brain_top_geq_0.200.txt','r')
    read_edge_file(edgeFile,G, nodeset)
    edgeFile.close()
    print(G.number_of_edges(), 'edges')
    print(G.number_of_nodes(), 'nodes')

    start=time.time()
    for t in range(0,timesteps):
        print('\n\n')
        print("t = " + str(t+1))
        iterativeMethod(G,t)
        inter_iterative_output(G)

        done=float(t)/float(timesteps)
        print('Time Elapsed:', time.time()-start)
        if done!=0:
            print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')

    print('Writing Output File')
    write_output(G, GeneMap)




    return




def read_gene_map_file(GeneMapfile, Graph):
    GeneMap=[]
    nodeset=set()
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
            # print(iname,ialiases,ientrez,iUID, iens)
            igene=[iname,ientrez]
            GeneMap.append(igene)
            # nodeset.add(ientrez)
            # Graph.add_node(ientrez, score=0.5, prev_score=0.5, label='Unlabeled', name=iname)
    print(len(GeneMap))
    return GeneMap, nodeset


def positiveReader(GeneFile, Graph, all_nodes):
    threshold=300
    x=0
    for line in GeneFile:
        x=x+1
        if x<threshold:
            entrezNumber=line.split('\t')
            entrezNumber=entrezNumber[0]
            if entrezNumber not in all_nodes:
                Graph.add_node(entrezNumber, prev_score=1.0, score=1.0, label='Positive', untouched=False)
            Graph.nodes[entrezNumber]['score']=1.0
            Graph.nodes[entrezNumber]['prev_score']=1.0
            Graph.nodes[entrezNumber]['label']='Positive'
            all_nodes.add(entrezNumber)

def negativeReader(GeneFile, Graph, all_nodes):
    for line in GeneFile:
        line=line.split(',')
        entrezNumber=line[0]
        if entrezNumber not in all_nodes:
            Graph.add_node(entrezNumber, prev_score=0.0, score=0.0, label='Negative', untouched=False)
        
        Graph.nodes[entrezNumber]['score']=0.0
        Graph.nodes[entrezNumber]['prev_score']=0.0
        Graph.nodes[entrezNumber]['label']='Negative'
        all_nodes.add(entrezNumber)


def read_edge_file(infile, Graph, all_nodes):

    start=time.time()

    x=0
    for line in infile:

        timepassedSinceStart=time.time()-start
        x=x+1
        done=x/3362057.0
        if x%100000==0:
            print() 
            print(done, 'percent completed')
            print('time remaining:', (1.0-done)*timepassedSinceStart/done, 'seconds')
        line=line.split('\t')
        line[2]=float(line[2][:len(line[2])-2])
        Graph.add_edge(line[0],line[1], weight=line[2])
        for i in range(0,2):
            node=line[i]
            if node not in all_nodes:
                Graph.add_node(node, prev_score=0.5, score=0.5, label='Unlabeled', untouched=True)
                all_nodes.add(node)
    return



#Takes as input a networkx graph
def iterativeMethod(Graph, t):
    positivechangesum=0
    sumofchanges=0
    changed=[]
    changedNegative=[]
    changedPositive=[]
    untouchedSet=set()
    #Note: Graph.nodes() is a list of all the nodes
    #Graph.nodes[node] is a dictionary of that node's attributes
    nodes = Graph.nodes()
    for node in nodes:
        #Want to keep positive scores at 1 and negative scores at 0 (or -1)
        if nodes[node]['label'] != 'Unlabeled':
            pass
        else:
            newConfidence = 0
            sumofWeights = 0
            #Note: Graph.adj[node].items() gives a list of tuples. Each tuple includes one of the
            #node's neighbors and a dictionary of the attributes that their shared edge has i.e. weight
            #neighbor is the node's neighbor, datadict is the dictionary of attributes
            for neighbor, datadict in Graph.adj[node].items(): 
                newConfidence = newConfidence + datadict['weight']*nodes[neighbor]['prev_score'] #use t-1 score
                sumofWeights = sumofWeights + datadict['weight']
            if sumofWeights>0:
                newScore = float(newConfidence)/float(sumofWeights)
                nodes[node]['score'] = newScore #update score to be the new score
        

    #After each time step is complete, the previous score is updated to be the current score
    #Prints the scores after each timestep is complete

    for node in nodes:
        if nodes[node]['prev_score'] != nodes[node]['score']:
            changed.append(node)
            nodes[node]['untouched']=False
            sumofchanges=sumofchanges+abs(nodes[node]['prev_score'] - nodes[node]['score'])
            if nodes[node]['prev_score'] > nodes[node]['score']:
                changedNegative.append(node)
            else:
                changedPositive.append(node)
                positivechangesum=positivechangesum+abs(nodes[node]['prev_score'] - nodes[node]['score'])


        if nodes[node]['untouched']==True:
            untouchedSet.add(node)





        nodes[node]['prev_score'] = nodes[node]['score']
    #     print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))
    print(len(changed), 'of', len(nodes), 'nodes changed')
    print(len(changedNegative), 'node scores decreased')
    print(len(changedPositive), 'node scores increased')
    print('Sum of absolute value of changes:', sumofchanges)
    print('Sum of positive changes:', positivechangesum)
    print('Untouched nodes:', len(untouchedSet))
    return



def write_output(Graph, GeneMap):
    start=time.time()
    nodes = Graph.nodes()
    x=1
    for node in nodes:
        x=x+1
        if x%1000==0:
            print() 
            print('assigning HUGO name to gene',x, 'of', len(nodes))
            done=x/41254.0
            print(done, 'percent completed')
            timepassedSinceStart=time.time()-start
            print('time remaining:', (1.0-done)*timepassedSinceStart/done, 'seconds')

        Found=False
        for gene in GeneMap:
            if gene[1]==node:
                Graph.nodes[node]['name']=gene[0]
                Found=True
        if not Found:
            Graph.nodes[node]['name']='No HUGO name'




    print('initializing list')
    nodeValues=[]
    for node in Graph.nodes:
        
        nodeValues.append([ node, Graph.nodes[node]['score'], Graph.nodes[node]['label'],Graph.nodes[node]['name']])

        

    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('gene_rankings.txt', 'w')
    y=open('gene_rankings_no_pos.txt', 'w')
    for node in nodeValues:
        # if G.nodes[node]['positive']==False:
        #     print(node)
        label=node[2]
        x.write(str(node)+'\n')
        if label=='Unlabeled':
            y.write(str(node)+'\n')
    return

def inter_iterative_output(Graph):
    start=time.time()
    nodes = Graph.nodes()






    nodeValues=[]
    for node in Graph.nodes:
        
        nodeValues.append([ node, Graph.nodes[node]['score'], Graph.nodes[node]['label']])

        

    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('mid_process_gene_rankings.txt', 'w')
    y=open('mide_process_gene_rankings_no_pos.txt', 'w')
    for node in nodeValues:
        # if G.nodes[node]['positive']==False:
        #     print(node)
        label=node[2]
        x.write(str(node)+'\n')
        if label=='Unlabeled':
            y.write(str(node)+'\n')
    return




main()










# GIANTbrain.close()


# print()

# for node in G.nodes():
#     if 'positive' not in G.nodes[node]:
#         G.nodes[node]['positive']=False
#     if 'name' not in G.nodes[node]:
#         G.nodes[node]['name']='no name'

# start=time.time()

# for iteration in range(iterations):
#     print('iteration', iteration+1)
#     dicv={}
#     nodes=G.nodes()


#     timepassedSinceStart=time.time()-start
#     x=x+1
#     done=x*100.0/float(iterations)

#     print() 
#     print(done, 'percent completed')
#     print('time remaining:', (100-done)*timepassedSinceStart/done, 'seconds')

    
#     for node in nodes:
        

#         try:
#             dicv[node]=G.nodes[node]['SZconfidence']
#         except:
#             dicv[node]=0.5
#             G.nodes[node]['SZconfidence']=0.5


#     for node in nodes:
#         newConfidence=0
#         sumofweights=0
#         for neighbor, datadict  in G.adj[node].items():
#             newConfidence=newConfidence+datadict['weight']*G.nodes[neighbor]['SZconfidence']
#             sumofweights=sumofweights+datadict['weight']
#         if sumofweights>0:
#             G.nodes[node]['SZconfidence']=float(newConfidence)/float(sumofweights)



















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



