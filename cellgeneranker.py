print('Initializing Modules')
import time
import networkx as nx
from operator import itemgetter


def main():
    timesteps = 150
    print('Opening Files')
    SZnegativeFile=open('SZnegatives.csv','r') #Uses SZ negatives for now
    CellpositiveFile=open('mastergenelist2.txt','r') 
    
    G = nx.Graph()
    nodeset = set()

    positiveReader(CellpositiveFile, G, nodeset) #adds positive nodes to graph
    negativeReader(SZnegativeFile, G, nodeset) #adds negative nodes to graph

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

        done=float(t)/float(timesteps)
        print('Time Elapsed:', time.time()-start)
        if done!=0:
            print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')

    print('Writing Output File')
    write_output(G)




    return

#Takes in positive gene file - two columns Gene, EntrezID
def positiveReader(GeneFile, Graph, all_nodes):
    for line in GeneFile:
        line=line.strip().split('\t')
        if line[0] == 'Gene':
            pass
        else:
            entrezNumber=line[1]
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



def write_output(Graph):
    print('initializing list')
    nodeValues=[]
    for node in Graph.nodes:
        
        nodeValues.append([ node, Graph.nodes[node]['score'], Graph.nodes[node]['label']])

        

    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('cellgene_rankings.txt', 'w')
    y=open('cellgene_rankings_no_pos.txt', 'w')
    for node in nodeValues:
        # if G.nodes[node]['positive']==False:
        #     print(node)
        label=node[2]
        x.write(str(node)+'\n')
        if label=='Unlabeled':
            y.write(str(node)+'\n')
    return






main()


