print('Initializing Modules')
import time
import networkx as nx
from operator import itemgetter
import random
import matplotlib.pyplot as plt

#Don't add nodes that aren't in the network 

#Gene ranker code, with a few changes 
def main():
    x = 1 #run the tests x times
    for i in range(x):
        #Runs the algorithm on the test positives (1/4)
        test_genes, hidden_genes, all_genes = random_positives('mastergenelist.txt')

        timesteps = 150
        print('Opening Files')
        SZnegativeFile=open('SZnegatives.csv','r') #Uses SZ negatives for now
        CelltestPos=test_genes     
        G = nx.Graph()
        nodeset = set()

        hidden_in_graph = positiveReader(CelltestPos, hidden_genes, G, nodeset) #adds positive nodes to graph
        negativeReader(SZnegativeFile, G, nodeset) #adds negative nodes to graph
        SZnegativeFile.close()

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
        hist_scores = get_output(G, hidden_in_graph,all_genes,timesteps,x) #writes all the output files and returns a list of 
        #hidden scores to make the histogram

        plot_histogram(hist_scores, timesteps, x)



    return




#Takes the master cell motility gene list, converts it into a set
#Uses random.sample to randomly choose 1/4 of the list 
def random_positives(positives):
    master_positives = set() #set of EntrezID master positives
    with open(positives, 'r') as pos:
        for line in pos:
            line = line.strip().split('\t')
            if line[0] == 'Gene': 
                pass
            else:
                master_positives.add(line[1])

    hidden_pos = random.sample(master_positives, (len(master_positives)//4)) #randomly choose 1/4
    #to remove from the set of positives
    hidden_pos = set(hidden_pos)

    test_pos = master_positives.difference(hidden_pos) #test positives: positives in 
    #master list but not in hidden positives
    return test_pos, hidden_pos, master_positives



#Takes in the test positives (set of positives/EntrezIDs), the Graph, and the nodes in the Graph
def positiveReader(GeneList, HiddenGeneList, Graph, all_nodes):
    GeneList = list(GeneList) #converts it from a set to a list in order to iterate through it
    HiddenGeneList = list(HiddenGeneList)
    notinGraph = 0 #number of known positives not in the graph
    HiddenNotInGraph = 0 #number of hidden positives not in the graph
    hiddenInGraph = set() #set of hidden positives not in the graph
    for gene in GeneList:
        if gene not in all_nodes: #if a positive isn't in the graph, do not add it. Hold on to it 
            #so that it's score is not checked later
            notinGraph += 1
            Graph.add_node(gene, prev_score=1.0, score=1.0, label='Positive', untouched=False)
        Graph.nodes[gene]['score']=1.0
        Graph.nodes[gene]['prev_score']=1.0
        Graph.nodes[gene]['label']='Positive'
        all_nodes.add(gene)

    for hid in HiddenGeneList:
        if hid not in all_nodes:
            HiddenNotInGraph += 1
        else:
            hiddenInGraph.add(hid)

    print('There are ' + str(notinGraph) + ' known positives not in the graph.')
    print('There are ' + str(HiddenNotInGraph) + ' hidden positives not in the graph.')  


    return hiddenInGraph

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
    changed=0
    changedNegative=0
    changedPositive=0
    untouchedSet=0
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
            changed += 1
            nodes[node]['untouched']=False
            sumofchanges=sumofchanges+abs(nodes[node]['prev_score'] - nodes[node]['score'])
            if nodes[node]['prev_score'] > nodes[node]['score']:
                changedNegative += 1
            else:
                changedPositive += 1
                positivechangesum=positivechangesum+abs(nodes[node]['prev_score'] - nodes[node]['score'])


        if nodes[node]['untouched']==True:
            untouchedSet += 1





        nodes[node]['prev_score'] = nodes[node]['score']
    #     print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))
    print(changed, 'of', len(nodes), 'nodes changed')
    print(changedNegative, 'node scores decreased')
    print(changedPositive, 'node scores increased')
    print('Sum of absolute value of changes:', sumofchanges)
    print('Sum of positive changes:', positivechangesum)
    print('Untouched nodes:', untouchedSet)
    return



def get_output(Graph, hidden_positives, all_pos, timesteps, x):
    print('initializing list')

    #Gets all the hidden positives scores for a histogram
    #and prints a separate file
    hid_pos_histogram=[] #gets a list of scores for the distribution histogram
    hidden_positives=list(hidden_positives) #converts from set to list
    with open('cell_testrankings_hidden_' + str(timesteps) + '_' + str(x) + '.txt', 'w') as file:
        for pos in hidden_positives:
            hid_pos_histogram.append(Graph.nodes[pos]['score'])
            file.write(str([pos, Graph.nodes[pos]['score'], Graph.nodes[pos]['label']]))


    #Gets all the node scores
    nodeValues=[]
    for node in Graph.nodes:
        nodeValues.append([ node, Graph.nodes[node]['score'], Graph.nodes[node]['label']])

        

    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('cellgene_testrankings_' + str(timesteps) + '_' + str(x) + '.txt', 'w')
    y=open('cellgene_testrankings_no_pos_' + str(timesteps) + '_' + str(x) + '.txt', 'w')
    for node in nodeValues:
        # if G.nodes[node]['positive']==False:
        #     print(node)
        label=node[2]
        x.write(str(node)+'\n')
        if label=='Unlabeled':
            y.write(str(node)+'\n')
    x.close()
    y.close()

    return hid_pos_histogram



#Takes in a list of the hidden gene scores and plots a distribution
def plot_histogram(hist_list, timesteps, x):
    fig = plt.figure(figsize=(6.5,4))
    plt.hist(hist_list, 20)
    plt.xlabel('Score')
    plt.ylabel('Frequency in Hidden Positives')
    plt.title('Histogram of Hidden Positive Rankings')
    plt.tight_layout()
    plt.savefig('cellgene_testhist_' + str(timesteps) + '_' + str(x) + '.png')
    print('Histogram done!')
    return




if __name__ == '__main__':
    main()

    

