print('Initializing Modules')
import time
import networkx as nx
from operator import itemgetter
import matplotlib.pyplot as plt


def main():
    timesteps = 50
    chartSteps = 500
    chartPrecision = timesteps/chartSteps
    chartMarker = 0

    print('Opening Files')
    SZnegativeFile=open('SZnegatives.csv','r') #Uses SZ negatives for now
    autismFile= open('autism.csv', 'r')
    E1_Positive_File=open('E1_positives.txt','r') 
    E2_Positive_File=open('E2_positives.txt','r') 
    E3_Positive_File=open('E3_positives.txt','r') 
    
    G = nx.Graph()
    nodeset = set()
    print("Initializing Graph")
    edgeFile=open('brain_top_geq_0.200.txt','r')
    read_edge_file(edgeFile,G, nodeset)
    edgeFile.close()

    print(G.number_of_edges(), 'edges')
    print(G.number_of_nodes(), 'nodes')

    autismReader(autismFile, G, nodeset)
    print(fbsbdfkhb)
    positiveReader(E1_Positive_File, E2_Positive_File, E3_Positive_File, G, nodeset) #adds positive nodes to graph
    negativeReader(SZnegativeFile, G, nodeset) #adds negative nodes to graph
    


    print(G.number_of_edges(), 'edges')
    print(G.number_of_nodes(), 'nodes')

    start=time.time()
    now=start
    iterationLogger=[]
    timelogger=[]
    for t in range(0,timesteps):
        print('\n\n')
        print("t = " + str(t+1))
        prev=now
        now=time.time()


        if t > chartMarker:
            iterationLogger.append(t)
            timelogger.append(now-prev)
            chartMarker=chartMarker+chartPrecision
            plt.plot(iterationLogger, timelogger)
            plt.ylabel('Time Per Iteration')
            plt.xlabel('Iteration')
            plt.savefig('timeCourse.png')



        for Evidence_Level in ['_E1','_E2','_E3']:
            iterativeMethod(G,t, Evidence_Level)
        prime_method(G, 0.5, ['_E1','_E2','_E3'])

        done=float(t)/float(timesteps)
        print('Time Elapsed:', now-start)
        if done!=0:
            print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')



    print('Writing Output File')
    prev=now
    now=time.time()
    iterationLogger.append(timesteps)
    timelogger.append(now-prev)
    chartMarker=chartMarker+chartPrecision
    plt.plot(iterationLogger, timelogger)
    plt.ylabel('Time Per Iteration')
    plt.xlabel('Iteration')
    plt.savefig('timeCourse.png')

    ranked_cand=write_output(G)

    plot_candidate_degrees(ranked_cand, G)





    return

#Takes in positive gene file - two columns Gene, EntrezID
def positiveReader(GeneFile1,GeneFile2, GeneFile3, Graph, all_nodes):
    for line in GeneFile1:
        line=line.strip().split('\t') 
        if line[0] == 'Gene':
            pass
        else:
            entrezNumber=line[0]
            if entrezNumber in all_nodes: #This will ignore nodes that are not in the graph.
                Graph.nodes[entrezNumber]['score_E1']=1.0
                Graph.nodes[entrezNumber]['prev_score_E1']=1.0
                Graph.nodes[entrezNumber]['label_E1']='Positive'
    for line in GeneFile2:
        line=line.strip().split('\t') 
        if line[0] == 'Gene':
            pass
        else:
            entrezNumber=line[0]
            if entrezNumber in all_nodes: #This will ignore nodes that are not in the graph.
                Graph.nodes[entrezNumber]['score_E2']=1.0
                Graph.nodes[entrezNumber]['prev_score_E2']=1.0
                Graph.nodes[entrezNumber]['label_E2']='Positive'
    for line in GeneFile3:
        line=line.strip().split('\t') 
        if line[0] == 'Gene':
            pass
        else:
            entrezNumber=line[0]
            if entrezNumber in all_nodes: #This will ignore nodes that are not in the graph.
                Graph.nodes[entrezNumber]['score_E3']=1.0
                Graph.nodes[entrezNumber]['prev_score_E3']=1.0
                Graph.nodes[entrezNumber]['label_E3']='Positive'

#Takes in positive gene file - two columns Gene, EntrezID
def negativeReader(GeneFile, Graph, all_nodes):
    for line in GeneFile:
        line=line.split(',')
        entrezNumber=line[0]
        if entrezNumber in all_nodes:
            Graph.nodes[entrezNumber]['score_E1']=0.0
            Graph.nodes[entrezNumber]['prev_score_E1']=0.0
            Graph.nodes[entrezNumber]['label_E1']='Negative'
            Graph.nodes[entrezNumber]['score_E2']=0.0
            Graph.nodes[entrezNumber]['prev_score_E2']=0.0
            Graph.nodes[entrezNumber]['label_E2']='Negative'
            Graph.nodes[entrezNumber]['score_E3']=0.0
            Graph.nodes[entrezNumber]['prev_score_E3']=0.0
            Graph.nodes[entrezNumber]['label_E3']='Negative'


            # all_nodes.add(entrezNumber)


def autismReader(GeneFile, Graph, all_nodes):
    for line in GeneFile:
        line=line.split(',')
        entrezNumber=line[1]
        if entrezNumber in all_nodes:
            Graph.nodes[entrezNumber]['score_prime']=line[5]
    plot_candidate_degrees2(write_output(Graph),Graph)


#Takes in edge list file from Humanbase - 3 columns gene 1, gene 2, functional interaction probability
def read_edge_file(infile, Graph, all_nodes):


    for line in infile:

        line=line.split('\t')
        line[2]=float(line[2][:len(line[2])-2]) #Wrestling with the formatting
        
        for i in range(0,2):
            node=line[i]
            if node not in all_nodes:
                Graph.add_node(node, prev_score_E1=0.5, score_E1=0.5, label_E1='Unlabeled',prev_score_E2=0.5, score_E2=0.5, label_E2='Unlabeled',prev_score_E3=0.5, score_E3=0.5, label_E3='Unlabeled', score_prime=0.5,untouched=True)
                all_nodes.add(node)
        Graph.add_edge(line[0],line[1], weight=line[2])
        Graph.add_edge(line[0],line[1], weight=line[2])
        Graph.add_edge(line[0],line[1], weight=line[2])
    return



#Takes as input a networkx graph
def iterativeMethod(Graph, t, level):
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
        if nodes[node]['label'+level] != 'Unlabeled':
            pass
        else:
            newConfidence = 0
            sumofWeights = 0
            #Note: Graph.adj[node].items() gives a list of tuples. Each tuple includes one of the
            #node's neighbors and a dictionary of the attributes that their shared edge has i.e. weight
            #neighbor is the node's neighbor, datadict is the dictionary of attributes
            for neighbor, datadict in Graph.adj[node].items(): 
                newConfidence = newConfidence + datadict['weight']*nodes[neighbor]['prev_score'+level] #use t-1 score
                sumofWeights = sumofWeights + datadict['weight']
            if sumofWeights>0:
                newScore = float(newConfidence)/float(sumofWeights)
                nodes[node]['score'+level] = newScore #update score to be the new score
        

    #After each time step is complete, the previous score is updated to be the current score
    #Prints the scores after each timestep is complete

    for node in nodes:
        if nodes[node]['prev_score'+level] != nodes[node]['score'+level]:
            changed += 1
            nodes[node]['untouched']=False
            sumofchanges=sumofchanges+abs(nodes[node]['prev_score'+level] - nodes[node]['score'+level])
            if nodes[node]['prev_score'+level] > nodes[node]['score'+level]:
                changedNegative += 1
            else:
                changedPositive += 1
                positivechangesum=positivechangesum+abs(nodes[node]['prev_score'+level] - nodes[node]['score'+level])


        if nodes[node]['untouched']==True:
            untouchedSet += 1





        nodes[node]['prev_score'+level] = nodes[node]['score'+level]
        # print(nodes[node]['prev_score'])
    #     print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))
    print(changed, 'of', len(nodes), 'nodes changed')
    print(changedNegative, 'node scores decreased')
    print(changedPositive, 'node scores increased')
    print('Sum of absolute value of changes:', sumofchanges)
    print('Sum of positive changes:', positivechangesum)
    print('Untouched nodes:', untouchedSet)
    return

def prime_method(Graph, alpha, Levels):
    nodes = Graph.nodes()
    for node in nodes:
        nodes[node]['score_prime']=(nodes[node]['score_E1']+0.75*nodes[node]['score_E2']+0.5*nodes[node]['score_E3'])/(1+0.75+0.5)
        for level in Levels:
            if nodes[node]['label'+level]=='Unlabeled':
                nodes[node]['score'+level]=(alpha*nodes[node]['score_prime']+nodes[node]['score'+level])/(alpha+1)



def write_output(Graph):
    print('initializing list')
    nodeValues=[]
    for node in Graph.nodes:
        nodeValues.append([ node, Graph.nodes[node]['score_prime'], Graph.nodes[node]['score_E1'],Graph.nodes[node]['score_E2'],Graph.nodes[node]['score_E3'],Graph.nodes[node]['label_E1'],Graph.nodes[node]['label_E2'],Graph.nodes[node]['label_E3'], Graph.degree(node)])

        

    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('gene_rankings_400_node.txt', 'w')
    y=open('gene_rankings_400_nopos_node.txt', 'w')
    for node in nodeValues:
        # if G.nodes[node]['positive']==False:
        #     print(node)
        x.write(str(node)+'\n')
        if Graph.nodes[node[0]]['label_E1']!= 'Positive' and Graph.nodes[node[0]]['label_E2']!= 'Positive' and Graph.nodes[node[0]]['label_E3']!= 'Positive':
            y.write(str(node)+'\n')
    return nodeValues


def plot_candidate_degrees(rank, Graph):
    print('Creating candidate plot...')
    fig = plt.figure(figsize=(4,4))
    y = []
    degrees=[]
    for cand in rank:
        node = cand[0]
        degrees.append(Graph.degree(node))
    print(sum(degrees)/len(degrees))

    for cand in rank: #iterates through list of 
        if cand[5]!= 'Positive' and cand[6]!= 'Positive' and cand[7]!= 'Positive':
            index=0
            node = cand[0]
            deg = Graph.degree(node) #looks up the degree of a candidate, returns [(cand, degree)]

            
            y.append(deg)

    plt.plot(y,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')

    plt.tight_layout()

    plt.savefig('candidate_degrees.png')
    print('Wrote to candidate_degrees.png')


    return

def plot_candidate_degrees2(rank, Graph):
    print('Creating candidate plot...')
    degreeList=[]
    valueList=[]
    fig = plt.figure(figsize=(4,4))
    y = []
    for cand in rank: #iterates through list of 
        index=0
        node = cand[0]
        deg = Graph.degree(node) #looks up the degree of a candidate, returns [(cand, degree)]
        degreeList.append(deg)
        val = Graph.nodes[node]['score_prime']
        valueList.append(val)

        

    print(len(degreeList))
    movingAverage=[]
    alpha=15
    for i in range(len(degreeList)-alpha):
        average=sum(degreeList[i:i+alpha])/(alpha+1)
        movingAverage.append(average)

    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('autism_degrees.png')
    

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('autism_scores.png')

    plt.figure()
    plt.plot(valueList,degreeList,'ob')
    plt.xlabel('Node Score')
    plt.ylabel('Node Degrees')
    plt.title('Candidate Nodes')
    plt.tight_layout()
    plt.savefig('autism_degrees_by_score.png')

    degreeList=degreeList[0:300]
    movingAverage=movingAverage[0:300]
    plt.figure()
    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_autism_degrees.png')

    valueList=valueList[0:300]

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_autism_scores.png')
    print(vdfvksjdnvlkjs)





    print('Wrote to candidate_degrees.png')




main()


