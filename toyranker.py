import networkx as nx
# from graphspace_python.api.client import GraphSpace
# from graphspace_python.graphs.classes.gsgraph import GSGraph
from datetime import datetime
import matplotlib.pyplot as plt
# graphspace = GraphSpace('mirbern@reed.edu', 'pumpkin')

def main():

    timesteps = 27

    G = nx.Graph()

    read_edge_file('toy-edges.txt', G)

    read_label_file('toy-labels.txt', G)

    print('t = 0')
    for node in G.nodes():
        print(str(node) + " Label: " + str(G.nodes[node]['label']) + ", Score: " + str(G.nodes[node]['score']))

    #Note: ranges from 0, 1, 2, ... so the current time step will be t+1 and the previous t 
    for t in range(0,timesteps):
        print('\n')
        print("t = " + str(t+1))
        iterativeMethod(G,t)

    # plot_network(graphspace, G)

    nx.draw(G, pos=nx.spring_layout(G), with_labels=True)
    plt.show()
    return


#Adds an edge to networkx graph with weight attribute
#Adds nodes to the networkx graph with attributes for their previous score (t-1) score (t) and label
#Initializes all nodes to be unlabeled and have scores of 0.5
def read_edge_file(infile, Graph):
    nodeset = set()
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            Graph.add_edge(line[0],line[1],weight=float(line[2]))
            for i in range(0,2): #two nodes in each line (edge)
                node = line[i]
                if node not in nodeset: #checks for duplicate nodes (the infile is an edge list)
                    Graph.add_node(node, prev_score=0.5, score=0.5, label='Unlabeled')
    return 

#Updates the labels and scores of the nodes in the graph
def read_label_file(infile, Graph):
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            node = line[0]
            if line[1] == 'pos':
                Graph.nodes[node]['label']='Positive'
                Graph.nodes[node]['prev_score']=1
                Graph.nodes[node]['score']=1
            elif line[1] == 'neg':
                Graph.nodes[node]['label']='Negative'
                Graph.nodes[node]['prev_score']=0
                Graph.nodes[node]['score']=0
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
            sumofchanges=sumofchanges+abs(nodes[node]['prev_score'] - nodes[node]['score'])
            if nodes[node]['prev_score'] > nodes[node]['score']:
                changedNegative += 1
            else:
                changedPositive += 1
                positivechangesum=positivechangesum+abs(nodes[node]['prev_score'] - nodes[node]['score'])







        nodes[node]['prev_score'] = nodes[node]['score']
        # print(nodes[node]['prev_score'])
        print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))
    print(changed, 'of', len(nodes), 'nodes changed')
    print(changedNegative, 'node scores decreased')
    print(changedPositive, 'node scores increased')
    print('Sum of absolute value of changes:', sumofchanges)
    print('Sum of positive changes:', positivechangesum)
    print('Untouched nodes:', untouchedSet)
    return

def rgb_to_hex(red,green,blue):
    maxHexValue = 255
    r = int(red*maxHexValue)
    g = int(green*maxHexValue)
    b = int(blue*maxHexValue)
    RR = format(r, '02x')
    GG = format(g, '02x')
    BB = format(b, '02x')
    return '#'+RR+GG+BB


#Posts iterative method on toy example network to graphspace
def plot_network(graphspace, GR):
    toy = GSGraph()
    toy.set_name('Toy Example Iterative Method' + ' ' + str(datetime.now()))
    toy.set_tags(['CREU'])

    #The code below colors the nodes according to their scores
    for node in GR.nodes():
        score = GR.nodes[node]['score']
        node_color = rgb_to_hex(1-score, 0, score)
        label = GR.nodes[node]['label']
        
        print(node)
        print(type(node))

        print(label)
        print(type(label))
        toy.add_node(node, label=node, popup=label + ', ' + score)
        toy.add_node_style(node, color=node_color)

    for edge in GR.edges():
        weight = GR.edges[edge]['weight']

        toy.add_edge(edge[0], edge[1])
        toy.add_edge_style(edge[0], edge[1], width=weight)
    

    graphspace.post_graph(toy)
    print("Graph posted!")

    return


if __name__ == '__main__':
    main()












