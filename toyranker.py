import networkx as nx
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
from datetime import datetime

graphspace = GraphSpace('mirbern@reed.edu', 'pumpkin')

def main():

    timesteps = 3

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


    return



#Adds an edge to networkx graph with weight attribute
#Adds nodes to the networkx graph with attributes for their score and label
#Initializes all nodes to be unlabeled and have scores of 0.5
#Initializes a probability dictionary that will keep track of every node's score at each time step
#Probability dictionary key:value = node:[scores], prob_dict[node][0] = score at t=0
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
    nodes = Graph.nodes()
    for node in nodes:
        #Want to keep positive scores at 1 and negative scores at 0 (or -1)
        if nodes[node]['label'] != 'Unlabeled':
            pass
        else:
            newConfidence = 0
            sumofWeights = 0
            for neighbor, datadict in Graph.adj[node].items(): 
                newConfidence = newConfidence + datadict['weight']*nodes[neighbor]['prev_score'] #use t-1 score
                sumofWeights = sumofWeights + datadict['weight']
            if sumofWeights>0:
                newScore = float(newConfidence)/float(sumofWeights)
                nodes[node]['score'] = newScore #update score to be the new score
        

    #After each time step is complete, the previous score is updated to be the current score
    #Prints the scores after each timestep is complete
    for node in nodes:
        nodes[node]['prev_score'] = nodes[node]['score']
        print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))

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


#Posts EGFR network to graphspace
def plot_network(graphspace, G):
    ITR = GSGraph()
    ITR.set_name('Toy Example Iterative Method' + ' ' + str(datetime.now()))
    ITR.set_tags(['CREU'])

    #The code below colors the nodes according to their scores
    for node in G.nodes():
        count = walk_dict[node]
        ratio = float(count)/float(count_max)
        node_color = rgb_to_hex(1-ratio, 0, ratio)
        
        EGFR.add_node(node, label=node)
        EGFR.add_node_style(node, color=node_color)

    for i in range(len(edge_list)):
        edge = edge_list[i]
        EGFR.add_edge(edge[0], edge[1], directed=True)

        edge_type = type_dict[i] #the index of an edge in the edge list is the key in the edge type dictionary 
    
        #Edges are colored based on the type of interaction 
        if edge_type == 'Interaction':
            edge_color = 'black'
        elif edge_type == 'protein_cleavage':
            edge_color = 'cyan'
        elif edge_type == 'Phosphorylation':
            edge_color = 'yellow'
        elif edge_type == 'Ubiquitination':
            edge_color = 'red'
        elif edge_type == 'Ligand_Binding':
            edge_color = 'blue'
        elif edge_type == 'Dephosphorylation':
            edge_color = 'green'

        EGFR.add_edge_style(edge[0], edge[1], directed=True, edge_style='dotted', color=edge_color)

    graphspace.post_graph(EGFR)

    return


if __name__ == '__main__':
    main()












