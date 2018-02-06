import networkx as nx

def main():

    timesteps = 27

    G = nx.Graph()

    prob_dict = read_edge_file('toy-edges.txt', G)

    read_label_file('toy-labels.txt', G, prob_dict)

    print('t = 0')
    for node in G.nodes():
        print(str(node) + " Label: " + str(G.nodes[node]['label']) + ", Score: " + str(G.nodes[node]['score']))

    #Note: ranges from 0, 1, 2, ... so the current time step will be t+1 and the previous t 
    for t in range(0,timesteps):
        print('\n')
        print("t = " + str(t+1))
        iterativeMethod(G,t,prob_dict)


    return



#Adds an edge to networkx graph with weight attribute
#Adds nodes to the networkx graph with attributes for their score and label
#Initializes all nodes to be unlabeled and have scores of 0.5
#Initializes a probability dictionary that will keep track of every node's score at each time step
#Probability dictionary key:value = node:[scores], prob_dict[node][0] = score at t=0
def read_edge_file(infile, Graph):
    prob_dict = {}
    nodeset = set()
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            Graph.add_edge(line[0],line[1],weight=float(line[2]))
            for i in range(0,2): #two nodes in each line (edge)
                node = line[i]
                if node not in nodeset: #checks for duplicate nodes (the infile is an edge list)
                    Graph.add_node(node, score=0.5, label='Unlabeled')
                    prob_dict[node] = [0.5]
    return prob_dict

#Updates the labels and scores of the nodes in the graph
def read_label_file(infile, Graph, prob_dict):
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            node = line[0]
            if line[1] == 'pos':
                Graph.nodes[node]['label']='Positive'
                Graph.nodes[node]['score']=1
                prob_dict[node] = [1]
            elif line[1] == 'neg':
                Graph.nodes[node]['label']='Negative'
                Graph.nodes[node]['score']=0
                prob_dict[node] = [0]
    return

#Takes as input a networkx graph
def iterativeMethod(Graph, t, prob_dict):
    nodes = Graph.nodes()
    for node in nodes:
        newConfidence = 0
        sumofWeights = 0
        for neighbor, datadict in Graph.adj[node].items(): #no clue what that is
            prob_dict[neighbor][t]
            newConfidence = newConfidence + datadict['weight']*prob_dict[neighbor][t]
            sumofWeights = sumofWeights + datadict['weight']
        if sumofWeights>0:
            newScore = float(newConfidence)/float(sumofWeights)
            nodes[node]['score'] = newScore
            prob_dict[node].append(newScore)

    #Prints the scores after each timestep is complete
    for node in nodes:
        print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))

    return



if __name__ == '__main__':
    main()












