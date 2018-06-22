    
import networkx as nx
import matplotlib.pyplot as plt
import math
import numpy
import glob
import os


def main():
    #Top edges of GIANT brain network downloaded from http://hb.flatironinstitute.org/download
    #Three column plain text file [entrez gene id 1][entrez gene id 2][posterior prob.]

    #gets network files #returns a list of all files with this path

    network_files = glob.glob('../source/infiles/networkfiles/brain_top_geq_*') 

    for file in network_files:
        print('READING FILE: ', file)
        threshold = os.path.splitext(file)[0][-5:] #gets the file name w/o the extension, then takes the last 5 characters
        #print(threshold)
        
        #Read the edge file - this function creates a networkx graph from our 3-column file of weighted edges 
        G = nx.read_weighted_edgelist(file)
        m = nx.number_of_edges(G)
        #print("Number of edges: ", m)

        node_list = list(G.nodes())
        n = nx.number_of_nodes(G)
        #print("Number of nodes: ", n)
        
        edge_list = read_edge_file(file)

        degree_dict = get_degrees(edge_list)
        #print('Degree dictionary: ', degree_dict)
        #print('Degree of gene ID 1:', degree_dict['1'])

        get_max_degree(degree_dict)

        avg_node = average_node_degree(degree_dict)
        #print('Average node degree: ', avg_node)

        histogram_dict = get_histogram(degree_dict)
        #print('Histogram dictionary: ', histogram_dict)

        plot_histogram(histogram_dict, threshold) 

        N_v = neighbor_dict(edge_list)
        #print("Neighbors dictionary (adjacency list): ", N_v)

        AND_dict = average_neighbor_degree(degree_dict, N_v)
        #print("Average node degree dictionary: ", AND_dict)

        avg_AND = degree_and_dictionary(degree_dict, AND_dict, histogram_dict)
        #print("Average AND degree:", avg_AND)

        plot_average_AND(avg_AND, threshold)

        path_list_100 = BFS(node_list, N_v, 100000)
        path_list_200 = BFS(node_list, N_v, 200000)

        plot_pathlength(path_list_100, threshold, 100000)
        plot_pathlength(path_list_200, threshold, 200000)
        
        N_v_weights = get_adj_list_with_weights(edge_list)
        #print("Adjacency list with weights: ", N_v_weights)
        neg_log_edges = log_probabilities(edge_list)
        SP_weights = dijkstra()


    return


def read_edge_file(edge_file):
    with open(edge_file, "r") as file:
        edges_list = []
        for line in file: 
            node_pair = line.strip() #strips whitespace characters
            edge = node_pair.split("\t") #splits at the tabs (there are tabs between each column)
            edges_list.append(edge)
    return edges_list


#Inputs: list of edges
#Outputs: list of edges where the probabilities are now the negative log
#This is done to compute shortest paths. Because highest weight edges (highest probabilities)
#are the "most important", but shortest paths depend on lowest weight edges, we convert
#those highest probabilities to lowest probabilites.
def log_probabilities(edge_list):
    for e in edge_list:
        e[2] = -numpy.log10(e[2])
    return edge_list

def get_max_degree(degree_dict):
    max_degree = 0
    max_node = []
    for node in degree_dict:
        if degree_dict[node] > max_degree:
            max_node = [node]
            max_degree = degree_dict[node]
        #We want to keep track of all high degree nodes so if 
        #a node has the same degree, it'll be added to the set of max degree nodes
        elif degree_dict[node] == max_degree:
            max_node.append(node)
            max_degree = degree_dict[node]
    print("These highest degree nodes ")
    print(max_node)
    print("have degree " + str(max_degree))
    print("\n")
    return max_node, max_degree



#Takes as input an edge list and outputs a dictionary of key:value pairs {node: degree}
def get_degrees(edge_list):
    degree = {} #initializes empty dictionary
    for edge in edge_list: #looks at every edge in the edge list
        for i in [0,1]:
            if edge[i] not in degree: #if the first node in the edge isn't already in the dictionary
                degree[edge[i]] = 1 #create a key value pair where the value is 1 because this is the first 
                #time it appears in the edge list (so far, it only has one connection and thus a degree of one)
            elif i == 1 and edge[i] == edge[0]: #If an edge is between a node and itself, you only want to 
            # increase the degree of that node by 1 and need to make sure to "discard" the second count because 
            # the degree is essentially the number of edges it is a part of 
                degree[edge[i]] += 0 
            elif edge[i] in degree: #if it's already in the dictionary
                degree[edge[i]] += 1 #add one to its degree
    return degree

#Takes as input a dictionary of the degrees of all the nodes
#and outputs the average node degree of the graph
def average_node_degree(degree_dict):
    average_node = 0
    for d in degree_dict:
        average_node += degree_dict[d]
    average_node = float(average_node)/float(len(degree_dict))

    return average_node


def get_histogram(degree_dict):
    histo = {}
    for key in degree_dict: #goes through all the keys (nodes) in the degree dictionary
        degree = degree_dict[key] #grabs the value of the key - the degree of the node 
        if degree not in histo: #checks if the degree is a key in the new dictionary
            histo[degree] = 1 #if not, add it and give it a value of 1
        else: #if it is, add one to its value
            histo[degree] += 1
    return histo


#Takes as input a histogram dictionary, which is a dictionary of the number of nodes with a certain degree
#{degree: # nodes with that degree}
#j is the probability threshold
#Output: a plot the histogram
def plot_histogram(histogram, j):
    fig = plt.figure(figsize=(6.5,4))
    x = [deg for deg in histogram] #x values are the degree
    y = [] #y values are the number of nodes with that degree
    for d in x:
        num_nodes = histogram[d]
        y.append(num_nodes)

    #Log plot values of the graph
    logx = [math.log(a) for a in x] 
    logy = [math.log(b) for b in y]

    plt.subplot(1,2,1)
    plt.plot(x,y,'or')
    plt.xlabel('Degree')
    plt.ylabel('# of Nodes')
    plt.title('Degree Histogram')

    plt.subplot(1,2,2)
    plt.plot(logx, logy, 'sb')
    plt.xlabel('Degree (log)')
    plt.ylabel('# of Nodes (log)')
    plt.title('Degree Histogram (log)')

    plt.tight_layout()

    plt.savefig('GIANT_degree_histogram_' + str(j) + '.png')
    print('Wrote to GIANT_degree_histogram_' + str(j) + '.png')
    return


def neighbor_dict(edge_list):
    #Takes as input an edge list and outputs a dictionary of nodes (key) and a list of their neighbor nodes (value)
    neighbors = {} #initializes empty dictionary
    for edge in edge_list: #looks at every edge in the edge list
        for i in [0,1]:
            neighbor_list = []
            index = -1*i + 1 #gives the index of the other node in the edge (0 if i=1, 1 if i=0)
            if edge[i] not in neighbors: #if the first node in the edge isn't already in the dictionary
                neighbor_list.append(edge[index])
                neighbors[edge[i]] = neighbor_list #create a key value pair where the value is a list with the other 
                #node in the edge 
                #time it appears in the edge list (so far, it only has one connection and thus a degree of one)
            elif i == 1 and edge[i] == edge[0]: #If an edge is between a node and itself, you don't want to add that node as 
            #a neighbor twice, so we do nothing 
                pass
            elif edge[i] in neighbors: #if it's already in the dictionary
                neighbors[edge[i]].append(edge[index]) #add one to its degree
    return neighbors


def average_neighbor_degree(degree_dict, neighbor_dict):
    #Takes as input the degree dictionary {node: degree} and dictionary of nodes' neighbors {node: [list of neighbors]}
    #And outputs a dictionary of {node: AND} key:value pairs
    AND_dict = {}
    for v in neighbor_dict: #iterates over the keys (proteins) in the neighbor dictionary
        Nv = neighbor_dict[v] #variable that holds the list of neighbors
        degree_sum = 0 #sum of the neighbors' degrees
        nb_sum = degree_dict[v] #number of neighbors (degree)
        for n in Nv: #iterates over the neighbors of the protein v
            degree = degree_dict[n] #finds the degree of that neighbor by looking in the degree dictionary
            degree_sum += degree #adds it to the sum of degrees
        AND = float(degree_sum)/float(nb_sum) #calculates average node degree (AND)
        AND_dict[v] = AND #adds the AND as the value whose key is the node
    return AND_dict


def degree_and_dictionary(degree_dict, and_dict, histogram_dict):
    #Takes as input the degree dictionary {node: degree} and 
    #average node degree dictionary {node: average neighbor degree} and
    #histogram dictionary {degree: number of nodes with that degree}
    #For all nodes with the same degree, it averages their AND values
    #and outputs a dictionary {degree: average AND value}
    average_AND = {}
    for protein in degree_dict:
        deg = degree_dict[protein] #protein's degree
        and_value = and_dict[protein] #protein's AND value
        if deg not in average_AND: #if that degree is not in the dictionary yet 
            average_AND[deg] = and_value #set a key with that protein's AND as the value
        else: #if the degree is already in the dictionary, we're going to add the AND value to it's current value
            average_AND[deg] += and_value

    #At this point, the values in the average_AND dictionary are just the sum of the AND values 
    #We need to divide by the number of nodes with that degree value in order to get the average AND value
    for key in average_AND:
        same_degree = histogram_dict[key] #looks up how many nodes have that degree by looking at histogram
        average_AND[key] = float(average_AND[key])/float(same_degree)
    return average_AND


#Input: a dictionary of the average AND values for each degree 
#keys are the x-values and the values (average AND) are the y-values
#j is the probability threshold
def plot_average_AND(average_AND, j):
    fig = plt.figure(figsize=(4.5,4.5))
    x = [deg for deg in average_AND] #x values are the degree
    y = [] #y values are the average AND values of each degreee
    for d in x:
        avg_and = average_AND[d]
        y.append(avg_and)


    plt.subplot(1,1,1)
    plt.plot(x,y,'*m')
    plt.xlim([0,100])
    plt.xlabel('Degree')
    plt.ylabel('Average AND')
    plt.title('Average AND of Nodes with Degree k') 

    plt.tight_layout()

    plt.savefig('GIANT_average_and_' + str(j) + '.png')
    print('Wrote to GIANT_average_and_' + str(j) + '.png')

    return


def create_node_list(edge_list):
    #Creates a node list from the edge_list
    nodes = []
    for edge in edge_list: #looks at every edge in the edge list
        for i in [0,1]:
            if edge[i] not in nodes: 
                nodes.append(edge[i])
    return nodes

#Takes as input a node list, adjacency list, and a number n which is the number of paths
#BFS should calculate up to and then stop. 
def BFS(node_list, adj_list, n):
    #Problem: Since every single node gets to be the source, the shortest path between
    #two nodes will be calculated twice, which will inflate the histogram values.
    #At the end, I just divided the pathlengths by 2 to eliminate this effect


    dist_histogram = [] #a list of all pathlengths to calculate the histogram
    p = 0 #p keeps track of the number of paths we've calculated 
    for s in node_list: #Modify the algorithm so that every node in the node list gets to be the source
        D = {} #dictionary of distances
        Q = [] #empty queue
        for node in node_list:
            D[node] = 1000000 #initialize all values with large integer
        D[s] = 0
        Q = [s]
        while len(Q) > 0:
            W = Q.pop(0) #dequeue 
            for neighbor in adj_list[W]: #look at neighbors of W given by adjacency list
                if D[neighbor] == 1000000: #if the distance is still the initialized value
                    D[neighbor] = D[W] + 1 #change the distance of the neighbor to be one more than its parent's distance
                    distance = D[neighbor]
                    dist_histogram.append(distance)
                    Q.append(neighbor) #add neighbor to end of queue
                    p += 1 #increment the number of paths found by 1 
                    if p > n: #if we've calculated n pathlengths, stop running 
                        return dist_histogram 

    return dist_histogram


#Takes as input a list of the first n pathlengths calculated on the graph 
#And the number n pathlengths that were calculated - this will be used to label the graph
#j is the probability threshold
def plot_pathlength(dist_histogram, j, n):   
    x = [d for d in dist_histogram]
    fig = plt.figure(figsize=(4,4))
    plt.subplot(1,1,1)
    plt.hist(x, 50, facecolor = 'green')
    plt.xlabel('Pathlength')
    plt.ylabel('# of Paths')
    plt.title('Pathlength Histogram of First ' + str(n) + ' Paths')
    plt.tight_layout()
    plt.savefig('GIANT_path_histogram_' + str(j) + '_' + str(n) + '.png')
    print('Wrote to GIANT_path_histogram_' + str(j) + '_' + str(n) + '.png')
    return


## Run Dijkstra's in the weighted, undirected graph.
## INPUT: set of nodes, 3-element list of edges [node1,node2,weight], source s
## OUTPUT: Dictionary of distances (D), Dictionary of predecessors (pi)
def dijkstra(nodes,edges,s):
    ## Build adjacency list that contains the weights of the edge.
    ## e.g., for edge (u,v), you can access the weight of that edge
    ## with adj_list[u][v] OR adj_list[v][u]
    adj_list = get_adj_list_with_weights(edges)

    LARGE_NUM = 1000000 ## like "infinity" here.

    ## initialize distances dictionary D.
    D = {n:LARGE_NUM for n in nodes}

    ## initialize predecessor dictionary pi.
    pi = {n:None for n in nodes}

    ## set distance to s to be 0
    D[s] = 0

    ## Queue is a dictionary (slow implementation)
    ## This could be sped up with a proper priority queue,
    ## but is fine for this homework.
    ## The queue values start as the distances for each node.
    Q = {n:D[n] for n in nodes} 
    
    while len(Q) > 0: ## While we haven't visited all the nodes...
        ## Find the node with the minimum weight.
        w = None 
        for n in Q: ## for every node in the Queue...
            if w == None or Q[n] < Q[w]: ## if we haven't set w yet or n is better...
                w = n ## set w to be this node.

        ## remove w from queue
        del Q[w] 
        
        ## Iterate through the neighbors of w
        for x in adj_list[w]:
            ## If the current distance to x is larger than coming from w, update
            if D[x] > D[w] + adj_list[w][x]:
                D[x] = D[w] + adj_list[w][x] ## update the distance
                pi[x] = w ## update the predecessor (we came from w)
                Q[x] = D[x] ## update the entry in the queue
                
    return D,pi

## Make an adjacency list that contains the weights of each edge.
## e.g., for edge (u,v), you can access the weight of that edge
## with adj_list[u][v] OR adj_list[v][u]
## INPUT: 3-element list of edges [node1,node2,weight]
## OUTPUT: dictionary of dictionaries
def get_adj_list_with_weights(edges):
    adj_list = {}
    ## loop over all edges.
    for u,v,w in edges: ## another way to specify elements of key

        ## We want to add the key-value pair (v,w) to adj_list[u].
        ## First see if u is a key in adj_list.
        if u not in adj_list:
            adj_list[u] = {}  ## add the key (value is a DICTIONARY)
        ## Add the key-value pair (v,w) to adj_list[u]
        adj_list[u][v] = w

        ## We want to add the key-value pair (u,w) to adj_list[v].
        ## First see if v is a key in adj_list.
        if v not in adj_list:
            adj_list[v] = {}  ## add the key (value is a DICTIONARY)
        ## Add the key-value pair (u,w) to adj_list[v]
        adj_list[v][u] = w
        
    return adj_list



if __name__ == '__main__':
    main()