    
import networkx as nx
import matplotlib.pyplot as plt
import math
import numpy
import glob
import os




def main():
    #Top edges of GIANT brain network downloaded from http://hb.flatironinstitute.org/download
    #Three column plain text file [entrez gene id 1][entrez gene id 2][posterior prob.]
    
    network_files = glob.glob('../source/infiles/networkfiles/brain_top_geq_*') 

    #List of lists containing the first 200,000 path lengths of each trimmed network
    #This will be used to plot the distributions on the same plot using plot_pathlength function 
    pathlength_list = []

    for file in network_files:

        print('READING FILE: ', file)

        #Read the edge file - this function creates a networkx graph from our 3-column file of 
        #weighted edges 
        G = nx.read_weighted_edgelist(file)
        
        node_list = list(G.nodes())

        edge_list = read_edge_file(file)

        N_v = neighbor_dict(edge_list)

        path_list_200 = BFS(node_list, N_v, 200000)

        pathlength_list.append(path_list_200)
   
    plot_pathlength(pathlength_list, 200000)


    return


def read_edge_file(edge_file):
    with open(edge_file, "r") as file:
        edges_list = []
        for line in file: 
            node_pair = line.strip() #strips whitespace characters
            edge = node_pair.split("\t") #splits at the tabs (there are tabs between each column)
            edges_list.append(edge)
    return edges_list



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


#Takes as input a list of lists of the first n pathlengths calculated on the trimmed networks
#And the number n pathlengths that were calculated - this will be used to label the graph
#Plots all the path length distributions of the networks on one plot

def plot_pathlength(dist_histogram, n):   
    fig = plt.figure(figsize=(6.5,4))
    face_color = ['green', 'red', 'cyan', 'blue', 'magenta', 'yellow', 'pink', 'black', 'olive', 'orange']
    
    plt.subplot(1,1,1)
    for i in range(len(dist_histogram)):
        x = dist_histogram[i]
        plt.hist(x, 50, facecolor = face_color[i])

    plt.xlabel('Pathlength')
    plt.ylabel('# of Paths')
    plt.title('Path length Histogram of First ' + str(n) + ' Paths')
    plt.tight_layout()
    plt.savefig('GIANT_path_all_networks' + '_' + str(n) + '.png')
    print('Wrote to GIANT_path_all_networks' + '_' + str(n) + '.png')
    return



if __name__ == '__main__':
    main()