#Miriam Bern
#HW4
#October 9, 2017

import random
import matplotlib.pyplot as plt
import math
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
from datetime import datetime

graphspace = GraphSpace('mirbern@reed.edu', 'pumpkin')

def main():
    random.seed(500)
    filename = ['test-edges.txt', 'EGFR1-edges.txt']
    q = [0.85, 0.95]
    t_sim = [1000, 100000] #timesteps for simulation
    t_prob = [100, 1000] #timesteps for probability table
    s = ['A', 'EGF']

    #Runs relevant code on test-edges.txt
    node_list, edge_list = read_edges(filename[0])

    out_dict = create_out_dict(edge_list, node_list)

    in_dict = create_in_dict(edge_list, node_list)

    out_degree = create_out_degree(out_dict)

    add_edges(node_list, out_dict, s[0]) 

    prob_table = rw_probs(node_list, in_dict, out_degree, s[0], q[0], t_prob[0])

    walker_dict = rw_simulate(node_list, out_dict, t_sim[0], s[0], q[0])

    #Runs relevant code on EGFR1-edges.txt
    node_list, edge_list, edge_type_dict = read_edges_with_edge_type(filename[1])

    out_dict = create_out_dict(edge_list, node_list)

    in_dict = create_in_dict(edge_list, node_list)

    out_degree = create_out_degree(out_dict)

    add_edges(node_list, out_dict, s[1]) 

    prob_table = rw_probs(node_list, in_dict, out_degree, s[1], q[1], t_prob[1])

    walker_dict = rw_simulate(node_list, out_dict, t_sim[1], s[1], q[1])

    plot_network(graphspace, node_list, edge_list, walker_dict, edge_type_dict)

    return

#Reads the contents of a 2-column edge file and outputs a list of nodes and list of edges
def read_edges(filename):
    nodes = set()
    edges = []
    with open(filename, 'r') as fin:
        for line in fin:
            row = line.strip().split('\t')
            nodes.add(row[0])
            nodes.add(row[1])
            edges.append([row[0],row[1]])
    print(len(nodes), 'nodes and', len(edges), 'edges')
    nodes = sorted(list(nodes))
    return nodes, edges

def read_edges_with_edge_type(filename):
    edge_type = {} #dictionary of each edge (key) and its edge type (value)
    nodes = set()
    edges = []
    index = 0 #index of an edge in the edge list
    with open(filename, 'r') as fin:
        for line in fin:
            row = line.strip().split('\t')
            if row[0] == '#Protein1': #Maybe not the best solution, but this skips over the first line
            #in the file that labels each column - prevents it from being added to edge list, etc.
                pass
            else:
                nodes.add(row[0])
                nodes.add(row[1])
                edges.append([row[0],row[1]])
                edge_type[index] = row[2]
                index += 1
    print(len(nodes), 'nodes and', len(edges), 'edges')
    nodes = sorted(list(nodes))
    print(edge_type)
    return nodes, edges, edge_type

#Creates a dictionary whose keys are nodes and whose values are a list of "out" neighbors
def create_out_dict(edge_list, node_list): 
    out_d = {}
    for edge in edge_list:
        node = edge[0]
        out_node = edge[1]
        if node not in out_d:
            out_d[node] = [out_node]
        else:
            out_d[node].append(out_node)

    #Because of the way the dictionary is initially constructed, only nodes that have out neighbors will be in the 
    #dictionary. This for loop ensures all nodes are in the dictionary
    for node in node_list:
        if node not in out_d:
            out_d[node] = []

    return out_d

#Creates a dictionary whose keys are nodes and whose values are a list of "in" neighbors
def create_in_dict(edge_list, node_list):
    in_d = {}
    for edge in edge_list:
        node = edge[1]
        in_node = edge[0]
        if node not in in_d:
            in_d[node] = [in_node]
        else:
            in_d[node].append(in_node)


    #Because of the way the dictionary is initially constructed, only nodes that have in neighbors will be in the 
    #dictionary. This for loop ensures all nodes are in the dictionary
    for node in node_list:
        if node not in in_d:
            in_d[node] = []

    return in_d

#Uses the out degree dictionary to create a dictionary of the number of out neighbors (degree) a node has
def create_out_degree(out_d):
    out_deg = {}
    for node in out_d:
        out_deg[node] = len(out_d[node])
    return out_deg


#Checks the dictionary of "out" neighbors and if a node has none, it adds the source as an out neighbor
#Inputs dictionary of outgoing neighbors and the source node and returns nothing
def add_edges(node_list, out_d, s):
    for node in node_list:
        out_list = out_d[node]
        if len(out_list) == 0:
            out_list.append(s) #add the source node as an outgoing neighbor
    return    



def rw_probs(node_list, in_d, out_deg, s, q, t_prob):
    node_list = node_list
    N_in = in_d
    out_deg = out_deg

    prob_dict = {}
    #Initializes a probability dictionary with the first probability in the list
    for node in node_list:
        if node == s:
            prob_dict[node] = [1]
        else:
            prob_dict[node] = [0]
    
    for i in range(1, t_prob):
        for node in node_list:
            n_sum = 0
            for u in N_in[node]: #for every neighbor coming into the node
                deg_out = out_deg[u] #the number of neighbors u points to (come out of u)
                prev_prob = prob_dict[u][i-1] #probability of u at t-1
                add_to_sum = float(prev_prob) / float(deg_out)
                n_sum += add_to_sum
            rand_surfer = (1-q) * prob_dict[node][0]
            prob = (q*n_sum) + rand_surfer #probability equation at time t
            prob_dict[node].append(prob)

    prob_table = []
    #Turns our probability dictionary into a list of lists
    for n in node_list:
        prob_table.append(prob_dict[n])

    return prob_table

#Inputs a node list, dictionary of outgoing neighbors for each node, time steps t, source s, probability q
#Outputs a dictionary of each node and number of times they've been visited
def rw_simulate(node_list, out_dict, t_sim, s, q):
    walk_dict = {node: 0 for node in node_list} #initializes a dictionary of key value pairs {node:# of times visited}
    walk_dict[s] += 1 #at time step 0, source node s has been visited so its visit count is incremented by 1
    current_node = s #keeps track of node we are currently at
    for time in range(1, t_sim): #at each time step
        out_neighbors = out_dict[current_node] #list of the current node's outgoing neighbors
        rand = random.random() #random number 0 to 1 generated
        if rand <= q: #if random number is less than or equal to probability that an outgoing neighbor will be chosen
            next_node = random.choice(out_neighbors) #the next node to visit is randomly chosen from current node's outgoing neighbors
        else: #if random number is greater than q
            next_node = s #go back to the source
        walk_dict[next_node] += 1 #increment the next node's visit count by 1
        current_node = next_node #set the current node to the next node and move on to the next time step

    return walk_dict #return the dictionary of nodes and number of times they've been visited


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
def plot_network(graphspace, node_list, edge_list, walk_dict, type_dict):
    EGFR = GSGraph()
    EGFR.set_name('Random Walk Along EGFR Network' + ' ' + str(datetime.now()))
    EGFR.set_tags(['HW4'])

    #The for loop below adds 1 to each visit count and then takes the log
    for node in walk_dict:
        walk_dict[node] += 1 #add one to the visit count to ensure log(0) doesn't happen
        walk_dict[node] = math.log(walk_dict[node])

    #The for loop below computes the max log(count) in order to normalize the counts
    count_max = 0
    for node in walk_dict:
        if walk_dict[node] > count_max:
            count_max = walk_dict[node]

    #The for loop below normalizes the counts
    for node in walk_dict:
        count = walk_dict[node]
        count = float(count) / float(count_max)

    #The code below colors the nodes according to how many times they've been visited
    for node in node_list:
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
            edge_color = 'green's

        EGFR.add_edge_style(edge[0], edge[1], directed=True, edge_style='dotted', color=edge_color)

    graphspace.post_graph(EGFR)

    return

if __name__ == '__main__':
    main()







