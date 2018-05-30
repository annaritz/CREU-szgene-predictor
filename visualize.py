import sys
from optparse import OptionParser, OptionGroup

import networkx as nx
from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph 
import matplotlib.pyplot as plt
graphspace = GraphSpace('mirbern@reed.edu', 'pumpkin')


def parse_arguments(argv):

    usage = 'visualize.py [options]'
    parser = OptionParser(usage=usage)
    
    #Input Files
    group = OptionGroup(parser,'Input Files')
    group.add_option('-n', '--network',\
        type='string',metavar='STR',default='brain_top_geq_0.200.txt',\
        help='Functional interaction network (default 0.200 threshold)')
    
    #Method Arguments
    group = OptionGroup(parser,'Method Argument.')
    group.add_option('-g', '--gene',\
        type='string',metavar='STR',default='56895',\
        help='Gene and neighbors to visualize (default is AGPAT4 entrez id 56895)')

    #parse command line arguments
    (opts, args) = parser.parse_args()

    return opts

def main(argv):
    opts = parse_arguments(argv)
    print(opts)

    if opts.gene and opts.network:
        print('Reading gene edge file %s' % (opts.network))

        G = nx.Graph()

        graph_nodes = read_edge_file(opts.network, G)

        create_neg_set('SZnegatives.csv', G, graph_nodes)
        create_SZ_set('SZpositives.txt', G, graph_nodes)
        create_CM_set('mastergenelist.txt', G, graph_nodes)

        

        print('Visualizing local graph...')

        visualize_graph(graphspace, G, opts.gene, opts.network[14:19])

    return

def read_edge_file(brain_network, Graph):
    #Add nodes and weighted edges to the networkx graph
    nodeset = set()
    with open(brain_network, 'r') as network:
        for line in network:
            line = line.strip().split('\t')
            Graph.add_edge(line[0], line[1], weight=float(line[2]))
            for i in range(0,2):
                node = line[i]
                if node not in nodeset: #checks for duplicates and adds nodes to the graph
                    Graph.add_node(node, SZ_pos=False, CM_pos=False, Neg=False)
                    nodeset.add(node)

    return nodeset

#Three different functions to handle the three different input files

def create_neg_set(file, Graph, Gnodes): #Reads negative file
    count = 0
    with open(file, 'r') as f:
        for line in f:
            line=line.strip().split(',')
            entrezID = line[0]
            if entrezID == 'gene id':
                pass
            else:
                if entrezID in Gnodes:
                    Graph.nodes[entrezID]['Neg'] = True
                    count += 1

    print('There are %d negatives in the network.' % count)

    return 

def create_SZ_set(file, Graph, Gnodes): #Reads SZ positive file
    count = 0
    with open(file, 'r') as f:
        for line in f:
            line=line.strip().split('\t')
            entrezID = line[0]
            if entrezID in Gnodes: #checks if in Graph before searching 
                Graph.nodes[entrezID]['SZ_pos'] = True
                count += 1

    print('There are %d SZ positives in the network.' % count)

    return 


def create_CM_set(file, Graph, Gnodes): #Reads CM positive file 
    count = 0
    with open(file, 'r') as f:
        for line in f:
            line=line.strip().split('\t')
            entrezID = line[1]
            if entrezID == 'Gene':
                pass
            else:
                if entrezID in Gnodes:
                    Graph.nodes[entrezID]['CM_pos'] = True
                    count += 1

    print('There are %s CM positives in the network.' % count)

    return 


def visualize_graph(graphspace, Graph, gene, threshold):
    if gene not in Graph.nodes():
        sys.exit('ERROR: %s is not in the network. It is possible that the probability threshold is too high.' % gene)

    weight_sum = 0
    neighbors = 0
    name = 'candidate_neighbors_' + str(gene) + '_' + str(threshold) + '.txt'
    with open(name, 'w') as c:
        c.write('EntrezID\tSZPos\tCMPos\tSZNeg\tEdgeWeight\n') #Title
        #For every neighbor, write a line of info to output file
        for neighbor,datadict in Graph.adj[gene].items():
            weight_sum += datadict['weight']
            neighbors += 1
            c.write(neighbor + '\t' + str(Graph.nodes[neighbor]['SZ_pos']) + '\t' + str(Graph.nodes[neighbor]['CM_pos']) + '\t' + str(Graph.nodes[neighbor]['Neg']) + '\t' + str(datadict['weight']) + '\n')
   
        average_weight = float(weight_sum/neighbors)

        #Summary line
        c.write('Entrez ID %s has %d neighbor(s) with an average edge weight score of %f' % (gene, neighbors, average_weight))


    """
    local_graph = GSGraph()
    local_graph.set_name('Local Graph of %s in %s' % (gene, threshold))
    local_graph.set_tags(['CREU'])
    local_graph.add_node(gene, label=gene)

    #Graph.adj[gene].items() gives a list of tuples. Each tuple includes one of the selected node's neighbors
    #and a dictionary of the attributes their shared edge has, i.e. weight
    #neighbor is the gene's neighbor, datadict is the dictionary of attributes
    for neighbor,datadict in Graph.adj[gene].items():
        #adds the gene and all of its neighbors to the graphspace graph
        local_graph.add_node(neighbor, label=neighbor)
        local_graph.add_edge(gene,neighbor)
        local_graph.add_edge_style(gene,neighbor, width=datadict['weight'])

    graphspace.post_graph(local_graph)
    print('Graph posted')
    """

    return






if __name__ == '__main__':
    main(sys.argv)



            




