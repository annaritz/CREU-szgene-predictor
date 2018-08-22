import sys
from optparse import OptionParser, OptionGroup

from datetime import datetime
from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph 
import matplotlib.pyplot as plt
#graphspace = GraphSpace('mirbern@reed.edu', 'pumpkin')
#graphspace = GraphSpace('aritz@reed.edu','platypus')
username = 'aritz@reed.edu'
password='platypus'
group='CREU'


def parse_arguments(argv):

    usage = 'visualize.py [options]'
    parser = OptionParser(usage=usage)
    
    #Input Files
    group = OptionGroup(parser,'Input Files')
    group.add_option('-n', '--network',\
        type='string',metavar='STR',default='../infiles/networkfiles/brain_top_geq_0.150.txt',\
        help='Functional interaction network (default 0.150 threshold)')
    
    #Method Arguments
    group = OptionGroup(parser,'Method Argument.')
    group.add_option('-g', '--gene',\
        type='string',metavar='STR',default='56895',\
        help='Gene and neighbors to visualize (default is AGPAT4 entrez id 56895)')

    group = OptionGroup(parser, 'Method Argument.')
    group.add_option('-m', '--motility_positive',\
        action='store_true',default=False,\
        help='Indicate if candidate is a cell motility positive (default=False)')

    group = OptionGroup(parser, 'Method Argument.')
    group.add_option('-s', '--SZ_positive',\
        action='store_true',default=False,\
        help='Indicate if candidate is a SZ positive (default=False')

    group = OptionGroup(parser, 'Method Argument.')
    group.add_option('-a', '--all_neighbors',\
        action='store_true',default=False,\
        help='Plot all neighbors instead of labeled neighbors (default=False)')

    #parse command line arguments
    (opts, args) = parser.parse_args()

    return opts

def main(argv):
    opts = parse_arguments(argv)
    print(opts)

    if opts.gene and opts.network:
        print('Reading gene edge file %s' % (opts.network))

        genemap = geneMapReader('../infiles/Homo_sapiens.txt')
        cand_neighbors = read_edge_file(opts.network, opts.gene)

        neg_neighbors = create_set('../infiles/SZ_negatives.txt', cand_neighbors, opts.gene)
        SZ_pos_neighbors = create_set('../infiles/SZ_positives.txt', cand_neighbors, opts.gene)
        CM_pos_neighbors = create_set('../infiles/motility_positives.txt', cand_neighbors, opts.gene)

        print('Visualizing local graph...')

        visualize_graph(opts.all_neighbors,opts.gene, opts.motility_positive, opts.SZ_positive, opts.network[-9:-4], cand_neighbors, neg_neighbors, SZ_pos_neighbors, CM_pos_neighbors, genemap)

        print('Done!')

    return

def read_edge_file(brain_network, candidate):
    adjacent_nodes = {} #initializes dictionary of neighbors:weights
    with open(brain_network, 'r') as network:
        for line in network:
            line = line.strip().split('\t')
            weight = line[2]
            if line[0] == candidate:
                adjacent_nodes[line[1]] = weight
            elif line[1] == candidate:
                adjacent_nodes[line[0]] = weight

    print(adjacent_nodes)
                
    return adjacent_nodes


#Reads in positive or negative file and returns set of neighbors that are in that set as well as two booleans
#Booleans indicate whether candidate is a positive
def create_set(file, neighbors, candidate): 
    SZ_positive = False
    CM_positive = False

    adj_set = set() 
    with open(file, 'r') as f:
        for line in f:
            line=line.strip().split('\t')
            entrezID = line[0]
            if entrezID in neighbors:
                adj_set.add(entrezID)

    print(file, 'Neighbors:\n', adj_set)

    return adj_set


def geneMapReader(infile):
    GeneMapfile=open(infile, 'r')
    genemap = {} #dictionary where each gene has two key:value entries {entrezID:symbol, symbol:entrezID}
    for i in GeneMapfile:
        i=i.split('\t')
        ientrez = i[1]
        iname=i[2]
        if iname!='Approved Symbol': #skips title 
            genemap[iname]=ientrez
            genemap[ientrez]=iname
    return genemap


def visualize_graph(all_neighbors,candidate, CM_positive, SZ_positive, threshold, neighbors, adj_neg, adj_SZ_pos, adj_CM_pos, genemap):
    #Red if negative, blue if CM positive, yellow if SZ positive, green if double positive, orange if unlabeled
    #TODO: make a color dictionary

    double_pos = adj_CM_pos.intersection(adj_SZ_pos) #get set of nodes in both positive sets
    print('Double pos:', double_pos)

    neg_pos = adj_neg.intersection(adj_SZ_pos).union(adj_neg.intersection(adj_CM_pos)) #get set of nodes that are both positives and negatives
    print('Neg and pos:', neg_pos)

    #Create the candidate node + node style
    symbol = genemap.get(candidate)

    local_graph = GSGraph()
    time = str(datetime.now())
    if all_neighbors:
        local_graph.set_name('All Neighbors %s (%s)' % (symbol,candidate))
    else:
        local_graph.set_name('Labeled Neighbors %s (%s)' % (symbol,candidate))
    local_graph.set_tags(['CREU'])

    
    if CM_positive:
        popup = '%s, %s, %s' % (candidate, 'Candidate', 'Cell Motility positive')
        local_graph.add_node(candidate, label=symbol, popup=popup)
        local_graph.add_node_style(candidate, color='#6699ff', height=40, width=100)
    elif SZ_positive:
        popup = '%s, %s, %s' % (candidate, 'Candidate', 'Schizophrenia positive')
        local_graph.add_node(candidate, label=symbol, popup=popup)
        local_graph.add_node_style(candidate, color='#fffc66', height=40, width=100)
    else:
        popup = '%s, %s' % (candidate, 'Candidate')
        local_graph.add_node(candidate, label=symbol, popup=popup)
        local_graph.add_node_style(candidate, color='#affef2', height=40, width=100)

    #Create the neighbor nodes + styles

    ## Add legend nodes
    local_graph.add_node('leg1', popup='SZ & CM Pos', label='SZ & CM Positive')
    local_graph.add_node_style('leg1', color='#5de83e', height=30, width=150)
    local_graph.add_node('leg2', popup='Unlabeled', label='Unlabeled')
    local_graph.add_node_style('leg2', color='#affef2', height=30, width=90)
    local_graph.add_node('leg3', popup='Negative', label='Negative')
    local_graph.add_node_style('leg3', color='#e56b52', height=30, width=90)     
    local_graph.add_node('leg4', popup='CM Positive', label='CM Positive')
    local_graph.add_node_style('leg4', color='#6699ff', height=30, width=100)
    local_graph.add_node('leg5', popup='SZ Positive', label='SZ Positive')
    local_graph.add_node_style('leg5', color='#fffc66', height=30, width=100)  

    for neighbor in neighbors:

        symbol = genemap.get(neighbor)
        
        if neighbor in neg_pos: #check if it's also a negative
            print('Neighbor %s found in both positive and negative set' % neighbor)
            continue

        elif neighbor in double_pos: #check if a neighbor is both a cell motility and SZ positive
            print('Neighbor %s found in both positive sets' % neighbor)
            popup = '%s, %s' % (neighbor, 'SZ and CM positive')
            local_graph.add_node(neighbor, popup=popup, label=symbol)
            local_graph.add_node_style(neighbor, color='#5de83e', height=30, width=90)

        elif neighbor in adj_neg:
            print('Negative')
            popup = '%s, %s' % (neighbor, 'Negative')
            local_graph.add_node(neighbor, popup=popup, label=symbol)
            local_graph.add_node_style(neighbor, color='#e56b52', height=30, width=90)              

        elif neighbor in adj_CM_pos:
            print('CM positive')
            popup = '%s, %s' % (neighbor, 'Cell Motility Positive')
            local_graph.add_node(neighbor, popup=popup, label=symbol)
            local_graph.add_node_style(neighbor, color='#6699ff', height=30, width=90)

        elif neighbor in adj_SZ_pos:
            print('SZ positive')
            popup = '%s, %s' % (neighbor, 'SZ Positive')
            local_graph.add_node(neighbor, popup=popup, label=symbol)
            local_graph.add_node_style(neighbor, color='#fffc66', height=30, width=90)  

        elif all_neighbors:
            local_graph.add_node(neighbor, popup=symbol, label='')
            local_graph.add_node_style(neighbor, color='#affef2', height=20, width=20)

        
        local_graph.add_edge(candidate,neighbor)
        local_graph.add_edge_style(candidate,neighbor, width=neighbors[neighbor]*50)

    gid = post(local_graph,username,password,group)
    print('Graph posted',gid)
    

    return



def post(G,username,password,group=None):
  '''
  Post a graph to graphspace.  If this graph is new, share it with
  the group.
  Inputs: Graph (GSGraph Object)
  Outputs: the graph ID (int)
  '''
  print('Posting graph...')
  
  # connect to GraphSpace with the username and password.
  gs = GraphSpace(username,password)
  try:
    # try updating the graph. If the graph does not exist, 
    # this will throw an error.
    graph = gs.update_graph(G)
  except:
    # catch the error and try posting a new graph.
    graph = gs.post_graph(G)
    if group:
        # share this graph with the group.
        gs.share_graph(graph_id=graph.id,group_name=group)
  return graph.id




if __name__ == '__main__':
    main(sys.argv)



            




