#Contains the functions that handle the edge files and write the output files
#Fixe partitionCurated

import random
import matplotlib.pyplot as plt
import networkx as nx 


def geneMapReader(infile):
    GeneMapfile=open(infile, 'r')
    genemap = {}
    for i in GeneMapfile:
        i=i.split('\t')
        ientrez = i[1]
        iname=i[2]
        if iname!='Approved Symbol':
            genemap[iname]=ientrez
            genemap[ientrez]=iname
    return genemap


def read_edge_file_multi(graph, layers):
    '''
    Inputs:
    Outputs:
    Note: Graph is modified in place (since it's a pointer to the networkx object) and isn't explicitly returned.
    '''
    multi_layer_dict = {} #dictionary that maps from original entrez ID to multi_layer/prime name {node:(node_prime, node_0, node_1)}
    seen_nodes = set() #keeps track of which nodes we've already seen in the edges

    print('Creating multi-layer network...')

    single_layer_edges = list(graph.edges()) #turns the single layer graph edges into a list 
    #this prevents it from being updated while trying to iterate the original edges
    weight_dict = nx.get_edge_attributes(graph,'weight') #dictionary of edges and weights

    for edge in single_layer_edges:
        for i in range(0,2):  ## iterates over the 2 nodes that are in each edge
            node = edge[i]
            if node not in seen_nodes: # if the node hasn't been seen yet...
                node_names = set() # node_names will contain all the copies of this node.
                # create and add the prime node
                node_names.add(node+'_prime')
                graph.add_node(node+'_prime', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=layers)
                for layer in range(layers): 
                    # for every layer, create and add the layer node.
                    node_names.add(node+'_'+str(layer))
                    graph.add_node(node+'_'+str(layer), prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=1)
                    # also add the edge from this layer's node to the prime node.
                    graph.add_edge(node+'_'+str(layer), node+'_prime', weight=1.0)

                # update seen_nodes and the node dictionary
                seen_nodes.add(node)
                multi_layer_dict[node] = node_names

        for layer in range(layers):
            # generate new node names for this layer
            node0 = edge[0]+'_'+str(layer)
            node1 = edge[1]+'_'+str(layer)
            # add the edge with the original edge weight
            graph.add_edge(node0, node1, weight=weight_dict[edge])

            # add the edge weights to the weighted degree attr.
            graph.nodes[node0]['weighted_degree'] += weight_dict[edge] #get edge weight 
            graph.nodes[node1]['weighted_degree'] += weight_dict[edge]
            

    #At the end, remove the old nodes, which have now been replaced with their layer duplicates + prime nodes
    #Edges are removed when nodes are removed
    ## TODO Suggestion: remove original nodes (from graph.nodes() before you add layer stuff)
    for old_node in seen_nodes:
        graph.remove_node(old_node)
    
    print('The multi-layer network contains %d edges and %d nodes' % (graph.number_of_edges(), graph.number_of_nodes()))

    return multi_layer_dict

def read_edge_file_single(filename, graph):
    all_nodes = set()
    with open(filename,'r') as fin:
        for line in fin:
            line=line.split('\t')
            line[2]=float(line[2][:len(line[2])-2]) #Wrestling with the formatting
            for i in range(0,2):
                node=line[i]
                if node not in all_nodes:
                    graph.add_node(node, prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=0)
                    all_nodes.add(node)
            graph.add_edge(line[0],line[1], weight=line[2])
            graph.nodes[line[0]]['weighted_degree'] += line[2]
            graph.nodes[line[1]]['weighted_degree'] += line[2]


# Handles reading curated files when there is one layer and the initial reading when there is more than one layer
# If # layers > 1, it'll read this, then check for overlap between positives and negatives
def curatedFileReader(filename,graph,verbose):
    #Issue: The multi graph is made so none of the original positives/negatives are in the graph 
    #Figure out a way to make the single graph, then make a multi graph and use that one instead
    #Note: Graph.nodes() is a list of all the nodes
    nodes = graph.nodes()
    tot = 0 #keeps track of number of nodes from curated set
    count = 0 #keeps track of number of nodes actually in the graph 
    curated = set()
    with open(filename,'r') as fin:
        for line in fin:
            entrezNumber = line.strip()
            tot+=1
            if entrezNumber in nodes:
                count+=1
                curated.add(entrezNumber)
            else:
                if verbose:
                    print('WARNING: EntrezID %s is not in graph.' % (entrezNumber))

    print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
    return curated



# The curated file has already been read and is taken as input here as a set 
# read_edge_file has already created multiple layers, this checks if the nodes are in the graph 
# partitions the positives/negatives and returns that set
def partitionCurated(original_curated,graph,verbose,layers):
    #Note: Graph.nodes() is a list of all the nodes
    nodes = graph.nodes()
    curated = set()

    print('Partitioning curated set...')

    #I'm questioning if we need this check - 
    #There is already a check in predict.py to make sure all original curated positives/negatives are in the graph,
    # which is used to create the multi-layer network
    #for entrezNumber in original_curated:
    #    for layer in range(layers):
    #        if entrezNumber+'_'+str(layer) in nodes:
    #            labeled_set.add(entrezNumber)
    #        else:
    #            if verbose:
    #                print('ERROR: EntrezID %s is not in graph.' % (entrezNumber))
    #            else:
    #                continue
                # TODO suggestion: if this id is not in the graph - quit.
                #sys.exit()

 
    #We need the positives to be randomly distributed throughout the layers. If we have multiple
    labeled_List=list(original_curated) #labeled nodes surrounding a prime node, that effectively blocks all score transmission

    ## "shuffle" the list randomly. random.sample() sample without replacement. TODO: confirm & put the link. 
    labeled_List=random.sample(labeled_List, k=len(labeled_List))

    ## Go through each curated positive
    for i in range(len(labeled_List)):
        ##Pick a layer.  Uses integer division (e.g., layer will be 0,1,2,0,1,2,etc. for 3 layers)
        layer=(i//(len(labeled_List)//layers))

        ## Once a layer is specified, append this to the curated positive node name.
        ## This "places" the positive in a particular layer.
        curated.add(labeled_List[i]+'_'+str(layer))

    #print('%d of %d nodes are in graph from file %s' % (len(labeled_set),len(curated_set),filename))
    return curated

#Formats results as a LaTeX table
#outfile variable is the name of the outfile, partially specified by --outprefix 
def formatCombinedResults(G,outfile,d_predictions,b_predictions,\
    disease_positives,biological_process_positives,negatives,blacklist,genemap, layers,normed=False):
    if normed: # if true, then normalize the d_preds and b_preds to be between 0 and 1.
        d_max = max([d_predictions[key] for key in d_predictions if key[-6:]=='_prime' or layers==1])
        d_predictions = {key:d_predictions[key]/float(d_max) for key in d_predictions if key[-6:]=='_prime'  or layers==1}
        b_max = max([b_predictions[key] for key in b_predictions if key[-6:]=='_prime' or layers==1])
        b_predictions = {key:b_predictions[key]/float(b_max) for key in b_predictions if key[-6:]=='_prime'  or layers==1}
    # write output
    out = open(outfile,'w')
    out.write('\\begin{table}[h]\n')
    out.write('\\centering\n')
    out.write('\\begin{tabular}{|llc|ccc|}\\hline\n')
    out.write('Gene Name & EntrezID & Deg & $f_{\mathcal{D}}$ & $f_{\mathcal{P}}$ & Score $g(v)$\\\\\\hline\n')
    if layers>1:
        nodeset = set([n for n in G.nodes() if n[-6:]=='_prime'])
    else:
        nodeset = G.nodes()

    for n in sorted(nodeset, key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):

        if layers==1 or n[-6:] == '_prime':
            score = d_predictions[n]*b_predictions[n]
            if score <= 0.25:
                continue
            if n[-6:] == '_prime':
                nname=n[:-6]
            else:
                nname = n
            name = genemap.get(nname,nname)
            entrezID = nname
            if layers > 1:
                deg = G.degree(nname+'_0')-1
            else:
                deg = G.degree(n)
            out.write('%s & %s & %d' % (name,entrezID,deg))
            names = [nname] + [nname+'_'+str(x) for x in range(layers)]
            #print(names)
            if any([x in disease_positives for x in names]):
                out.write(' & \\textit{%.2f}' % (d_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (d_predictions[n]))
            if any([x in biological_process_positives for x in names]):
                out.write(' & \\textit{%.2f}' % (b_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (b_predictions[n]))
            out.write(' & %.2f\\\\\n' % (score))
    out.write('\\end{tabular}\n')
    out.write('\\end{table}\n')
    print('Wrote to %s' % (outfile))
    return

#Formats results as a LaTeX table
#outfile variable is the name of the outfile, partially specified by --outprefix 
def formatCombinedResults_unlabeled(G,outfile,d_predictions,b_predictions,\
    disease_positives,biological_process_positives,negatives,blacklist,genemap,layers,normed=False,union_predictions=False):
    if normed: # if true, then normalize the d_preds and b_preds to be between 0 and 1.
        d_max = max([d_predictions[key] for key in d_predictions if key[-6:]=='_prime' or layers==1])
        d_predictions = {key:d_predictions[key]/float(d_max) for key in d_predictions if key[-6:]=='_prime' or layers==1}
        b_max = max([b_predictions[key] for key in b_predictions if key[-6:]=='_prime' or layers==1])
        b_predictions = {key:b_predictions[key]/float(b_max) for key in b_predictions if key[-6:]=='_prime' or layers==1}
        if union_predictions:
            u_max = max([union_predictions[key] for key in union_predictions if key[-6:]=='_prime' or layers==1])
            union_predictions = {key:union_predictions[key]/float(u_max) for key in union_predictions if key[-6:]=='_prime' or layers==1}

    # write output
    out = open(outfile,'w')
    out.write('\\begin{table}[h]\n')
    out.write('\\centering\n')
    if union_predictions:
        out.write('\\begin{tabular}{|llcc|ccc|c|}\\hline\n')
        out.write('Gene Name & EntrezID & Rank & Deg & $f_{\mathcal{D}}$ & $f_{\mathcal{P}}$ & $f_{\mathcal{D}\cup\mathcal{P}}$ & $g(v)$\\\\\\hline\n')
    else:
        out.write('\\begin{tabular}{|llcc|ccc|}\\hline\n')
        out.write('Gene Name & EntrezID & Rank & Deg & $f_{\mathcal{D}}$ & $f_{\mathcal{P}}$ & Score $g(v)$\\\\\\hline\n')
    if layers>1:
        nodeset = set([n for n in G.nodes() if n[-6:]=='_prime'])
    else:
        nodeset = G.nodes()
    rank = 0
    for n in sorted(nodeset, key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        rank+=1
        if layers==1 or n[-6:] == '_prime':
            score = d_predictions[n]*b_predictions[n]
            if score <= 0.25:
                continue
            if n[-6:] == '_prime':
                nname=n[:-6]
            else:
                nname = n
            name = genemap.get(nname,nname)
            entrezID = nname
            if layers > 1:
                deg = G.degree(nname+'_0')-1
            else:
                deg = G.degree(n)

            ## NEW: ignore if node is BOTH a disease and a process positive.
            names = [nname] + [nname+'_'+str(x) for x in range(layers)]
            if any([x in disease_positives for x in names]) and any([x in biological_process_positives for x in names]):
                #print('SKIPPING',nname,name,any([x in disease_positives for x in names]),any([x in biological_process_positives for x in names]))
                continue
            out.write('%s & %s & %d & %d' % (name,entrezID,rank,deg))
            
            #print(names)
            if any([x in disease_positives for x in names]):
                out.write(' & \\textit{%.2f}' % (d_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (d_predictions[n]))

            if any([x in biological_process_positives for x in names]):
                out.write(' & \\textit{%.2f}' % (b_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (b_predictions[n]))
            if union_predictions:
                out.write(' & %.2f & %.2f\\\\\n' % (union_predictions[n],score))
            else:
                out.write(' & %.2f\\\\\n' % (score))
    out.write('\\end{tabular}\n')
    out.write('\\end{table}\n')
    print('Wrote to %s' % (outfile))
    return

# Write output file of nodes and process score * disease score 
def writeCombinedResults(G,outfile,d_predictions,b_predictions,\
    disease_positives,biological_process_positives,negatives,blacklist,genemap, layers, normed=False,union_predictions=False):
    if normed: # if true, then normalize the d_preds and b_preds to be between 0 and 1.
        d_max = max([d_predictions[key] for key in d_predictions if key[-6:]=='_prime' or layers==1])
        d_predictions = {key:d_predictions[key]/float(d_max) for key in d_predictions if key[-6:]=='_prime'  or layers==1}
        b_max = max([b_predictions[key] for key in b_predictions if key[-6:]=='_prime' or layers==1])
        b_predictions = {key:b_predictions[key]/float(b_max) for key in b_predictions if key[-6:]=='_prime'  or layers==1}
        if union_predictions:
            u_max = max([union_predictions[key] for key in union_predictions if key[-6:]=='_prime' or layers==1])
            union_predictions = {key:union_predictions[key]/float(u_max) for key in union_predictions if key[-6:]=='_prime' or layers==1}
    # write output
    out = open(outfile,'w')
    fig = plt.figure(figsize=(4,4))
    if union_predictions:
        out.write('#EntrezID\tName\tDisLabel\tDisScore\tProcLabel\tProcScore\tCombined\tUNION\tConflict?\tDegree\n')
    else:
        out.write('#EntrezID\tName\tDisLabel\tDisScore\tProcLabel\tProcScore\tCombined\tConflict?\tDegree\n')
    degreeList=[]
    if layers>1:
        nodeset = set([n for n in G.nodes() if n[-6:]=='_prime'])
    else:
        nodeset = G.nodes()
    for n in sorted(nodeset, key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        #print('writing out ',n)
        disLabel='Unlabeled'
        procLabel='Unlabeled'
        if layers > 1:
            for layer in range(layers):
                layer_n = n[:-6]+'_'+str(layer)
                if layer_n in negatives:
                    disLabel = 'Negative'
                    procLabel = 'Negative'
                if layer_n in disease_positives:
                    disLabel='Positive'
                if layer_n in biological_process_positives:
                    procLabel='Positive'
        else:
            if n in negatives:
                disLabel = 'Negative'
                procLabel = 'Negative'
            if n in disease_positives:
                disLabel='Positive'
            if n in biological_process_positives:
                procLabel='Positive'
        final_score = d_predictions[n]*b_predictions[n]
        if n in blacklist: ## "blacklist" means it was both a pos and a neg...it's now unlabeled. Mark it if it had this
            bl = 'YES'
        else:
            bl = 'NO'
        #if layers>1:
        #    if n[-6:] == '_prime':
        #        out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (n,genemap.get(n,n[:-6]),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, G.degree(n[:-6]+'0')))
        #        degreeList.append(G.degree(n[:-6]+'0'))
        #else:
        #    out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (n,genemap.get(n,n),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, G.degree(n)))
        if layers > 1:
            name = n[:-6]
            deg = G.degree(name+'_0')-1
        else:
            name = n
            deg = G.degree(n)
        if union_predictions:
            out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%f\t%s\t%s\n' % (name,genemap.get(name,name),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,union_predictions[n],bl, deg))
        else:
            out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (name,genemap.get(name,name),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, deg))
    out.close()
    print('Wrote to %s' % (outfile))

    return

# Output results for multi-layer method
# Called once each for CM and SZ
def writeResultsMulti(statsfile,outprefix,outfile,times,changes,predictions,genemap,G,layers,pos,neg,sinksource_constant,name,sinksource_method):
    #### TODO -- if ever plotting the scores (instead of the ranks), consider the normed=False option above
    ### and normalize predictions so the max is 1.0.

    degreePos=[] #just positives ranked by where they appear in the prime rankings
    scorePos = []
    rankPos = []
    degreeNeg = [] #just negatives 
    scoreNeg = []
    rankNeg = []
    degreeUnlabeled = [] #just unlabeled
    scoreUnlabeled = []
    rankUnlabeled = []
    relativerankUnlabeled = []
    degreeList = [] #all prime nodes
    degree_file=open(outprefix+'_degree_lists_%s_layer%sconstant_degrees.txt' % (layers, sinksource_constant), 'w')
    out = open(statsfile,'w')
    out.write('#Iter\tTime\tChange\n')
    for i in range(len(times)):
        out.write('%d\t%.2f\t%e\n' % (i,times[i],changes[i]))
    out.close()
    print('Wrote to %s' % (statsfile))

    out = open(outfile,'w') 
    out.write('#EntrezID\tName\tScore\tDegree\n') #writes non-prime degree of nodes corresponding to prime node
    sorted_predict = sorted(predictions, key=lambda x:predictions[x], reverse=True) #list of nodes sorted by high to low rank
    rank = 1 #relative rank of prime nodes - second prime node found has rank 1, third 2, etc. (updated when a prime node seen)
    unlabeled_rank = 1
    for n in sorted_predict:
        if n[-6:] == '_prime':
            entrez = n[:-6]
            node = entrez+'_0' #look at the nodes in the first layer to check the degree because corresponding nodes in 
            #layers all have the same degree
            degreeList.append(G.degree(node))
            out.write('%s\t%s\t%f\t%s\n' % (entrez,genemap.get(entrez,entrez),predictions[n],G.degree(node)))
            is_pos = False
            is_neg = False
            for layer in range(layers): #search through all k nodes prime node is connected to to see if one is a pos or neg
                if entrez+'_'+str(layer) in pos:
                    is_pos = True
                
                elif entrez+'_'+str(layer) in neg:
                    is_neg = True
            if is_pos == True:
                degreePos.append(G.degree(node))
                scorePos.append(predictions[n])
                rankPos.append(rank)
            elif is_neg == True:
                degreeNeg.append(G.degree(node))
                scoreNeg.append(predictions[n])
                rankNeg.append(rank)
            else: #Unlabeled
                degreeUnlabeled.append(G.degree(node))
                scoreUnlabeled.append(predictions[n])
                rankUnlabeled.append(rank)
                relativerankUnlabeled.append(unlabeled_rank)
                unlabeled_rank += 1
            rank += 1
            

    out.close()
    degree_file.close()
    print('Wrote to %s' % (outfile))
    movingAverage=[]
    alpha=15

    #goes through list of sorted candidates and gets the average prediction scores of every 15 scores
    for i in range(len(degreeList)-alpha): 
        average = sum(degreeList[i:i+alpha])/(alpha+1)
        movingAverage.append(average)

    print('all deg', len(degreeList))
    print('moving avg', len(movingAverage))

    movingAvgUnlabeled=[]
    window=15
    for i in range(len(degreeUnlabeled)-window):
        avg = sum(degreeUnlabeled[i:i+alpha])/(window+1)
        movingAvgUnlabeled.append(avg)

    print('unlabeled deg',len(degreeUnlabeled))
    print('moving avg',len(movingAvgUnlabeled))

    fig = plt.figure()
    plt.clf()

    #why does moving average not work??
    #implement moving average for just unlabeled

    #color=rgb array
    #Scatter plot + moving average line for all prime rankings
    plt.scatter(rankUnlabeled, degreeUnlabeled, alpha=0.2, color=[.7, .7, .7], label='Unlabeled')
    plt.scatter(rankPos, degreePos, alpha=0.4, color='b', label='Positives')
    plt.scatter(rankNeg, degreeNeg, alpha=0.4, color='r', label='Negatives')

    plt.plot(movingAverage, color='black')
    plt.legend(loc='upper right')
    plt.xlabel('Candidate Rank')
    plt.ylabel('Candidate Degrees')
    
    if sinksource_method: 
        plt.savefig(outprefix+'_%s_%slayers_SS+%s_rankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to %s_%s_%slayers_SS+%s_rankxDegree.png' % (outprefix,name, layers, sinksource_constant))
    else:
        plt.savefig(outprefix+'_%s_%slayers_rankxDegree.png' % (name, layers))
        print('Wrote to %s_%s_%slayers_rankxDegree.png' % (outprefix,name, layers))

    #Scatter plot + moving average line for just unlabeled prime rankings
    fig = plt.figure()
    plt.scatter(relativerankUnlabeled, degreeUnlabeled, alpha=0.2, color=[0.3, 0, 0.8])

    plt.plot(movingAvgUnlabeled, color='black')
    plt.xlabel('Unlabeled Candidate Rank')
    plt.ylabel('Unlabeled Candidate Degree')

    if sinksource_method:
        plt.savefig(outprefix+'_%s_%slayers_SS+%s_UnlabeledrankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to %s_%s_%slayers_SS+%s_UnlabeledrankxDegree.png' % (outprefix,name, layers, sinksource_constant))
    else:
        plt.savefig(outprefix+'_%s_%slayers_UnlabeledrankxDegree.png' % (name, layers))
        print('Wrote to %s_%s_%slayers_UnlabeledrankxDegree.png' % (outprefix,name, layers))

    #Truncated plots
    fig = plt.figure()
    plt.clf()
    plt.scatter(rankUnlabeled, degreeUnlabeled, alpha=0.2, color=[.7, .7, .7], label='Unlabeled')
    plt.scatter(rankPos, degreePos, alpha=0.4, color='b', label='Positives')
    plt.scatter(rankNeg, degreeNeg, alpha=0.4, color='r', label='Negatives')

    plt.plot(movingAverage, color='black')
    plt.xlim(0,1000)
    plt.legend(loc='upper right')
    plt.xlabel('Candidate Rank')
    plt.ylabel('Candidate Degrees')

    if sinksource_method:
        plt.savefig(outprefix+'_%s_%slayers_SS+%s_Truncated_rankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to %s_%s_%slayers_SS+%s_Truncated_rankxDegree.png' % (outprefix,name, layers, sinksource_constant))
    else:
        plt.savefig(outprefix+'_%s_%slayers_Truncated_rankxDegree.png' % (name, layers))
        print('Wrote to %s_%s_%slayers_Truncated_rankxDegree.png' % (outprefix,name, layers))

    return


# Output results for single-layer method 
# Called once for each set of positives (SZ and then CM)
# predictions is the output from the method
def writeResultsSingle(statsfile,outprefix,outfile,times,changes,predictions,genemap, G,pos,neg, sinksource_constant,name,sinksource_method):
    print('Writing %s single-layer output files...' % name)
    degreePos=[]
    scorePos = []
    rankPos = []
    degreeNeg = []
    scoreNeg = []
    rankNeg = []
    degreeUnlabeled = []
    scoreUnlabeled = []
    rankUnlabeled = [] #unlabeled ranked nodes relative to all ranked nodes
    relativerankUnlabeled = [] #unlabeled ranked nodes relative to each other
    degreeList = []
    rank = 1
    Unlabeled_rank = 1 #keeps track of rank relative to unlabeled nodes
    degree_file=open(outprefix+'_degree_lists_1-layer%sconstant_degrees.txt' % sinksource_constant, 'w')

    out = open(statsfile,'w')
    out.write('#Iter\tTime\tChange\n')
    for i in range(len(times)):
        out.write('%d\t%.2f\t%e\n' % (i,times[i],changes[i]))
    out.close()
    print('Wrote to %s' % (statsfile))

    out = open(outfile,'w')
    out.write('#EntrezID\tName\tScore\n')
    #predictions dictionary keys (nodes) are sorted by score from highest to lowest - this list is iterated through
    for n in sorted(predictions, key=lambda x:predictions[x], reverse=True):
        out.write('%s\t%s\t%f\t%s\n' % (n,genemap.get(n,n),predictions[n], G.degree(n)))
        degreeList.append(G.degree(n))
        if n in pos:
            degreePos.append(G.degree(n))
            scorePos.append(predictions[n])
            rankPos.append(rank)
        elif n in neg: 
            degreeNeg.append(G.degree(n))
            scoreNeg.append(predictions[n])
            rankNeg.append(rank)
        else:
            degreeUnlabeled.append(G.degree(n))
            scoreUnlabeled.append(predictions[n])
            rankUnlabeled.append(rank)
            relativerankUnlabeled.append(Unlabeled_rank)
            Unlabeled_rank +=1
        rank += 1

    out.close()
    print('Wrote to %s' % (outfile))

    movingAverage=[]
    alpha = 15
    for i in range(len(degreeList)-alpha):
        average=sum(degreeList[i:i+alpha])/(alpha+1)
        movingAverage.append(average)


    movingAvgUnlabeled=[]
    window=15
    for i in range(len(degreeUnlabeled)-window):
        avg = sum(degreeUnlabeled[i:i+alpha])/(window+1)
        movingAvgUnlabeled.append(avg)

    fig = plt.figure()
    plt.clf()

    #color=rgb array
    plt.scatter(rankUnlabeled, degreeUnlabeled, alpha=0.2, color=[.7, .7, .7], label='Unlabeled')
    plt.scatter(rankPos, degreePos, alpha=0.4, color='b', label='Positives')
    plt.scatter(rankNeg, degreeNeg, alpha=0.4, color='r', label='Negatives')
    
    plt.plot(movingAverage, color='black')
    plt.legend(loc='upper right')
    plt.xlabel('Candidate Rank')
    plt.ylabel('Candidate Degrees')
    
    if sinksource_method:
        plt.savefig(outprefix+'_%s_1layer_SS+%s_rankxDegree.png' % (name, sinksource_constant))
        print('Wrote to %s_%s_1layer_SS+%s_rankxDegree.png' % (outprefix,name, sinksource_constant))
    else:
        plt.savefig('%s_%s_1layer_rankxDegree.png' % (outprefix,name))
        print('Wrote to %s_%s_1layer_rankxDegree.png' % (outprefix,name))

    fig = plt.figure()
    plt.scatter(relativerankUnlabeled, degreeUnlabeled, alpha=0.2, color=[0.3, 0, 0.8])

    plt.plot(movingAvgUnlabeled, color='black')
    plt.xlabel('Unlabeled Candidate Rank')
    plt.ylabel('Unlabeled Candidate Degree')

    if sinksource_method:
        plt.savefig('%s_%s_1layer_SS+%s_UnlabeledrankxDegree.png' % (outprefix,name, sinksource_constant))
        print('Wrote to %s_%s_1layer_SS+%s_UnlabeledrankxDegree.png' % (outprefix,name, sinksource_constant))
    else:
        plt.savefig('%s_%s_1layer_UnlabeledrankxDegree.png' % (outprefix,name))
        print('Wrote to %s_%s_1layer_UnlabeledrankxDegree.png' % (outprefix,name))

    #Truncated
    fig = plt.figure()
    plt.clf()

    #color=rgb array
    plt.scatter(rankUnlabeled, degreeUnlabeled, alpha=0.2, color=[.7, .7, .7], label='Unlabeled')
    plt.scatter(rankPos, degreePos, alpha=0.4, color='b', label='Positives')
    plt.scatter(rankNeg, degreeNeg, alpha=0.4, color='r', label='Negatives')
    
    plt.plot(movingAverage, color='black')
    plt.xlim(0,1000)
    plt.legend(loc='upper right')
    plt.xlabel('Candidate Rank')
    plt.ylabel('Candidate Degrees')

    if sinksource_method:
        plt.savefig('%s_%s_1layer_SS+%s_TruncatedrankxDegree.png' % (outprefix,name, sinksource_constant))
        print('Wrote to %s_%s_1layer_SS+%s_TruncaterankxDegree.png' % (outprefix,name, sinksource_constant))
    else:
        plt.savefig('%s_%s_1layer_TruncatedrankxDegree.png' % (outprefix,name))
        print('Wrote to %s_%s_1layer_TruncatedrankxDegree.png' % (outprefix,name))

    return


def readResults(statsfile,outfile,layers):
    times = []
    changes = []
    predictions = {}
    with open(statsfile,'r') as fin:
        for line in fin:
            if line[0]=='#' or 'Iter' in line:
                continue
            row = line.strip().split('\t')
            times.append(float(row[1]))
            changes.append(float(row[2]))
    with open(outfile) as fin:
        for line in fin:
            if line[0] == '#' or 'EntrezID' in line:
                continue
            row = line.strip().split('\t')
            if layers > 1:
                row[0] = row[0]+'_prime'
            predictions[row[0]] = float(row[2])

    return times,changes,predictions


