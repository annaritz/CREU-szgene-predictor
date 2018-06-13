#Contains the functions that handle the edge files and write the output files, and 

import random
import matplotlib.pyplot as plt


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


#Takes in edge list file from Humanbase - 3 columns gene 1, gene 2, functional interaction probability

def read_edge_file(filename, graph, layers):
    if layers > 1: #creates new edges for multi-layered method (default)
        all_nodes = set()
        with open(filename,'r') as fin:
            for line in fin:
                line=line.split('\t')
                line[2]=float(line[2][:len(line[2])-2]) #Wrestling with the formatting
                for i in range(0,2):
                    node=line[i]
                    if node not in all_nodes:
                        graph.add_node(node+'_prime', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=layers)
                        for layer in range(layers):
                            graph.add_node(node+'_'+str(layer), prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=1)
                            graph.add_edge(node+'_'+str(layer),node+'_prime', weight=1.0)
                        all_nodes.add(node)
                for layer in range(layers):
                    node0 = line[0]+'_'+str(layer)
                    node1 = line[1]+'_'+str(layer)
                    graph.nodes[node0]['weighted_degree'] += line[2]
                    graph.nodes[node1]['weighted_degree'] += line[2]
                    graph.add_edge(node0,node1, weight=line[2])

    else: #Single-layered method 
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
        # nodes_to_remove=set()
        # for node in graph.nodes():
        #     adj_dict=graph.adj[node]
        #     if len(adj_dict)<2:
        #         nodes_to_remove.add(node)
        # for node in nodes_to_remove:
        #     graph.remove_node(node)

# Handles reading curated files when there is one layer 
def curatedFileReader(filename,graph,verbose, layers):
    #Note: Graph.nodes() is a list of all the nodes
    if layers==1:
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
    else:
        return curatedFileReaderMulti(filename, graph, verbose, layers)


# This handles reading the curated files when there is more than one layer 
# read_edge_file has already created multiple layers, this checks if the nodes are in the graph
def curatedFileReaderMulti(filename,graph,verbose, layers):
    #Note: Graph.nodes() is a list of all the nodes
    nodes = graph.nodes()
    tot = 0 #keeps track of number of nodes from curated positives
    count = 0 #keeps track of number of nodes actually in the graph
    curated = set()
    labeled_set = set()

    with open(filename,'r') as fin:
        for line in fin:
            entrezNumber=line.strip()
            tot += 1
            for layer in range(layers):
                if entrezNumber+'_'+str(layer) in nodes:
                    labeled_set.add(entrezNumber)
                else:
                    if verbose:
                        print('WARNING: EntrezID %s is not in graph.' % (entrezNumber))
                    else:
                        continue
        count=len(labeled_set)

        random.seed('Alexander_King') #We need the positives to be randomly distributed throughout the layers. If we have multiple
        labeled_List=list(labeled_set) #labeled nodes surrounding a prime node, that effectively blocks all score transmission
        labeled_List=random.sample(labeled_List, k=len(labeled_List))

        for i in range(len(labeled_List)):
            layer=(i//(len(labeled_List)//layers))
            curated.add(labeled_List[i]+'_'+str(layer))


    print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
    return curated

#Formats results as a LaTeX table
#outfile variable is the name of the outfile, partially specified by --outprefix 
def formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap, layers):
    # write output
    out = open(outfile,'w')
    out.write('\\begin{table}[h]\n')
    out.write('\\centering\n')
    out.write('\\begin{tabular}{|ll|ccc|}\\hline\n')
    out.write('Gene Name & EntrezID & $f_{\mathcal{D}}$ & $f_{\mathcal{P}}$ & Score $g(v)$\\\\\\hline\n')
    for n in sorted(G.nodes(), key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        if layers==1 or n[-6:] == '_prime':
            score = d_predictions[n]*b_predictions[n]
            if d_predictions[n] < 0.5 or b_predictions[n] < 0.5 or score == 0.25:
                continue
            if n[-6:] == '_prime':
                n=n[:-6]
            name = genemap.get(n,n)
            entrezID = n
            out.write('%s & %s ' % (name,entrezID))
            if n in disease_positives:
                out.write(' & \\textit{%.2f}' % (d_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (d_predictions[n]))
            if n in biological_process_positives:
                out.write(' & \\textit{%.2f}' % (b_predictions[n]))
            else:
                out.write(' & \\textbf{%.2f}' % (b_predictions[n]))
            out.write(' & %.2f\\\\\n' % (score))
    out.write('\\end{tabular}\n')
    out.write('\\end{table}\n')
    print('Wrote to %s' % (outfile))
    return

# Write output file of nodes and process score * disease score 
def writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap, layers):
    # write output
    out = open(outfile,'w')
    fig = plt.figure(figsize=(4,4))
    out.write('#EntrezID\tName\tDisLabel\tDisScore\tProcLabel\tProcScore\tCombined\tConflict?\n')
    degreeList=[]
    for n in sorted(G.nodes(), key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        disLabel='Unlabeled'
        procLabel='Unlabeled'
        if n in negatives:
            disLabel = 'Negative'
            procLabel = 'Negative'
        if n in disease_positives:
            disLabel='Positive'
        if n in biological_process_positives:
            procLabel='Positive'
        final_score = d_predictions[n]*b_predictions[n]
        if n in blacklist:
            bl = 'YES'
        else:
            bl = 'NO'
        if layers>1:
            if n[-6:] == '_prime':
                out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (n[:-6],genemap.get(n[:-6],n[:-6]),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, G.degree(n[:-6]+'0')))
                degreeList.append(G.degree(n[:-6]+'0'))
        else:
            out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (n,genemap.get(n,n),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, G.degree(n)))
    out.close()
    print('Wrote to %s' % (outfile))

    return

# Output results for multi-layer method
# Called once each for CM and SZ
def writeResultsMulti(statsfile,outfile,times,changes,predictions,genemap,G,layers,pos,neg,sinksource_constant,name,sinksource_method):
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
    degree_file=open('degree_lists/%s_layer%sconstant_degrees.txt' % (layers, sinksource_constant), 'w')
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
        plt.savefig('outfiles/%s_%slayers_SS+%s_rankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to outfiles/%s_%slayers_SS+%s_rankxDegree.png' % (name, layers, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_%slayers_rankxDegree.png' % (name, layers))
        print('Wrote to outfiles/%s_%slayers_rankxDegree.png' % (name, layers))

    #Scatter plot + moving average line for just unlabeled prime rankings
    fig = plt.figure()
    plt.scatter(relativerankUnlabeled, degreeUnlabeled, alpha=0.2, color=[0.3, 0, 0.8])

    plt.plot(movingAvgUnlabeled, color='black')
    plt.xlabel('Unlabeled Candidate Rank')
    plt.ylabel('Unlabeled Candidate Degree')

    if sinksource_method:
        plt.savefig('outfiles/%s_%slayers_SS+%s_UnlabeledrankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to outfiles/%s_%slayers_SS+%s_UnlabeledrankxDegree.png' % (name, layers, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_%slayers_UnlabeledrankxDegree.png' % (name, layers))
        print('Wrote to outfiles/%s_%slayers_UnlabeledrankxDegree.png' % (name, layers))

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
        plt.savefig('outfiles/%s_%slayers_SS+%s_Truncated_rankxDegree.png' % (name, layers, sinksource_constant))
        print('Wrote to outfiles/%s_%slayers_SS+%s_Truncated_rankxDegree.png' % (name, layers, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_%slayers_Truncated_rankxDegree.png' % (name, layers))
        print('Wrote to outfiles/%s_%slayers_Truncated_rankxDegree.png' % (name, layers))

    return


# Output results for single-layer method 
# Called once for each set of positives (SZ and then CM)
# predictions is the output from the method
def writeResultsSingle(statsfile,outfile,times,changes,predictions,genemap, G,pos,neg, sinksource_constant,name,sinksource_method):
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
    degree_file=open('degree_lists/1-layer%sconstant_degrees.txt' % sinksource_constant, 'w')

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
        plt.savefig('outfiles/%s_1layer_SS+%s_rankxDegree.png' % (name, sinksource_constant))
        print('Wrote to outfiles/%s_1layer_SS+%s_rankxDegree.png' % (name, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_1layer_rankxDegree.png' % name)
        print('Wrote to outfiles/%s_1layer_rankxDegree.png' % name)

    fig = plt.figure()
    plt.scatter(relativerankUnlabeled, degreeUnlabeled, alpha=0.2, color=[0.3, 0, 0.8])

    plt.plot(movingAvgUnlabeled, color='black')
    plt.xlabel('Unlabeled Candidate Rank')
    plt.ylabel('Unlabeled Candidate Degree')

    if sinksource_method:
        plt.savefig('outfiles/%s_1layer_SS+%s_UnlabeledrankxDegree.png' % (name, sinksource_constant))
        print('Wrote to outfiles/%s_1layer_SS+%s_UnlabeledrankxDegree.png' % (name, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_1layer_UnlabeledrankxDegree.png' % name)
        print('Wrote to outfiles/%s_1layer_UnlabeledrankxDegree.png' % name)

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
        plt.savefig('outfiles/%s_1layer_SS+%s_TruncatedrankxDegree.png' % (name, sinksource_constant))
        print('Wrote to outfiles/%s_1layer_SS+%s_TruncaterankxDegree.png' % (name, sinksource_constant))
    else:
        plt.savefig('outfiles/%s_1layer_TruncatedrankxDegree.png' % name)
        print('Wrote to outfiles/%s_1layer_TruncatedrankxDegree.png' % name)

    return


def readResults(statsfile,outfile):
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
            predictions[row[0]] = float(row[2])

    return times,changes,predictions


