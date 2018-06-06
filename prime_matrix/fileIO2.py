#Contains the functions that handle the edge files, write the output files, and 

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

def read_edge_file(filename, graph, single):
    if not single: #creates new edges for multi-layered method (default)
        Levels=['_E1', '_E2', '_E3']
        all_nodes = set()
        with open(filename,'r') as fin:
            for line in fin:
                line=line.split('\t')
                line[2]=float(line[2][:len(line[2])-2]) #Wrestling with the formatting
                for i in range(0,2):
                    node=line[i]
                    if node not in all_nodes:
                        graph.add_node(node+'_E1', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=1)
                        graph.add_node(node+'_E2', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=1)
                        graph.add_node(node+'_E3', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=1)
                        graph.add_node(node+'_prime', prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=3)
                        graph.add_edge(node+'_E1',node+'_prime', weight=1.0)
                        graph.add_edge(node+'_E2',node+'_prime', weight=1.0)
                        graph.add_edge(node+'_E3',node+'_prime', weight=1.0)
                        all_nodes.add(node)

                graph.nodes[line[0]+'_E1']['weighted_degree'] += line[2]
                graph.nodes[line[0]+'_E2']['weighted_degree'] += line[2]
                graph.nodes[line[0]+'_E3']['weighted_degree'] += line[2]
                graph.nodes[line[1]+'_E1']['weighted_degree'] += line[2]
                graph.nodes[line[1]+'_E2']['weighted_degree'] += line[2]
                graph.nodes[line[1]+'_E3']['weighted_degree'] += line[2]
                graph.add_edge(line[0]+'_E1',line[1]+'_E1', weight=line[2])
                graph.add_edge(line[0]+'_E2',line[1]+'_E2', weight=line[2])
                graph.add_edge(line[0]+'_E3',line[1]+'_E3', weight=line[2])

        # nodes_to_remove=set()
        # for node in graph.nodes():
        #     adj_dict=graph.adj[node]
        #     if len(adj_dict)<=2 and (node[-3:]=='_E1' or node[-3:]=='_E2' or node[-3:]=='_E3'):
        #         nodes_to_remove.add(node[:-3]+'_E1')
        #         nodes_to_remove.add(node[:-3]+'_E2')
        #         nodes_to_remove.add(node[:-3]+'_E3')
        #         nodes_to_remove.add(node[:-3]+'_prime')

        # for node in nodes_to_remove:
        #     graph.remove_node(node)

        # edges_to_add=[]
        # x=0
        # pairs=set()


        # for node in graph.nodes():
        #     if x%1==0:
        #         print(x, len(edges_to_add), len(pairs))

        #     x=x+1
        #     if node[-6:] != '_prime':
        #         adj_dict=graph.adj[node]
        #         if len(adj_dict)<10:
        #             for neighbor in adj_dict:
        #                 for double_neighbor in graph.adj[neighbor]:
        #                     if double_neighbor not in adj_dict and (node, double_neighbor) not in pairs:
        #                         edges_to_add.append([node,double_neighbor])
        #                         pairs.add((node, double_neighbor))

        # for edge in edges_to_add:
        #     graph.add_edge(edge[0],edge[1],weight=0.05)



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



#Formats results as a LaTeX table
#outfile variable is the name of the outfile, partially specified by --outprefix 
def formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap):
    # write output
    out = open(outfile,'w')
    out.write('\\begin{table}[h]\n')
    out.write('\\centering\n')
    out.write('\\begin{tabular}{|ll|ccc|}\\hline\n')
    out.write('Gene Name & EntrezID & $f_{\mathcal{D}}$ & $f_{\mathcal{P}}$ & Score $g(v)$\\\\\\hline\n')
    for n in sorted(G.nodes(), key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        score = d_predictions[n]*b_predictions[n]
        if d_predictions[n] < 0.5 or b_predictions[n] < 0.5 or score == 0.25:
            continue
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

def writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap, single_layer):
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
        if single_layer:
            if n[-6:] == '_prime':
                out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\t%s\n' % (n[:-6],genemap.get(n[:-6],n[:-6]),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl, G.degree(n[:-6]+'_E1')))

                degreeList.append(G.degree(n[:-6]+'_E1'))
        else:
            out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\n' % (n,genemap.get(n,n),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl))
    out.close()
    print('Wrote to %s' % (outfile))
    degreeList=degreeList

    return

# Output results for multi-layer method
def writeResults(statsfile,outfile,times,changes,predictions,genemap, G):
    degreeList=[]
    valueList=[]
    fig = plt.figure(figsize=(4,4))
    out = open(statsfile,'w')
    out.write('#Iter\tTime\tChange\n')
    for i in range(len(times)):
        out.write('%d\t%.2f\t%e\n' % (i,times[i],changes[i]))
    out.close()
    print('Wrote to %s' % (statsfile))

    out = open(outfile,'w')
    out.write('#EntrezID\tName\tScore\n')
    for n in sorted(predictions, key=lambda x:predictions[x], reverse=True):
        if n[-6:] == '_prime':
            out.write('%s\t%s\t%f\t%s\n' % (n[:-6],genemap.get(n[:-6],n[:-6]),predictions[n],G.degree(n[:-6]+'_E1')))
            if G.node[n[:-6]+'_E1']!= 'Positive' and G.node[n[:-6]+'_E2']!= 'Positive' and G.node[n[:-6]+'_E3']!= 'Positive':
                degreeList.append(G.degree(n[:-6]+'_E1'))
                valueList.append(predictions[n])
    out.close()
    print('Wrote to %s' % (outfile))
    # degreeList=degreeList[0:300]
    movingAverage=[]
    alpha=15
    for i in range(len(degreeList)-alpha):
        average=sum(degreeList[i:i+alpha])/(alpha+1)


        movingAverage.append(average)

    
    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('candidate_degrees_multilayer.png')
    

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('candidate_scores_multilayer.png')

    plt.figure()
    plt.plot(valueList,degreeList,'ob')
    plt.xlabel('Node Score')
    plt.ylabel('Node Degrees')
    plt.title('Candidate Nodes')
    plt.tight_layout()
    plt.savefig('candidate_degrees_by_score_multilayer.png')

    degreeList=degreeList[0:300]
    movingAverage=movingAverage[0:300]
    plt.figure()
    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_candidate_degrees_multilayer.png')

    valueList=valueList[0:300]

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_candidate_scores_multilayer.png')




    #print(fnjansfkns)

    return


# Output results for single-layer method 
# Called once for each set of positives (SZ and then CM)
def writeResultsSingle(statsfile,outfile,times,changes,predictions,genemap, G):
    degreeList=[]
    valueList=[]
    out = open(statsfile,'w')
    out.write('#Iter\tTime\tChange\n')
    for i in range(len(times)):
        out.write('%d\t%.2f\t%e\n' % (i,times[i],changes[i]))
    out.close()
    print('Wrote to %s' % (statsfile))

    out = open(outfile,'w')
    out.write('#EntrezID\tName\tScore\n')
    for n in sorted(predictions, key=lambda x:predictions[x], reverse=True):
        out.write('%s\t%s\t%f\n' % (n,genemap.get(n,n),predictions[n]))
        degreeList.append(G.degree(n))
        valueList.append(predictions[n])
    out.close()
    print('Wrote to %s' % (outfile))


    movingAverage=[]
    alpha=15
    for i in range(len(degreeList)-alpha):
        average=sum(degreeList[i:i+alpha])/(alpha+1)


        movingAverage.append(average)

    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('candidate_degrees_singlelayer.png')
    

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('candidate_scores_singlelayer.png')

    plt.figure()
    plt.plot(valueList,degreeList,'ob')
    plt.xlabel('Node Score')
    plt.ylabel('Node Degrees')
    plt.title('Candidate Nodes')
    plt.tight_layout()
    plt.savefig('candidate_degrees_by_score_singlelayer.png')

    degreeList=degreeList[0:300]
    movingAverage=movingAverage[0:300]
    plt.figure()
    plt.plot(degreeList,'ob')
    plt.plot(movingAverage)
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_candidate_degrees_singlelayer.png')

    valueList=valueList[0:300]

    plt.figure()
    plt.plot(valueList,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Score')
    plt.title('Candidate Degrees')
    plt.tight_layout()
    plt.savefig('top_300_candidate_scores_singlelayer.png')

    #print(vdfvksjdnvlkjs)

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


def curatedFileReader(filename,graph,verbose, single, minimum_labeled):
    #Note: Graph.nodes() is a list of all the nodes
    if single:
        nodes = graph.nodes()
        tot = 0
        count = 0
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


        if len(curated)<minimum_labeled:
            minimum_labeled=len(curated)
        else:
            count=minimum_labeled
        print(minimum_labeled)
        print(len(curated))

        curated=set(random.sample(curated, minimum_labeled))

        print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
        return curated, minimum_labeled
    else:
        return curatedFileReaderMulti(filename, graph, verbose, minimum_labeled)

def curatedFileReader1(filename,graph,verbose):
    #Note: Graph.nodes() is a list of all the nodes
    nodes = graph.nodes()
    tot = 0
    count = 0
    curated = set()
    with open(filename,'r') as fin:
        for line in fin:
            entrezNumber = line.strip()
            tot+=1
            if entrezNumber+'_E1' in nodes or entrezNumber+'_E2' in nodes or entrezNumber+'_E3' in nodes:
                count+=1
                curated.add(entrezNumber+'_E1')
                curated.add(entrezNumber+'_E2')
                curated.add(entrezNumber+'_E3')
            else:
                if verbose:
                    print('WARNING: EntrezID %s is not in graph.' % (entrezNumber))

    print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
    return curated, minimum_labeled


def curatedFileReaderMulti(filename,graph,verbose, minimum_labeled):
    #Note: Graph.nodes() is a list of all the nodes
    nodes = graph.nodes()
    tot = 0
    count = 0
    curated = set()
    labeled_set = set()

    with open(filename,'r') as fin:
        for line in fin:
            entrezNumber=line.strip()
            tot+=1
            if entrezNumber+'_E1' in nodes or entrezNumber+'_E2' in nodes or entrezNumber+'_E3' in nodes:
                count+=1
                labeled_set.add(entrezNumber)
            else:
                if verbose:
                    print('WARNING: EntrezID %s is not in graph.' % (entrezNumber))
                else:
                    continue
        if count<minimum_labeled:
            minimum_labeled=count
        else:
            count=minimum_labeled
        # labeled_List=random.shuffle(labeled_List)

        # labeled_set=random.shuffle(labeled_set)
        E1=set(random.sample(labeled_set, minimum_labeled//3))
        labeled_set=labeled_set-E1
        E2=set(random.sample(labeled_set, len(E1)))
        labeled_set=labeled_set-E1
        E3=set(random.sample(labeled_set, len(E1)))

        
        for entrezNumber in E1:
            curated.add(entrezNumber+'_E1')

        for entrezNumber in E2:
            curated.add(entrezNumber+'_E2')

        for entrezNumber in E3:
            curated.add(entrezNumber+'_E3')


    print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
    return curated, minimum_labeled


def plot_candidate_degrees(rank, Graph):
    print('Creating candidate plot...')
    fig = plt.figure(figsize=(4,4))
    y = []
    degrees=[]
    for cand in rank:
        node = cand[0]
        degrees.append(Graph.degree(node))
    print(sum(degrees)/len(degrees))

    for cand in rank: #iterates through list of 
        if cand[5]!= 'Positive' and cand[6]!= 'Positive' and cand[7]!= 'Positive':
            index=0
            node = cand[0]
            deg = Graph.degree(node) #looks up the degree of a candidate, returns [(cand, degree)]

            
            y.append(deg)

    plt.plot(y,'ob')
    plt.xlabel('Node Rank')
    plt.ylabel('Node Degree')
    plt.title('Candidate Degrees')

    plt.tight_layout()

    plt.savefig('candidate_degrees.png')
    print('Wrote to candidate_degrees.png')


    return




def write_output(Graph, timesteps):
    print('initializing list')
    nodeValues=[]
    for node in Graph.nodes:

        nodeValues.append([ node, Graph.nodes[node]['score'], Graph.nodes[node]['label']])



    nodeValues=sorted(nodeValues, key=itemgetter(1), reverse=True)


    x=open('cellgene_rankings_.txt', 'w')
    y=open('cellgene_rankings_no_pos.txt', 'w')
    for node in nodeValues:

        label=node[2]
        x.write(str(node)+'\n')
        if label=='Unlabeled':
            y.write(str(node)+'\n')
    return


