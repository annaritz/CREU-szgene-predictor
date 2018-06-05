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

def writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap):
    # write output
    out = open(outfile,'w')
    out.write('#EntrezID\tName\tDisLabel\tDisScore\tProcLabel\tProcScore\tCombined\tConflict?\n')
    

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
        out.write('%s\t%s\t%s\t%f\t%s\t%f\t%f\t%s\n' % (n,genemap.get(n,n),disLabel,d_predictions[n],procLabel,b_predictions[n],final_score,bl))
    out.close()
    print('Wrote to %s' % (outfile))
    return

def writeResults(statsfile,outfile,times,changes,predictions,genemap):
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
    out.close()
    print('Wrote to %s' % (outfile))
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

def curatedFileReader(filename,graph,verbose):
    #Note: Graph.nodes() is a list of all the nodes
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

    print('%d of %d nodes are in graph from file %s' % (count,tot,filename))
    return curated

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
def read_edge_file(filename, graph):
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