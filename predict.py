## Full program that calls all parts of the analysis.
import sys
import os.path
import time
import math
import networkx as nx
from operator import itemgetter
import matplotlib.pyplot as plt
from optparse import OptionParser

def parse_arguments(argv):
    usage = 'predict.py [options]\n'
    parser = OptionParser(usage=usage)

    ## FILES
    parser.add_option('-g','--interaction_graph',\
        type='string',metavar='STR',default='brain_top_geq_0.200.txt',\
        help='Functional interaction network (default="brain_top_geq_0.200.txt").')
    parser.add_option('-b','--biological_process_positives',\
        type='string',metavar='STR',default='infiles/motility_positives.txt',\
        help='File of positives for the biological process (default="infiles/motility_positives.txt")')
    parser.add_option('-d','--disease_positives',\
        type='string',metavar='STR',default='infiles/SZ_positives.txt',\
        help='File of positives for the disease (default="infiles/SZ_positives.txt")')
    parser.add_option('-n','--negatives',\
        type='string',metavar='STR',default='infiles/SZ_negatives.txt',\
        help='File of negatives (default="infiles/SZ_negatives.txt")')
    parser.add_option('-m','--gene_map_file',\
        type='string',metavar='STR',default='Homo_sapiens.txt',\
        help='File of EntrezID to Common Names (default="Homo_sapiens.txt")')
    parser.add_option('-o','--outprefix',\
        type='string',metavar='STR',default='outfiles/out',\
        help='output file prefixes (default="outfiles/out").')

    ## ARGUMENTS
    parser.add_option('','--matrix',\
        action='store_true',default=False,\
        help='Run the matrix version of the iterative method (default=False).')
    parser.add_option('-e','--epsilon',\
        type='float',metavar='FLOAT',default=0.001,\
        help='Epsilon to terminate the iterative process (default = 0.001).')
    parser.add_option('-t','--timesteps',\
        type='int',metavar='INT',default=150,\
        help='Maximum number of timesteps to consider (default=150).')

    ## OUTPUTS
    parser.add_option('','--stats',\
        action='store_true',default=False,\
        help='Plot statistics about the network. Default=False.')
    parser.add_option('','--plot',\
        action='store_true',default=False,\
        help='Plot outputs about method iterations. Default=False.')
    parser.add_option('','--format',\
        action='store_true',default=False,\
        help='Format LaTeX tables. Default=False.')
    parser.add_option('','--force',\
        action='store_true',default=False,\
        help='Run the learner on both positive sets, overwriting files if necessary.  Default=False.')
    parser.add_option('','--verbose',\
        action='store_true',default=False,\
        help='Print extra statements to the screen. Default=False.')
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()
    return opts

def main(argv):
    opts = parse_arguments(argv)
    print(opts)
    
    genemap = geneMapReader(opts.gene_map_file)
    G = nx.Graph()

    print('\nLoading Data...')

    print(" reading edge file...")
    edgeFile=open('brain_top_geq_0.200.txt','r')
    read_edge_file(opts.interaction_graph,G)
    print(' ',G.number_of_edges(), 'edges')
    print(' ',G.number_of_nodes(), 'nodes')

    print(' reading positive and negative files...')
    disease_positives = curatedFileReader(opts.disease_positives,G,opts.verbose)
    biological_process_positives = curatedFileReader(opts.biological_process_positives,G,opts.verbose)
    negatives = curatedFileReader(opts.negatives,G,opts.verbose)
    blacklist = set()
    if disease_positives.intersection(negatives):
        overlap_set = disease_positives.intersection(negatives)
        print('WARNING: %d nodes are disease positives and negatives. Ignoring.' % (len(overlap_set)))
        negatives = negatives.difference(overlap_set)
        disease_positives = disease_positives.difference(overlap_set)
        blacklist.update(overlap_set)

    if biological_process_positives.intersection(negatives):
        overlap_set = biological_process_positives.intersection(negatives)
        print('WARNING: %d nodes are biological process positives and negatives. Ignoring.' % (len(overlap_set)))
        negatives = negatives.difference(overlap_set)
        biological_process_positives = biological_process_positives.difference(overlap_set)
        blacklist.update(overlap_set)
        
    print('Final Curated Sets: %d Disease Positives, %d Biological Process Positives, and %d Negatives.' % \
        (len(disease_positives),len(biological_process_positives),len(negatives)))
    print('%d nodes have been blacklisted because they were in both positive and negative sets.' % (len(blacklist)))

    print('\nRunning Learning Algorithms...')

    print(' disease predictions...')
    statsfile = opts.outprefix + '_disease_stats.txt'
    outfile = opts.outprefix+'_disease_output.txt'
    if not opts.force and os.path.isfile(outfile):
        print('  File %s exists. Not running (use --force to override)' % (outfile))
        d_times,d_changes,d_predictions = readResults(statsfile,outfile)
    else:
        if opts.matrix:
            pass
        else:
            setGraphAttrs(G,disease_positives,negatives)
            d_times,d_changes,d_predictions = iterativeUpdate(G,opts.epsilon,opts.timesteps,opts.verbose)
            writeResults(statsfile,outfile,d_times,d_changes,d_predictions,genemap)

    print(' biological process predictions...')
    statsfile = opts.outprefix + '_biological_process_stats.txt'
    outfile = opts.outprefix+'_biological_process_output.txt'
    if not opts.force and os.path.isfile(outfile):
        print('  File %s exists. Not running (use --force to override)' % (outfile))
        b_times,b_changes,b_predictions = readResults(statsfile,outfile)
    else:
        if opts.matrix:
            pass
        else:
            setGraphAttrs(G,biological_process_positives,negatives)
            b_times,b_changes,b_predictions = iterativeUpdate(G,opts.epsilon,opts.timesteps,opts.verbose)
            writeResults(statsfile,outfile,b_times,b_changes,b_predictions,genemap)
    
    print('\nWriting Output and Post-Processing...')
    outfile = opts.outprefix+'_combined_output.txt'
    writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap)

    if opts.plot:
        plt.clf()
        plt.plot(range(len(d_times)),d_times,'-r',label='Disease Predictions')
        plt.plot(range(len(b_times)),b_times,'-b',label='Biological Process Predictions')
        plt.legend()
        plt.ylabel('Time Per Iteration')
        plt.xlabel('Iteration')
        plt.savefig(opts.outprefix+'_timeCourse.png')

        plt.clf()
        plt.plot(range(len(d_changes)),[math.log(x) for x in d_changes],'-r',label='Disease Predictions')
        plt.plot(range(len(b_changes)),[math.log(x) for x in b_changes],'-b',label='Biological Process Predictions')
        plt.legend()
        plt.ylabel('Log Absolute Value Change Per Iteration')
        plt.xlabel('Iteration')
        plt.savefig(opts.outprefix+'_scoreChange.png')

    if opts.format:
        outfile = opts.outprefix+'_formatted.txt'
        formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap)

    print('Done.')

    return

def formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap):
    # write output
    out = open(outfile,'w')
    out.write('\\begin{table*}[h]\n')
    out.write('\\centering\n')
    out.write('\\begin{tabular}{ll|c|cc|cc|}\n')
    out.write(' & & & \\multicolumn{2}{|c|}{Schizophrenia ($\mathcal{D}$)} & \\multicolumn{2}{|c|}{Cell Motility ($\mathcal{P})$} \\\\\n')
    out.write('Name & EntrezID & Score $g(v)$ & In $C_{\mathcal{D}}$? & $f_{\mathcal{D}}$ & In $C_{\mathcal{P}}$? & $f_{\mathcal{P}}$\\\\\\hline\n')
    for n in sorted(G.nodes(), key=lambda x:d_predictions[x]*b_predictions[x], reverse=True):
        score = d_predictions[n]*b_predictions[n]
        if d_predictions[n] < 0.5 or b_predictions[n] < 0.5 or score == 0.25:
            continue
        name = genemap.get(n,n)
        entrezID = n
        if n in disease_positives:
            d_pos = '$\\checkmark$'
        else:
            d_pos = ''
        if n in biological_process_positives:
            b_pos = '$\\checkmark$'
        else:
            b_pos = ''
        
        out.write('%s & %s & %.2f & %s & %.2f & %s & %.2f\\\\\n' % (name,entrezID,score,d_pos,d_predictions[n],b_pos,b_predictions[n]))
    out.write('\\end{tabular}\n')
    out.write('\\end{table*}\n')
    print('Wrote to %s' % (outfile))
    return

def writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap):
    # write output
    out = open(outfile,'w')
    out.write('#EntrezID\tName\tDisLabel\tDisScore\tProcLabel\tProcScore\tCombined\tBlackList?\n')
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

def setGraphAttrs(graph,pos,neg):
    for node in graph.nodes():
        if node in neg and node in pos:
            sys.exit('ERROR: Node %s is both a negative and a positive.' % (node))
        if node in pos:
            graph.nodes[node]['score']=1.0
            graph.nodes[node]['prev_score']=1.0
            graph.nodes[node]['label']='Positive'
        elif node in neg:
            graph.nodes[node]['score']=0.0
            graph.nodes[node]['prev_score']=0.0
            graph.nodes[node]['label']='Negative'
        else:
            graph.nodes[node]['score']=0.5
            graph.nodes[node]['prev_score']=0.5
            graph.nodes[node]['label']='Unlabeled'
    return

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
                    graph.add_node(node, prev_score=0.5, score=0.5, label='Unlabeled', untouched=True)
                    all_nodes.add(node)
            graph.add_edge(line[0],line[1], weight=line[2])
    return

def iterativeUpdate(G,epsilon,timesteps,verbose):
    changeLogger=[]
    timeLogger=[]
    for t in range(0,timesteps):
        
        start = time.time()
        
        changes = iterativeMethod(G,t,verbose)
        end = time.time()
        timeLogger.append(end-start)
        changeLogger.append(changes)

        print("t = %d: change = %.4f" % (t,changes))

        if changes < epsilon:
            print('BELOW THRESHOLD OF %.2e! Breaking out of loop.' % (epsilon))
            break
            
        if verbose:    
            done=float(t)/float(timesteps)
            print('Time Elapsed:', end-start)
            if done!=0:
                print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')

    # predictions is a dictionary of nodes to values.
    predictions = {}
    nodes = G.nodes()
    for n in nodes:
        predictions[n] = nodes[n]['score']

    return timeLogger,changeLogger, predictions

#Takes as input a networkx graph
def iterativeMethod(Graph, t,verbose):
    positivechangesum=0
    sumofchanges=0
    changed=0
    changedNegative=0
    changedPositive=0
    untouchedSet=0
    #Note: Graph.nodes() is a list of all the nodes
    #Graph.nodes[node] is a dictionary of that node's attributes
    nodes = Graph.nodes()
    for node in nodes:
        #Want to keep positive scores at 1 and negative scores at 0 (or -1)
        if nodes[node]['label'] != 'Unlabeled':
            pass
        else:
            newConfidence = 0
            sumofWeights = 0
            #Note: Graph.adj[node].items() gives a list of tuples. Each tuple includes one of the
            #node's neighbors and a dictionary of the attributes that their shared edge has i.e. weight
            #neighbor is the node's neighbor, datadict is the dictionary of attributes
            for neighbor, datadict in Graph.adj[node].items(): 
                newConfidence = newConfidence + datadict['weight']*nodes[neighbor]['prev_score'] #use t-1 score
                sumofWeights = sumofWeights + datadict['weight']
            if sumofWeights>0:
                newScore = float(newConfidence)/float(sumofWeights)
                nodes[node]['score'] = newScore #update score to be the new score
        

    #After each time step is complete, the previous score is updated to be the current score
    #Prints the scores after each timestep is complete

    for node in nodes:
        if nodes[node]['prev_score'] != nodes[node]['score']:
            changed += 1
            nodes[node]['untouched']=False
            sumofchanges=sumofchanges+abs(nodes[node]['prev_score'] - nodes[node]['score'])
            if nodes[node]['prev_score'] > nodes[node]['score']:
                changedNegative += 1
            else:
                changedPositive += 1
                positivechangesum=positivechangesum+abs(nodes[node]['prev_score'] - nodes[node]['score'])


        if nodes[node]['untouched']==True:
            untouchedSet += 1

        nodes[node]['prev_score'] = nodes[node]['score']

    if verbose:
        #     print(str(node) + " Label: " + str(nodes[node]['label']) + ", Score: " + str(nodes[node]['score']))
        print(changed, 'of', len(nodes), 'nodes changed')
        print(changedNegative, 'node scores decreased')
        print(changedPositive, 'node scores increased')
        print('Sum of absolute value of changes:', sumofchanges)
        print('Sum of positive changes:', positivechangesum)
        print('Untouched nodes:', untouchedSet)
    return sumofchanges



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






if __name__ == '__main__':
    main(sys.argv)


