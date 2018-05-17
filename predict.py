## Full program that calls all parts of the analysis.
import sys
import os.path
import time
import math
import networkx as nx
from operator import itemgetter
import matplotlib.pyplot as plt
import random
from optparse import OptionParser
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import scipy.stats as stats

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
    parser.add_option('-e','--epsilon',\
        type='float',metavar='FLOAT',default=0.001,\
        help='Epsilon to terminate the iterative process (default = 0.001).')
    parser.add_option('-t','--timesteps',\
        type='int',metavar='INT',default=500,\
        help='Maximum number of timesteps to consider (default=500).')
    parser.add_option('-k','--k_fold',\
        type='int',metavar='INT',default=5,\
        help='k for k-fold cross validation (default=5).')
    parser.add_option('-a','--auc_samples',\
        type='int',metavar='INT',default=50,\
        help='number of cross validation iterations to compute AUC (default=50).')
    parser.add_option('','--matrix',\
        action='store_true',default=False,\
        help='Run the matrix version of the iterative method (default=False).')

    ## OUTPUTS
    parser.add_option('','--stats',\
        action='store_true',default=False,\
        help='Plot statistics about the network. Default=False.')
    parser.add_option('','--plot',\
        action='store_true',default=False,\
        help='Plot outputs about method iterations. Default=False.')
    parser.add_option('','--auc',\
        action='store_true',default=False,\
        help='Run k-fold cross validation and compute AUC values.')
    parser.add_option('','--roc',\
        action='store_true',default=False,\
        help='Compute ROC cuves after holding out the overlap set.')
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
    print('')
    print('%d nodes have been blacklisted because they were in both positive and negative sets.' % (len(blacklist)))

    print('\nRunning Learning Algorithms...')

    print(' disease predictions...')
    statsfile = opts.outprefix + '_disease_stats.txt'
    outfile = opts.outprefix+'_disease_output.txt'
    d_predictions = learn(outfile,statsfile,genemap,G,disease_positives,negs,opts.epsilon,opts.timesteps,opts.verbose,opts.force,write=True)

    print(' biological process predictions...')
    statsfile = opts.outprefix + '_biological_process_stats.txt'
    outfile = opts.outprefix+'_biological_process_output.txt'
    b_predictions = learn(outfile,statsfile,genemap,G,biological_process_positives,negs,opts.epsilon,opts.timesteps,opts.verbose,opts.force,write=True)
    
    outfile = opts.outprefix+'_combined_output.txt'
    writeCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap)

    print(' union of positives...')
    statsfile = opts.outprefix + '_union_stats.txt'
    outfile = opts.outprefix+'_union_output.txt'
    union_positives = disease_positives.union(biological_process_positives)
    u_predictions = learn(outfile,statsfile,genemap,G,union_positives,negs,opts.epsilon,opts.timesteps,opts.verbose,opts.force,write=True)
    
    if opts.auc:
        print('\nCalculating AUC w/ matrix method...')
        d_AUCs = []
        for i in range(opts.auc_samples):
            print('#%d of %d' % (i,opts.auc_samples))
            # subsaample 1/k of the positives...
            hidden_genes = random.sample(disease_positives,int(len(disease_positives)/opts.k_fold))
            test_positives = disease_positives.difference(hidden_genes)
            print('%d hidden %d test genes' % (len(hidden_genes),len(test_positives)))
            ignore,ignore,d_predictions = matrixLearn(G,test_positives,negatives,opts.epsilon,opts.timesteps,opts.verbose)
            d_AUCs.append(Mann_Whitney_U_test(d_predictions, hidden_genes,negatives))

        b_AUCs = []
        for i in range(opts.auc_samples):
            print('#%d of %d' % (i,opts.auc_samples))
            # subsaample 1/k of the positives...
            hidden_genes = random.sample(biological_process_positives,int(len(biological_process_positives)/opts.k_fold))
            test_positives = biological_process_positives.difference(hidden_genes)
            print('%d hidden %d test genes' % (len(hidden_genes),len(test_positives)))
            ignore,ignore,d_predictions = matrixLearn(G,test_positives,negatives,opts.epsilon,opts.timesteps,opts.verbose)
            b_AUCs.append(Mann_Whitney_U_test(b_predictions, hidden_genes,negatives))
        plt.clf()
        plt.boxplot([d_AUCs,b_AUCs])
        plt.savefig(opts.outprefix+'_auc.png')

        out = open(opts.outprefix+'_auc.txt','w')
        out.write('#DiseaseAUCs\tBiologicalProcessAUCs\n')
        for i in range(len(d_AUCs)):
            out.write('%f\t%f\n' % (d_AUCs[i],b_AUCs[i]))
        out.close()

    if opts.roc:
        print('\nHolding out overlap set and running procedure.')
        

    print('\nWriting Output and Post-Processing...')
    
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

        plt.clf()
        names = ['Disease Predictor $f_{\mathcal{D}}$','Biological Process Predictor $f_{\mathcal{P}}$','Union Predictor','Score $g$']
        colors =['g','b','k','r']
        n = G.number_of_nodes()
        preds = [d_predictions,b_predictions,u_predictions,{x:d_predictions[x]*u_predictions[x] for x in G.nodes()}]
        for i in range(len(names)):
            yvals =sorted(preds[i].values(),reverse=True)
            plt.plot(range(n),yvals,color=colors[i],label=names[i])
        plt.plot([0,n-1],[0.5,0.5],':k',label='_nolabel_')
        plt.legend()
        plt.xlabel('Node ($n=%s$)' % (n))
        plt.ylabel('Ranking ($f$ or $g$)')
        plt.xlim(0,n-1)
        plt.ylim(0.01,1.01)
        plt.savefig(opts.outprefix+'_nodeRankings.png')

    if opts.format:
        outfile = opts.outprefix+'_formatted.txt'
        formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap)

    print('Done.')

    return

def learn(outfile,statsfile,genemap,G,pos,negs,epsilon,timesteps,verbose,force,write=False):
    if not force and os.path.isfile(outfile):
        print('  File %s exists. Not running (use --force to override)' % (outfile))
        times,changes,predictions = readResults(statsfile,outfile)
    else:
        if opts.matrix:
            times,changes,predictions = matrixLearn(G,pos,negs,epsilon,timesteps,verbose)
        else:
            setGraphAttrs(G,disease_positives,negatives)
            times,changes,predictions = iterativeLearn(G,epsilon,timesteps,verbose)
        if write:
            writeResults(statsfile,outfile,times,changes,predictions,genemap)
    return predictions

def Mann_Whitney_U_test(predictions, hidden_nodes, negatives):
    #Runs a Mann-Whitney U test on the lists
    kfold_ranks=[] #This will be the results we have computed without the positives
    test_ranks=[] #This will be the results we have computed with all positives

    nodeValues=[] #This is the vehicle by which we extract graph information
    hiddenNodeValues=[]
    notPositiveNodeValues=[]

    for node in predictions:
        if node in hidden_nodes:
            hiddenNodeValues.append(predictions[node])
        else:
            notPositiveNodeValues.append(predictions[node])

    U, p=stats.mannwhitneyu(hiddenNodeValues, notPositiveNodeValues, alternative="two-sided")
    #print(U,p)
    AUC=U/(len(hiddenNodeValues)*len(notPositiveNodeValues))
    #print(AUC)
    return AUC

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
                    graph.add_node(node, prev_score=0.5, score=0.5, label='Unlabeled', untouched=True, weighted_degree=-1)
                    all_nodes.add(node)
            graph.add_edge(line[0],line[1], weight=line[2])
    return

def matrixLearn(G,pos,neg,epsilon,timesteps,verbose):

    ## Takes the form of f = M * f + c. 

    ## sort unlabeled nodes.
    unlabeled = set(G.nodes()).difference(pos).difference(neg)
    unlabeled_list = sorted(unlabeled)
    unlabeled_inds = {unlabeled_list[i]:i for i in range(len(unlabeled_list))}
    print('%d unlabeled nodes.' % (len(unlabeled)))

    print('Preparing matrix data')
    #Add summed degree to graph.
    #Note: Graph.adj[node].items() gives a list of tuples. Each tuple includes one of the
    #node's neighbors and a dictionary of the attributes that their shared edge has i.e. weight
    #neighbor is the node's neighbor, datadict is the dictionary of attributes
    print(' computing weighted degree...')
    for node in G.nodes():
        G.nodes[node]['weighted_degree'] = 0
        for neighbor, datadict in G.adj[node].items():
            G.nodes[node]['weighted_degree'] += datadict['weight']

    #Make sparse M matrix.
    print(' making M matrix...')
    M = lil_matrix((len(unlabeled),len(unlabeled)))
    for u,v in G.edges():
        if u in unlabeled and v in unlabeled:
            i = unlabeled_inds[u]
            j = unlabeled_inds[v]
            M[i,j] = float(G.edges[u,v]['weight'])/G.nodes[u]['weighted_degree']
            M[j,i] = float(G.edges[u,v]['weight'])/G.nodes[v]['weighted_degree']

    #Make c vector
    print(' making c vector...')
    c = [0]*len(unlabeled_list)
    for i in range(len(unlabeled_list)):
        v = unlabeled_list[i]
        for neighbor, datadict in G.adj[v].items():
            if neighbor not in unlabeled: # it is labeled
                if neighbor in pos:
                    c[i] += datadict['weight']
                else: # neighbor in negative; label is 0 so nothing is added.
                    pass

        c[i] = float(c[i])/G.nodes[v]['weighted_degree']

    #Make initial f vector.  This is a value of 0.5 for all unlabeled nodes.
    f = [0.5]*len(unlabeled_list)

    changeLogger=[]
    timeLogger=[]
    for t in range(0,timesteps):
        
        start = time.time()
        
        ## conduct sparse matrix operation.
        f_prev = f
        f = M.dot(f)+c

        ## sum changes
        changes = sum([abs(f[i]-f_prev[i]) for i in range(len(f))])

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
    for n in pos:
        predictions[n] = 1
    for n in neg:
        predictions[n] = 0
    for i in range(len(unlabeled_list)):
        predictions[unlabeled_list[i]] = f[i]
    
    return timeLogger,changeLogger, predictions

def iterativeLearn(G,epsilon,timesteps,verbose):
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


