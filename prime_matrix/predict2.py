## Full program that calls all parts of the analysis.

## this needs to be at the top so matplotlib doesn't
## throw an error if being called remotely.
import matplotlib
matplotlib.use('Agg')

## system and utility packages
import sys
import os.path
import random
import time
from operator import itemgetter

## argument parsing
from optparse import OptionParser, OptionGroup

## graph and plotting packages
import networkx as nx
import matplotlib.pyplot as plt

## math/stats packages
import math
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import scipy.stats as stats

## Functions from other local python files
import fileIO2 as fileIO
import learners2 as learners


def parse_arguments(argv):
    '''
    Argument parser.  Run with predict.py -h to see all the options.
    Inputs: arguments passed in a command line
    Output: object that contains all parsed option.
    TODO: OptionParser is deprecatedl should replaced with argparse module
    (https://docs.python.org/2/library/optparse.html)
    '''

    usage = 'predict.py [options]\nRun experiments for submitted manuscript.  At least one of the following methods/experiments must be specified:\n\t--stats: compute statistics\n\t--single: run single experiment\n\t--auc: run k-fold cross validation\n\t--roc: compute ROC of positives that appear in both curated sets.\n'
    parser = OptionParser(usage=usage)

    ## Common Options.
    parser.add_option('-o','--outprefix',\
        type='string',metavar='STR',default='outfiles/out',\
        help='output file prefixes (default="outfiles/out").')
    parser.add_option('','--force',\
        action='store_true',default=False,\
        help='Run the learner on both positive sets, overwriting files if necessary.  Default=False.')
    parser.add_option('','--verbose',\
        action='store_true',default=False,\
        help='Print extra statements to the screen. Default=False.')

    ## Algorithms/Methods
    group = OptionGroup(parser,'Algorithms/Methods.')
    group.add_option('','--stats',\
        action='store_true',default=False,\
        help='Plot statistics about the network. Default=False.')
    group.add_option('','--single',\
        action='store_true',default=False,\
        help='Run single experiments for disease and biological process and computes score. Also generates figures of results and formats a LaTeX table.  Default=False.')
    group.add_option('','--union',\
        action='store_true',default=False,\
        help='Run single experiments for union of positives and compares to single-set positives. --single must be specified.  Default=False.')
    group.add_option('','--auc',\
        action='store_true',default=False,\
        help='Run k-fold cross validation and compute AUC values and make a figure. Default = False.')
    group.add_option('','--roc',\
        action='store_true',default=False,\
        help='Compute ROC cuves after holding out the overlap set and makes a figure. Default=False.')
    parser.add_option_group(group)

    ## Input Files
    group = OptionGroup(parser,'Input Files (all have default values)')
    group.add_option('-g','--interaction_graph',\
        type='string',metavar='STR',default='networkfiles/brain_top_geq_0.200.txt',\
        help='Functional interaction network (default="networkfiles/brain_top_geq_0.200.txt").')
    group.add_option('-b','--biological_process_positives',\
        type='string',metavar='STR',default='infiles/motility_positives.txt',\
        help='File of positives for the biological process (default="infiles/motility_positives.txt")')
    group.add_option('-d','--disease_positives',\
        type='string',metavar='STR',default='infiles/SZ_positives.txt',\
        help='File of positives for the disease (default="infiles/SZ_positives.txt")')
    group.add_option('-n','--negatives',\
        type='string',metavar='STR',default='infiles/SZ_negatives.txt',\
        help='File of negatives (default="infiles/SZ_negatives.txt")')
    group.add_option('-m','--gene_map_file',\
        type='string',metavar='STR',default='infiles/Homo_sapiens.txt',\
        help='File of EntrezID to Common Names (default="infiles/Homo_sapiens.txt")')
    parser.add_option_group(group)

    ## Algorithms/Method Arguments
    group = OptionGroup(parser,'Method Arguments (all have default values)')
    group.add_option('-e','--epsilon',\
        type='float',metavar='FLOAT',default=0.001,\
        help='Epsilon to terminate the iterative process (default = 0.001).')
    group.add_option('-t','--timesteps',\
        type='int',metavar='INT',default=500,\
        help='Maximum number of timesteps to consider (default=500).')
    group.add_option('','--iterative_update',\
        action='store_true',default=False,\
        help='Run the update (original) version of the iterative method intsead of the matrix method (default=False).')
    group.add_option('-k','--k_fold',\
        type='int',metavar='INT',default=5,\
        help='k for k-fold cross validation (default=5).')
    group.add_option('-s','--auc_samples',\
        type='int',metavar='INT',default=50,\
        help='number of cross validation iterations to compute AUC (default=50).')
    group.add_option('-l', '--layers',\
        type='int', default=3,\
        help='Run the experiments for disease and biological process with n nodes per gene. Each of the gene\'s nodes are connected to a node that represents the score for the gene. Positives are distributed among these layers. Reducing the nodes to 1 eliminates the process')
    group.add_option('-c', '--sinksource_constant',\
        type='float',metavar='FLOAT',default=1,\
        help='How much to add to the denominator of the node value (default=1)')
    parser.add_option_group(group)

    group = OptionGroup(parser,'Aggregate Analysis (combines runs into one figure)')
    group.add_option('','--aggregate',\
        action='append',default=[],\
        type='string',metavar='STR1,STR2,',
        help='plot aggregate information (takes a list of prefixes)')
    group.add_option('','--aggregate_names',\
        action='append',default=[],\
        type='string',metavar='STR1,STR2,',
        help='plot aggregate information (takes a list of dataset names)')
    parser.add_option_group(group)

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    # checks consistency of arguments.
    if opts.union and not opts.single:
        sys.exit('ERROR: if --union is specified, --single must also be specified. Exiting.')
    if (opts.aggregate_names and not opts.aggregate) or (opts.aggregate and not opts.aggregate_names):
        sys.exit('ERROR: --aggregate requires --aggregate_names to be specified (and vice versa). Exiting.')
    if not (opts.stats or opts.single or opts.auc or opts.roc or opts.aggregate):
        sys.exit('ERROR: method does not specify one or more experiments to run.  At least one of the following methods/experiments must be specified:\n\t--stats: compute statistics\n\t--single: run single experiment\n\t--auc: run k-fold cross validation\n\t--roc: compute ROC of positives that appear in both curated sets.\nExiting.')
    return opts

def main(argv):
    opts = parse_arguments(argv)
    print(opts)

    if opts.stats:
        sys.exit('ERROR: --stats option not available yet.')

    ## if opts.aggregate, don't need to load input datasets. 
    ## for all other experiments, we need them.
    if not opts.aggregate or opts.stats or opts.single or opts.auc or opts.roc:
        print('\nLoading Data...')

        print(" reading gene map file %s..." % (opts.gene_map_file))
        genemap = fileIO.geneMapReader(opts.gene_map_file)

        print(" reading edge file %s..." % (opts.interaction_graph))
        G = nx.Graph()
        layers=3
        fileIO.read_edge_file(opts.interaction_graph,G, opts.layers)
        print(' ',G.number_of_edges(), 'edges')
        print(' ',G.number_of_nodes(), 'nodes')

        print(' reading positive and negative files %s %s %s...' % (opts.disease_positives, opts.biological_process_positives, opts.negatives))
        minimum_labeled=25825
        disease_positives, minimum_labeled = fileIO.curatedFileReader(opts.disease_positives,G,opts.verbose, minimum_labeled, opts.layers)
        biological_process_positives, minimum_labeled = fileIO.curatedFileReader(opts.biological_process_positives,G,opts.verbose, minimum_labeled, opts.layers)
        # negatives, minimum_labeled = fileIO.curatedFileReader(opts.negatives,G,opts.verbose, opts.single_layer, minimum_labeled, opts.gene_mania)
        negatives=set()

        ## some nodes appear in both positive and negative sets; identify these and remove
        ## them from the curated set.
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

        print('Final Curated Sets: %d Disease Positives, %d Biological Process Positives, and %d Negatives.\n' % \
            (len(disease_positives),len(biological_process_positives),len(negatives)))
        print('%d nodes have been blacklisted because they were in both positive and negative sets.' % (len(blacklist)))

    ##########################
    ## opts.single: run three individual learning experiments; one with disease positives, one with
    ## biological process positives, and one with the union of the two positive sets.  
    ##
    ## Note that the learner functions won't overwrite the files if they exist unless the --force option is used.
    if opts.single:
        print('\nRunning Learning Algorithms...')

        print('Disease predictions...')
        statsfile = opts.outprefix + '_disease_stats.txt'
        outfile = opts.outprefix+'_disease_output.txt'
        name = 'disease'
        d_times,d_changes,d_predictions = learners.learn(outfile,statsfile,genemap,G,disease_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant, opts.layers,name,write=True)

        print('Biological process predictions...')
        statsfile = opts.outprefix + '_biological_process_stats.txt'
        outfile = opts.outprefix+'_biological_process_output.txt'
        name = 'process'
        b_times,b_changes,b_predictions = learners.learn(outfile,statsfile,genemap,G,biological_process_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,name,write=True)

        ## write combined results for disease and biological process predictions, including the final score 
        ## which is the product of the two sets of predictions.
        outfile = opts.outprefix+'_combined_output.txt'
        fileIO.writeCombinedResults(G,outfile,d_predictions,b_predictions,\
            disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers)

        ## Formats the results as a LaTeX table.
        outfile = opts.outprefix+'_formatted.txt'
        fileIO.formatCombinedResults(G,outfile,d_predictions,b_predictions,disease_positives,biological_process_positives,negatives,blacklist,genemap, opts.layers)

        ## Plot the time course of the runs per iteration.
        plt.clf()
        plt.plot(range(len(d_times)),d_times,'-r',label='Disease Predictions')
        plt.plot(range(len(b_times)),b_times,'-b',label='Biological Process Predictions')
        plt.legend()
        plt.ylabel('Time Per Iteration')
        plt.xlabel('Iteration')
        plt.savefig(opts.outprefix+'_timeCourse.png')
        print('wrote to '+opts.outprefix+'_timeCourse.png')

        ## Plot the change per iteration
        plt.clf()
        plt.plot(range(len(d_changes)),[math.log(x) for x in d_changes],'-r',label='Disease Predictions')
        plt.plot(range(len(b_changes)),[math.log(x) for x in b_changes],'-b',label='Biological Process Predictions')
        plt.legend()
        plt.ylabel('Log Absolute Value Change Per Iteration')
        plt.xlabel('Iteration')
        plt.savefig(opts.outprefix+'_scoreChange.png')
        print('wrote to '+opts.outprefix+'_scoreChange.png')

        ##########################
        ## opts.union: run the experiment for the union of disease and biological process sets. 
        ## Can only run when opts.single is specified.
        if opts.union:
            print(' union of positives...')
            statsfile = opts.outprefix + '_union_stats.txt'
            outfile = opts.outprefix+'_union_output.txt'
            union_positives = disease_positives.union(biological_process_positives)
            u_times,u_changes,u_predictions = learners.learn(outfile,statsfile,genemap,G,union_positives,negatives,\
                opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,write=True)

            ## plot node rankings.
            plt.clf()
            names = ['Disease Predictor $f_{\mathcal{D}}$','Biological Process Predictor $f_{\mathcal{P}}$','Union Predictor','Score $g$']
            colors =['g','b','k','r']
            n = G.number_of_nodes()
            preds = [d_predictions,b_predictions,u_predictions,{x:d_predictions[x]*u_predictions[x] for x in G.nodes()}]
            for i in range(len(names)):
                yvals =sorted(preds[i].values(),reverse=True)
                plt.plot(range(n),yvals,color=colors[i],label=names[i],lw=2)
            plt.plot([0,n-1],[0.5,0.5],':k',label='_nolabel_')
            plt.legend()
            plt.xlabel('Node ($n=%s$)' % (n))
            plt.ylabel('Ranking ($f$ or $g$)')
            plt.xlim(0,n-1)
            plt.ylim(0.01,1.01)
            plt.savefig(opts.outprefix+'_nodeRankings.png')
            pritn('wrote to '+opts.outprefix+'_nodeRankings.png')

    ##########################
    ## opts.auc: subsample 1/opts.k_fold of the positives for disease, run the matrix method, calculate
    ## the AUC value, and append it to a list.  Plot the distribution of these AUCs.
    if opts.auc:
        print('\nCalculating AUC w/ matrix method...')
        outfile = opts.outprefix+'_auc.txt'
        if not opts.force and os.path.isfile(outfile):
            print('  File %s exists. Not running (use --force to override)' % (outfile))

            ## read AUC lists from the existing file.
            d_AUCs = []
            b_AUCs = []
            with open(outfile,'r') as fin:
                for line in fin:
                    if line[0]=='#': # skip header
                        continue
                    row = line.strip().split()
                    d_AUCs.append(float(row[0]))
                    b_AUCs.append(float(row[1]))
        else:
            ## disease k-fold validation
            d_AUCs = []
            for i in range(opts.auc_samples):
                print('#%d of %d' % (i,opts.auc_samples))
                # subsample 1/k of the positives...
                hidden_genes = random.sample(disease_positives,int(len(disease_positives)/opts.k_fold))
                test_positives = disease_positives.difference(hidden_genes)
                print('%d hidden %d test genes' % (len(hidden_genes),len(test_positives)))
                ignore,ignore,d_predictions = learners.matrixLearn(G,test_positives,negatives,\
                    opts.epsilon,opts.timesteps,opts.verbose, sinksource_constant)
                d_AUCs.append(Mann_Whitney_U_test(d_predictions, hidden_genes,negatives))

            ## biological process k-fold validation
            b_AUCs = []
            for i in range(opts.auc_samples):
                print('#%d of %d' % (i,opts.auc_samples))
                # subsample 1/k of the positives...
                hidden_genes = random.sample(biological_process_positives,int(len(biological_process_positives)/opts.k_fold))
                test_positives = biological_process_positives.difference(hidden_genes)
                print('%d hidden %d test genes' % (len(hidden_genes),len(test_positives)))
                ignore,ignore,b_predictions = learners.matrixLearn(G,test_positives,negatives,\
                    opts.epsilon,opts.timesteps,opts.verbose)
                b_AUCs.append(Mann_Whitney_U_test(b_predictions, hidden_genes,negatives))

            ## write the output file.
            out = open(outfile,'w')
            out.write('#DiseaseAUCs\tBiologicalProcessAUCs\n')
            for i in range(len(d_AUCs)):
                out.write('%f\t%f\n' % (d_AUCs[i],b_AUCs[i]))
            out.close()

        ## plot the AUC distribution.
        plt.clf()
        plt.boxplot([d_AUCs,b_AUCs])
        plt.xticks([1,2],['SZ','Cell Motility'])
        plt.ylabel('AUC')
        plt.ylim([0,1])
        plt.title('5-Fold Cross Validation (AUC of 50 Iterations)')
        plt.savefig(opts.outprefix+'_auc.png')
        print('wrote to '+opts.outprefix+'_auc.png')

    ##########################
    ## opts.roc: "hide" the positives that appear in both disease and biological process
    ## positives and see how well the methods recover the hidden positives.  Unlike the AUC
    ## experiment, there's only one run of each method (we're not subsampling).
    if opts.roc:
        print('\nHolding out overlap set and running procedure.')
        hidden_genes = disease_positives.intersection(biological_process_positives)
        test_biological_process_positives = biological_process_positives.difference(hidden_genes)
        test_disease_positives = disease_positives.difference(hidden_genes)

        print('%d hidden genes, %d test disease genes, and %d test biological process genes' % \
            (len(hidden_genes),len(test_disease_positives),len(test_biological_process_positives)))
        print(' disease predictions...')
        statsfile = opts.outprefix + '_holdout_disease_stats.txt'
        outfile = opts.outprefix+'_holdout_disease_output.txt'
        ignore,ignore,holdout_d_predictions = learners.learn(outfile,statsfile,genemap,G,test_disease_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,write=True)

        print(' biological process predictions...')
        statsfile = opts.outprefix + '_holdout_biological_process_stats.txt'
        outfile = opts.outprefix+'_holdout_biological_process_output.txt'
        ignore,ignore,holdout_b_predictions = learners.learn(outfile,statsfile,genemap,G,test_biological_process_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,write=True)

        ## write combined results for disease and biological process predictions, including the final score 
        ## which is the product of the two sets of predictions.
        outfile = opts.outprefix+'_holdout_combined_output.txt'
        fileIO.writeCombinedResults(G,outfile,holdout_d_predictions,holdout_b_predictions,\
            test_disease_positives,test_biological_process_positives,negatives,hidden_genes,genemap)

        ## plot ROC.
        names = ['Disease Predictor $f_{\mathcal{D}}$','Biological Process Predictor $f_{\mathcal{P}}$','Score $g$']
        colors =['g','b','r']
        preds = [holdout_d_predictions,holdout_b_predictions,{x:holdout_d_predictions[x]*holdout_b_predictions[x] for x in G.nodes()}]
        test_union_positives=test_disease_positives.union(test_biological_process_positives)
        pos = [test_disease_positives,test_biological_process_positives,test_union_positives]
        plt.clf()
        for i in range(len(names)):
            x,y = getROCvalues(preds[i],pos[i],hidden_genes)
            plt.plot(x,y,color=colors[i],label=names[i])
            AUC = Mann_Whitney_U_test(preds[i], hidden_genes,negatives)
            print(names[i],AUC)
        plt.xlabel('# False Positives')
        plt.ylabel('# True Positives')
        plt.title('Receiver Operator Characteristic (ROC)')
        plt.legend(loc='lower right')
        plt.savefig(opts.outprefix+'_ROC.png')
        print('wrote to '+opts.outprefix+'_ROC.png')


    ##########################
    ## opts.aggregate: take all AUC files for the list of netorks and plot the result.
    if opts.aggregate:
        print('\nAggregating information')
        files = opts.aggregate
        names = opts.aggregate_names
        # aggregate AUCs
        plt.clf()
        plt.figure(figsize=(8,4))
        AUCs = []
        AUC_names = []
        for i in range(len(files)):
            d = []
            b = []
            with open(files[i]+'_auc.txt') as fin:
                for line in fin:
                    if line[0]=='#':
                        continue
                    row = line.strip().split()
                    d.append(float(row[0]))
                    b.append(float(row[1]))
            print(names[i]+' Disease Average:',sum(d)/len(d))
            print(names[i]+' Biological Process Average:',sum(b)/len(b))
            AUCs = AUCs + [d,b]
            AUC_names = AUC_names + [names[i]+'\n$\mathcal{D}$',names[i]+'\n$\mathcal{P}$']

        bplot = plt.boxplot(AUCs,patch_artist=True)
        for patch in bplot['boxes']:
            patch.set_facecolor('lightblue')
        plt.xticks(range(1,len(files)*2+1),AUC_names)
        plt.ylabel('AUC')
        plt.ylim([0,1])
        plt.title('5-Fold Cross Validation (AUC of 50 Iterations)')
        plt.savefig(opts.outprefix+'_aggregate_auc.png')
        print('wrote to '+opts.outprefix+'_aggregate_auc.png')


        plt.clf()
        vals = []
        colors =['g','b','k','r','y','c','m']
        for i in range(len(files)):
            vals.append([])
            with open(files[i]+'_combined_output.txt') as fin:
                for line in fin:
                    if line[0]=='#':
                        continue
                    vals[i].append(float(line.strip().split()[6]))
        for i in range(len(vals)):
            yvals =sorted(vals[i],reverse=True)
            plt.plot(range(len(vals[i])),yvals,color=colors[i],label=names[i],lw=2)
        plt.plot([0,len(yvals)-1],[0.5,0.5],':k',label='_nolabel_')
        plt.legend()
        plt.xlabel('Node')
        plt.ylabel('Combined Score $g$')
        plt.xlim(0,500)
        plt.ylim(0.01,1.01)
        plt.savefig(opts.outprefix+'_aggregate_nodeRankings.png')
        print('wrote to '+opts.outprefix+'_aggregate_nodeRankings.png')
    print('Done.')

    return

def getROCvalues(preds,pos,hidden):
    x = [0]
    y = [0]
    sorted_preds = sorted(preds, key=lambda x:preds[x], reverse=True)
    runningx = 0
    runningy = 0
    for node in sorted_preds:
        if node in hidden:
            runningy += 1
        elif node not in pos: # ignore positives
            runningx += 1
        x.append(runningx)
        y.append(runningy)
        if runningy == len(hidden):
            x.append(1)
            y.append(runningy)
            break
    return x,y

def Mann_Whitney_U_test(predictions, hidden_nodes, negatives):
    #Runs a Mann-Whitney U test on the lists
    kfold_ranks=[] #This will be the results we have computed without the positives
    test_ranks=[] #This will be the results we have computed with all positives

    nodeValues=[] #This is the vehicle by which we extract graph information
    hiddenNodeValues=[]
    notPositiveNodeValues=[]

    for node in predictions:
        # if node[-6:] == '_prime':
        if node in hidden_nodes:
            hiddenNodeValues.append(predictions[node])
        else:
            notPositiveNodeValues.append(predictions[node])

    U, p=stats.mannwhitneyu(hiddenNodeValues, notPositiveNodeValues, alternative="two-sided")
    #print(U,p)
    AUC=U/(len(hiddenNodeValues)*len(notPositiveNodeValues))
    #print(AUC)
    return AUC












if __name__ == '__main__':
    main(sys.argv)
