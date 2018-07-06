## Full program that calls all parts of the analysis.

## this needs to be at the top so matplotlib doesn't
## throw an error if being called remotely.
import matplotlib
matplotlib.use('Agg')

## system and utility packages
import sys
import os.path
import random
random.seed('Alexander_King') #TODO: every once in a while comment this out and re-run a couple times. 

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
import fileIO
import learners


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
        type='string',metavar='STR',default='../outfiles/out',\
        help='output file prefixes (default="../outfiles/out").')
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
        type='string',metavar='STR',default='../infiles/networkfiles/brain_top_geq_0.400.txt',\
        help='Functional interaction network (default="../infiles/networkfiles/brain_top_geq_0.400.txt").')
    group.add_option('-b','--biological_process_positives',\
        type='string',metavar='STR',default='../infiles/motility_positives.txt',\
        help='File of positives for the biological process (default="../infiles/motility_positives.txt")')
    group.add_option('-d','--disease_positives',\
        type='string',metavar='STR',default='../infiles/SZ_positives.txt',\
        help='File of positives for the disease (default="../infiles/SZ_positives.txt")')
    group.add_option('-n','--negatives',\
        type='string',metavar='STR',default='../infiles/SZ_negatives.txt',\
        help='File of negatives (default="../infiles/SZ_negatives.txt")')
    group.add_option('-m','--gene_map_file',\
        type='string',metavar='STR',default='../infiles/Homo_sapiens.txt',\
        help='File of EntrezID to Common Names (default="../infiles/Homo_sapiens.txt")')
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
    group.add_option('','--sinksource_method',\
        action='store_true',default=False,\
        help='Run the sinksource+ algorithm (default=False).')
    group.add_option('-k','--k_fold',\
        type='int',metavar='INT',default=5,\
        help='k for k-fold cross validation (default=5).')
    group.add_option('-s','--auc_samples',\
        type='int',metavar='INT',default=50,\
        help='number of cross validation iterations to compute AUC (default=50).')
    group.add_option('-l', '--layers',\
        type='int', default=4,\
        help='Run the experiments for disease and biological process with n nodes per gene. Each of the gene\'s nodes are connected to a node that represents the score for the gene. Positives are distributed among these layers. Reducing the nodes to 1 eliminates the process')
    group.add_option('-c', '--sinksource_constant',\
        type='float',metavar='FLOAT',default=None,\
        help='How much to add to the denominator of the node value (must set)')
    group.add_option('-w', '--with_negatives',\
        action='store_false',default=True,\
        help='Include negative set (default=True).')
    group.add_option('-E', '--evidence_levels',\
        action='store_false', default=True,\
        help='Factor evidence levels for gold standards (default=True)')
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
    if opts.sinksource_constant != None and opts.sinksource_method == False:
        sys.exit('ERROR: sinksource method must be specified to use constant.')
    if opts.sinksource_method==True and opts.sinksource_constant==None:
        sys.exit('ERROR: must specify constant for sinksource method.')
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
        fileIO.read_edge_file_single(opts.interaction_graph,G)
    
        print(' The network contains %d edges and %d nodes' % (G.number_of_edges(), G.number_of_nodes()))

        print(' reading positive and negative files %s %s %s...' % (opts.disease_positives, opts.biological_process_positives, opts.negatives))
        if opts.evidence_levels:
            disease_positives, disease_evidence_dictionary = fileIO.curatedEvidenceFileReader(opts.disease_positives,G,opts.verbose)
            autism_positives, autism_evidence_dictionary = fileIO.curatedEvidenceFileReader('../infiles/ASD_positives_E.txt',G,opts.verbose)
            biological_process_positives, biological_process_evidence_dictionary = fileIO.curatedEvidenceFileReader(opts.biological_process_positives,G,opts.verbose)
        else:
            disease_positives= fileIO.curatedFileReader(opts.disease_positives,G,opts.verbose)
            autism_positives= fileIO.curatedFileReader('../infiles/ASD_positives.txt', G, opts.verbose)
            biological_process_positives= fileIO.curatedFileReader(opts.biological_process_positives,G,opts.verbose)

        if opts.with_negatives: #if True (default) pass in negative set and remove overlapping positives and negatives
            if opts.evidence_levels:
                negatives, negatives_evidence_dictionary= fileIO.curatedEvidenceFileReader(opts.negatives,G,opts.verbose)
            else:
                negatives= fileIO.curatedFileReader(opts.negatives,G,opts.verbose)
            ## some nodes appear in both positive and negative sets; identify these and remove
            ## them from the curated set.
            blacklist = set()
            if disease_positives.intersection(negatives):
                overlap_set = disease_positives.intersection(negatives)
                print('WARNING: %d genes are disease positives and negatives. Ignoring.' % (len(overlap_set)))
                negatives = negatives.difference(overlap_set)
                disease_positives = disease_positives.difference(overlap_set)

                blacklist.update(overlap_set)

            if autism_positives.intersection(negatives):
                overlap_set = autism_positives.intersection(negatives)
                print('WARNING: %d genes are autism positives and negatives. Ignoring.' % (len(overlap_set)))
                negatives = negatives.difference(overlap_set)
                autism_positives = autism_positives.difference(overlap_set)

                blacklist.update(overlap_set)

            if biological_process_positives.intersection(negatives):
                overlap_set = biological_process_positives.intersection(negatives)
                print('WARNING: %d genes are biological process positives and negatives. Ignoring.' % (len(overlap_set)))
                negatives = negatives.difference(overlap_set)
                biological_process_positives = biological_process_positives.difference(overlap_set)
                blacklist.update(overlap_set)

            print('%d genes have had their labels removed because they were in both positive and negative sets.' % (len(blacklist)))
            # for node in blacklist:
            #     G.remove_node(node)
       
        else: #if opts.with_negatives is False, it'll be an empty set (no negatives)
            negatives = set()

        #Rename the variables of the original positives and partitioned positives in order to use code for single and multi layer
        #and keep track of original and partitioned positives 

        #After checking the graph, we need to modify it to include multiple layers + primes
        multi_node_dict = fileIO.read_edge_file_multi(G, opts.layers)

        orig_disease_positives = disease_positives
        orig_autism_positives = autism_positives
        orig_biological_process_positives = biological_process_positives
        orig_negatives = negatives


        if opts.evidence_levels:
            disease_positives = fileIO.partitionEvidenceCurated(orig_disease_positives,G,opts.verbose,opts.layers, disease_evidence_dictionary)
            autism_positives = fileIO.partitionEvidenceCurated(orig_autism_positives,G,opts.verbose,opts.layers, autism_evidence_dictionary)
            biological_process_positives = fileIO.partitionEvidenceCurated(orig_biological_process_positives,G,opts.verbose,opts.layers, biological_process_evidence_dictionary)
            negatives = fileIO.partitionEvidenceCurated(orig_negatives,G,opts.verbose,opts.layers, negatives_evidence_dictionary)

        else:
            disease_positives = fileIO.partitionCurated(orig_disease_positives,G,opts.verbose,opts.layers)
            autism_positives = fileIO.partitionCurated(orig_autism_positives,G,opts.verbose,opts.layers)
            biological_process_positives = fileIO.partitionCurated(orig_biological_process_positives,G,opts.verbose,opts.layers)
            negatives = fileIO.partitionCurated(orig_negatives,G,opts.verbose,opts.layers)

            


        print('Final Curated Sets: %d Labeled Disease Nodes, %d Labeled Autism Nodes,%d Labeled Biological Process Nodes, and %d Labeled Negative Nodes.\n' % \
        (len(disease_positives),len(autism_positives),len(biological_process_positives),len(negatives)))
        



    ##########################
    ## opts.single: run three individual learning experiments; one with disease positives, one with
    ## biological process positives, and one with the union of the two positive sets.  
    ##
    ## Note that the learner functions won't overwrite the files if they exist unless the --force option is used.
    ## TODO: single is a bad name -- individual? 
    if opts.single:
        print('\nRunning Learning Algorithms...')
       
        print('Disease predictions...')
        statsfile = opts.outprefix + str(opts.layers) + 'layers_disease_stats.txt'
        outfile = opts.outprefix + str(opts.layers) + 'layers_disease_output.txt'
        name = 'disease'
        d_times,d_changes,d_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,disease_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,name,opts.sinksource_method,write=True)

        print('Biological process predictions...')
        statsfile = opts.outprefix + str(opts.layers) + 'layers_biological_process_stats.txt'
        outfile = opts.outprefix + str(opts.layers) + 'layers_biological_process_output.txt'
        name = 'process'
        b_times,b_changes,b_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,biological_process_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,name,opts.sinksource_method,write=True)

        ## NEW 6/14 by Anna: normed=True means that both predictions are normalized so the maximum is 1.0.
        normed=True
        ## write combined results for disease and biological process predictions, including the final score 
        ## which is the product of the two sets of predictions.
        outfile = opts.outprefix+'_combined_output.txt'
        fileIO.writeCombinedResults(G,outfile,d_predictions,b_predictions,\
            disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed)

        ## Formats the results as a LaTeX table.
        outfile = opts.outprefix+'_formatted.txt'
        fileIO.formatCombinedResults(G,outfile,d_predictions,b_predictions,\
            disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed)
        outfile = opts.outprefix+'_formatted_unlabeled.txt'
        fileIO.formatCombinedResults_unlabeled(G,outfile,d_predictions,b_predictions,\
            disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed)

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
            u_times,u_changes,u_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,union_positives,negatives,\
                opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,\
                opts.layers,name,opts.sinksource_method,write=True)

            ## plot node rankings.
            names = ['Disease Predictor $f_{\mathcal{D}}$','Biological Process Predictor $f_{\mathcal{P}}$','Union Predictor $f_{\mathcal{D} \cup \mathcal{P}}$','Combined Score $g$']
            colors =['g','b','k','r']
            styles = ['--',':','-','-']
            if opts.layers > 1:
                nodeset = set([n for n in G.nodes() if n[-6:]=='_prime'])
            else:
                nodeset = G.nodes()

            outfile = opts.outprefix+'_combined_output_withunion.txt'
            fileIO.writeCombinedResults(G,outfile,d_predictions,b_predictions,\
                disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed,union_predictions=u_predictions)
            outfile = opts.outprefix+'_formatted_unlabeled_withunion.txt'
            fileIO.formatCombinedResults_unlabeled(G,outfile,d_predictions,b_predictions,\
                disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed,union_predictions=u_predictions)

            if normed:
                d_max = max([d_predictions[key] for key in nodeset])
                d_predictions = {key:d_predictions[key]/float(d_max) for key in nodeset}
                b_max = max([b_predictions[key] for key in nodeset])
                b_predictions = {key:b_predictions[key]/float(b_max) for key in nodeset}
                u_max = max([u_predictions[key] for key in nodeset])
                u_predictions = {key:u_predictions[key]/float(u_max) for key in nodeset}
            n = len(nodeset)
            preds = [{x:d_predictions[x] for x in nodeset},\
                {x:b_predictions[x] for x in nodeset},\
                {x:u_predictions[x] for x in nodeset},\
                {x:d_predictions[x]*b_predictions[x] for x in nodeset}]

            ## checking max value (should all be 1 since now normalizing)
            print('----')
            for i in range(len(preds)):
                print(names[i],max(preds[i].values()))
            print('----')

            plt.figure(figsize=(5,3))
            plt.clf()

            for i in range(len(names)):
                yvals = sorted(preds[i].values(),reverse=True)
                plt.plot(range(len(yvals)),yvals,color=colors[i],ls=styles[i],label=names[i],lw=2)

            ##for i in range(len(names)):
            ##    yvals =sorted(preds[i].values(),reverse=True)
            ##    ax.plot(range(n),yvals,color=colors[i],ls=styles[i],label=names[i],lw=2)
            #plt.plot([0,n-1],[0.5,0.5],':k',label='_nolabel_')
            plt.legend()
            plt.xlabel('Candidate Rank')
            plt.ylabel('Score ($f$ or $g$)')
            plt.title('Ranked Genes by Different Predictors')
            #plt.xlim(0,n-1)
            #plt.ylim(0.01,1.01)

            plt.tight_layout()
            plt.savefig(opts.outprefix+'_nodeRankings.png')
            print('wrote to '+opts.outprefix+'_nodeRankings.png')

            ## plot union scatter.
            plt.figure(figsize=(5,4))
            plt.clf()
            ## get nodes to plot
            top_num = 1000
            sorted_by_combo = sorted(preds[3],key=lambda x:preds[3][x],reverse=True)
            sorted_by_union = sorted(preds[2],key=lambda x:preds[2][x],reverse=True)
            nodes_to_plot = set(sorted_by_union[:top_num]).union(sorted_by_combo[:top_num])

            xpos=[]
            ypos=[]
            xbothpos=[]
            ybothpos=[]
            xunl=[]
            yunl=[]
            ob_union = []
            ob_combo = []
            ob_union2 = []
            ob_combo2 = []

            ## TODO: double-check we're using correct positive and negative sets for multiple layers.
            ## TODO: pull out all figures and put in chartmaker.py?

            for n in nodes_to_plot:
                names = [n] + [n[:-6]+'_'+str(x) for x in range(opts.layers)]
                #print(names)
                if any([x in disease_positives for x in names]):
                    d_pos = True
                else:
                    d_pos = False
                if any([x in biological_process_positives for x in names]):
                    b_pos = True
                else:
                    b_pos = False
                #print(n,b_pos,d_pos)
                if d_pos and b_pos:
                    xbothpos.append(sorted_by_combo.index(n))
                    ybothpos.append(sorted_by_union.index(n))
                elif d_pos or b_pos:
                    xpos.append(sorted_by_combo.index(n))
                    ypos.append(sorted_by_union.index(n))
                else:
                    xunl.append(sorted_by_combo.index(n))
                    yunl.append(sorted_by_union.index(n))
                if sorted_by_union.index(n) <= top_num and sorted_by_combo.index(n) > top_num:
                    ob_union.append(sorted_by_union.index(n))
                    ob_combo.append(sorted_by_combo.index(n))
                if sorted_by_union.index(n) > top_num and sorted_by_combo.index(n) <= top_num:
                    ob_union2.append(sorted_by_union.index(n))
                    ob_combo2.append(sorted_by_combo.index(n))
            #print(xunl,xbothpos,xpos)

            print('Out of Bounds for COMBO:')
            print('UNION (%d total): %d-%d avg %.2f; COMBO ranks from %d-%d avg %.2f' % (len(ob_union),min(ob_union),max(ob_union),sum(ob_union)/len(ob_union),min(ob_combo),max(ob_combo),sum(ob_combo)/len(ob_combo)))
            print('Out of Bounds for UNION:')
            print('COMBO (%d total): %d-%d avg %.2f; UNION ranks from %d-%d avg %.2f' % (len(ob_union2),min(ob_combo2),max(ob_combo2),sum(ob_combo2)/len(ob_combo2),min(ob_union2),max(ob_union2),sum(ob_union2)/len(ob_union2)))
            
            plt.scatter(xpos, ypos, alpha=.5, color=[.7,.7,1], label='SZ or CM Positive')
            plt.scatter(xbothpos, ybothpos, alpha=.5, color=[0,0,1], label='SZ and CM Positive')
            if len(xunl)>0:
                print(len(xunl))
                plt.scatter(xunl, yunl, alpha=.5, color=[.6,.6,.6], label='Unlabeled')
            
            plt.plot([0,top_num],[0,top_num],':k',label='_nolabel_')
            plt.plot([0,0,top_num,top_num,0,0],[top_num,top_num,top_num,0,0,top_num],':k',label='_nolabel_')
            ##for i in range(len(names)):
            ##    yvals =sorted(preds[i].values(),reverse=True)
            ##    ax.plot(range(n),yvals,color=colors[i],ls=styles[i],label=names[i],lw=2)
            #plt.plot([0,n-1],[0.5,0.5],':k',label='_nolabel_')
            plt.legend()

            
            plt.xlabel('Candidate Rank by Combined Predictor $g$')
            plt.ylabel('Candidate Rank by Union Predictor $f_{\mathcal{D} \cup \mathcal{P}}$')
            plt.title('Ranked Genes by Combined vs. Union Predictors')
            #plt.xlim(0,n-1)
            #plt.ylim(0.01,1.01)
            plt.tight_layout()
            plt.savefig(opts.outprefix+'_scatter_union_full.png')
            print('wrote to '+opts.outprefix+'_scatter_union_full.png')
            plt.xlim(0,top_num*1.5)
            plt.ylim(0,top_num*1.5)
            plt.tight_layout()
            plt.savefig(opts.outprefix+'_scatter_union.png')
            print('wrote to '+opts.outprefix+'_scatter_union.png')

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
            a_AUCs = []
            b_AUCs = []

            with open(outfile,'r') as fin:
                for line in fin:
                    if line[0]=='#': # skip header
                        continue
                    row = line.strip().split()
                    d_AUCs.append(float(row[0]))
                    a_AUCs.append(float(row[1]))
                    b_AUCs.append(float(row[12]))
        else:
            ## disease k-fold validation
            print('\nDisease Tests')
            d_AUCs = []
            start=time.time()
            
            for i in range(opts.auc_samples):
                done=float(i)/float(opts.auc_samples)
                if done!=0:
                    time_remaining=(3.0-done)*(time.time()-start)/done
                    print('Estimated Time Remaining:', time_remaining, 'seconds')
                    print('Estimated Time of Completion:', time.strftime('%I:%M:%S %p',time.localtime(time.time()+time_remaining)))
                print('#%d of %d' % (i+1,opts.auc_samples))
                # subsample 1/k of the positives...
                hidden_genes=random.sample(orig_disease_positives,int(len(orig_disease_positives)/opts.k_fold))
                hidden_nodes = set(node for gene in hidden_genes for node in multi_node_dict[gene] if node in disease_positives)
                test_positives = disease_positives.difference(hidden_nodes)
                print('%d hidden %d test nodes' % (len(hidden_nodes),len(test_positives)))
                if not opts.sinksource_method:
                    ignore,ignore,d_predictions = learners.matrixLearn(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose)
                else:
                    ignore,ignore,d_predictions = learners.matrixLearnSinkSource(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose, opts.sinksource_constant)
                MWU = Mann_Whitney_U_test(d_predictions, hidden_nodes, negatives, test_positives, multi_node_dict)
                print('Disease AUC = ', MWU)
                d_AUCs.append(MWU)

    
            # autism k-fold validation
            a_AUCs=[]
            start=time.time()
            done=0
            print('\nASD Tests')
            for i in range(opts.auc_samples):
                done=float(i)/float(opts.auc_samples)
                if done!=0:
                    time_remaining=(2.0-done)*(time.time()-start)/done
                    print('Estimated Time Remaining:', time_remaining, 'seconds')
                    print('Estimated Time of Completion:', time.strftime('%I:%M:%S %p',time.localtime(time.time()+time_remaining)))
                print('#%d of %d' % (i+1,opts.auc_samples))
                # subsample 1/k of the positives...
                hidden_genes = random.sample(orig_autism_positives,int(len(orig_autism_positives)/opts.k_fold))
                hidden_nodes = set(node for gene in hidden_genes for node in multi_node_dict[gene] if node in autism_positives)
                test_positives = autism_positives.difference(hidden_nodes)
                print('%d hidden %d test nodes' % (len(hidden_nodes),len(test_positives)))
                if not opts.sinksource_method:
                    ignore,ignore,a_predictions = learners.matrixLearn(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose)
                else:
                    ignore,ignore,a_predictions = learners.matrixLearnSinkSource(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose, opts.sinksource_constant)
                MWU = Mann_Whitney_U_test(a_predictions, hidden_nodes, negatives, test_positives, multi_node_dict)
                print('Autism AUC = ', MWU)
                a_AUCs.append(MWU)

            ## biological process k-fold validation
            print('\nBiological Process Tests')
            b_AUCs = []
            start=time.time()
            for i in range(opts.auc_samples):
                done=float(i)/float(opts.auc_samples)
                if done!=0:
                    time_remaining=(1.0-done)*(time.time()-start)/done
                    print('Estimated Time Remaining:', time_remaining, 'seconds')
                    print('Estimated Time of Completion:', time.strftime('%I:%M:%S %p',time.localtime(time.time()+time_remaining)))
                print('#%d of %d' % (i+1,opts.auc_samples))
                # subsample 1/k of the positives...
                hidden_genes = random.sample(orig_biological_process_positives,int(len(orig_biological_process_positives)/opts.k_fold))
                hidden_nodes = set(node for gene in hidden_genes for node in multi_node_dict[gene] if node in biological_process_positives)
                test_positives = biological_process_positives.difference(hidden_nodes)
                print('%d hidden %d test nodes' % (len(hidden_nodes),len(test_positives)))
                if not opts.sinksource_method:
                    ignore,ignore,b_predictions = learners.matrixLearn(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose)
                else:
                    ignore,ignore,b_predictions = learners.matrixLearnSinkSource(G,test_positives,negatives,\
                        opts.epsilon,opts.timesteps,opts.verbose, opts.sinksource_constant)
                MWU = Mann_Whitney_U_test(b_predictions, hidden_nodes, negatives, test_positives, multi_node_dict)
                print('Biological Process AUC = ', MWU)
                b_AUCs.append(MWU)
            
            ## write the output file.
            out = open(outfile,'w')
            out.write('#DiseaseAUCs\t#AutismAUCs\tBiologicalProcessAUCs\n')
            for i in range(len(d_AUCs)):
                out.write('%f\t%f\t%f\n' % (d_AUCs[i],a_AUCs[i],b_AUCs[i]))
    
            out.close()

        ## plot the AUC distribution.
        plt.clf()
        plt.boxplot([d_AUCs,a_AUCs,b_AUCs])
        plt.xticks([1,2],['SZ','ASD','Cell Motility'])
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
        ignore,ignore,holdout_d_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_disease_positives,negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,write=True)

        print(' biological process predictions...')
        statsfile = opts.outprefix + '_holdout_biological_process_stats.txt'
        outfile = opts.outprefix+'_holdout_biological_process_output.txt'
        ignore,ignore,holdout_b_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_biological_process_positives,negatives,\
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
            a = []
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

def Mann_Whitney_U_test(predictions, hidden_nodes, negatives, test_positives, layer_dict):
    #predictions is a dictionary of nodes:scores
    #Runs a Mann-Whitney U test on the lists
    kfold_ranks=[] #This will be the results we have computed without the positives
    test_ranks=[] #This will be the results we have computed with all positives

    nodeValues=[] #This is the vehicle by which we extract graph information
    hiddenNodeValues=[]
    notPositiveNodeValues=[] #Newest version: holds value of unlabeled prime nodes (positives and negative excluded)
    negative_count = 0
    positive_count = 0 
    sorted_preds = sorted(predictions, key=lambda x:predictions[x], reverse=True)
    #Iterate through layer_dict instead?
    i=1
    for node in sorted_preds:
        if node[-6:] == '_prime': #only want to look at prime nodes
            entrez = node[:-6] 
            names = layer_dict[entrez] #gives set of duplicate + prime names for a given entrez ID
            if bool(names.intersection(hidden_nodes)): #bool() is True if the prime node is attached to a hidden node, False if not
                hiddenNodeValues.append(predictions[node])
                if i < 100:
                    i=i+1
            else: #if it's not a hidden node, check if it's unlabeled
                if bool(names.intersection(test_positives)):
                    positive_count += 1
                    continue ## I think just check this one.
                if bool(names.intersection(negatives)):
                    negative_count += 1
                notPositiveNodeValues.append(predictions[node])
                if i < 100:

            

    print('Negative count: ', negative_count)
    print('Positive count: ', positive_count)
    print('# hidden nodes: ', len(hiddenNodeValues))
    print('# unlabeled prime nodes: ', len(notPositiveNodeValues))

    ## TODO suggestion: sys.exit() with an error if these numbesr aren't what we expect.  

    U, p=stats.mannwhitneyu(hiddenNodeValues, notPositiveNodeValues, alternative="two-sided")
    AUC=U/(len(hiddenNodeValues)*len(notPositiveNodeValues))
    #print(AUC)
    return AUC












if __name__ == '__main__':
    main(sys.argv)
