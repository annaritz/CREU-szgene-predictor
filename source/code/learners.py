#Contains all the algorithms for calculating scores - iterative, matrix, matrix with sinksource+ constant, etc.
import os.path
import time
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from numpy import diagflat
import numpy as np

import fileIO

def learn(outprefix,outfile,statsfile,genemap,G,pos,negs,\
    epsilon,timesteps,iterative_update,verbose,force,sinksource_constant,layers,name,sinksource_method,write=False):
    if not force and os.path.isfile(outfile):
        print('  File %s exists. Not running (use --force to override)' % (outfile))
        times,changes,predictions = fileIO.readResults(statsfile,outfile,layers)
    else:
        if not iterative_update:
            if not sinksource_method:
                times,changes,predictions = matrixLearn(G,pos,negs,epsilon,timesteps,verbose)
            else:
                times,changes,predictions = matrixLearnSinkSource(G,pos,negs,epsilon,timesteps,verbose,sinksource_constant)
        else:
            setGraphAttrs(G,pos,negs) #intitializes the scores
            times,changes,predictions = iterativeLearn(G,epsilon,timesteps,verbose)
        if write:
            fileIO.writeResults(statsfile,outprefix,outfile,times,changes,predictions,genemap, G, layers,pos, negs,sinksource_constant, name, sinksource_method)
    return times,changes,predictions

# Utilizes sinksource constant 
def matrixLearnSinkSource(G,pos,neg,epsilon,timesteps,verbose, sinksource_constant):

    ## Takes the form of f = M * f + c where each entry for M 
    ## is calculated by dividing the weight by the weighted degree + constant from sinksource+ algorithm

    ## sort unlabeled nodes.
    unlabeled = set(G.nodes()).difference(pos).difference(neg)
    unlabeled_list = sorted(unlabeled)
    unlabeled_inds = {unlabeled_list[i]:i for i in range(len(unlabeled_list))}
    print('%d unlabeled nodes.' % (len(unlabeled)))

    print('Preparing matrix data')

    #Make sparse M matrix.
    start = time.time()
    print(' making M matrix for sinksource+...')
    # from https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    data_for_M = {}
    for u,v in G.edges():
        if u in unlabeled and v in unlabeled:
            i = unlabeled_inds[u]
            j = unlabeled_inds[v]

            data_for_M[(i,j)] = float(G.edges[u,v]['weight'])/(G.nodes[u]['weighted_degree']+sinksource_constant)
            data_for_M[(j,i)] = float(G.edges[u,v]['weight'])/(G.nodes[v]['weighted_degree']+sinksource_constant)
    keys = data_for_M.keys()
    data = [data_for_M[key] for key in keys]
    row = [key[0] for key in keys]
    col = [key[1] for key in keys]
    n = len(unlabeled)
    M = coo_matrix((data, (row,col)), shape=(n,n))
    M = M.tocsr()
    end = time.time()
    print(' %f seconds' % (end-start))
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

        c[i] = float(c[i])/(G.nodes[v]['weighted_degree']+sinksource_constant)

    #Make initial f vector.  This is a value of 0.5 for all unlabeled nodes.
    f = [0.5]*len(unlabeled_list)

    changeLogger=[]
    timeLogger=[]
    for t in range(0,timesteps):

        ## conduct sparse matrix operation.
        f_prev = f
        start = time.time()
        f = M.dot(f)+c
        end = time.time()

        ## sum changes
        changes = sum([abs(f[i]-f_prev[i]) for i in range(len(f))])

        timeLogger.append(end-start)
        changeLogger.append(changes)
        if t % 10 == 0:
            print("t = %d: change = %.4f" % (t,changes))

        if changes < epsilon:
            print('BELOW THRESHOLD OF %.2e! Breaking out of loop.' % (epsilon))
            break

        if verbose:
            done=float(t)/float(timesteps)
            print('Time Elapsed:', end-start)
            if done!=0:
                print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')
    print('Average Mult. Time: %f seconds' % (sum(timeLogger)/len(timeLogger)))

    # predictions is a dictionary of nodes to values.
    predictions = {}
    for n in pos:
        predictions[n] = 1
    for n in neg:
        predictions[n] = 0
    for i in range(len(unlabeled_list)):
        predictions[unlabeled_list[i]] = f[i]
    return timeLogger,changeLogger, predictions


# Does not utilize sinksource constant (uses original formula)
def matrixLearn(G,pos,neg,epsilon,timesteps,verbose):

    ## Takes the form of f = M * f + c.

    ## sort unlabeled nodes.
    unlabeled = set(G.nodes()).difference(pos).difference(neg)
    unlabeled_list = sorted(unlabeled)
    unlabeled_inds = {unlabeled_list[i]:i for i in range(len(unlabeled_list))}
    print('%d unlabeled nodes.' % (len(unlabeled)))

    print('Preparing matrix data')

    #Make sparse M matrix.
    start = time.time()
    print(' making M matrix...')
    # from https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    data_for_M = {}
    for u,v in G.edges():
        if u in unlabeled and v in unlabeled:
            i = unlabeled_inds[u]
            j = unlabeled_inds[v]
            data_for_M[(i,j)] = float(G.edges[u,v]['weight'])/G.nodes[u]['weighted_degree']
            data_for_M[(j,i)] = float(G.edges[u,v]['weight'])/G.nodes[v]['weighted_degree']
    keys = data_for_M.keys()
    data = [data_for_M[key] for key in keys]
    row = [key[0] for key in keys]
    col = [key[1] for key in keys]
    n = len(unlabeled)
    M = coo_matrix((data, (row,col)), shape=(n,n))
    M = M.tocsr()
    end = time.time()
    print(' %f seconds' % (end-start))
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

        ## conduct sparse matrix operation.
        f_prev = f
        start = time.time()
        f = M.dot(f)+c
        end = time.time()

        ## sum changes
        changes = sum([abs(f[i]-f_prev[i]) for i in range(len(f))])

        timeLogger.append(end-start)
        changeLogger.append(changes)
        if t % 10 == 0:
            print("t = %d: change = %.4f" % (t,changes))

        if changes < epsilon:
            print('BELOW THRESHOLD OF %.2e! Breaking out of loop.' % (epsilon))
            break

        if verbose:
            done=float(t)/float(timesteps)
            print('Time Elapsed:', end-start)
            if done!=0:
                print('Estimated Time Remaining:', (1.0-done)*(time.time()-start)/done, 'seconds')
    print('Average Mult. Time: %f seconds' % (sum(timeLogger)/len(timeLogger)))



    # predictions is a dictionary of nodes to values.
    predictions = {}
    for n in pos:
        predictions[n] = 1
    for n in neg:
        predictions[n] = 0
    for i in range(len(unlabeled_list)):
        predictions[unlabeled_list[i]] = f[i]

    return timeLogger,changeLogger, predictions



# Utilizes iterative method instead of matrix method
# Called when --iterative_update is specified (matrix method is default)
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
    for n in pos:
        predictions[n] = 1
    for n in neg:
        predictions[n] = 0
    for i in range(len(unlabeled_list)):

        predictions[unlabeled_list[i]] = f[i]


    return timeLogger,changeLogger, predictions

#Takes as input a networkx graph
#Verbose set to False by default, --verbose True will print more descriptive statements of nodes and scores changed in each iteration
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


