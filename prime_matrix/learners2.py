import os.path
import time
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

import fileIO2 as fileIO

def learn(outfile,statsfile,genemap,G,pos,negs,epsilon,timesteps,iterative_update,verbose,force,single_layer, write=False):
    if not force and os.path.isfile(outfile):
        print('  File %s exists. Not running (use --force to override)' % (outfile))
        times,changes,predictions = fileIO.readResults(statsfile,outfile)
    else:
        if not iterative_update:
            times,changes,predictions = matrixLearn(G,pos,negs,epsilon,timesteps,verbose)

        else:
            setGraphAttrs(G,pos,negs)
            times,changes,predictions = iterativeLearn(G,epsilon,timesteps,verbose)
        if write:
            if single_layer:
                fileIO.writeResultsSingle(statsfile,outfile,times,changes,predictions,genemap, G)
            else:
                fileIO.writeResults(statsfile,outfile,times,changes,predictions,genemap, G)
    return times,changes,predictions

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

def matrixLearn2(G,pos,neg,epsilon,timesteps,verbose):

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
    for i in range(len(unlabeled_list)):
        if unlabeled_list[i][-6:] == '_prime':
            predictions[unlabeled_list[i][:-6]] = f[i]

    return timeLogger,changeLogger, predictions


def matrixLearn3(G,pos,neg,epsilon,timesteps,verbose):

    pos1, pos2, pos3 = pos[0], pos[1], pos[2]

    ## Takes the form of f = M * f + c.

    ## sort unlabeled nodes.
    unlabeled1, unlabeled2, unlabeled3 = set(G.nodes()).difference(pos[0]).difference(neg), set(G.nodes()).difference(pos[1]).difference(neg), set(G.nodes()).difference(pos[2]).difference(neg)
    unlabeled_list1, unlabeled_list2, unlabeled_list3  = sorted(unlabeled1), sorted(unlabeled2), sorted(unlabeled3)
    unlabeled_inds1, unlabeled_inds2, unlabeled_inds3 = {unlabeled_list1[i]:i for i in range(len(unlabeled_list1))}, {unlabeled_list2[i]:i for i in range(len(unlabeled_list2))}, {unlabeled_list3[i]:i for i in range(len(unlabeled_list3))}
    print('%d unlabeled nodes' % (len(unlabeled1)+len(unlabeled2)+len(unlabeled3)))


    print('Preparing matrix data')

    #Make sparse M matrix.
    start = time.time()
    print(' making M matrix...')
    # from https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    data_for_M1 = {}
    data_for_M2 = {}
    data_for_M3 = {}
    for u,v in G.edges():
        if u in unlabeled1 and v in unlabeled1:
            i = unlabeled_inds1[u]
            j = unlabeled_inds1[v]
            data_for_M1[(i,j)] = float(G.edges[u,v]['weight'])/G.nodes[u]['weighted_degree']
            data_for_M1[(j,i)] = float(G.edges[u,v]['weight'])/G.nodes[v]['weighted_degree']
        if u in unlabeled2 and v in unlabeled2:
            i = unlabeled_inds2[u]
            j = unlabeled_inds2[v]
            data_for_M2[(i,j)] = float(G.edges[u,v]['weight'])/G.nodes[u]['weighted_degree']
            data_for_M2[(j,i)] = float(G.edges[u,v]['weight'])/G.nodes[v]['weighted_degree']
        if u in unlabeled3 and v in unlabeled3:
            i = unlabeled_inds3[u]
            j = unlabeled_inds3[v]
            data_for_M3[(i,j)] = float(G.edges[u,v]['weight'])/G.nodes[u]['weighted_degree']
            data_for_M3[(j,i)] = float(G.edges[u,v]['weight'])/G.nodes[v]['weighted_degree']


    keys1 = data_for_M1.keys()
    data1 = [data_for_M1[key] for key in keys1]
    row1 = [key[0] for key in keys1]
    col1 = [key[1] for key in keys1]

    keys1 = data_for_M1.keys()
    data1 = [data_for_M1[key] for key in keys1]
    row1 = [key[0] for key in keys1]
    col1 = [key[1] for key in keys1]



    n = len(unlabeled)
    M = coo_matrix((data, (row,col)), shape=(n,n))
    M = M.tocsr()
    end = time.time()
    print(' %f seconds' % (end-start))
    #Make c vector
    print(' making c vectors...')
    cs = [[0]*len(unlabeled_list), [0]*len(unlabeled_list), [0]*len(unlabeled_list)]
    for impact in cs:
        for i in range(len(unlabeled_list)):
            v = unlabeled_list[i]
            for neighbor, datadict in G.adj[v].items():
                if neighbor not in unlabeled: # it is labeled
                    if neighbor in pos1:
                        impact[i] += datadict['weight']
                    else: # neighbor in negative; label is 0 so nothing is added.
                        pass

            impact[i] = float(c1[i])/G.nodes[v]['weighted_degree']

    #Make initial f vector.  This is a value of 0.5 for all unlabeled nodes.
    f1, f2, f3 = [0.5]*len(unlabeled_list), [0.5]*len(unlabeled_list), [0.5]*len(unlabeled_list)

    changeLogger=[]
    timeLogger=[]
    for t in range(0,timesteps):

        ## conduct sparse matrix operation.
        f_prev = f
        start = time.time()
        f1 = M.dot(f1)+c1
        f2 = M.dot(f2)+c1
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