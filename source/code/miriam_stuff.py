'''
    ##########################
    ## New version (Miriam 07/10): Does the same thing as opts.auc (k-fold cross validation)
    ## but it plots the ROC curve for each run to see how it changes

    if opts.roc:
        #print('\nHolding out overlap set and running method.')

        #disease_predictions_list = [] #collect the disease predictions
        #disease_hidden_list = [] #collect all the hidden disease sets
        #disease_test_list = [] #collect all the test disease positive sets
        
        #process_predictions_list = [] #collect the process predictions
        #process_hidden_list = [] #collect all the hiddden process sets
        #process_test_list = [] #collect all the rest process positive sets

        #for j in range(opts.auc_samples): 
        #    print()
        #    print('Run #', j+1)
            #if opts.layers == 1:
            #    hidden_genes = disease_positives.intersection(biological_process_positives)
            #    test_biological_process_positives = biological_process_positives.difference(hidden_genes)
            #    test_disease_positives = disease_positives.difference(hidden_genes)
            #else: # adjust for layers
            #    d_pos = set([x.split('_')[0] for x in disease_positives])
            #    b_pos = set([x.split('_')[0] for x in biological_process_positives])
            #    h_genes = d_pos.intersection(b_pos)
            #    hidden_genes = set([x for x in disease_positives if x.split('_')[0] in h_genes])
            #    test_biological_process_positives = set([x for x in biological_process_positives if x.split('_')[0] not in h_genes])
            #    test_disease_positives = set([x for x in disease_positives if x.split('_')[0] not in h_genes])


            #hidden_disease_genes = random.sample(disease_positives,int(len(disease_positives)/opts.k_fold))
            #disease_hidden_list.append(hidden_disease_genes)
            #test_disease_positives = disease_positives.difference(hidden_disease_genes)
            #disease_test_list.append(test_disease_positives)

            #hidden_process_genes = random.sample(biological_process_positives,int(len(biological_process_positives)/opts.k_fold))
            #process_hidden_list.append(hidden_process_genes)
            #test_biological_process_positives = biological_process_positives.difference(hidden_process_genes)
            #process_test_list.append(test_biological_process_positives)

            #print('ROC CURVE: %d hidden disease genes, %d test disease genes, %d hidden biological process genes, and %d test biological process genes' % \
            #   (len(hidden_disease_genes),len(test_disease_positives),len(hidden_process_genes),len(test_biological_process_positives)))
            #print('Hidden Genes:',sorted([x for x in hidden_genes]))
            #print(' disease predictions...')
            #statsfile = opts.outprefix + '_holdout_disease_stats.txt'
            #outfile = opts.outprefix+'_holdout_disease_output.txt'
            #name = 'holdout_disease'
            #ignore,ignore,holdout_d_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_disease_positives,negatives,\
            #    opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,\
            #    name,opts.sinksource_method,write=False)
            #disease_predictions_list.append(holdout_d_predictions)

            #print(' biological process predictions...')
            #statsfile = opts.outprefix + '_holdout_biological_process_stats.txt'
            #outfile = opts.outprefix+'_holdout_biological_process_output.txt'
            #name = 'holdout_biological_process'
            #ignore,ignore,holdout_b_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_biological_process_positives,negatives,\
            #    opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,\
            #    name,opts.sinksource_method,write=False)
            #process_predictions_list.append(holdout_b_predictions)

            ## NEW 7/2 by Anna: normed=True means that both predictions are normalized so the maximum is 1.0.
            #normed=True
            ## write combined results for disease and biological process predictions, including the final score 
            ## which is the product of the two sets of predictions.
            
            #outfile = opts.outprefix+'_holdout_combined_output.txt'
            #fileIO.writeCombinedResults(G,outfile,holdout_d_predictions,holdout_b_predictions,\
            #    disease_positives,biological_process_positives,negatives,blacklist,genemap,opts.layers,normed=normed)
            
            ## NEW 7/2: adjust predictions to ONLY be prime nodes 
            #This checks for primes and removes suffixes before MWU test, which does the same thing - need to decide which is better
            #if opts.layers > 1:
            #    holdout_d_predictions = {x[:-6]:holdout_d_predictions[x] for x in holdout_d_predictions if '_prime' in x}
            #    holdout_b_predictions = {x[:-6]:holdout_b_predictions[x] for x in holdout_b_predictions if '_prime' in x}
            #    test_disease_positives = set([x[:-2] for x in test_disease_positives])
            #    test_biological_process_positives = set([x[:-2] for x in test_biological_process_positives])
            #   hidden_genes = set([x[:-2] for x in hidden_genes])

        print('\nHolding out overlap set and running procedure.')
        hidden_positive_genes = orig_disease_positives.intersection(orig_biological_process_positives)
        hidden_positive_nodes = set(node for gene in hidden_positive_genes for node in multi_node_dict[gene] if node in disease_positives and node in biological_process_positives)

        test_disease_positives = disease_positives.difference(hidden_positive_nodes)
        test_biological_process_positives = biological_process_positives.difference(hidden_positive_nodes)

        if opts.with_negatives:
            hidden_negative_genes = random.sample(orig_negatives,int(len(orig_negatives)/opts.k_fold))
            hidden_negative_nodes = set(node for gene in hidden_negative_genes for node in multi_node_dict[gene] if node in negatives)
            test_negatives = negatives.difference(hidden_negative_nodes)
        else:
            test_negatives=negatives
        

        print('ROC CURVE: %d hidden genes, %d test disease genes, and %d test biological process genes' % \
            (len(hidden_positive_nodes),len(test_disease_positives),len(test_biological_process_positives)))
        print('Hidden Genes:',sorted([genemap[x] for x in hidden_positive_genes]))
        print(' disease predictions...')
        statsfile = opts.outprefix + '_holdout_disease_stats.txt'
        outfile = opts.outprefix+'_holdout_disease_output.txt'
        name='holdout_disease'
        ignore,ignore,holdout_d_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_disease_positives,test_negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,\
            name,opts.sinksource_method,write=True)

        print(' biological process predictions...')
        statsfile = opts.outprefix + '_holdout_biological_process_stats.txt'
        outfile = opts.outprefix+'_holdout_biological_process_output.txt'
        name = 'holdout_biological_process'
        ignore,ignore,holdout_b_predictions = learners.learn(opts.outprefix,outfile,statsfile,genemap,G,test_biological_process_positives,test_negatives,\
            opts.epsilon,opts.timesteps,opts.iterative_update,opts.verbose,opts.force,opts.sinksource_constant,opts.layers,\
            name,opts.sinksource_method,write=True)


        ## NEW 7/2 by Anna: normed=True means that both predictions are normalized so the maximum is 1.0.
        normed=True
        ## write combined results for disease and biological process predictions, including the final score 
        ## which is the product of the two sets of predictions.
        fileIO.writeCombinedResults(G,outfile,holdout_d_predictions,holdout_b_predictions,\
            disease_positives,biological_process_positives,test_negatives,blacklist,genemap,opts.layers,normed=normed)


      ## plot ROC.
        names = ['SZ $f_{\mathcal{D}}$','CM $f_{\mathcal{P}}$','Combined $g$']
        colors =['g','b','r']
        preds = [holdout_d_predictions,holdout_b_predictions,{x:holdout_d_predictions[x]*holdout_b_predictions[x] for x in holdout_d_predictions}]
        test_union_positives=test_disease_positives.union(test_biological_process_positives)
        pos = [test_disease_positives,test_biological_process_positives,test_union_positives]
        plt.clf()
        for i in range(len(names)):
            if opts.with_negatives:
                AUC = Mann_Whitney_U_test2(preds[i], multi_node_dict, pos[i],hidden_positive_nodes, test_negatives, hidden_negative_nodes)
                x,y = getROCvalues2(preds[i],hidden_positive_nodes,pos[i], multi_node_dict, hidden_negative_nodes)

            else:
                AUC = Mann_Whitney_U_test(preds[i], multi_node_dict, pos[i],hidden_positive_nodes, negatives)
                x,y = getROCvalues(preds[i],hidden_positive_nodes,pos[i], multi_node_dict)
            plt.plot(x,y,color=colors[i],label=names[i]+' (AUC=%.2f)' % AUC)
            # plt.xlim(0,len(preds)/(opts.layers+1))
            print(names[i],AUC)
        plt.xlabel('# False Positives')
        plt.ylabel('# True Positives')
        plt.title('ROC (%d layers, $\lambda$=%.2f)' % (opts.layers,opts.sinksource_constant))
        plt.legend(loc='lower right')
        plt.savefig(opts.outprefix+'_ROC.png')
        print('wrote to '+opts.outprefix+'_ROC.png')
        print('Done.')
        return

            ## plot ROC.
            
            #names = ['SZ $f_{\mathcal{D}}$','CM $f_{\mathcal{P}}$','Combined $g$']
            #colors =['g','b','r']
            #preds = [holdout_d_predictions,holdout_b_predictions,{x:holdout_d_predictions[x]*holdout_b_predictions[x] for x in holdout_d_predictions}]
            #test_union_positives=test_disease_positives.union(test_biological_process_positives)
            #pos = [test_disease_positives,test_biological_process_positives,test_union_positives]

        plt.figure()

       
        disease_AUC_sum = 0
        process_AUC_sum = 0

        for k in range(opts.auc_samples): #go through all the runs of k-fold cross validation, collect the sets used and predictions, plot ROC
            names = ['Disease', 'Process']
            colors = ['g','b']
            preds = [disease_predictions_list[k], process_predictions_list[k]]
            hidden_genes = [disease_hidden_list[k], process_hidden_list[k]]
            pos = [disease_test_list[k], process_test_list[k]]
            
            for i in range(len(names)):
                AUC = Mann_Whitney_U_test(preds[i], hidden_genes[i], None, pos[i], multi_node_dict)
                if names[i] == 'Disease':
                    disease_AUC_sum += AUC
                    avg = float(disease_AUC_sum)/float(opts.auc_samples)
                else:
                    process_AUC_sum += AUC
                    avg = float(process_AUC_sum)/float(opts.auc_samples)
                x,y = getROCvalues(preds[i],hidden_genes[i],pos[i])
                print(names[i],AUC)
                plt.xlabel('# False Positives')
                plt.ylabel('# True Positives')
                if k == (opts.auc_samples-1): #Just want one line of each color to have a label - use the last run to get average! 
                    plt.plot(x,y,color=colors[i],label=names[i]+' (Avg AUC=%.2f)' % avg, alpha=0.5)
                else:
                    plt.plot(x,y,color=colors[i], alpha=0.7)
        plt.title('ROC (%d layers, $\lambda$=%.2f)' % (opts.layers,opts.sinksource_constant))
        plt.legend(loc='lower right')
        plt.savefig(opts.outprefix+'_ROC.png')
        print('wrote to '+opts.outprefix+'_ROC.png')

    print('Done.')

    return

'''
def getROCvalues(preds, hidden, pos):
    '''
    Return two lists, which contain coordiantes (x,y) representing
    the number of false positives (x) and the number of true positives (y) 
    as we walk down the list of predictions.
    '''
    
    # x and y are lists of the same length.
    x = [0] # this will be a list of false positives
    y = [0] # this will be a list of true positives


def getROCvalues(preds, hidden, pos, layer_dict):
    '''
    Return two lists, which contain coordiantes (x,y) representing
    the number of false positives (x) and the number of true positives (y) 
    as we walk down the list of predictions.
    '''
    
    # x and y are lists of the same length.
    x = [0] # this will be a list of false positives
    y = [0] # this will be a list of truw positives

    # sort the predictions by their value, largest to smallest
    # this will be a list of nodes
    sorted_preds = sorted(preds.keys(), key=lambda x:preds[x], reverse=True)

    runningx = 0 # current FP counter
    runningy = 0 # current TP counter
    for node in sorted_preds:
        
        ## update running y value (increment if a true positive)
        if node in hidden:
            runningy += 1
        ## update running x value (increment if a false positive)
        elif node not in pos: # ignore positives
            runningx += 1

        ## append (runningx,runningy) as a coordinate
        ## if it's a new coordinate (one of x or y was incremented)
        if runningx != x[-1] or runningy != y[-1]:
            x.append(runningx)
            y.append(runningy)
        if node[-6:] == '_prime': 
            entrez = node[:-6] 
            names = layer_dict[entrez]
        
            ## update running y value (increment if a true positive)
            if bool(names.intersection(hidden)):
                runningy += 1
            ## update running x value (increment if a false positive)
            elif not bool(names.intersection(pos)): # ignore positives
                runningx += 1

            ## append (runningx,runningy) as a coordinate
            ## if it's a new coordinate (one of x or y was incremented)
            if runningx != x[-1] or runningy != y[-1]:
                x.append(runningx)
                y.append(runningy)

    return x,y


def getROCvalues2(preds, hidden_pos, pos, layer_dict, hidden_neg):
    '''
    Return two lists, which contain coordiantes (x,y) representing
    the number of false positives (x) and the number of true positives (y) 
    as we walk down the list of predictions.
    '''
    
    # x and y are lists of the same length.
    x = [0] # this will be a list of false positives
    y = [0] # this will be a list of truw positives

    # sort the predictions by their value, largest to smallest
    # this will be a list of nodes
    sorted_preds = sorted(preds.keys(), key=lambda x:preds[x], reverse=True)

    runningx = 0 # current FP counter
    runningy = 0 # current TP counter
    for node in sorted_preds:
        if node[-6:] == '_prime': 
            entrez = node[:-6] 
            names = layer_dict[entrez]
        
            ## update running y value (increment if a true positive)
            if bool(names.intersection(hidden_pos)):
                runningy += 1
            ## update running x value (increment if a false positive)
            elif bool(names.intersection(hidden_neg)): # ignore positives
                runningx += 1

            ## append (runningx,runningy) as a coordinate
            ## if it's a new coordinate (one of x or y was incremented)
            if runningx != x[-1] or runningy != y[-1]:
                x.append(runningx)
                y.append(runningy)

    return x,y

