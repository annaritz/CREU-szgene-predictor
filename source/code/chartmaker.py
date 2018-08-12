import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 
import sys

#SET USE_SD=False to get Error is IRQ, not standard deviation
USE_SD=True
OUTFILE_DIR = '../outfiles/'
LAMBDAS = ['0','0.01','0.1','1','10','50']
LAYERS = [1,2,3]
EXPERIMENTS = ['SZ','CM_SZ','ASD','CM_ASD'] 
NAMES = {'SZ':'Schizophrenia (SZ)',
        'CM_SZ': 'Cell Motility-SZ',
        'ASD':'Autism (ASD)',
        'CM_ASD':'Cell Motility-ASD'}

def main():
    print('Generating Results from Directory %s' % (OUTFILE_DIR))

    ## read all the files at the beginning.
    data = read_data()

    probplot(data)  ## to make sure that t-test is OK
    
    figure_1(data) ## Layer=1, varying lambda
    
    figure_2(data) ## Negs. vs. NoNegs, Layer=1, varying lambda
    figure_2_full(data) ## Negs vs. NoNegs vs. RandNegs vs. RandNegsPreserveDegree, Layer=1, varying lambda
    
    figure_3_full(data) ## All layers, all lambdas

    ## get best inds
    best_lambda_inds = {} #indices of highest-average lambda values for each experiment.
    for name in EXPERIMENTS:
        best_lambda_inds[name] = {}
        for layer in LAYERS:
            to_check = [mean(data[name][layer][i]) for i in range(len(LAMBDAS))]
            for i in range(len(LAMBDAS)):
                if i == 0 or mean(data[name][layer][best_lambda_inds[name][layer]]) < to_check[i]:
                    best_lambda_inds[name][layer] = i
    
    figure_3(data,best_lambda_inds) ## All layers, best lambda for each layer/experiment.
    

    return

def read_data():
    ## READ DATA
    print('Processing files...')
    data = {}

    ## initialize names
    for disease in EXPERIMENTS:
        data[disease] = {l:[] for l in LAYERS}
        data[disease+'_no_neg'] = {l:[] for l in LAYERS}
        data[disease+'_rand_neg'] = {l:[] for l in LAYERS}
        data[disease+'_rand_neg_deg'] = {l:[] for l in LAYERS}

    for layer in LAYERS:
        for l in LAMBDAS:
            # build list of (file,name1,name2) tuples
            to_process = []
            ## general experiments
            to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource_auc.txt' % (layer,l),'SZ','CM_SZ'))
            to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource_auc.txt' % (layer,l),'ASD','CM_ASD'))
            
            if layer == 1:
                ## no negatives (layer=1 only)
                to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource_no_neg_auc.txt' % (layer,l),'SZ_no_neg','CM_SZ_no_neg'))
                to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource_no_neg_auc.txt' % (layer,l),'ASD_no_neg','CM_ASD_no_neg'))
                ## random negatives (layer 1 only)
                to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource-random_negatives_auc.txt' % (layer,l),'SZ_rand_neg','CM_SZ_rand_neg'))
                to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource-random_negatives_auc.txt' % (layer,l),'ASD_rand_neg','CM_ASD_rand_neg'))
                ## degree-preserving random negatives (layer 1 only)
                to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource-random_negatives_degree_auc.txt' % (layer,l),'SZ_rand_neg_deg','CM_SZ_rand_neg_deg'))
                to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource-random_negatives_degree_auc.txt' % (layer,l),'ASD_rand_neg_deg','CM_ASD_rand_neg_deg'))

            ## process each file in to_proces
            for infile,diseasename,cellmotilityname in to_process:
                #print(diseasename,cellmotilityname,'reading from file',infile)
                disease,process = file_parser(infile)
                data[diseasename][layer].append(disease)
                data[cellmotilityname][layer].append(process)
    return data

#Appends AUC values of each positive set to a list and returns the 3 lists 
def file_parser(auc_file):
    file=open(auc_file,'r')
    x=0
    SZ=[]
    ASD=[]
    CM=[]
    for line in file:
        if x>0:
            line=line.strip('\n').split('\t')
            SZ.append(float(line[0].strip('\''))) #SZ AUCs are always in first column
            CM.append(float(line[1].strip('\''))) #CM AUCs always in third column
        x+=1

    means_IRQs = [np.mean(SZ),np.mean(ASD),np.mean(CM),IQR(SZ),IQR(CM)]

    #return SZ, ASD, CM, means_IRQs
    return SZ, CM


def IQR(dist):
    return [np.percentile(dist, 75) - np.mean(dist), np.mean(dist)-np.percentile(dist, 25)]

def probplot(data):
    fig1, grid = plt.subplots(ncols=len(EXPERIMENTS), nrows=len(LAMBDAS), figsize=(12,14))
    for i in range(len(LAMBDAS)):
        for j in range(len(EXPERIMENTS)):
            name = EXPERIMENTS[j]
            stats.probplot(data[name][1][i],plot=grid[i][j])
            grid[i][j].set_title('QQ Plot: %s $\lambda=%s$' % (name,LAMBDAS[i]))

    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'probplot.png')
    print('Created '+OUTFILE_DIR+'probplot.png')
    return

def figure_1(data):

    with open(OUTFILE_DIR+'sinksource+_constant_statistical_tests.txt', 'w') as out:
        sig = {} # asterisks
        for name in ['SZ','CM_SZ','ASD','CM_ASD']:
            sig[name] = []

            out.write('--- %s ---\n' % (name))

            ## Write Means of AUC Distribution
            out.write('AUC Means:\n')
            for i in range(len(LAMBDAS)):
                avg = mean(data[name][1][i])
                out.write('\tLambda %s: %.4f\n' % (LAMBDAS[i],avg))

            ## 1-way ANOVA to assess means
            f_value, p_value = stats.f_oneway(data[name][1][0], data[name][1][1], data[name][1][2], data[name][1][3], data[name][1][4], data[name][1][5])
            out.write('1-way ANOVA\n\tf-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (f_value,p_value,p_value))
            
            ## Welch's t-test to assess 0 vs. a weighted lambda.
            for i in range(1,len(LAMBDAS)):
                t_value,p_value = stats.ttest_ind(data[name][1][0], data[name][1][i], equal_var=False)
                if t_value < 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                    sig[name].append(i+1)
                else: 
                    res = ''
                out.write('Welch\'s t-test 0 vs. %s: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[i],res,t_value,p_value,p_value))
            out.write('\n')


    #fig1, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, nrows=1, sharey=True, figsize=(7,4))
    fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8,6))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        bp = ax.boxplot(data[name][1], notch=True, widths=.7, sym='', \
                patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        if len(sig[name]) > 0:
            ax.plot(sig[name],[0.88]*len(sig[name]),'*k')
        # set title, axis limits, and labels
        ax.set_title(NAMES[name])
        ax.set_xticks([j+1 for j in range(len(data[name][1]))])
        ax.set_xticklabels(['0','0.01','0.1','1','10','50'])
        
        ax.set_ylim(0.5,0.9)
        
        ax.set_ylabel('AUC')
        ax.set_xlabel('$\lambda$')
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_lambda.png')

    print('Created '+OUTFILE_DIR+'compare_lambda.png')
    
    return


def figure_2(data):
    
    with open(OUTFILE_DIR+'Neg_vs_NoNeg_statistical_tests.txt', 'w') as out:
        
        for disease in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s ---\n' % (disease))
            for j in range(len(LAMBDAS)): 
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(data[disease][1][j], data[disease+'_no_neg'][1][j], equal_var=False)
                if t_value > 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                else: 
                    res = ''
                out.write('Welch\'s t-test Lambda=%s (negs vs. no-negs): %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[j],res,t_value,p_value,p_value))
            out.write('\n')

        out.write('Two-way ANOVA\n')
        organized_data_1=[data[exp][1] for exp in EXPERIMENTS]
        organized_data_2=[data[exp+'_no_neg'][1] for exp in EXPERIMENTS]
        for i in range(len(organized_data_1)):
            if i == 0:
                out.write('Schizophrenia\t')
            if i == 1:
                out.write('Cell Motility-SZ\t')
            if i == 2:
                out.write('ASD\t')
            if i == 3:
                out.write('Cell Motility-ASD\t')
        
            out.write('\tssq\tdf\tF\tPR(>F)\n')

            two_way_ANOVA_results=two_way_ANOVA(organized_data_1[i], organized_data_2[i])
            for j in range(len(two_way_ANOVA_results)):
                if j == 0:
                    name = 'Including Negatives\t'
                elif j == 1:
                    name = 'sinksource constant\t'
                else:
                    name = 'Interaction\t'
                for k in range(len(two_way_ANOVA_results[j])):
                    name=name+str(two_way_ANOVA_results[j][k])+'\t'
                name=name+'\n'
                out.write(name)
            out.write('\n')

    fig2, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8,6))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        bp1 = ax.boxplot(data[name][1], notch=True, positions=[2,7,12,17,22,27], widths=1.5, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        bp2 = ax.boxplot(data[name+'_no_neg'][1], notch=True, positions=[4,9,14,19,24,29], widths=1.5, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        
        ax.set_title(NAMES[name])
        ax.set_xticks([3,8,13,18,23,28])
        ax.set_xticklabels(['0','0.01','0.1','1','10','50'])
        ax.set_xlim(0,31)
        ax.set_ylim(0.5,0.9)
        ax.set_ylabel('AUC')
        ax.set_xlabel('$\lambda$')
        
        if i == 0:
            ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize='x-small')
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_neg_noneg.png')

    print('Created '+OUTFILE_DIR+'compare_neg_noneg.png')

    return

def figure_2_full(data):

    #compare negatives vs. no negatives vs. random negatives.
    
    with open(OUTFILE_DIR+'RandNeg_statistical_tests.txt', 'w') as out:
        
        for disease in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s RANDOM NEGATIVES ---\n' % (disease))
            for j in range(len(LAMBDAS)): 
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(data[disease][1][j], data[disease+'_rand_neg'][1][j], equal_var=False)
                if t_value > 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                else: 
                    res = ''
                out.write('Welch\'s t-test Lambda=%s (negs vs. rand-negs): %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[j],res,t_value,p_value,p_value))
            out.write('\n')

        out.write('Two-way ANOVA RANDOM NEGATIVES\n')
        organized_data_1=[data[exp][1] for exp in EXPERIMENTS]
        organized_data_2=[data[exp+'_no_neg'][1] for exp in EXPERIMENTS]
        for i in range(len(organized_data_1)):
            if i == 0:
                out.write('Schizophrenia\t')
            if i == 1:
                out.write('Cell Motility-SZ\t')
            if i == 2:
                out.write('ASD\t')
            if i == 3:
                out.write('Cell Motility-ASD\t')
        
            out.write('\tssq\tdf\tF\tPR(>F)\n')

            two_way_ANOVA_results=two_way_ANOVA(organized_data_1[i], organized_data_2[i])
            for j in range(len(two_way_ANOVA_results)):
                if j == 0:
                    name = 'Including Random Negatives\t'
                elif j == 1:
                    name = 'sinksource constant\t'
                else:
                    name = 'Interaction\t'
                for k in range(len(two_way_ANOVA_results[j])):
                    name=name+str(two_way_ANOVA_results[j][k])+'\t'
                name=name+'\n'
                out.write(name)
            out.write('\n')

        for disease in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s RANDOM NEGATIVES (PRESERVING DEGREE) ---\n' % (disease))
            for j in range(len(LAMBDAS)): 
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(data[disease][1][j], data[disease+'_rand_neg_deg'][1][j], equal_var=False)
                if t_value > 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                else: 
                    res = ''
                out.write('Welch\'s t-test Lambda=%s (negs vs. degree-preserving rand_negs): %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[j],res,t_value,p_value,p_value))
            out.write('\n')

        out.write('Two-way ANOVA RANDOM NEGATIVES (PRESERVING DEGREE)\n')
        organized_data_1=[data[exp][1] for exp in EXPERIMENTS]
        organized_data_2=[data[exp+'_no_neg'][1] for exp in EXPERIMENTS]
        for i in range(len(organized_data_1)):
            if i == 0:
                out.write('Schizophrenia\t')
            if i == 1:
                out.write('Cell Motility-SZ\t')
            if i == 2:
                out.write('ASD\t')
            if i == 3:
                out.write('Cell Motility-ASD\t')
        
            out.write('\tssq\tdf\tF\tPR(>F)\n')

            two_way_ANOVA_results=two_way_ANOVA(organized_data_1[i], organized_data_2[i])
            for j in range(len(two_way_ANOVA_results)):
                if j == 0:
                    name = 'Including Degree-Preserving Random Negatives\t'
                elif j == 1:
                    name = 'sinksource constant\t'
                else:
                    name = 'Interaction\t'
                for k in range(len(two_way_ANOVA_results[j])):
                    name=name+str(two_way_ANOVA_results[j][k])+'\t'
                name=name+'\n'
                out.write(name)
            out.write('\n')

    fig2, ((ax1),(ax2),(ax3),(ax4)) = plt.subplots(ncols=1, nrows=4, figsize=(8,12))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        bp1 = ax.boxplot(data[name][1], notch=True, positions=[1,6,11,16,21,26], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        bp2 = ax.boxplot(data[name+'_no_neg'][1], notch=True, positions=[2,7,12,17,22,27], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        bp3 = ax.boxplot(data[name+'_rand_neg'][1], notch=True, positions=[3,8,13,18,23,28], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#eeb5ff',color='#A152B8'))
        bp4 = ax.boxplot(data[name+'_rand_neg_deg'][1], notch=True, positions=[4,9,14,19,24,29], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#FFE888',color='#E1C23E'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        for x in [5,10,15,20,25]:
            ax.plot([x,x],[0.5,0.9],color='lightgrey',alpha=0.5)
        ax.set_title(NAMES[name])
        ax.set_xticks([2.5,7.5,12.5,17.5,22.5,27.5])
        ax.set_xticklabels(['0','0.01','0.1','1','10','50'])
        ax.set_xlim(0,30)
        ax.set_ylim(0.5,0.9)
        ax.set_ylabel('AUC')
        ax.set_xlabel('$\lambda$')
        
        ax.legend([bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0], bp4['boxes'][0]], ['With Negatives', 'Without Negatives', 'Rand Negatives','Deg. Preserving Rand Negatives'], loc='best', fontsize='x-small')
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_rand_negs.png')

    print('Created '+OUTFILE_DIR+'compare_rand_negs.png')

    return
    


def figure_3(data,best_lambda_inds):
    best_lambda_data = {}
    best_lambdas = {}
    for name in EXPERIMENTS:
        best_lambdas[name] = {}
        best_lambda_data[name] = {}
        for layer in LAYERS:
            i = best_lambda_inds[name][layer]
            best_lambdas[name][layer] = LAMBDAS[i]
            best_lambda_data[name][layer] = data[name][layer][i]

    with open(OUTFILE_DIR+'Best_Lambdas_statistical_tests.txt', 'w') as out:

        for name in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s ---\n' % (name))

            ## Write Means of AUC Distribution
            out.write('  AUC Mean (best lambda):\n')
            for layer in LAYERS:
                out.write('\tLayer %d Lambda %s: %.4f\n' % (layer,best_lambdas[name][layer],mean(best_lambda_data[name][layer])))
            out.write('\n')

            ## Welch's t-test to assess pairs of layers.
            # Layer1 vs. Layer2
            t_value,p_value = stats.ttest_ind(best_lambda_data[name][1], best_lambda_data[name][2], equal_var=False)
            res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
            out.write('  Welch\'s t-test Layer1 vs. Layer2: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))

            # Layer1 vs. Layer3
            t_value,p_value = stats.ttest_ind(best_lambda_data[name][1], best_lambda_data[name][3], equal_var=False)
            res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
            out.write('  Welch\'s t-test Layer1 vs. Layer3: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))

            # Layer2 vs. Layer3
            t_value,p_value = stats.ttest_ind(best_lambda_data[name][2], best_lambda_data[name][3], equal_var=False)
            res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
            out.write('  Welch\'s t-test Layer2 vs. Layer3: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))
            out.write('\n')

    '''
    ## TODO check QQ plot for t-test vs. MWU test.
    Mann_test_lists = [SZ_sets, ASD_sets, CM_sets]
    with open(OUTFILE_DIR+'MWU_PVALUES_123layers.txt', 'w') as out:
        for i in range(len(Mann_test_lists)):
            if i == 0:
                out.write('Schizophrenia\t')
            elif i == 1:
                out.write('Autism\t')
            else:
                out.write('Cell Motility\t')
            out.write('\n')
            U1, p1 = stats.mannwhitneyu(Mann_test_lists[i][0], Mann_test_lists[i][1], alternative='two-sided') 
            out.write('2 layer vs. 1 layer\t')
            out.write(str(p1) + '\n')
            U3, p3 = stats.mannwhitneyu(Mann_test_lists[i][0], Mann_test_lists[i][2], alternative='two-sided')
            out.write('2 layer vs. 3 layer\t')
            out.write(str(p3) + '\n') 
            out.write('\n')
    '''

    fig3, axes = plt.subplots(ncols=4, nrows=1, sharey=True, figsize=(7,3))
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        bp1 = ax.boxplot(best_lambda_data[name][1], notch=True, positions=[2], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        bp2 = ax.boxplot(best_lambda_data[name][2], notch=True, positions=[4], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        bp3 = ax.boxplot(best_lambda_data[name][3], notch=True, positions=[6], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#eeb5ff',color='#A152B8'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        ax.set_ylim(0.5,0.9)
        ax.set_xlim(0,8)
        ax.set_xticks([2,4,6])
        ax.set_xticklabels(['1','2','3'])
        if i==0:
            ax.set_ylabel('AUC')
        ax.set_xlabel('Layers')
        ax.set_title(NAMES[name])
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'comparing_layers_best_lambdas.png')
    print('Created '+OUTFILE_DIR+'comparing_best_lambdas.png')

    return


def figure_3_full(data):

    with open(OUTFILE_DIR+'Layers_statistical_tests.txt', 'w') as out:

        for name in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s ---\n' % (name))

            for i in range(len(LAMBDAS)):
                out.write(' +++ LAMBDA %s +++\n' % (LAMBDAS[i]))

                ## Write Means of AUC Distribution
                out.write('  AUC Means:\n')
                for layer in LAYERS:
                    out.write('\tLayer %d: %.4f\n' % (layer,mean(data[name][layer][i])))
                out.write('\n')

                ## Welch's t-test to assess pairs of layers.
                # Layer1 vs. Layer2
                t_value,p_value = stats.ttest_ind(data[name][1][i], data[name][2][i], equal_var=False)
                res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
                out.write('  Welch\'s t-test Layer1 vs. Layer2: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))

                # Layer1 vs. Layer3
                t_value,p_value = stats.ttest_ind(data[name][1][i], data[name][3][i], equal_var=False)
                res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
                out.write('  Welch\'s t-test Layer1 vs. Layer3: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))

                # Layer2 vs. Layer3
                t_value,p_value = stats.ttest_ind(data[name][2][i], data[name][3][i], equal_var=False)
                res = 'SIGNIFICANT (two-tailed, 0.01)' if p_value < 0.01 else ''
                out.write('  Welch\'s t-test Layer2 vs. Layer3: %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (res,t_value,p_value,p_value))
            out.write('\n')


        # Unclear to me why to do a 2-way ANOVA here...
        '''
        out.write('Two-way ANOVA\n')
        for i in range(len(one_list)):
            if i == 0:
                out.write('Schizophrenia\t')
            if i == 1:
                out.write('Cell Motility-SZ\t')
            if i == 2:
                out.write('ASD\t')
            if i == 3:
                out.write('Cell Motility-ASD\t')
        
            out.write('\tssq\tdf\tF\tPR(>F)\n')

            two_way_ANOVA_results=two_way_ANOVA(one_list[i], two_list[i])
            for j in range(len(two_way_ANOVA_results)):
                if j == 0:
                    name = 'Layers\t'
                elif j == 1:
                    name = 'sinksource constant\t'
                else:
                    name = 'Interaction\t'
                for k in range(len(two_way_ANOVA_results[j])):
                    name=name+str(two_way_ANOVA_results[j][k])+'\t'
                name=name+'\n'
                out.write(name)
            out.write('\n')
        '''

    fig3, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8,6))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        bp1 = ax.boxplot(data[name][1], notch=True, positions=[2,9,16,23,30,37], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        bp2 = ax.boxplot(data[name][2], notch=True, positions=[4,11,18,25,32,39], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        bp3 = ax.boxplot(data[name][3], notch=True, positions=[6,13,20,27,34,41], \
            widths=1.5, sym='', patch_artist=True, boxprops=dict(facecolor='#eeb5ff',color='#A152B8'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

        ax.set_xticks([4,11,18,25,32,39])
        ax.set_xticklabels(['0','0.01','0.1','1','10','50'])
        ax.set_ylim(0.5,0.9)
        ax.set_xlim(0,43)
        ax.set_ylabel('AUC')
        ax.set_xlabel('$\lambda$')
        ax.set_title(NAMES[name])
        if i==0:
            ax.legend([bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0]], ['1 Layer', '2 Layers', '3 Layers'], loc='best', fontsize='x-small')
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'comparing_layers.png')

    print('Created ../outfiles/comparing_layers.png')

    return
    
    

def two_way_ANOVA(List_of_lists_A, List_of_lists_B):
    number_of_settings=2
    number_of_constants=6

    N=50*number_of_settings*number_of_constants
    df_settings=number_of_settings-1
    df_constants=number_of_constants-1
    df_settingsxconstants = df_settings*df_constants
    df_w = (50-1)* number_of_settings*number_of_constants
    grandList=[]

    A_setting_List=[]
    B_setting_List=[]

    zero_list=[]
    point_zero_one_list=[]
    point_one_list=[]
    one_list=[]
    ten_list=[]
    fifty_list=[]

    constants_list_of_lists=[[],[],[],[],[],[]]



    for i in range(len(List_of_lists_A)):
        for j in range(len(List_of_lists_A[i])):
            grandList.append(List_of_lists_A[i][j])
            A_setting_List.append(List_of_lists_A[i][j])
            constants_list_of_lists[i].append(List_of_lists_A[i][j])

    for i in range(len(List_of_lists_B)):
        for j in range(len(List_of_lists_B[i])):
            grandList.append(List_of_lists_B[i][j])
            B_setting_List.append(List_of_lists_B[i][j])
            constants_list_of_lists[i].append(List_of_lists_B[i][j])


    ssq_settings=50*number_of_constants*((mean(A_setting_List)-mean(grandList))**2 + (mean(B_setting_List)-mean(grandList))**2) 
    MS_settings=ssq_settings/df_settings

    ssq_constants=50*number_of_settings* sum([(mean(constant_values)-mean(grandList))**2 for constant_values in constants_list_of_lists])
    MS_constants=ssq_constants/df_constants

    ssq_settingsxconstants=0

    ssq_within=0

    for i in range(len(List_of_lists_A)):
        within_mean=mean(List_of_lists_A[i])
        for j in range(len(List_of_lists_A[i])):
            ssq_within = ssq_within + (List_of_lists_A[i][j] - within_mean)**2

            ssq_settingsxconstants = ssq_settingsxconstants + (within_mean - mean(A_setting_List) - mean(constants_list_of_lists[i]) + mean(grandList))**2

    for i in range(len(List_of_lists_B)):
        within_mean=mean(List_of_lists_B[i])
        for j in range(len(List_of_lists_B[i])):
            ssq_within = ssq_within + (List_of_lists_B[i][j] - within_mean)**2

            ssq_settingsxconstants = ssq_settingsxconstants + (within_mean - mean(B_setting_List) - mean(constants_list_of_lists[i]) + mean(grandList))**2


    ssq_settingsxconstants = 50*ssq_settingsxconstants

    # ssq_settingsxconstants = ssq_total-ssq_settings-ssq_constants-ssq_within

    MS_settingsxconstants = ssq_settingsxconstants/df_settingsxconstants

    MS_within = ssq_within/df_w


    f_settings = MS_settings/MS_within
    f_constants = MS_constants/MS_within
    f_settingsxconstants = MS_settingsxconstants/MS_within

    p_a = stats.f.sf(f_settings, df_settings, df_w)
    p_b = stats.f.sf(f_constants, df_constants, df_w)
    p_axb = stats.f.sf(f_settingsxconstants, df_settingsxconstants, df_w)

    results = [[ssq_settings, df_settings, f_settings, p_a],[ssq_constants, df_constants, f_constants, p_b],[ssq_settingsxconstants, df_settingsxconstants, f_settingsxconstants, p_axb]]

    return results

def two_way_ANOVA(List_of_lists_A, List_of_lists_B):
    number_of_settings=2
    number_of_constants=6
    replications=50

    N=replications*number_of_settings*number_of_constants
    df_settings=number_of_settings-1
    df_constants=number_of_constants-1
    df_settingsxconstants = df_settings*df_constants
    df_w = (replications-1)* number_of_settings*number_of_constants
    grandList=[]

    A_setting_List=[]
    B_setting_List=[]

    zero_list=[]
    point_zero_one_list=[]
    point_one_list=[]
    one_list=[]
    ten_list=[]
    fifty_list=[]

    constants_list_of_lists=[[],[],[],[],[],[]]



    for i in range(len(List_of_lists_A)):
        for j in range(len(List_of_lists_A[i])):
            grandList.append(List_of_lists_A[i][j])
            A_setting_List.append(List_of_lists_A[i][j])
            constants_list_of_lists[i].append(List_of_lists_A[i][j])

    for i in range(len(List_of_lists_B)):
        for j in range(len(List_of_lists_B[i])):
            grandList.append(List_of_lists_B[i][j])
            B_setting_List.append(List_of_lists_B[i][j])
            constants_list_of_lists[i].append(List_of_lists_B[i][j])


    ssq_settings=replications*number_of_constants*((mean(A_setting_List)-mean(grandList))**2 + (mean(B_setting_List)-mean(grandList))**2) 
    MS_settings=ssq_settings/df_settings

    ssq_constants=replications*number_of_settings* sum([(mean(constant_values)-mean(grandList))**2 for constant_values in constants_list_of_lists])
    MS_constants=ssq_constants/df_constants



    ssq_within=0

    for i in range(len(List_of_lists_A)):
        within_mean=mean(List_of_lists_A[i])
        for j in range(len(List_of_lists_A[i])):
            ssq_within = ssq_within + (List_of_lists_A[i][j] - within_mean)**2


    for i in range(len(List_of_lists_B)):
        within_mean=mean(List_of_lists_B[i])
        for j in range(len(List_of_lists_B[i])):
            ssq_within = ssq_within + (List_of_lists_B[i][j] - within_mean)**2


    ssq_total=sum([(value-mean(grandList))**2 for value in grandList])

    ssq_settingsxconstants = ssq_total-ssq_settings-ssq_constants-ssq_within

    MS_settingsxconstants = ssq_settingsxconstants/df_settingsxconstants

    MS_within = ssq_within/df_w


    f_settings = MS_settings/MS_within
    f_constants = MS_constants/MS_within
    f_settingsxconstants = MS_settingsxconstants/MS_within

    p_a = stats.f.sf(f_settings, df_settings, df_w)
    p_b = stats.f.sf(f_constants, df_constants, df_w)
    p_axb = stats.f.sf(f_settingsxconstants, df_settingsxconstants, df_w)

    results = [[ssq_settings, df_settings, f_settings, p_a],[ssq_constants, df_constants, f_constants, p_b],[ssq_settingsxconstants, df_settingsxconstants, f_settingsxconstants, p_axb]]

    return results

def two_way_ANOVA_test(List_of_lists_A, List_of_lists_B):
    number_of_settings=2
    number_of_constants=2
    replications=5

    N=replications*number_of_settings*number_of_constants
    df_settings=number_of_settings-1
    df_constants=number_of_constants-1
    df_settingsxconstants = df_settings*df_constants
    df_w = (replications-1)* number_of_settings*number_of_constants
    grandList=[]

    A_setting_List=[]
    B_setting_List=[]

    zero_list=[]
    point_zero_one_list=[]
    point_one_list=[]
    one_list=[]
    ten_list=[]
    fifty_list=[]

    constants_list_of_lists=[[],[]]



    for i in range(len(List_of_lists_A)):
        for j in range(len(List_of_lists_A[i])):
            grandList.append(List_of_lists_A[i][j])
            A_setting_List.append(List_of_lists_A[i][j])
            constants_list_of_lists[i].append(List_of_lists_A[i][j])

    for i in range(len(List_of_lists_B)):
        for j in range(len(List_of_lists_B[i])):
            grandList.append(List_of_lists_B[i][j])
            B_setting_List.append(List_of_lists_B[i][j])
            constants_list_of_lists[i].append(List_of_lists_B[i][j])


    ssq_settings=replications*number_of_constants*((mean(A_setting_List)-mean(grandList))**2 + (mean(B_setting_List)-mean(grandList))**2) 
    MS_settings=ssq_settings/df_settings

    ssq_constants=replications*number_of_settings* sum([(mean(constant_values)-mean(grandList))**2 for constant_values in constants_list_of_lists])
    MS_constants=ssq_constants/df_constants

    ssq_total=sum([(value-mean(grandList))**2 for value in grandList])
    ssq_within=0

    for i in range(len(List_of_lists_A)):
        within_mean=mean(List_of_lists_A[i])
        for j in range(len(List_of_lists_A[i])):
            ssq_within = ssq_within + (List_of_lists_A[i][j] - within_mean)**2


    for i in range(len(List_of_lists_B)):
        within_mean=mean(List_of_lists_B[i])
        for j in range(len(List_of_lists_B[i])):
            ssq_within = ssq_within + (List_of_lists_B[i][j] - within_mean)**2



    ssq_total=sum([(value-mean(grandList))**2 for value in grandList])
    ssq_settingsxconstants = ssq_total-ssq_settings-ssq_constants-ssq_within

    MS_settingsxconstants = ssq_settingsxconstants/df_settingsxconstants

    MS_within = ssq_within/df_w


    f_settings = MS_settings/MS_within
    f_constants = MS_constants/MS_within
    f_settingsxconstants = MS_settingsxconstants/MS_within

    p_a = stats.f.sf(f_settings, df_settings, df_w)
    p_b = stats.f.sf(f_constants, df_constants, df_w)
    p_axb = stats.f.sf(f_settingsxconstants, df_settingsxconstants, df_w)

    results = [[ssq_settings, df_settings, f_settings, p_a],[ssq_constants, df_constants, f_constants, p_b],[ssq_settingsxconstants, df_settingsxconstants, f_settingsxconstants, p_axb]]

    return results


def mean(bunch_of_numbers):
    return np.mean(bunch_of_numbers)


main()

# manure=[[13.7, 16.8, 14.9, 17.6, 16.5],[17.4, 13.5, 15.1, 15.4, 13.2]]
# no_manure=[[16.0, 16.1, 13.0, 16.7, 13.2],[13.4, 11.6, 14.7, 9.7, 11.9]]


# x= two_way_ANOVA_test(manure, no_manure)

# for i in x:
#     print(i)



























