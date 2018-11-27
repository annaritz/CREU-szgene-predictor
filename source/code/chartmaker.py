import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 
import sys
import os 
import itertools

#SET USE_SD=False to get Error is IRQ, not standard deviation
USE_SD=True
OUTFILE_DIR = '../outfiles/'
LAMBDAS = ['0','0.01','0.1','1','10','50']
LAYERS = [1,2,3]

EXPERIMENTS = ['SZ','CM_SZ','ASD','CM_ASD'] 
NAMES = {'SZ':'Schizophrenia (SZ)',
        'CM_SZ': 'Cell Motility-SZ',
        'ASD':'Autism (ASD)',
        'CM_ASD':'Cell Motility-ASD',
        'CM': 'Cell Motility'}

EXP_COLORS = {'SZ':'#68C8C3',
        'CM_SZ': '#FFC300',
        'ASD':'#A4B97D',
        'CM_ASD':'#C70039'}
TITLE_SIZE=14
LABEL_SIZE=12
TICK_SIZE=10

def main():
    print('Generating Results from Directory %s' % (OUTFILE_DIR))

    ## read all the files at the beginning.
    data = read_data()

    

    '''
    OLD FIGURES (before removing the multilayer stuff...)
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
    '''




    deg_dist_fig('0',ymax=500)
    deg_dist_fig('0.01',ymax=500)
    deg_dist_fig('0.1',ymax=500)
    deg_dist_fig('1',ymax=500)
    deg_dist_fig('10',ymax=500)
    deg_dist_fig('50',ymax=500)
    
    figure_2_full(data) ## Negs vs. NoNegs vs. RandNegs vs. RandNegsPreserveDegree vs. KrishnanNegs, Layer=1, varying lambda
    figure_2(data)

    probplot(data)  ## to make sure that t-test is OK
    figure_3_full(data) ## All layers, all lambdas
    sinksource_fig(data)
    roc_fig(data)


    full_ss_scatter_fig(None,'0')
    full_ss_scatter_fig(None,'1')

    ss_scatter_fig(500,'0')    
    ss_scatter_fig(500,'0',2000)
    ss_scatter_fig(500,'0.01',2000)
    ss_scatter_fig(500,'0.1',2000)
    ss_scatter_fig(500,'1',2000)
    ss_scatter_fig(500,'10',2000)
    ss_scatter_fig(500,'50',2000)
    

    return


def read_data():
    ## READ DATA
    print('Processing files...')
    data = {}

    ## initialize names
    for disease in EXPERIMENTS:
        data[disease] = {l:[] for l in LAYERS}
        data[disease+'_no_neg'] = {l:[] for l in LAYERS}
        data[disease+'_old_neg'] = {l:[] for l in LAYERS}
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

                ## old (krishnan et al.) negatives (layer 1 only)
                to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource-old_negatives_auc.txt' % (layer,l),'SZ_old_neg','CM_SZ_old_neg'))
                to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource-old_negatives_auc.txt' % (layer,l),'ASD_old_neg','CM_ASD_old_neg'))

            ## process each file in to_proces
            for infile,diseasename,cellmotilityname in to_process:
                #print(diseasename,cellmotilityname,'reading from file',infile)
                disease,process = file_parser(infile)
                data[diseasename][layer].append(disease)
                data[cellmotilityname][layer].append(process)
    return data

def read_rocs():
    ## READ DATA
    print('Processing ROC files...')
    data = {}

    ## initialize names
    for disease in EXPERIMENTS:
        data[disease] = {l:[] for l in LAYERS}

    for layer in LAYERS:
        for l in ['0']:
            # build list of (file,name1,name2) tuples
            to_process = []
            ## general experiments
            to_process.append((OUTFILE_DIR+'SZ_%d-layer_%s-sinksource_roccurves.txt' % (layer,l),'SZ','CM_SZ'))
            to_process.append((OUTFILE_DIR+'ASD_%d-layer_%s-sinksource_roccurves.txt' % (layer,l),'ASD','CM_ASD'))
            
            ## process each file in to_proces
            for infile,diseasename,cellmotilityname in to_process:
                if not os.path.isfile(infile):
                    print('ERROR: file %s does not exist. Skipping.' % (infile))
                    continue

                with open(infile) as fin:
                    for line in fin:
                        if line[0] == '#':  # skip header
                            continue
                        row = line.strip().split()
                        d_rec = [int(a) for a in row[1].split(',')]
                        d_prec = [int(a) for a in row[2].split(',')]
                        b_rec = [int(a) for a in row[3].split(',')]
                        b_prec = [int(a) for a in row[4].split(',')]

                        ## normalize to be between 0 and 1.
                        d_rec = [a/d_rec[-1] for a in d_rec]
                        d_prec = [a/d_prec[-1] for a in d_prec]
                        b_rec = [a/b_rec[-1] for a in b_rec]
                        b_prec = [a/b_prec[-1] for a in b_prec]

                        data[diseasename][layer].append([d_rec,d_prec])
                        data[cellmotilityname][layer].append([b_rec,b_prec])
                    print(infile,':',len(data[diseasename][layer]),len(data[cellmotilityname][layer]))
    return data

#Appends AUC values of each positive set to a list and returns the 3 lists 
def file_parser(auc_file):
    
    if not os.path.isfile(auc_file):
        print('ERROR: File %s does not exist.' % (auc_file))
        return [],[]

    disease=[]
    CM=[]
    with open(auc_file,'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            line=line.strip().split('\t')
            disease.append(float(line[0])) #SZ AUCs are always in first column
            CM.append(float(line[1])) #CM AUCs always in third column

    #means_IRQs = [np.mean(SZ),np.mean(ASD),np.mean(CM),IQR(SZ),IQR(CM)]

    #return SZ, ASD, CM, means_IRQs
    return disease, CM


def IQR(dist):
    return [np.percentile(dist, 75) - np.mean(dist), np.mean(dist)-np.percentile(dist, 25)]

####### NEW FIGURES

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

    plt.savefig(OUTFILE_DIR+'probplot.pdf')
    print('Created '+OUTFILE_DIR+'probplot.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'probplot.pdf '+OUTFILE_DIR+'probplot.pdf')
    return


def deg_dist_fig(l,ymax=None):

    pos = {ex:[] for ex in ['SZ','ASD','CM']}
    data = {ex:[] for ex in ['SZ','ASD','CM']}
    infile = OUTFILE_DIR+'SZ_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        process = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            if row[2] == 'Unlabeled':
                disease.append([float(row[3]),int(row[8])])
            elif row[2] == 'Positive':
                pos['SZ'].append(int(row[8]))
            if row[4] == 'Unlabeled':
                process.append([float(row[5]),int(row[8])])
            elif row[4] == 'Positive':
                pos['CM'].append(int(row[8]))
        disease.sort(reverse=True,key=lambda x: x[0])
        process.sort(reverse=True,key=lambda x: x[0])
        data['SZ'] = [disease[i][1] for i in range(len(pos['SZ']))]
        data['CM'] = [process[i][1] for i in range(len(pos['CM']))]

    infile = OUTFILE_DIR+'ASD_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            if row[2] == 'Unlabeled':
                disease.append([float(row[3]),int(row[8])])
            elif row[2] == 'Positive':
                pos['ASD'].append(int(row[8]))
        disease.sort(reverse=True,key=lambda x: x[0])
        data['ASD'] = [disease[i][1] for i in range(len(pos['ASD']))]

    fig2, ((ax1,ax2,ax3)) = plt.subplots(ncols=3, nrows=1, figsize=(10,3))
    axes = [ax1,ax2,ax3]
    for i in range(len(['SZ','ASD','CM'])):
        name = ['SZ','ASD','CM'][i]
        print(name,len(data[name]),len(pos[name]),sum(data[name]),sum(pos[name]))
        ax = axes[i]
        ax.hist(pos[name],bins=range(0,3500,100), color='k',edgecolor='k', linewidth=1.2,label='Positives (avg=%.2f)' % (sum(pos[name])/len(pos[name])))#, alpha=0.2)#, color=[.7, .7, .7], label=None)
        ax.hist(data[name],bins=range(0,3500,100), color='#A6DBF7', alpha=0.8, rwidth=0.7,edgecolor='k', linewidth=1.2, label='Pseudo SS (avg=%.2f)' % (sum(data[name])/len(data[name])))
        ax.set_title(NAMES[name] +'\n$\lambda=%s$' % (l),fontsize=TITLE_SIZE)
        #ax.set_xlim(.5,6.5)
        ax.set_ylabel('Count',fontsize=LABEL_SIZE)
        ax.set_xlabel('Degree',fontsize=LABEL_SIZE)
        
        if ymax:
            ax.set_ylim(0,ymax)
        ax.legend(loc='best',fontsize=TICK_SIZE)
        #if i == 0:
        #    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    prefix = 'deg_fig_%s' % (l.replace('.','pt'))
    plt.savefig(OUTFILE_DIR+prefix+'.png')
    print('Created '+OUTFILE_DIR+prefix+'.png')

    plt.savefig(OUTFILE_DIR+prefix+'.pdf')
    print('Created '+OUTFILE_DIR+prefix+'.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+prefix+'.pdf '+OUTFILE_DIR+prefix+'.pdf')
    return


def ss_scatter_fig(num,l,ymax=None):
    data = {ex:[] for ex in EXPERIMENTS}
    avg = {ex:[] for ex in EXPERIMENTS}
    alpha=15

    infile = OUTFILE_DIR+'SZ_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        process = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            if row[2] == 'Unlabeled':
                disease.append([float(row[3]),int(row[8])])
            if row[4] == 'Unlabeled':
                process.append([float(row[5]),int(row[8])])
        disease.sort(reverse=True,key=lambda x: x[0])
        process.sort(reverse=True,key=lambda x: x[0])
        print(infile,disease[1:10])
        data['SZ'] = [disease[i][1] for i in range(num)]
        data['CM_SZ'] = [process[i][1] for i in range(num)]
        avg['SZ'] = [0]*(len(data['SZ'])-alpha)
        avg['CM_SZ'] = [0]*(len(data['CM_SZ'])-alpha)
        for i in range(len(data['SZ'])-alpha):
            avg['SZ'][i] = sum(data['SZ'][i:i+alpha])/alpha
        for i in range(len(data['SZ'])-alpha):
            avg['CM_SZ'][i] = sum(data['CM_SZ'][i:i+alpha])/alpha

    infile = OUTFILE_DIR+'ASD_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        process = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            if row[2] == 'Unlabeled':
                disease.append([float(row[3]),int(row[8])])
            if row[4] == 'Unlabeled':
                process.append([float(row[5]),int(row[8])])
        disease.sort(reverse=True,key=lambda x: x[0])
        process.sort(reverse=True,key=lambda x: x[0])
        
        data['ASD'] = [disease[i][1] for i in range(num)]
        data['CM_ASD'] = [process[i][1] for i in range(num)]
        avg['ASD'] = [0]*(len(data['ASD'])-alpha)
        avg['CM_ASD'] = [0]*(len(data['CM_ASD'])-alpha)
        for i in range(len(data['ASD'])-alpha):
            avg['ASD'][i] = sum(data['ASD'][i:i+alpha])/alpha
        for i in range(len(data['ASD'])-alpha):
            avg['CM_ASD'][i] = sum(data['CM_ASD'][i:i+alpha])/alpha

    fig2, ((ax1,ax2,ax3,ax4)) = plt.subplots(ncols=4, nrows=1, figsize=(12,3))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        ax.scatter(range(len(data[name])),data[name], alpha=0.2, color=[.7, .7, .7], label=None)
        ax.plot(range(len(avg[name])),avg[name], color='k', label=None)
        ax.set_xlim(0,num)
        if ymax:
            ax.set_ylim(0,ymax)
            ax.set_title(NAMES[name]+'\n$\lambda=%s$' % (l),fontsize=TITLE_SIZE)
        else:
            ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_ylabel('Degree',fontsize=LABEL_SIZE)
        ax.set_xlabel('Rank',fontsize=LABEL_SIZE)
        if num > 200:
            print(name,'Average over first 100:',sum(data[name][:100])/100)
        
        #if i == 0:
        #    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    if ymax:
        prefix = 'ss_scatter_fig_%s' % (l.replace('.','pt'))
    else:
        prefix = 'ss_scatter_fig_%s_noymax' % (l.replace('.','pt'))

    plt.savefig(OUTFILE_DIR+prefix+'.png')
    print('Created '+OUTFILE_DIR+prefix+'.png')

    plt.savefig(OUTFILE_DIR+prefix+'.pdf')
    print('Created '+OUTFILE_DIR+prefix+'.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+prefix+'.pdf '+OUTFILE_DIR+prefix+'.pdf')
    return

def full_ss_scatter_fig(num,l,ymax=None):
    data = {ex:[] for ex in EXPERIMENTS}
    avg = {ex:[] for ex in EXPERIMENTS}
    label = {ex:[] for ex in EXPERIMENTS}
    alpha=15

    infile = OUTFILE_DIR+'SZ_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        process = []
        label_d = []
        label_p = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            disease.append([float(row[3]),int(row[8])])
            process.append([float(row[5]),int(row[8])])
            label_d.append([float(row[3]),row[2]])
            label_p.append([float(row[5]),row[4]])
        disease.sort(reverse=True,key=lambda x: x[0])
        process.sort(reverse=True,key=lambda x: x[0])
        label_d.sort(reverse=True,key=lambda x: x[0])
        label_p.sort(reverse=True,key=lambda x: x[0])
        print(infile,disease[1:10])
        if not num:
            num = len(disease)
        data['SZ'] = [disease[i][1] for i in range(num)]
        data['CM_SZ'] = [process[i][1] for i in range(num)]
        label['SZ'] = [label_d[i][1] for i in range(num)]
        label['CM_SZ'] = [label_p[i][1] for i in range(num)]
        avg['SZ'] = [0]*(len(data['SZ'])-alpha)
        avg['CM_SZ'] = [0]*(len(data['CM_SZ'])-alpha)
        for i in range(len(data['SZ'])-alpha):
            avg['SZ'][i] = sum(data['SZ'][i:i+alpha])/alpha
        for i in range(len(data['SZ'])-alpha):
            avg['CM_SZ'][i] = sum(data['CM_SZ'][i:i+alpha])/alpha

    infile = OUTFILE_DIR+'ASD_1-layer_%s-sinksource_combined_output.txt' % (l)
    with open(infile) as fin:
        disease = []
        process = []
        label_d = []
        label_p = []
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            disease.append([float(row[3]),int(row[8])])
            process.append([float(row[5]),int(row[8])])
            label_d.append([float(row[3]),row[2]])
            label_p.append([float(row[5]),row[4]])
        disease.sort(reverse=True,key=lambda x: x[0])
        process.sort(reverse=True,key=lambda x: x[0])
        label_d.sort(reverse=True,key=lambda x: x[0])
        label_p.sort(reverse=True,key=lambda x: x[0])
        if not num:
            num = len(disease)
        data['ASD'] = [disease[i][1] for i in range(num)]
        data['CM_ASD'] = [process[i][1] for i in range(num)]
        label['ASD'] = [label_d[i][1] for i in range(num)]
        label['CM_ASD'] = [label_p[i][1] for i in range(num)]
        avg['ASD'] = [0]*(len(data['ASD'])-alpha)
        avg['CM_ASD'] = [0]*(len(data['CM_ASD'])-alpha)
        for i in range(len(data['ASD'])-alpha):
            avg['ASD'][i] = sum(data['ASD'][i:i+alpha])/alpha
        for i in range(len(data['ASD'])-alpha):
            avg['CM_ASD'][i] = sum(data['CM_ASD'][i:i+alpha])/alpha

    fig2, ((ax1,ax2,ax3,ax4)) = plt.subplots(ncols=4, nrows=1, figsize=(12,3))
    axes = [ax1,ax2,ax3,ax4]
    scatter_colors = {'Unlabeled':[.7,.7,.7],'Positive':[0,0,.7],'Negative':[.7,0,0]}
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        for this_label in ['Positive','Negative','Unlabeled']:
            x = [j for j in range(len(data[name])) if label[name][j]==this_label]
            y = [data[name][j] for j in range(len(data[name])) if label[name][j]==this_label]
            if len(x)>0:
                ax.scatter(x,y, alpha=0.2, color=scatter_colors[this_label], label=this_label)
        ax.plot(range(len(avg[name])),avg[name], color='k', label=None)
        ax.set_xlim(0,num)
        if ymax:
            ax.set_ylim(0,ymax)
            ax.set_title(NAMES[name]+'\n$\lambda=%s$' % (l),fontsize=TITLE_SIZE)
        else:
            ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_ylabel('Degree',fontsize=LABEL_SIZE)
        ax.set_xlabel('Rank',fontsize=LABEL_SIZE)

        if i==0:
            ax.legend(loc='best',fontsize='small')
        #if i == 0:
        #    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    if ymax:
        prefix = 'full_ss_scatter_fig_%s' % (l.replace('.','pt'))
    else:
        prefix = 'full_ss_scatter_fig_%s_noymax' % (l.replace('.','pt'))

    plt.savefig(OUTFILE_DIR+prefix+'.png')
    print('Created '+OUTFILE_DIR+prefix+'.png')

    plt.savefig(OUTFILE_DIR+prefix+'.pdf')
    print('Created '+OUTFILE_DIR+prefix+'.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+prefix+'.pdf '+OUTFILE_DIR+prefix+'.pdf')


    fig2, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,4))
    scatter_colors = {'Unlabeled':[.7,.7,.7],'Positive':[0,0,.7],'Negative':[.7,0,0]}
    for i in [0]:
        name = 'SZ'
        for this_label in ['Positive','Negative','Unlabeled']:
            x = [j for j in range(len(data[name])) if label[name][j]==this_label]
            y = [data[name][j] for j in range(len(data[name])) if label[name][j]==this_label]
            if len(x)>0:
                ax.scatter(x,y, alpha=0.2, color=scatter_colors[this_label], label=this_label)
        ax.plot(range(len(avg[name])),avg[name], color='k', label=None)
        ax.set_ylim(0,3500)
        ax.set_xlim(0,num)
        if ymax:
            ax.set_title(NAMES[name]+'\n$\lambda=%s$' % (l),fontsize=TITLE_SIZE)
        else:
            ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_ylabel('Degree',fontsize=LABEL_SIZE)
        ax.set_xlabel('Rank',fontsize=LABEL_SIZE)

        if i==0:
            ax.legend(loc='best',fontsize='small')
        #if i == 0:
        #    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    if ymax:
        prefix = 'full_SZ_ss_scatter_fig_%s' % (l.replace('.','pt'))
    else:
        prefix = 'full_SZ_ss_scatter_fig_%s_noymax' % (l.replace('.','pt'))

    plt.savefig(OUTFILE_DIR+prefix+'.png')
    print('Created '+OUTFILE_DIR+prefix+'.png')

    plt.savefig(OUTFILE_DIR+prefix+'.pdf')
    print('Created '+OUTFILE_DIR+prefix+'.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+prefix+'.pdf '+OUTFILE_DIR+prefix+'.pdf')
    return


def sinksource_fig(data):
    sig = {} # asterisks
    with open(OUTFILE_DIR+'sinksource_fig.txt', 'w') as out:
        for disease in ['SZ','CM_SZ','ASD','CM_ASD']:
            sig[disease] = []
            out.write('--- %s ---\n' % (disease))
            for j in range(1,len(LAMBDAS)): 
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(data[disease][1][0], data[disease+'_no_neg'][1][j], equal_var=False)
                if t_value > 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                    sig[disease].append(j+1)
                else: 
                    res = ''
                out.write('Welch\'s t-test Lambda=%s (negs vs. no-negs): %s\n\tno-neg mean: %.4f\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[j],res,mean(data[disease+'_no_neg'][1][j]),t_value,p_value,p_value))
            out.write('\n')
        out.write('----- pairs of SS experiments -----\n')
        for e1,e2 in itertools.combinations(['SZ','CM_SZ','ASD','CM_ASD'],2):
            print(e1,e2)
            t_value,p_value = stats.ttest_ind(data[e1][1][0],data[e2][1][0],equal_var=False)
            if p_value < 0.01:
                res = 'SIGNIFICANT (two-tailed, 0.01)'
            else:
                res = ''
            out.write('Welch\'s t-test SS (%s vs %s): %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (e1,e2,res,t_value,p_value,p_value))

    fig2, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8,6))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        if len(data[name][1][0]) > 0:
            bp1 = ax.boxplot(data[name][1][0], notch=True, positions=[1], widths=.75, sym='', \
                patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        if len(data[name+'_no_neg'][1][1]) > 0:
            bp2 = ax.boxplot(data[name+'_no_neg'][1][1:], notch=True, positions=[2,3,4,5,6], widths=.75, sym='', \
                patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        if len(sig[name]) > 0:
            ax.plot(sig[name],[0.85]*len(sig[name]),'*k')
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_xticks([1,2,3,4,5,6])
        ax.set_xticklabels(['SS','SS+\n$\lambda\!=\!0.01$','SS+\n$\lambda\!=\!0.1$','SS+\n$\lambda\!=\!1$','SS+\n$\lambda\!=\!10$','SS+\n$\lambda\!=\!50$'])
        ax.set_xlim(.5,6.5)
        ax.set_ylim(0.5,0.9)
        ax.set_ylabel('AUC',fontsize=LABEL_SIZE)
        #ax.set_xlabel('Method',fontsize=LABEL_SIZE)
        
        #if i == 0:
        #    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'sinksource_fig.png')
    print('Created '+OUTFILE_DIR+'sinksource_fig.png')

    plt.savefig(OUTFILE_DIR+'sinksource_fig.pdf')
    print('Created '+OUTFILE_DIR+'sinksource_fig.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'sinksource_fig.pdf '+OUTFILE_DIR+'sinksource_fig.pdf')

    return

def roc_fig(data):
    roc_data = read_rocs()
    fig2, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(8,6))
    axes = [ax1,ax2,ax3,ax4]
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]
        label=False
        print('plotting %d' % (len(roc_data[name][1])))
        for i in range(len(roc_data[name][1])):
            rec = roc_data[name][1][i][0]
            prec = roc_data[name][1][i][1]
            if not label:
                ax.plot(rec,prec,color=EXP_COLORS[name],alpha=0.4,label='Mean AUC= %.2f' % (mean(data[name][1][0])))
                label=True
            else:    
                ax.plot(rec,prec,color=EXP_COLORS[name],alpha=0.4,label=None)
        ax.plot([0,1],[0,1],'--',color='k',label=None)
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        #ax.set_xlim(.5,6.5)
        #ax.set_ylim(0.5,0.9)
        ax.set_ylabel('True Positive Rate',fontsize=LABEL_SIZE)
        ax.set_xlabel('False Positive Rate',fontsize=LABEL_SIZE)
        ax.legend(loc='best',fontsize=LABEL_SIZE)

    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'roc_fig.png')
    print('Created '+OUTFILE_DIR+'roc_fig.png')

    plt.savefig(OUTFILE_DIR+'roc_fig.pdf')
    print('Created '+OUTFILE_DIR+'roc_fig.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'roc_fig.pdf '+OUTFILE_DIR+'roc_fig.pdf')

#def read_roc_data():
    ## read ROC curve data for all Layer 1 


####### OLD FIGURES
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
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_xticks([j+1 for j in range(len(data[name][1]))])
        ax.set_xticklabels(['0','0.01','0.1','1','10','50'])
        
        ax.set_ylim(0.5,0.9)
        
        ax.set_ylabel('AUC',fontsize=LABEL_SIZE)
        ax.set_xlabel('$\lambda$',fontsize=LABEL_SIZE)
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_lambda.png')
    print('Created '+OUTFILE_DIR+'compare_lambda.png')

    plt.savefig(OUTFILE_DIR+'compare_lambda.pdf')
    print('Created '+OUTFILE_DIR+'compare_lambda.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'compare_lambda.pdf '+OUTFILE_DIR+'compare_lambda.pdf')
    
    return

def figure_2(data):

    ## SS negative analysis

    fig2, axes = plt.subplots(ncols=4, nrows=1, sharey=True, figsize=(10,4))
    for i in range(len(EXPERIMENTS)):
        name = EXPERIMENTS[i]
        ax = axes[i]

        bp1 = ax.boxplot(data[name][1][0], notch=True, positions=[1], widths=.9, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#8193ef'))
        bp2 = ax.boxplot(data[name+'_old_neg'][1][0], notch=True, positions=[2], widths=.9, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        bp3 = ax.boxplot(data[name+'_rand_neg'][1][0], notch=True, positions=[3,], widths=.9, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#eeb5ff',color='#A152B8'))
        bp4 = ax.boxplot(data[name+'_rand_neg_deg'][1][0], notch=True, positions=[4], widths=.9, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#FFE888',color='#E1C23E'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

        sig_x = []
        sig_y = []
        negs_to_check = [name+'_old_neg',name+'_rand_neg',name+'_rand_neg_deg']
        for j in range(len(negs_to_check)):
            t_value,p_value = stats.ttest_ind(data[name][1][0], data[negs_to_check[j]][1][0], equal_var=False)
            if t_value > 0 and p_value/2 < 0.01:
                sig_x.append(j+2)
                sig_y.append(0.85)
        ax.plot(sig_x,sig_y,'*k')

        for x in [5,10,15,20,25]:
            ax.plot([x,x],[0.5,0.9],color='lightgrey',alpha=0.5)
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_xticks([])
        #ax.set_xticklabels(['Curated','Krishnan et al.','Random','Random\n(Deg. Pres.)'])
        ax.set_xlim(0,5)
        ax.set_ylim(0.5,0.9)
        if i==0:
            ax.set_ylabel('AUC',fontsize=LABEL_SIZE)
        #ax.set_xlabel('Negative Set',fontsize=LABEL_SIZE)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        
        if i==0:
            handles = [bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0], bp4['boxes'][0]]
            leg_names = ['Curated Negatives', 'Krishnan et al. Negatives', 'Random Negatives', 'Degree-Preserving\nRandom Negatives']
            lgd = ax.legend(handles, leg_names, ncol=4, fontsize='small',bbox_to_anchor=(2.3,-.1),loc='center')
    
    #plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_rand_negs_SS.png', bbox_extra_artists=(lgd,axes[0],axes[1],axes[2],axes[3]), bbox_inches='tight',bbox_pad=1)
    print('Created '+OUTFILE_DIR+'compare_rand_negs_SS.png')

    plt.savefig(OUTFILE_DIR+'compare_rand_negs_SS.pdf', bbox_extra_artists=(lgd,axes[0],axes[1],axes[2],axes[3]), bbox_inches='tight',bbox_pad=1)
    print('Created '+OUTFILE_DIR+'compare_rand_negs_SS.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'compare_rand_negs_SS.pdf '+OUTFILE_DIR+'compare_rand_negs_SS.pdf')

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

        for disease in ['SZ','CM_SZ','ASD','CM_ASD']:
            out.write('--- %s KRISHNAN ET AL NEGATIVES ---\n' % (disease))
            for j in range(len(LAMBDAS)): 
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(data[disease][1][j], data[disease+'_old_neg'][1][j], equal_var=False)
                if t_value > 0 and p_value/2 < 0.01:
                    res = 'SIGNIFICANT (one-tailed, 0.01)'
                else: 
                    res = ''
                out.write('Welch\'s t-test Lambda=%s (negs vs. old negs): %s\n\tt-value: %.4f\n\tp-value: %.2e (%.4f)\n' % (LAMBDAS[j],res,t_value,p_value,p_value))
            out.write('\n')

        out.write('Two-way ANOVA KRISHNAN ET AL NEGATIVES\n')
        organized_data_1=[data[exp][1] for exp in EXPERIMENTS]
        organized_data_2=[data[exp+'_old_neg'][1] for exp in EXPERIMENTS]
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
                    name = 'Including Old NEgatives\t'
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
        bp2 = ax.boxplot(data[name+'_old_neg'][1], notch=True, positions=[2,7,12,17,22,27], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#b0f2c2',color='#3A9A54'))
        bp3 = ax.boxplot(data[name+'_rand_neg'][1], notch=True, positions=[3,8,13,18,23,28], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#eeb5ff',color='#A152B8'))
        bp4 = ax.boxplot(data[name+'_rand_neg_deg'][1], notch=True, positions=[4,9,14,19,24,29], widths=1, sym='', \
            patch_artist=True, boxprops=dict(facecolor='#FFE888',color='#E1C23E'))
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        for x in [5,10,15,20,25]:
            ax.plot([x,x],[0.5,0.9],color='lightgrey',alpha=0.5)
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        ax.set_xticks([2.5,7.5,12.5,17.5,22.5,27.5])
        ax.set_xticklabels(['SS\n($\lambda$=0)','$\lambda$=0.01','$\lambda$=0.1','$\lambda$=1','$\lambda$=10','$\lambda$=50'])
        ax.set_xlim(0,30)
        ax.set_ylim(0.5,0.9)
        ax.set_ylabel('AUC',fontsize=LABEL_SIZE)
        ax.set_xlabel('$\lambda$',fontsize=LABEL_SIZE)
        
        ax.legend([bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0], bp4['boxes'][0]], ['Curated Negatives', 'Krishnan et al. Negatives', 'Random Negatives','Deg. Preserving Random Negatives'], ncol=2, loc='best', fontsize='small')
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'compare_rand_negs.png')
    print('Created '+OUTFILE_DIR+'compare_rand_negs.png')

    plt.savefig(OUTFILE_DIR+'compare_rand_negs.pdf')
    print('Created '+OUTFILE_DIR+'compare_rand_negs.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'compare_rand_negs.pdf '+OUTFILE_DIR+'compare_rand_negs.pdf')

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
            ax.set_ylabel('AUC',fontsize=TICK_SIZE)
        ax.set_xlabel('Layers',fontsize=TICK_SIZE)
        ax.set_title(NAMES[name],fontsize=LABEL_SIZE)
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'comparing_layers_best_lambdas.png')
    print('Created '+OUTFILE_DIR+'comparing_layers_best_lambdas.png')

    plt.savefig(OUTFILE_DIR+'comparing_layers_best_lambdas.pdf')
    print('Created '+OUTFILE_DIR+'comparing_layers_best_lambdas.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'comparing_layers_best_lambdas.pdf '+OUTFILE_DIR+'comparing_layers_best_lambdas.pdf')

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
        ax.set_ylabel('AUC',fontsize=LABEL_SIZE)
        ax.set_xlabel('$\lambda$',fontsize=LABEL_SIZE)
        ax.set_title(NAMES[name],fontsize=TITLE_SIZE)
        if i==0:
            ax.legend([bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0]], ['1 Layer', '2 Layers', '3 Layers'], loc='best', fontsize=TICK_SIZE)
    
    plt.tight_layout()
    plt.savefig(OUTFILE_DIR+'comparing_layers.png')
    print('Created '+OUTFILE_DIR+'comparing_layers.png')

    plt.savefig(OUTFILE_DIR+'comparing_layers.pdf')
    print('Created '+OUTFILE_DIR+'comparing_layers.pdf')
    os.system('pdfcrop '+OUTFILE_DIR+'comparing_layers.pdf '+OUTFILE_DIR+'comparing_layers.pdf')

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

def median(bunch_of_numbers):
    return np.median(bunch_of_numbers)


main()

# manure=[[13.7, 16.8, 14.9, 17.6, 16.5],[17.4, 13.5, 15.1, 15.4, 13.2]]
# no_manure=[[16.0, 16.1, 13.0, 16.7, 13.2],[13.4, 11.6, 14.7, 9.7, 11.9]]


# x= two_way_ANOVA_test(manure, no_manure)

# for i in x:
#     print(i)



























