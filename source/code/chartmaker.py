import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 


#SET USE_SD=False to get Error is IRQ, not standard deviation
USE_SD=True

def main():
    
    figure_1()
    figure_2()
    figure_3()
    
    return

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

            ASD.append(float(line[1].strip('\''))) #ASD AUCs always in second column
            CM.append(float(line[2].strip('\''))) #CM AUCs always in third column
        x+=1

    means_IRQs = [np.mean(SZ),np.mean(ASD),np.mean(CM),IQR(SZ),IQR(ASD),IQR(CM)]

    #return SZ, ASD, CM, means_IRQs
    return SZ, ASD, CM


def IQR(dist):
    return [np.percentile(dist, 75) - np.mean(dist), np.mean(dist)-np.percentile(dist, 25)]


def figure_1():
    
    zero_file='outfiles/SZ_1-layer_0.150_auc.txt'
    SZ_zero, ASD_zero, CM_zero = file_parser(zero_file)

    point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    SZ_pt_zero_one, ASD_pt_zero_one, CM_pt_zero_one = file_parser(point_zero_one_file)

    point_one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    SZ_pt_one, ASD_pt_one, CM_pt_one = file_parser(point_one_file)

    ten_file='outfiles/SZ_1-layer_10-sinksource_0.150_auc.txt'
    SZ_ten, ASD_ten, CM_ten = file_parser(ten_file)

    fifty_file='outfiles/SZ_1-layer_50-sinksource_0.150_auc.txt'
    SZ_fifty, ASD_fifty, CM_fifty = file_parser(fifty_file)

    SZ_data = [SZ_zero, SZ_pt_zero_one, SZ_pt_one, SZ_ten, SZ_fifty]
    ASD_data = [ASD_zero, ASD_pt_zero_one, ASD_pt_one, ASD_ten, ASD_fifty]
    CM_data = [CM_zero, CM_pt_zero_one, CM_pt_one, CM_ten, CM_fifty]

    fig1, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, sharey=True, figsize=(7,4))

    bp1 = ax1.boxplot(SZ_data, notch=True, positions=[2,4,6,8,10], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','10','50'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,12)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia')
    

    bp3 = ax2.boxplot(ASD_data, notch=True, positions=[2,4,6,8,10], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','10','50'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,12)
    ax2.set_xlabel('$\lambda$')
    ax2.set_title('Autism')


    bp5 = ax3.boxplot(CM_data, notch=True, positions=[2,4,6,8,10], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','10','50'])
    ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax3.set_xlim(0,12)
    ax3.set_title('Cell Motility')

    plt.savefig('outfiles/compare_lambda.png')

    print('Created outfiles/compare_lambda.png')
    
    return


def figure_2():
    
    #Get negative data
    zero_file='outfiles/SZ_1-layer_0.150_auc.txt'
    SZ_zero, ASD_zero, CM_zero = file_parser(zero_file)

    point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    SZ_pt_zo, ASD_pt_zo, CM_pt_zo = file_parser(point_zero_one_file)

    point_one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    SZ_pt_o, ASD_pt_o, CM_pt_o = file_parser(point_one_file)

    ten_file='outfiles/SZ_1-layer_10-sinksource_0.150_auc.txt'
    SZ_ten, ASD_ten, CM_ten = file_parser(ten_file)

    fifty_file='outfiles/SZ_1-layer_50-sinksource_0.150_auc.txt'
    SZ_fifty, ASD_fifty, CM_fifty = file_parser(fifty_file)

    #Get no negative data
    no_neg_zero_file='outfiles/SZ_1-layer_no_neg_0.150_auc.txt'
    SZ_no_neg_zero, ASD_no_neg_zero, CM_no_neg_zero = file_parser(no_neg_zero_file)

    no_neg_point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_no_neg_0.150_auc.txt'
    SZ_no_neg_pt_zo, ASD_no_neg_pt_zo, CM_no_neg_pt_zo = file_parser(no_neg_point_zero_one_file)

    no_neg_point_one_file='outfiles/SZ_1-layer_0.1-sinksource_no_neg_0.150_auc.txt'
    SZ_no_neg_pt_o, ASD_no_neg_pt_o, CM_no_neg_pt_o = file_parser(no_neg_point_one_file)

    no_neg_ten_file='outfiles/SZ_1-layer_10-sinksource_no_neg_0.150_auc.txt'
    SZ_no_neg_ten, ASD_no_neg_ten, CM_no_neg_ten = file_parser(no_neg_ten_file)

    no_neg_fifty_file='outfiles/SZ_1-layer_50-sinksource_no_neg_0.150_auc.txt'
    SZ_no_neg_fifty, ASD_no_neg_fifty, CM_no_neg_fifty = file_parser(no_neg_fifty_file)

    #For each lambda value, generate three p-values: SZ neg vs no neg, ASD neg vs no neg, CM neg vs no neg
    #Create lists to iterate through for each set 

    zero_neg = [SZ_zero, ASD_zero, CM_zero]
    zero_no_neg = [SZ_no_neg_zero, ASD_no_neg_zero, CM_no_neg_zero]

    pt_zo_neg = [SZ_pt_zo, ASD_pt_zo, CM_pt_zo]
    pt_zo_no_neg = [SZ_no_neg_pt_zo, ASD_no_neg_pt_zo, CM_no_neg_pt_zo]

    pt_one_neg = [SZ_pt_o, ASD_pt_o, CM_pt_o]
    pt_one_no_neg = [SZ_no_neg_pt_o, ASD_no_neg_pt_o, CM_no_neg_pt_o]

    pt_ten_neg = [SZ_ten, ASD_ten, CM_ten]
    pt_ten_no_neg = [SZ_no_neg_ten, ASD_no_neg_ten, CM_no_neg_ten]

    pt_fifty_neg = [SZ_fifty, ASD_fifty, CM_fifty] 
    pt_fifty_no_neg = [SZ_no_neg_fifty, ASD_no_neg_fifty, CM_no_neg_fifty]

    neg_lists = [zero_neg, pt_zo_neg, pt_one_neg, pt_ten_neg, pt_fifty_neg]
    no_neg_lists = [zero_no_neg, pt_zo_no_neg, pt_one_no_neg, pt_ten_no_neg, pt_fifty_no_neg]

    #Compare no negatives to negatives and print p-values to an outfile

    with open('outfiles/MWU_PVALUES_Neg_vs_NoNeg.txt', 'w') as out:
        #first iterate through each constant, then iterate through each positive set for that constant to get the list of AUCs
        for i in range(len(neg_lists)): #the two data set lists are always the same length 
            if i == 0:
                out.write('Lambda=0\n\n')
            elif i == 1:
                out.write('Lambda=0.01\n\n')
            elif i == 2:
                out.write('Lambda=0.1\n\n')
            elif i == 3:
                out.write('Lambda=10\n\n')
            else:
                out.write('Lambda=50\n\n')
            for j in range(len(zero_neg)):
                if j == 0:
                    out.write('Schizophrenia\t')
                elif j == 1:
                    out.write('Autism\t')
                else:
                    out.write('Cell Motility\t')
                U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                out.write(str(p_value)+'\n')
            out.write('\n')


    
    #SZ_data = [SZ_zero, SZ_no_neg_zero, SZ_pt_zo, SZ_no_neg_pt_zo, SZ_pt_o, SZ_no_neg_pt_o, SZ_ten, SZ_no_neg_ten, SZ_fifty, SZ_no_neg_fifty]
    SZ_neg_data = [SZ_zero, SZ_pt_zo, SZ_pt_o, SZ_ten, SZ_fifty]
    SZ_no_neg_data = [SZ_no_neg_zero, SZ_no_neg_pt_zo, SZ_no_neg_pt_o, SZ_no_neg_ten, SZ_no_neg_fifty]

    #ASD_data = [ASD_zero, ASD_no_neg_zero, ASD_pt_zo, ASD_no_neg_pt_zo, ASD_pt_o, ASD_no_neg_pt_o, ASD_ten, ASD_no_neg_ten, ASD_fifty, ASD_no_neg_fifty]
    ASD_neg_data = [ASD_zero, ASD_pt_zo, ASD_pt_o, ASD_ten, ASD_fifty]
    ASD_no_neg_data = [ASD_no_neg_zero, ASD_no_neg_pt_zo, ASD_no_neg_pt_o, ASD_no_neg_ten, ASD_no_neg_fifty]

    #CM_data = [CM_zero, CM_no_neg_zero, CM_pt_zo, CM_no_neg_pt_zo, CM_pt_o, CM_no_neg_pt_o, CM_ten, CM_no_neg_ten, CM_fifty, CM_no_neg_fifty]
    CM_neg_data = [CM_zero, CM_pt_zo, CM_pt_o, CM_ten, CM_fifty]
    CM_no_neg_data = [CM_no_neg_zero, CM_no_neg_pt_zo, CM_no_neg_pt_o, CM_no_neg_ten, CM_no_neg_fifty]


    fig2, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, sharey=True, figsize=(7,4))
    
    bp1 = ax1.boxplot(SZ_neg_data, notch=True, positions=[2,7,12,17,22], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp2 = ax1.boxplot(SZ_no_neg_data, notch=True, positions=[4,9,14,19,24], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','10','50'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,26)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia')
    ax1.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='upper left', fontsize='x-small')
    

    bp3 = ax2.boxplot(ASD_neg_data, notch=True, positions=[2,7,12,17,22], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp4 = ax2.boxplot(ASD_no_neg_data, notch=True, positions=[4,9,14,19,24], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','10','50'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,26)
    ax2.set_xlabel('$\lambda$')
    ax2.set_title('Autism')


    bp5 = ax3.boxplot(CM_neg_data, notch=True, positions=[2,7,12,17,22], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp6 = ax3.boxplot(CM_no_neg_data, notch=True, positions=[4,9,14,19,24], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','10','50'])
    ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax3.set_xlim(0,26)
    ax3.set_title('Cell Motility')

    plt.savefig('outfiles/compare_neg_noneg.png')

    print('Created outfiles/compare_neg_noneg.png')

    return
    


def figure_3():

    one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    SZ_first_one, ASD_first_one, CM_first_one = file_parser(one_file)

    one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    SZ_other_one, ASD_other_one, CM_other_one = file_parser(one_file)

    two_file='outfiles/SZ_2-layer_10-sinksource_0.150_auc.txt'
    SZ_two, ASD_two, CM_two = file_parser(two_file)

    three_file='outfiles/SZ_3-layer_10-sinksource_0.150_auc.txt'
    SZ_three, ASD_three, CM_three = file_parser(three_file)

    #Compare 2 layer to 1 layer and 3 layer for each positive set
    #Create lists for each positive set of what data to compare [2 layer AUCs, 1 layer, 3 layer]

    SZ_sets = [SZ_two, SZ_first_one, SZ_three]
    ASD_sets = [ASD_two, ASD_other_one, ASD_three]
    CM_sets = [CM_two, CM_other_one, CM_three]

    Mann_test_lists = [SZ_sets, ASD_sets, CM_sets]

    with open('outfiles/MWU_PVALUES_123layers.txt', 'w') as out:
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

    SZ_data = [SZ_first_one, SZ_two, SZ_three]
    ASD_data = [ASD_other_one, ASD_two, ASD_three]
    CM_data = [CM_other_one, CM_two, CM_three]

    fig3, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, sharey=True, figsize=(7,4))

    bp1 = ax1.boxplot(SZ_data, notch=True, positions=[2,4,6], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['1','2','3'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,8)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia')
    

    bp3 = ax2.boxplot(ASD_data, notch=True, positions=[2,4,6], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['1','2','3'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,8)
    ax2.set_xlabel('Layers')
    ax2.set_title('Autism')


    bp5 = ax3.boxplot(CM_data, notch=True, positions=[2,4,6], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['1','2','3'])
    ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax3.set_xlim(0,8)
    ax3.set_title('Cell Motility')

    plt.savefig('outfiles/comparing_layers.png')

    print('Created outfiles/comparing_layers.png')


    return
    





main()