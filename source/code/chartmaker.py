import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 


#SET USE_SD=False to get Error is IRQ, not standard deviation
USE_SD=True

def main():
    
    figure_1()
    figure_2()
    figure_3_full()
    
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

            CM.append(float(line[1].strip('\''))) #CM AUCs always in third column
        x+=1

    means_IRQs = [np.mean(SZ),np.mean(ASD),np.mean(CM),IQR(SZ),IQR(CM)]

    #return SZ, ASD, CM, means_IRQs
    return SZ, CM


def IQR(dist):
    return [np.percentile(dist, 75) - np.mean(dist), np.mean(dist)-np.percentile(dist, 25)]


def figure_1():
    
    zero_file='../outfiles/SZ_1-layer_0-sinksource_auc.txt'
    SZ_zero, CM_SZ_zero = file_parser(zero_file)

    point_zero_one_file='../outfiles/SZ_1-layer_0.01-sinksource_auc.txt'
    SZ_pt_zero_one, CM_SZ_pt_zero_one = file_parser(point_zero_one_file)

    point_one_file='../outfiles/SZ_1-layer_0.1-sinksource_auc.txt'
    SZ_pt_one, CM_SZ_pt_one = file_parser(point_one_file)

    one_file='../outfiles/SZ_1-layer_1-sinksource_auc.txt'
    SZ_one, CM_SZ_one = file_parser(one_file)

    ten_file='../outfiles/SZ_1-layer_10-sinksource_auc.txt'
    SZ_ten, CM_SZ_ten = file_parser(ten_file)

    fifty_file='../outfiles/SZ_1-layer_50-sinksource_auc.txt'
    SZ_fifty, CM_SZ_fifty = file_parser(fifty_file)

    SZ_data = [SZ_zero, SZ_pt_zero_one, SZ_pt_one, SZ_one, SZ_ten, SZ_fifty]
    CM_SZ_data = [CM_SZ_zero, CM_SZ_pt_zero_one, CM_SZ_pt_one,CM_SZ_one, CM_SZ_ten, CM_SZ_fifty]


    zero_file='../outfiles/ASD_1-layer_0-sinksource_auc.txt'
    ASD_zero, CM_ASD_zero = file_parser(zero_file)

    point_zero_one_file='../outfiles/ASD_1-layer_0.01-sinksource_auc.txt'
    ASD_pt_zero_one, CM_ASD_pt_zero_one = file_parser(point_zero_one_file)

    point_one_file='../outfiles/ASD_1-layer_0.1-sinksource_auc.txt'
    ASD_pt_one, CM_ASD_pt_one = file_parser(point_one_file)

    one_file='../outfiles/ASD_1-layer_1-sinksource_auc.txt'
    ASD_one, CM_ASD_one = file_parser(one_file)

    ten_file='../outfiles/ASD_1-layer_10-sinksource_auc.txt'
    ASD_ten, CM_ASD_ten = file_parser(ten_file)

    fifty_file='../outfiles/ASD_1-layer_50-sinksource_auc.txt'
    ASD_fifty, CM_ASD_fifty = file_parser(fifty_file)

    ASD_data = [ASD_zero, ASD_pt_zero_one, ASD_pt_one, ASD_one, ASD_ten, ASD_fifty]
    CM_ASD_data = [CM_ASD_zero, CM_ASD_pt_zero_one, CM_ASD_pt_one,CM_ASD_one, CM_ASD_ten, CM_ASD_fifty]


    all_Data=[SZ_data, CM_SZ_data, ASD_data, CM_ASD_data]
    with open('../outfiles/sinksource+_constant_statistical_tests.txt', 'w') as out:

        f_value, p_value = stats.f_oneway(SZ_zero, SZ_pt_zero_one, SZ_pt_one, SZ_one, SZ_ten, SZ_fifty)
        out.write('Schizophrenia\t'+str(p_value)+'\n')

        f_value, p_value = stats.f_oneway(CM_SZ_zero, CM_SZ_pt_zero_one, CM_SZ_pt_one,CM_SZ_one, CM_SZ_ten, CM_SZ_fifty)
        out.write('Cell Motility-Schizophrenia\t'+str(p_value)+'\n')

        f_value, p_value = stats.f_oneway(ASD_zero, ASD_pt_zero_one, ASD_pt_one, ASD_one, ASD_ten, ASD_fifty)
        out.write('ASD\t'+str(p_value)+'\n')

        f_value, p_value = stats.f_oneway(CM_ASD_zero, CM_ASD_pt_zero_one, CM_ASD_pt_one,CM_ASD_one, CM_ASD_ten, CM_ASD_fifty)
        out.write('Cell Motility-ASD\t'+str(p_value)+'\n')

    fig1, (ax1, ax2,ax3,ax4) = plt.subplots(ncols=4, nrows=1, sharey=True, figsize=(7,4))

    bp1 = ax1.boxplot(SZ_data, notch=True, positions=[2,4,6,8,10,12], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','1','10','50'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,14)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia', fontsize=8)
    


    bp3 = ax2.boxplot(CM_SZ_data, notch=True, positions=[2,4,6,8,10,12], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','1','10','50'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,14)
    ax2.set_title('Cell Motility-Schizophrenia', fontsize=8)



    bp5 = ax3.boxplot(ASD_data, notch=True, positions=[2,4,6,8,10,12], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','1','10','50'])
    ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax3.set_xlim(0,14)
    ax3.set_ylabel('AUC')
    ax3.set_title('ASD', fontsize=8)
    


    bp7 = ax4.boxplot(CM_ASD_data, notch=True, positions=[2,4,6,8,10,12], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['0','0.01','0.1','1','10','50'])
    ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax4.set_xlim(0,14)
    ax4.set_title('Cell Motility-ASD', fontsize=8)

    plt.savefig('../outfiles/compare_lambda.png')

    print('Created ../outfiles/compare_lambda.png')
    
    return


def figure_2():
    
    #Get negative data
    zero_file='../outfiles/SZ_1-layer_0-sinksource_auc.txt'
    SZ_zero, CM_SZ_zero = file_parser(zero_file)

    point_zero_one_file='../outfiles/SZ_1-layer_0.01-sinksource_auc.txt'
    SZ_pt_zo, CM_SZ_pt_zo = file_parser(point_zero_one_file)

    point_one_file='../outfiles/SZ_1-layer_0.1-sinksource_auc.txt'
    SZ_pt_o, CM_SZ_pt_o = file_parser(point_one_file)

    one_file='../outfiles/SZ_1-layer_1-sinksource_auc.txt'
    SZ_o, CM_SZ_o = file_parser(one_file)

    ten_file='../outfiles/SZ_1-layer_10-sinksource_auc.txt'
    SZ_ten, CM_SZ_ten = file_parser(ten_file)

    fifty_file='../outfiles/SZ_1-layer_50-sinksource_auc.txt'
    SZ_fifty, CM_SZ_fifty = file_parser(fifty_file)

    #Get no negative data
    no_neg_zero_file='../outfiles/SZ_1-layer_0-sinksource_no_neg_auc.txt'
    SZ_no_neg_zero, CM_SZ_no_neg_zero = file_parser(no_neg_zero_file)

    no_neg_point_zero_one_file='../outfiles/SZ_1-layer_0.01-sinksource_no_neg_auc.txt'
    SZ_no_neg_pt_zo, CM_SZ_no_neg_pt_zo = file_parser(no_neg_point_zero_one_file)

    no_neg_point_one_file='../outfiles/SZ_1-layer_0.1-sinksource_no_neg_auc.txt'
    SZ_no_neg_pt_o, CM_SZ_no_neg_pt_o = file_parser(no_neg_point_one_file)

    no_neg_one_file='../outfiles/SZ_1-layer_1-sinksource_no_neg_auc.txt'
    SZ_no_neg_o, CM_SZ_no_neg_o = file_parser(no_neg_one_file)

    no_neg_ten_file='../outfiles/SZ_1-layer_10-sinksource_no_neg_auc.txt'
    SZ_no_neg_ten, CM_SZ_no_neg_ten = file_parser(no_neg_ten_file)

    no_neg_fifty_file='../outfiles/SZ_1-layer_50-sinksource_no_neg_auc.txt'
    SZ_no_neg_fifty, CM_SZ_no_neg_fifty = file_parser(no_neg_fifty_file)


    #Get negative data
    zero_file='../outfiles/ASD_1-layer_0-sinksource_auc.txt'
    ASD_zero, CM_ASD_zero = file_parser(zero_file)

    point_zero_one_file='../outfiles/ASD_1-layer_0.01-sinksource_auc.txt'
    ASD_pt_zo, CM_ASD_pt_zo = file_parser(point_zero_one_file)

    point_one_file='../outfiles/ASD_1-layer_0.1-sinksource_auc.txt'
    ASD_pt_o, CM_ASD_pt_o = file_parser(point_one_file)

    one_file='../outfiles/ASD_1-layer_1-sinksource_auc.txt'
    ASD_o, CM_ASD_o = file_parser(one_file)

    ten_file='../outfiles/ASD_1-layer_10-sinksource_auc.txt'
    ASD_ten, CM_ASD_ten = file_parser(ten_file)

    fifty_file='../outfiles/ASD_1-layer_50-sinksource_auc.txt'
    ASD_fifty, CM_ASD_fifty = file_parser(fifty_file)

    #Get no negative data
    no_neg_zero_file='../outfiles/ASD_1-layer_0-sinksource_no_neg_auc.txt'
    ASD_no_neg_zero, CM_ASD_no_neg_zero = file_parser(no_neg_zero_file)

    no_neg_point_zero_one_file='../outfiles/ASD_1-layer_0.01-sinksource_no_neg_auc.txt'
    ASD_no_neg_pt_zo, CM_ASD_no_neg_pt_zo = file_parser(no_neg_point_zero_one_file)

    no_neg_point_one_file='../outfiles/ASD_1-layer_0.1-sinksource_no_neg_auc.txt'
    ASD_no_neg_pt_o, CM_ASD_no_neg_pt_o = file_parser(no_neg_point_one_file)

    no_neg_one_file='../outfiles/ASD_1-layer_1-sinksource_no_neg_auc.txt'
    ASD_no_neg_o, CM_ASD_no_neg_o = file_parser(no_neg_one_file)

    no_neg_ten_file='../outfiles/ASD_1-layer_10-sinksource_no_neg_auc.txt'
    ASD_no_neg_ten, CM_ASD_no_neg_ten = file_parser(no_neg_ten_file)

    no_neg_fifty_file='../outfiles/ASD_1-layer_50-sinksource_no_neg_auc.txt'
    ASD_no_neg_fifty, CM_ASD_no_neg_fifty = file_parser(no_neg_fifty_file)

    #For each lambda value, generate three p-values: SZ neg vs no neg, ASD neg vs no neg, CM_SZ neg vs no neg
    #Create lists to iterate through for each set 

    zero_neg = [SZ_zero, CM_SZ_zero, ASD_zero, CM_ASD_zero]
    zero_no_neg = [SZ_no_neg_zero, CM_SZ_no_neg_zero, ASD_no_neg_zero, CM_ASD_no_neg_zero]

    pt_zo_neg = [SZ_pt_zo, CM_SZ_pt_zo, ASD_pt_zo, CM_ASD_pt_zo]
    pt_zo_no_neg = [SZ_no_neg_pt_zo, CM_SZ_no_neg_pt_zo, ASD_no_neg_pt_zo, CM_ASD_no_neg_pt_zo]

    pt_one_neg = [SZ_pt_o, CM_SZ_pt_o, ASD_pt_o, CM_ASD_pt_o]
    pt_one_no_neg = [SZ_no_neg_pt_o, CM_SZ_no_neg_pt_o, ASD_no_neg_pt_o, CM_ASD_no_neg_pt_o]

    one_neg = [SZ_o, CM_SZ_o, ASD_o, CM_ASD_o]
    one_no_neg = [SZ_no_neg_o, CM_SZ_no_neg_o, ASD_no_neg_o, CM_ASD_no_neg_o]

    pt_ten_neg = [SZ_ten, CM_SZ_ten, ASD_ten, CM_ASD_ten]
    pt_ten_no_neg = [SZ_no_neg_ten, CM_SZ_no_neg_ten, ASD_no_neg_ten, CM_ASD_no_neg_ten]

    pt_fifty_neg = [SZ_fifty,  CM_SZ_fifty, ASD_fifty,  CM_ASD_fifty] 
    pt_fifty_no_neg = [SZ_no_neg_fifty, CM_SZ_no_neg_fifty, ASD_no_neg_fifty, CM_ASD_no_neg_fifty]

    neg_lists = [zero_neg, pt_zo_neg, pt_one_neg, one_neg, pt_ten_neg, pt_fifty_neg]
    no_neg_lists = [zero_no_neg, pt_zo_no_neg, pt_one_no_neg, one_no_neg, pt_ten_no_neg, pt_fifty_no_neg]

    #Compare no negatives to negatives and print p-values to an outfile



    #SZ_data = [SZ_zero, SZ_no_neg_zero, SZ_pt_zo, SZ_no_neg_pt_zo, SZ_pt_o, SZ_no_neg_pt_o, SZ_ten, SZ_no_neg_ten, SZ_fifty, SZ_no_neg_fifty]
    SZ_neg_data = [SZ_zero, SZ_pt_zo, SZ_pt_o, SZ_o,SZ_ten, SZ_fifty]
    SZ_no_neg_data = [SZ_no_neg_zero, SZ_no_neg_pt_zo, SZ_no_neg_pt_o, SZ_no_neg_o, SZ_no_neg_ten, SZ_no_neg_fifty]

    #CM_SZ_data = [CM_SZ_zero, CM_SZ_no_neg_zero, CM_SZ_pt_zo, CM_SZ_no_neg_pt_zo, CM_SZ_pt_o, CM_SZ_no_neg_pt_o, CM_SZ_ten, CM_SZ_no_neg_ten, CM_SZ_fifty, CM_SZ_no_neg_fifty]
    CM_SZ_neg_data = [CM_SZ_zero, CM_SZ_pt_zo, CM_SZ_pt_o, CM_SZ_o, CM_SZ_ten, CM_SZ_fifty]
    CM_SZ_no_neg_data = [CM_SZ_no_neg_zero, CM_SZ_no_neg_pt_zo, CM_SZ_no_neg_pt_o, CM_SZ_no_neg_o, CM_SZ_no_neg_ten, CM_SZ_no_neg_fifty]

    # #ASD_data = [ASD_zero, ASD_no_neg_zero, ASD_pt_zo, ASD_no_neg_pt_zo, ASD_pt_o, ASD_no_neg_pt_o, ASD_ten, ASD_no_neg_ten, ASD_fifty, ASD_no_neg_fifty]
    ASD_neg_data = [ASD_zero, ASD_pt_zo, ASD_pt_o, ASD_o,ASD_ten, ASD_fifty]
    ASD_no_neg_data = [ASD_no_neg_zero, ASD_no_neg_pt_zo, ASD_no_neg_pt_o, ASD_no_neg_o,ASD_no_neg_ten, ASD_no_neg_fifty]

    CM_ASD_neg_data = [CM_ASD_zero, CM_ASD_pt_zo, CM_ASD_pt_o, CM_ASD_o, CM_ASD_ten, CM_ASD_fifty]
    CM_ASD_no_neg_data = [CM_ASD_no_neg_zero, CM_ASD_no_neg_pt_zo, CM_ASD_no_neg_pt_o, CM_ASD_no_neg_o, CM_ASD_no_neg_ten, CM_ASD_no_neg_fifty]

    with open('../outfiles/Neg_vs_NoNeg_statistical_tests.txt', 'w') as out:
        #first iterate through each constant, then iterate through each positive set for that constant to get the list of AUCs
        out.write('T-Tests\n')
        for i in range(len(neg_lists)): #the two data set lists are always the same length 
            if i == 0:
                out.write('Lambda=0\n\n')
            elif i == 1:
                out.write('Lambda=0.01\n\n')
            elif i == 2:
                out.write('Lambda=0.1\n\n')
            elif i == 3:
                out.write('Lambda=1\n\n')
            elif i == 4:
                out.write('Lambda=10\n\n')
            else:
                out.write('Lambda=50\n\n')
            for j in range(len(zero_neg)):
                if j == 0:
                    out.write('Schizophrenia\t')
                if j == 1:
                    out.write('Cell Motility-SZ\t')
                if j == 2:
                    out.write('ASD\t')
                if j == 3:
                    out.write('Cell Motility-ASD\t')
                # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')
                t_value,p_value = stats.ttest_ind(neg_lists[i][j], no_neg_lists[i][j], equal_var=False)
                out.write(str(p_value)+'\n')
            out.write('\n')

        out.write('Two-way ANOVA\n')
        organized_data_1=[SZ_neg_data, CM_SZ_neg_data, ASD_neg_data, CM_ASD_neg_data]
        organized_data_2=[SZ_no_neg_data, CM_SZ_no_neg_data, ASD_no_neg_data, CM_ASD_no_neg_data]
        

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

        


    


    fig2, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, sharey=True, figsize=(7,4))
    
    bp1 = ax1.boxplot(SZ_neg_data, notch=True, positions=[2,7,12,17,22,27], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp2 = ax1.boxplot(SZ_no_neg_data, notch=True, positions=[4,9,14,19,24,29], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','1','10','50'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,32)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia')
    ax1.legend([bp1['boxes'][0], bp2['boxes'][0]], ['With Negatives', 'Without Negatives'], loc='upper left', fontsize='x-small')
    

    bp3 = ax2.boxplot(CM_SZ_neg_data, notch=True, positions=[2,7,12,17,22,27], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp4 = ax2.boxplot(CM_SZ_no_neg_data, notch=True, positions=[4,9,14,19,24,29], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','1','10','50'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,32)
    ax2.set_xlabel('$\lambda$')
    ax2.set_title('Cell Motility-Schizophrenia')


    bp5 = ax3.boxplot(ASD_neg_data, notch=True, positions=[2,7,12,17,22,27], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp6 = ax3.boxplot(ASD_no_neg_data, notch=True, positions=[4,9,14,19,24,29], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','1','10','50'])
    ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax3.set_xlim(0,32)
    ax3.set_title('ASD')

    bp7 = ax4.boxplot(CM_ASD_neg_data, notch=True, positions=[2,7,12,17,22,27], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp7 = ax4.boxplot(CM_ASD_no_neg_data, notch=True, positions=[4,9,14,19,24,29], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'), labels=['0','0.01','0.1','1','10','50'])
    ax4.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax4.set_xlim(0,32)
    ax4.set_title('CM-ASD')


    plt.savefig('../outfiles/compare_neg_noneg.png')

    print('Created ../outfiles/compare_neg_noneg.png')

    return
    


def figure_3():

    one_file='../outfiles/SZ_1-layer_0.01-sinksource_auc.txt'
    SZ_first_one, ASD_first_one, CM_first_one = file_parser(one_file)

    one_file='../outfiles/SZ_1-layer_0.1-sinksource_auc.txt'
    SZ_other_one, ASD_other_one, CM_other_one = file_parser(one_file)

    two_file='../outfiles/SZ_2-layer_10-sinksource_auc.txt'
    SZ_two, ASD_two, CM_two = file_parser(two_file)

    three_file='../outfiles/SZ_3-layer_10-sinksource_auc.txt'
    SZ_three, ASD_three, CM_three = file_parser(three_file)

    #Compare 2 layer to 1 layer and 3 layer for each positive set
    #Create lists for each positive set of what data to compare [2 layer AUCs, 1 layer, 3 layer]

    SZ_sets = [SZ_two, SZ_first_one, SZ_three]
    ASD_sets = [ASD_two, ASD_other_one, ASD_three]
    CM_sets = [CM_two, CM_other_one, CM_three]

    Mann_test_lists = [SZ_sets, ASD_sets, CM_sets]

    with open('../outfiles/MWU_PVALUES_123layers.txt', 'w') as out:
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

    plt.savefig('../outfiles/comparing_layers.png')

    print('Created ../outfiles/comparing_layers.png')


    return


def figure_3_full():

    zero_file='../outfiles/SZ_1-layer_0-sinksource_auc.txt'
    SZ_zero, CM_SZ_zero = file_parser(zero_file)

    point_zero_one_file='../outfiles/SZ_1-layer_0.01-sinksource_auc.txt'
    SZ_pt_zo, CM_SZ_pt_zo = file_parser(point_zero_one_file)

    point_one_file='../outfiles/SZ_1-layer_0.1-sinksource_auc.txt'
    SZ_pt_o, CM_SZ_pt_o = file_parser(point_one_file)

    one_file='../outfiles/SZ_1-layer_1-sinksource_auc.txt'
    SZ_o, CM_SZ_o = file_parser(one_file)

    ten_file='../outfiles/SZ_1-layer_10-sinksource_auc.txt'
    SZ_ten, CM_SZ_ten = file_parser(ten_file)

    fifty_file='../outfiles/SZ_1-layer_50-sinksource_auc.txt'
    SZ_fifty, CM_SZ_fifty = file_parser(fifty_file)

    SZ_1_data = [SZ_zero, SZ_pt_zo, SZ_pt_o, SZ_o,SZ_ten, SZ_fifty]
    CM_SZ_1_data = [CM_SZ_zero, CM_SZ_pt_zo, CM_SZ_pt_o, CM_SZ_o,CM_SZ_ten, CM_SZ_fifty]


    #Get negative data
    zero_file='../outfiles/SZ_2-layer_0-sinksource_auc.txt'
    SZ_zero_2, CM_SZ_zero_2 = file_parser(zero_file)

    point_zero_one_file='../outfiles/SZ_2-layer_0.01-sinksource_auc.txt'
    SZ_pt_zo_2, CM_SZ_pt_zo_2 = file_parser(point_zero_one_file)

    point_one_file='../outfiles/SZ_2-layer_0.1-sinksource_auc.txt'
    SZ_pt_o_2, CM_SZ_pt_o_2 = file_parser(point_one_file)

    one_file='../outfiles/SZ_2-layer_1-sinksource_auc.txt'
    SZ_o_2, CM_SZ_o_2 = file_parser(one_file)

    ten_file='../outfiles/SZ_2-layer_10-sinksource_auc.txt'
    SZ_ten_2, CM_SZ_ten_2 = file_parser(ten_file)

    fifty_file='../outfiles/SZ_2-layer_50-sinksource_auc.txt'
    SZ_fifty_2, CM_SZ_fifty_2 = file_parser(fifty_file)

    SZ_2_data = [SZ_zero_2, SZ_pt_zo_2, SZ_pt_o_2, SZ_o_2, SZ_ten_2, SZ_fifty_2]
    CM_SZ_2_data = [CM_SZ_zero_2, CM_SZ_pt_zo_2, CM_SZ_pt_o_2,CM_SZ_o_2, CM_SZ_ten_2, CM_SZ_fifty_2]


    zero_file='../outfiles/SZ_3-layer_0-sinksource_auc.txt'
    SZ_zero_3, CM_SZ_zero_3 = file_parser(zero_file)

    point_zero_one_file='../outfiles/SZ_3-layer_0.01-sinksource_auc.txt'
    SZ_pt_zo_3, CM_SZ_pt_zo_3 = file_parser(point_zero_one_file)

    point_one_file='../outfiles/SZ_3-layer_0.1-sinksource_auc.txt'
    SZ_pt_o_3, CM_SZ_pt_o_3 = file_parser(point_one_file)

    one_file='../outfiles/SZ_3-layer_1-sinksource_auc.txt'
    SZ_o_3, CM_SZ_o_3 = file_parser(one_file)

    ten_file='../outfiles/SZ_3-layer_10-sinksource_auc.txt'
    SZ_ten_3, CM_SZ_ten_3 = file_parser(ten_file)

    fifty_file='../outfiles/SZ_3-layer_50-sinksource_auc.txt'
    SZ_fifty_3, CM_SZ_fifty_3 = file_parser(fifty_file)

    SZ_3_data = [SZ_zero_3, SZ_pt_zo_3, SZ_pt_o_3, SZ_o_3, SZ_ten_3, SZ_fifty_3]
    CM_SZ_3_data = [CM_SZ_zero_3, CM_SZ_pt_zo_3, CM_SZ_pt_o_3,CM_SZ_o_3, CM_SZ_ten_3, CM_SZ_fifty_3]

    #Compare 2 layer to 1 layer and 3 layer for each positive set
    #Create lists for each positive set of what data to compare [2 layer AUCs, 1 layer, 3 layer]


    one_list = [SZ_1_data, CM_SZ_1_data]
    two_list = [SZ_2_data, CM_SZ_2_data]
    three_list = [SZ_3_data, CM_SZ_3_data]







    # Mann_test_lists = [SZ_sets, ASD_sets, CM_sets]

    # with open('../outfiles/MWU_PVALUES_123layers.txt', 'w') as out:
    #     for i in range(len(Mann_test_lists)):
    #         if i == 0:
    #             out.write('Schizophrenia\t')
    #         elif i == 1:
    #             out.write('Autism\t')
    #         else:
    #             out.write('Cell Motility\t')
    #         out.write('\n')
    #         U1, p1 = stats.mannwhitneyu(Mann_test_lists[i][0], Mann_test_lists[i][1], alternative='two-sided') 
    #         out.write('2 layer vs. 1 layer\t')
    #         out.write(str(p1) + '\n')
    #         U3, p3 = stats.mannwhitneyu(Mann_test_lists[i][0], Mann_test_lists[i][2], alternative='two-sided')
    #         out.write('2 layer vs. 3 layer\t')
    #         out.write(str(p3) + '\n') 
    #         out.write('\n')

    with open('../outfiles/Layers_statistical_tests.txt', 'w') as out:
        #first iterate through each constant, then iterate through each positive set for that constant to get the list of AUCs
        out.write('T-Tests\n')

        for i in range(len(one_list)): #the two data set lists are always the same length 
            if i == 0:
                out.write('Schizophrenia\n\n')
            if i == 1:
                out.write('Cell Motility-SZ\n\n')
            # if j == 2:
            #     out.write('ASD\t')
            # if j == 3:
            #     out.write('Cell Motility-ASD\t')
            for j in range(len(one_list[i])):
                if j == 0:
                    out.write('Lambda=0\n')
                elif j == 1:
                    out.write('Lambda=0.01\n')
                elif j == 2:
                    out.write('Lambda=0.1\n')
                elif j == 3:
                    out.write('Lambda=1\n')
                elif j == 4:
                    out.write('Lambda=10\n')
                else:
                    out.write('Lambda=50\n')
                    # U, p_value = stats.mannwhitneyu(neg_lists[i][j], no_neg_lists[i][j], alternative='two-sided')

                t_value,p_value = stats.ttest_ind(one_list[i][j], two_list[i][j])
                out.write('1 vs 2: '+str(p_value)+'\n')

                t_value,p_value = stats.ttest_ind(one_list[i][j], three_list[i][j])
                out.write('1 vs 3: '+str(p_value)+'\n')

                t_value,p_value = stats.ttest_ind(two_list[i][j], three_list[i][j])
                out.write('2 vs 3: '+str(p_value)+'\n\n')

            out.write('\n')

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




    SZ_data = [SZ_1_data, SZ_2_data]
    CM_SZ_data = [CM_SZ_1_data,CM_SZ_2_data]
    # ASD_data = [ASD_other_one, ASD_two, ASD_three]
    

    fig3, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, sharey=True, figsize=(7,4))

    bp1 = ax1.boxplot(SZ_1_data, notch=True, positions=[2,9,16,23,30,37], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp2 = ax1.boxplot(SZ_2_data, notch=True, positions=[4,11,18,25,32,39], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'))
    bp3 = ax1.boxplot(SZ_3_data, notch=True, positions=[6,13,20,27,34,41], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#eeb5ff'), labels=['0','0.01','0.1','1','10','50'])
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_xlim(0,43)
    ax1.set_ylabel('AUC')
    ax1.set_title('Schizophrenia')
    ax1.legend([bp1['boxes'][0], bp2['boxes'][0], bp3['boxes'][0]], ['1 layer', '2 layer', '3 layer'], loc='upper left', fontsize='x-small')
    

    bp4 = ax2.boxplot(CM_SZ_1_data, notch=True, positions=[2,9,16,23,30,37], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'))
    bp5 = ax2.boxplot(CM_SZ_2_data, notch=True, positions=[4,11,18,25,32,39], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#b0f2c2'))
    bp6 = ax2.boxplot(CM_SZ_3_data, notch=True, positions=[6,13,20,27,34,41], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#eeb5ff'), labels=['0','0.01','0.1','1','10','50'])
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax2.set_xlim(0,43)
    ax2.set_xlabel('$\lambda$')
    ax2.set_title('Cell Motility-Schizophrenia')


    # bp5 = ax3.boxplot(CM_data, notch=True, positions=[2,4,6], widths=1.5, sym='k+', patch_artist=True, boxprops=dict(facecolor='#8193ef'), labels=['1','2','3'])
    # ax3.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    # ax3.set_xlim(0,8)
    # ax3.set_title('Cell Motility')

    plt.savefig('../outfiles/comparing_layers.png')

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



























