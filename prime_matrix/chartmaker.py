import matplotlib.pyplot as plt
import numpy as np

#SET USE_SD=False to get Error is IQR, not standard deviation
USE_SD=True

def main():
    figure_1()
    figure_2()
    figure_3()
    return

def figure_1():

    

    #SET SD=False to get Error is IQR, not standard deviation
    #SD=False
    
    zero_file='outfiles/SZ_1-layer_0.150_auc.txt'
    zero=file_parser(zero_file)

    point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    point_zero_one=file_parser(point_zero_one_file)

    point_one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    point_one=file_parser(point_one_file)

    ten_file='outfiles/SZ_1-layer_10-sinksource_0.150_auc.txt'
    ten=file_parser(ten_file)

    fifty_file='outfiles/SZ_1-layer_50-sinksource_0.150_auc.txt'
    fifty=file_parser(fifty_file)

    values=[zero,point_zero_one,point_one, ten, fifty]
    
    sz_means=[]
    sz_plus_err=[]
    sz_minus_err=[]
    
    asd_means=[]
    asd_plus_err=[]
    asd_minus_err=[]

    cm_means=[]
    cm_plus_err=[]
    cm_minus_err=[]


    for test in values:
        sz_means.append(test[0])
        sz_plus_err.append(test[3][0])
        sz_minus_err.append(test[3][1])

        asd_means.append(test[1])
        asd_plus_err.append(test[4][0])
        asd_minus_err.append(test[4][1])

        cm_means.append(test[2])
        cm_plus_err.append(test[5][0])
        cm_minus_err.append(test[5][1])

    sz_err=[sz_plus_err]+[sz_minus_err]
    asd_err=[asd_plus_err]+[asd_minus_err]
    cm_err=[cm_plus_err]+[cm_minus_err]



    



    fig1 = plt.figure(figsize=(7,4))
    ax1 = plt.subplot(1,3,1)
    ax2 = plt.subplot(1,3,2)
    ax3 = plt.subplot(1,3,3)
    #fig1, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1)

    x_axis=np.arange(len(values))
    ax1.bar(x_axis,sz_means, tick_label=['0','0.01','0.1','10','50'], yerr=sz_err, capsize=7, color='b', alpha=0.8)
    ax1.set_title('SZ')
    ax1.set_ylabel('AUC')
    ax1.set_ylim(0.5,0.85)
    ax1.tick_params(labelsize=9)



    ax2.bar(x_axis,asd_means,tick_label=['0','0.01','0.1','10','50'], yerr=asd_err, capsize=7,color='b', alpha=0.8) 
    ax2.set_title('ASD')
    ax2.set_xlabel('$\lambda$') #Try to see if there is a way to put the set name above the chart
    ax2.set_ylim(0.5,0.85)
    ax2.tick_params(labelsize=9)

    ax3.bar(x_axis,cm_means,tick_label=['0','0.01','0.1','10','50'], yerr=cm_err, capsize=7,color='b', alpha=0.8) 
    ax3.set_title('Cell Motility')
    ax3.set_ylim(0.5,0.85)
    ax3.tick_params(labelsize=9)






    fig1.tight_layout()
    plt.savefig('Figure_1.png')
    #plt.show()
    return


def figure_2():

    #Error is IQR, not standard deviation

    
    zero_file='outfiles/SZ_1-layer_0.150_auc.txt'
    zero=file_parser(zero_file)

    point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    point_zero_one=file_parser(point_zero_one_file)

    point_one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    point_one=file_parser(point_one_file)

    ten_file='outfiles/SZ_1-layer_10-sinksource_0.150_auc.txt'
    ten=file_parser(ten_file)

    fifty_file='outfiles/SZ_1-layer_50-sinksource_0.150_auc.txt'
    fifty=file_parser(fifty_file)


    no_neg_zero_file='outfiles/SZ_1-layer_no_neg_0.150_auc.txt'
    no_neg_zero=file_parser(no_neg_zero_file)




    no_neg_point_zero_one_file='outfiles/SZ_1-layer_0.01-sinksource_no_neg_0.150_auc.txt'
    no_neg_point_zero_one=file_parser(no_neg_point_zero_one_file)

    no_neg_point_one_file='outfiles/SZ_1-layer_0.1-sinksource_no_neg_0.150_auc.txt'
    no_neg_point_one=file_parser(no_neg_point_one_file)

    no_neg_ten_file='outfiles/SZ_1-layer_10-sinksource_no_neg_0.150_auc.txt'
    no_neg_ten=file_parser(no_neg_ten_file)

    no_neg_fifty_file='outfiles/SZ_1-layer_50-sinksource_no_neg_0.150_auc.txt'
    no_neg_fifty=file_parser(no_neg_fifty_file)

    # values=[zero,point_zero_one, ten, fifty]
    # no_neg_values=[no_neg_zero,no_neg_point_zero_one,no_neg_ten,no_neg_fifty]

    values=[zero, point_zero_one,point_one,ten, fifty]
    no_neg_values=[no_neg_zero,no_neg_point_zero_one,no_neg_point_one, no_neg_ten,no_neg_fifty]

    SZ_neg_means=[]
    SZ_neg_plus_err=[]
    SZ_neg_minus_err=[]

    SZ_no_neg_means=[]
    SZ_no_neg_plus_err=[]
    SZ_no_neg_minus_err=[]

    ASD_neg_means=[]
    ASD_neg_plus_err=[]
    ASD_neg_minus_err=[]

    ASD_no_neg_means=[]
    ASD_no_neg_plus_err=[]
    ASD_no_neg_minus_err=[]

    CM_neg_means=[]
    CM_neg_plus_err=[]
    CM_neg_minus_err=[]

    CM_no_neg_means=[]
    CM_no_neg_plus_err=[]
    CM_no_neg_minus_err=[]
    
    x_axis=np.arange(len(values))

    for test in values:
        SZ_neg_means.append(test[0])
        SZ_neg_plus_err.append(test[3][0])
        SZ_neg_minus_err.append(test[3][1])

        ASD_neg_means.append(test[1])
        ASD_neg_plus_err.append(test[4][0])
        ASD_neg_minus_err.append(test[4][1])

        CM_neg_means.append(test[2])
        CM_neg_plus_err.append(test[5][0])
        CM_neg_minus_err.append(test[5][1])

    for test in no_neg_values:
        SZ_no_neg_means.append(test[0])
        SZ_no_neg_plus_err.append(test[3][0])
        SZ_no_neg_minus_err.append(test[3][1])

        ASD_no_neg_means.append(test[1])
        ASD_no_neg_plus_err.append(test[4][0])
        ASD_no_neg_minus_err.append(test[4][1])

        CM_no_neg_means.append(test[2])
        CM_no_neg_plus_err.append(test[5][0])
        CM_no_neg_minus_err.append(test[5][1])


    SZ_neg_err=[SZ_neg_plus_err]+[SZ_neg_minus_err]
    SZ_no_neg_err=[SZ_no_neg_plus_err]+[SZ_no_neg_minus_err]

    ASD_neg_err=[ASD_neg_plus_err]+[ASD_neg_minus_err]
    ASD_no_neg_err=[ASD_no_neg_plus_err]+[ASD_no_neg_minus_err]

    CM_neg_err=[CM_neg_plus_err]+[CM_neg_minus_err]
    CM_no_neg_err=[CM_no_neg_plus_err]+[CM_no_neg_minus_err]




    


    fig2 = plt.figure(figsize=(7,4))
    ax1 = plt.subplot(1,3,1)
    ax2 = plt.subplot(1,3,2)
    ax3 = plt.subplot(1,3,3)
    #fig2, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1)
    bar_width = 0.3
    error_config = {'ecolor': '0.3'}
    x_axis=np.arange(len(values))

    # rects1=ax1.bar(x_axis,neg_means, bar_width, tick_label=['0','0.1','10','50'], yerr=neg_err, capsize=7)
    # rects2=ax1.bar(x_axis,no_neg_means,bar_width, tick_label=['0','0.1','10','50'], yerr=no_neg_err, capsize=7, color='r')
    rects1=ax1.bar(x_axis,SZ_neg_means, bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=SZ_neg_err, capsize=2,\
        error_kw=error_config, label='Pseudo-SS+', color='b', alpha=0.8)
    rects2=ax1.bar(x_axis+bar_width,SZ_no_neg_means,bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=SZ_no_neg_err, capsize=7, \
            color='#67CB8A',edgecolor='k',error_kw=error_config, label='SS+', alpha=0.8)
    ax1.set_ylabel('AUC')
    ax1.set_xticks(x_axis + bar_width / 2)
    ax1.set_ylim(0.5,0.85)
    ax1.legend()
    ax1.set_title('SZ')



    rects1=ax2.bar(x_axis,ASD_neg_means, bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=ASD_neg_err, capsize=2,error_kw=error_config, label='with Negatives', color='b', alpha=0.8)
    rects2=ax2.bar(x_axis+bar_width,ASD_no_neg_means,bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=ASD_no_neg_err, capsize=7, color='#67CB8A',edgecolor='k',error_kw=error_config, label='without Negatives', alpha=0.8)
    ax2.set_xlabel('$\lambda$')
    ax2.set_xticks(x_axis + bar_width / 2)
    ax2.set_ylim(0.5,0.85)
    ax2.set_title('ASD')


    rects1=ax3.bar(x_axis,CM_neg_means, bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=CM_neg_err, capsize=2,error_kw=error_config, label='Pseudo-SS+', color='b', alpha=0.8)
    rects2=ax3.bar(x_axis+bar_width,CM_no_neg_means,bar_width, tick_label=['0','0.01','0.1','10','50'], yerr=CM_no_neg_err, capsize=7, color='#67CB8A',edgecolor='k',error_kw=error_config, label='without Negatives', alpha=0.8)
    ax3.set_xticks(x_axis + bar_width / 2)
    ax3.set_ylim(0.5,0.85)
    ax3.set_title('Cell Motility')



    fig2.tight_layout()
    plt.savefig('Figure_2.png')
    #plt.show()


def figure_3():

    #Error is IQR, not standard deviation
    
    

    one_file='outfiles/SZ_1-layer_0.01-sinksource_0.150_auc.txt'
    first_one=file_parser(one_file)

    one_file='outfiles/SZ_1-layer_0.1-sinksource_0.150_auc.txt'
    other_one=file_parser(one_file)

    two_file='outfiles/SZ_2-layer_10-sinksource_0.150_auc.txt'
    two=file_parser(two_file)

    three_file='outfiles/SZ_3-layer_10-sinksource_0.150_auc.txt'
    three=file_parser(three_file)


    one=[first_one[0],other_one[1], other_one[2], first_one[3], other_one[4], other_one[5]]



    

    values=[one,two,three]
    
    sz_means=[]
    sz_plus_err=[]
    sz_minus_err=[]
    
    asd_means=[]
    asd_plus_err=[]
    asd_minus_err=[]

    cm_means=[]
    cm_plus_err=[]
    cm_minus_err=[]


    for test in values:
        sz_means.append(test[0])
        sz_plus_err.append(test[3][0])
        sz_minus_err.append(test[3][1])

        asd_means.append(test[1])
        asd_plus_err.append(test[4][0])
        asd_minus_err.append(test[4][1])

        cm_means.append(test[2])
        cm_plus_err.append(test[5][0])
        cm_minus_err.append(test[5][1])

    sz_err=[sz_plus_err]+[sz_minus_err]
    asd_err=[asd_plus_err]+[asd_minus_err]
    cm_err=[cm_plus_err]+[cm_minus_err]

    



    



    fig3 = plt.figure(figsize=(7,4))
    ax1 = plt.subplot(1,3,1)
    ax2 = plt.subplot(1,3,2)
    ax3 = plt.subplot(1,3,3)
    #fig3, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1)

    x_axis=np.arange(len(values))
    ax1.bar(x_axis,sz_means, tick_label=['1','2','3'], yerr=sz_err, capsize=7,color='b', alpha=0.8)
    ax1.set_title('SZ')
    ax1.set_ylabel('AUC')
    ax1.set_ylim(0.5,0.85)



    ax2.bar(x_axis,asd_means,tick_label=['1','2','3'], yerr=asd_err, capsize=7,color='b', alpha=0.8) 
    ax2.set_title('ASD')
    ax2.set_xlabel('Layers') #Try to see if there is a way to put the set name above the chart
    ax2.set_ylim(0.5,0.85)

    ax3.bar(x_axis,cm_means,tick_label=['1','2','3'], yerr=cm_err, capsize=7,color='b', alpha=0.8) 
    ax3.set_title('Cell Motility')
    ax3.set_ylim(0.5,0.85)




    print('SZ means:',sz_means)
    print('ASD means:',asd_means)
    print('motility means:',cm_means)

    plt.tight_layout()
    plt.savefig('Figure_3.png')
    #plt.show()
    return
    



def file_parser(auc_file):
    file=open(auc_file,'r')
    x=0
    SZ=[]
    ASD=[]
    CM=[]
    for line in file:
        if x>0:
            line=line.strip('\n').split('\t')
            SZ.append(float(line[0].strip('\'')))

            ASD.append(float(line[1].strip('\'')))
            CM.append(float(line[2].strip('\'')))
        x+=1

    if USE_SD:
        return [np.mean(SZ),np.mean(ASD),np.mean(CM),SD(SZ),SD(ASD),SD(CM)]
    else:    
        return [np.mean(SZ),np.mean(ASD),np.mean(CM),IQR(SZ),IQR(ASD),IQR(CM)]

def SD(dist):
    return [np.std(dist),np.std(dist)]

def IQR(dist):
    return [np.percentile(dist, 75) - np.mean(dist), np.mean(dist)-np.percentile(dist, 25)]


main()