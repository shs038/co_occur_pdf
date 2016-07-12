import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
%matplotlib inline
def co_occur(n):
    '''
    input is a matix contains all data of the genomic position of motif in chr1.
    output is a matix contains the frequency of co-occurance of every two motifs. 
    '''
    (row,col)=n.shape #find the dimension of the data
    count_matrix=np.zeros((col,col),dtype=np.int) #a matrix to store future data
    col_vector=np.zeros((row,1),dtype=np.int) #a column vector use in loop
    for i in range (col):#use the column vector to store each column one by one
        col_vector=n[:,i]!= -1# Find if the motifs exist
        for j in range (col): 
            log_v=n[:,j]!= -1#Find if the motifs exist
            vector_sum = 1*col_vector + 1*log_v
            log_input = vector_sum == 2
            input=np.sum(1*log_input)#convert logical ou
            count_matrix[i,j]=count_matrix[i,j]+input
            np.fill_diagonal(count_matrix, 0)# change diagonal back to zero
    return count_matrix
# read in data
df = pd.read_excel('/Users/sun/Desktop/Lab/chr1new.xlsx',parse_cols="B:GO")#load part of file contain data
data_list=df.as_matrix()#conver data into matrix form
# find co-occurence
co=co_occur(data_list)
motifs = df.columns.values
def Find_pdf(n):
    '''
    For all pairs of motifs - is there a pair that co-occurs more or less often than you would expect.
    input: a matirc contains frequency of co-occurence of every pair of motifs.
    output:an array contains z score of each pair of motifs.
    '''
    p_matrix = np.zeros((n.shape[0],n.shape[1]-1),dtype=np.float)
    for i in range (n.shape[0]):
        #find meand and stv for each motif
        co_motif = co[i,:]
        co_motif = np.delete(co_motif,i)#remove data of the motif co-occur with itself
        mean=np.mean(co_motif)
        std=np.std(co_motif)
        pdf=stats.norm.pdf(co_motif,loc=mean, scale=std)# find z socre 
        p_matrix[i,:]=pdf
    return p_matrix
p_value=Find_pdf(co)
p_value
#assign co-occurance z socre of one motif with itself as 100 to make dataframe 
p_for_dataframe=np.zeros((p_value.shape[0],p_value.shape[0]),dtype=np.float)
for i in range (p_value.shape[0]):
        p_motif_self=p_value[i,:]
        p_motif_self=np.insert(p_motif_self,i,100)
        p_for_dataframe[i,:]=p_motif_self
print p_for_dataframe
p_frame = pd.DataFrame(p_for_dataframe, columns=motifs)
p_frame.index = motifs
p_frame
def Find_sigpair(siglevel):
    threshold = siglevel # threshold of z score
    significantPairCount = 0
    sigpair=[]
    for i in range(co.shape[0]):
        for j in range(co.shape[0]):
            current_p = p_for_dataframe[i][j]
            if current_p <= threshold:
                # display motif pair
                motif1 = motifs[i]
                motif2 = motifs[j]
                sigpair.append((motif1,motif2))
                significantPairCount += 1
    print significantPairCount
    return sigpair
#pairs with p<0.025/195
sigpair=Find_sigpair(0.025/195)
sigpair_frame = pd.DataFrame(sigpair,columns=['motif1','motif2'])
#remove repeating sigpair
for i in range (len(sigpair_frame)-1):
    for j in range (i+1,len(sigpair_frame)):
        if sigpair_frame.loc[i]['motif1']==sigpair_frame.loc[j]['motif2'] and sigpair_frame.loc[i]['motif2']==sigpair_frame.loc[j]['motif1']:
            sigpair_frame.loc[j]['motif1']=0
sigpair_frame=sigpair_frame[sigpair_frame['motif1']!=0]
sigpair_frame = sigpair_frame.reset_index(drop=True)
sigpair_frame.to_csv('/Users/sun/Desktop/sigpair.tsv', sep='\t')
#pairs with p<0.05/195
sigpair_more=Find_sigpair(0.05/195)
sigpair_more_frame = pd.DataFrame(sigpair_more,columns=['motif1','motif2'])
#remove repeating sigpair
for i in range (len(sigpair_more_frame)-1):
    for j in range (i+1,len(sigpair_more_frame)):
        if sigpair_more_frame.loc[i]['motif1']==sigpair_more_frame.loc[j]['motif2'] and sigpair_more_frame.loc[i]['motif2']==sigpair_more_frame.loc[j]['motif1']:
            sigpair_more_frame.loc[j]['motif1']=0
sigpair_more_frame=sigpair_more_frame[sigpair_more_frame['motif1']!=0]
sigpair_more_frame = sigpair_more_frame.reset_index(drop=True)
sigpair_more_frame.to_csv('/Users/sun/Desktop/sigpair_more.tsv', sep='\t')
