import pandas as pd
from subprocess import Popen, PIPE
from random import randint, choice
import numpy as np

def BPA_Aptamer_design(n_pool,n_gen,n_candidates,desired_length=False):
    print('Starting program')
    print('-------------------------------')
    aptamer1_seq='CGGTGGGTGGTCAGGTGGGATA'.lower()
    aptamer1_str,_=Structure_Aptamer('CGGTGGGTGGTCAGGTGGGATAGCGTTCCGCGTATGGCCCAGCG')
    aptamer2_seq='GGATA'.lower()
    aptamer2_str,_=Structure_Aptamer('GGATAGCGGGTTCC')
    HRP_DNAzyme='GTGGGGCATTGTGGGTGGGTGTGG'
    pool=initial_pool_gen(n_pool,desired_length)
    for generation in range(n_gen):
        print('Current generation: '+str(generation+1))
        candidates=pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
        pool=new_pool_gen(candidates,n_pool,n_candidates)
    final_candidates=pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
    final_candidates.to_excel('Final_candidates.xlsx')
    print('------------------------------------')
    print('Succesfully executed!')
    
    

def MFE_Hybridization(seq_aptamer,seq_target):
    '''Return the MFE (kcal/mol) of the hybridization of two given RNA sequences'''
    MFE_calc=Popen('C:\Program Files (x86)\ViennaRNA Package\RNAcofold.exe', stdin=PIPE, stdout=PIPE,shell=True)
    Input='>Seq_toehold\n'+seq_aptamer+'\n>Seq_target\n'+seq_target
    Result=MFE_calc.communicate(Input.encode())
    return float(Result[0][-9:-3])


def Structure_Aptamer(seq):
    '''Return the MFE (kcal/mol) of the RNA secondary strucure of a given sequence'''
    MFE_calc=Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE,shell=True)
    Result=MFE_calc.communicate(seq.encode())
    Structure=Result[0][len(seq):-10]
    return Structure.decode('UTF-8'),float(Result[0][-9:-3])



def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])
        
def initial_pool_gen(n_pool,desired_length):
    seqs=[]
    for n in range(n_pool):
        seq=''
        if desired_length!=False:
            for m in range(desired_length):
                seq=seq+choice('atgc')
        else:
            length=randint(12,100)
            for m in range(length):
                seq=seq+choice('atgc')
        seqs.append(seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def new_pool_gen(candidates,n_pool,n_candidates):
    div=int((n_pool/n_candidates)-1)
    seqs=[]
    for seq in list(candidates.sequences):
        seqs.append(seq)
        for n in range(div):
            seqs.append(mutation(seq))
    return pd.DataFrame(seqs,columns=['sequences'])

def mutation(seq):
    seq_list=list(seq)
    mutation_type=randint(0,2)
    n_to_mutate=randint(0,len(seq)-1)
    if mutation_type==0:
        del seq_list[n_to_mutate]
    elif mutation_type==1:
        seq_list.insert(n_to_mutate,choice('atgc'))
    else:
        if seq_list[n_to_mutate]=='a':
            seq_list[n_to_mutate]=choice('tgc')
        elif seq_list[n_to_mutate]=='t':
            seq_list[n_to_mutate]=choice('agc')
        elif seq_list[n_to_mutate]=='g':
            seq_list[n_to_mutate]=choice('atc')
        else:
            seq_list[n_to_mutate]=choice('atg')
    return ''.join(seq_list)
    
def pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme):
    EditDis_apt1=[]
    EditDis_apt2=[]
    EditDis_str1=[]
    EditDis_str2=[]
    dMFE=[]
    for ind in list(pool.index.values):
        EditDis_apt1.append(editDistance(pool.loc[ind,'sequences'],aptamer1_seq))
        EditDis_apt2.append(editDistance(pool.loc[ind,'sequences'],aptamer2_seq))
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        EditDis_str1.append(editDistance(candidate_str,aptamer1_str))
        EditDis_str2.append(editDistance(candidate_str,aptamer2_str))
        Hybrid_MFE=MFE_Hybridization(pool.loc[ind,'sequences'],HRP_DNAzyme)
        if Hybrid_MFE<0:
            if Hybrid_MFE<candidate_str_MFE:
                dMFE.append(abs(abs(Hybrid_MFE)-abs(candidate_str_MFE)))
            else:
                dMFE.append(0.0)
        else:
            dMFE.append(0.0)
    pool['edit_apt1']=EditDis_apt1
    pool['edit_apt2']=EditDis_apt2
    pool['edit_apt1_str']=EditDis_str1
    pool['edit_apt2_str']=EditDis_str2
    pool['dMFE']=dMFE
    pool['edit_apt1_norm']=(1-(pool['edit_apt1']/max(pool['edit_apt1'])))
    pool['edit_apt2_norm']=(1-(pool['edit_apt2']/max(pool['edit_apt2'])))
    pool['edit_apt1_str_norm']=(1-(pool['edit_apt1_str']/max(pool['edit_apt1_str'])))
    pool['edit_apt2_str_norm']=(1-(pool['edit_apt2_str']/max(pool['edit_apt2_str'])))
    if max(pool['dMFE'])==0.0:
        pool['dMFE_norm']=pool['dMFE']
    else:
        pool['dMFE_norm']=pool['dMFE']/max(pool['dMFE'])
    pool['score']=np.average([pool['edit_apt1_norm'],pool['edit_apt2_norm'],pool['edit_apt1_str_norm'],pool['edit_apt2_str_norm'],pool['dMFE_norm']])
    return pool.nlargest(n_candidates,['score'])

BPA_Aptamer_design(100,1000,10)                    
        
    
                
        
