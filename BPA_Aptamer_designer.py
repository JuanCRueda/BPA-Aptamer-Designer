import pandas as pd
from subprocess import Popen, PIPE
from random import randint, choice
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def BPA_Aptamer_design(n_pool,n_gen,n_candidates,desired_length=False,starting_sequence=False,visualize=True):
    print('Starting program')
    print('-------------------------------')
    aptamer1_seq='CGGTGGGTGGTCAGGTGGGATA'.lower()
    aptamer1_str,_=Structure_Aptamer('CGGTGGGTGGTCAGGTGGGATAGCGTTCCGCGTATGGCCCAGCG')
    aptamer2_seq='GGATA'.lower()
    aptamer2_str,_=Structure_Aptamer('GGATAGCGGGTTCC')
    HRP_DNAzyme='GTGGGGCATTGTGGGTGGGTGTGG'.lower()
    if starting_sequence==False:
        pool=initial_pool_gen(n_pool,desired_length)
    else:
        pool=initial_pool_with_seq(n_pool,starting_sequence,desired_length)
    candidates=pd.DataFrame({'sequences':['','','','','','','','','',''],'score':[0,0,0,0,0,0,0,0,0,0]})
    if visualize:
        generations=[]
        scores=[]
        line=[]
        line2=[]
    for generation in range(n_gen):
        print('Current generation: '+str(generation+1))
        candidates_past=candidates
        candidates=pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
        if visualize:
            generations.append(generation+1)
            scores.append(np.max(candidates.score))
            line,line2=plotter(generations,scores,line,line2,pause_time=0.1)
        if np.average(candidates.score)>0.7:
            print('Average score treshold reached')
            break
        if candidate_similarities(candidates,candidates_past)>0.7:
            pool=new_pool_explosion(candidates,n_pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
            candidates=pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
            pool=new_pool_gen(candidates,n_pool,n_candidates)
        else:
            pool=new_pool_gen(candidates,n_pool,n_candidates)
    final_candidates=pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme)
    final_candidates.to_excel('Final_candidates.xlsx')
    print('------------------------------------')
    print('Succesfully executed!')

def plotter(generations,scores,line,line2,pause_time=30):
    if line==[]:
        plt.ion()
        fig=plt.figure(figsize=[6.4,6.0])
        ax1=fig.add_subplot(2,1,1)
        ax2=fig.add_subplot(2,1,2)
        line,=ax1.plot(generations,scores)
        line2,=ax2.plot(generations[-10:],scores[-10:])
        ax1.set_ylabel('Max. score')
        ax2.set_ylabel('Max. score')
        ax2.set_xlabel('Generations')
        plt.show()
    line.set_data(generations,scores)
    line2.set_data(generations,scores)
    fig_current=plt.gcf()
    ax_list=fig_current.get_axes()
    if np.min(scores)<=line.axes.get_ylim()[0] or np.max(scores)>=line.axes.get_ylim()[1]:
        ax_list[0].set_ylim([np.min(scores)-0.05,np.max(scores)+0.05])
    if np.max(generations)>=line.axes.get_xlim()[1]:
        ax_list[0].set_xlim([0,np.max(generations)+0.05*np.max(generations)])
    ax_list[1].set_xlim([np.min(generations[-10:])-0.1,np.max(generations[-10:])+0.1])
    ax_list[1].set_ylim([np.min(scores[-10:])-0.001,np.max(scores[-10:])+0.001])
    plt.pause(pause_time)
    return line,line2
    
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
    #try:
    MFE=float(Result[0][-9:-3])
    #except ValueError:
        #print(Result)
        #MFE=float(Result[0][-8:-3])
    return Structure.decode('UTF-8'),MFE



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

def initial_pool_with_seq(n_pool,seq,desired_length):
    seq=seq.lower()
    seqs=[seq]
    for n in range(n_pool-1):
        new_seq=mutation(seq)
        while new_seq in seqs:
            new_seq=mutation(new_seq)
        seqs.append(new_seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def new_pool_gen(candidates,n_pool,n_candidates):
    div=int((n_pool/n_candidates)-1)
    seqs=[]
    for seq in list(candidates.sequences):
        if seq not in seqs:
            seqs.append(seq)
            div_s=div
        else:
            div_s=div+1
        for n in range(div_s):
            new_seq=mutation(seq)
            while len(new_seq)<10:
                new_seq=n_addition_mut(new_seq)
            while new_seq in seqs:
                new_seq=mutation(new_seq)
                while len(new_seq)<10:
                    new_seq=n_addition_mut(new_seq)
            seqs.append(new_seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def new_pool_explosion(candidates,n_pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme):
    print('------------------------------------------')
    print('Starting hyperdiverse subgenerations to liberate plateau')
    div=int((n_pool/n_candidates)-1)
    seqs=list(candidates.sequences)
    seqs_len=[len(x) for x in seqs]
    #n_new_gens=int(min(seqs_len)/2)
    #print('Number of hyperdiverse subgenerations to generate: '+str(n_new_gens))
    print('------------------------------------------')
    for n in range(3):
        seqs_2=[]
        print('Current subgeneration: '+str(n+1))
        for seq in seqs:
            seqs_2.append(seq)
            for n in range(div):
                new_seq=mutation(seq)
                while len(new_seq)<10:
                    new_seq=n_addition_mut(new_seq)
                while (new_seq in seqs) or (new_seq in seqs_2):
                    new_seq=mutation(new_seq)
                    while len(new_seq)<10:
                        new_seq=n_addition_mut(new_seq)
                seqs_2.append(new_seq)
        seqs=seqs_2[:]
    print('------------------------------------------')
    print('End of hyperdiverse period')
    print('------------------------------------------')
    big_gen=pd.DataFrame(seqs,columns=['sequences'])
    print('------------------------------------------')
    print('Selecting subset from hyperdiverse period')
    print('------------------------------------------')
    EditDis_apt1=[]
    EditDis_apt2=[]
    dMFE=[]
    for ind in list(big_gen.index.values):
        EditDis_apt1.append(editDistance(big_gen.loc[ind,'sequences'],aptamer1_seq))
        EditDis_apt2.append(editDistance(big_gen.loc[ind,'sequences'],aptamer2_seq))
        Hybrid_MFE=MFE_Hybridization(big_gen.loc[ind,'sequences'],HRP_DNAzyme)
        if Hybrid_MFE<0:
            dMFE.append(abs(Hybrid_MFE))
        else:
            dMFE.append(0.0)

    big_gen['edit_apt1']=EditDis_apt1
    big_gen['edit_apt2']=EditDis_apt2
    big_gen['dMFE']=dMFE
    big_gen['edit_apt1_norm']=1/(big_gen['edit_apt1']+1)
    big_gen['edit_apt2_norm']=1/(big_gen['edit_apt2']+1)
    if max(big_gen['dMFE'])==0.0:
        big_gen['dMFE_norm']=big_gen['dMFE']
    else:
        big_gen['dMFE_norm']=big_gen['dMFE']/max(big_gen['dMFE'])
    big_gen['score']=np.average(big_gen.loc[:,'edit_apt1_norm':'dMFE_norm'])
    print('------------------------------------------')
    print('Subset created, return to main algorithm')
    print('------------------------------------------')
    return big_gen.nlargest(n_pool,['score'])


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

def n_addition_mut(seq):
    seq_list=list(seq)
    n_to_mutate=randint(0,len(seq)-1)
    seq_list.insert(n_to_mutate,choice('atgc'))
    return ''.join(seq_list)

def pool_evaluation(pool,n_candidates,aptamer1_seq,aptamer1_str,aptamer2_seq,aptamer2_str,HRP_DNAzyme):
    EditDis_apt1=[]
    EditDis_apt2=[]
    EditDis_str1=[]
    EditDis_str2=[]
    dMFE=[]
    for ind in list(pool.index.values):
        EditDis_apt1.append(editDistance(aptamer1_seq,pool.loc[ind,'sequences']))
        EditDis_apt2.append(editDistance(aptamer2_seq,pool.loc[ind,'sequences']))
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        EditDis_str1.append(editDistance(aptamer1_str,candidate_str))
        EditDis_str2.append(editDistance(aptamer2_str,candidate_str))
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
    #pool['edit_apt1_norm']=(1-(pool['edit_apt1']/max(pool['edit_apt1'])))
    #pool['edit_apt2_norm']=(1-(pool['edit_apt2']/max(pool['edit_apt2'])))
    #pool['edit_apt1_str_norm']=(1-(pool['edit_apt1_str']/max(pool['edit_apt1_str'])))
    #pool['edit_apt2_str_norm']=(1-(pool['edit_apt2_str']/max(pool['edit_apt2_str'])))
    pool['edit_apt1_norm']=1/(pool['edit_apt1']+1)
    pool['edit_apt2_norm']=1/(pool['edit_apt2']+1)
    pool['edit_apt1_str_norm']=1/(pool['edit_apt1_str']+1)
    pool['edit_apt2_str_norm']=1/(pool['edit_apt2_str']+1)
    if max(pool['dMFE'])==0.0:
        pool['dMFE_norm']=pool['dMFE']
    else:
        pool['dMFE_norm']=pool['dMFE']/max(pool['dMFE'])
    pool['score']=0.175*pool['edit_apt1_norm']+0.175*pool['edit_apt2_norm']+0.25*pool['edit_apt1_str_norm']+0.25*pool['edit_apt2_str_norm']+0.15*pool['dMFE_norm']
    return pool.nlargest(n_candidates,['score'])

def candidate_similarities(candidates,candidates_past):
    counter=0
    past_list=list(candidates_past.sequences)
    for seq in list(candidates.sequences):
        if seq in past_list:
            counter+=1
    return float(counter)/float(len(past_list))

BPA_Aptamer_design(100,1000,10,starting_sequence="CCACACCCACCCACAATGCCCCAC")                    
        
    
                
        
