#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np 
import pandas as pd
import itertools

def open_file(filename):
    with open(filename,'r') as seq:
        fasta=seq.read().split('\n')[1:]
    fasta_seq=''.join(fasta)
    return fasta_seq

def get_scoring_matrix(scoring_mat):
    if scoring_mat=='blosum62.txt':
        with open(scoring_mat,'r') as f:
            blosum=f.read().split('\n')
        blosum=blosum[:-2]
        blosum=[blosum[i].split(' ') for i in range(len(blosum))]
        blosum62=[]
        for i in range(len(blosum)):
            el=[]
            for j in range(1,len(blosum[i])):
                if blosum[i][j]!='':
                    el.append(blosum[i][j])
            blosum62.append(el)
        for i in range(1,len(blosum62)):
            for j in range(len(blosum62[i])):
                blosum62[i][j]=int(blosum62[i][j])
        blo=pd.DataFrame(blosum62[1:],index=blosum62[0], columns=blosum62[0])
    elif scoring_mat=='blosum50.txt':
        with open('blosum50.txt','r') as f:
            blosum=f.read().split('\n')
        blosum=blosum[:-2]
        blosum=[blosum[i].split('\t') for i in range(len(blosum))]
        for i in range(1,len(blosum)):
            blosum[i]=blosum[i][1:]
        for i in range(len(blosum)):
            for j in range(len(blosum[i])):
                blosum[i][j]=blosum[i][j].replace(' ','')
                if i>0:
                    blosum[i][j]=int(blosum[i][j])
        blo=pd.DataFrame(blosum[1:],index=blosum[0], columns=blosum[0])
    return blo

def local_alignment(seqA,seqB,w,blosum):    
    if w>0:
        w=(-w)
    al=np.zeros((len(seqA)+1,len(seqB)+1),dtype=int)
    tb=np.zeros((len(seqA)+1,len(seqB)+1),dtype=list)
    for i in range(len(seqB)+1):
        al[0,i]=i*w
        if i>0:
            tb[0,i]=[0,i-1]
        elif i==0:
            tb[0,i]=[0,0]
    for i in range(len(seqA)+1):
        al[i,0]=i*w
        if i>0:
            tb[i,0]=[i-1,0]
    for i in range(1,len(seqB)+1):
        tb[1,i]=[1,i-1]
    for i in range(1,len(seqA)+1):
        tb[i,1]=[i-1,1]
    for i in range(1,len(seqA)+1):
        for j in range(1,len(seqB)+1):
            gap_1=al[i-1,j]+w
            gap_2=al[i,j-1]+w
            match=al[i-1,j-1]+blosum[seqA[i-1]][seqB[j-1]] #if multiple elements have the same value, the first one is picked
            al[i,j]=max(match,0,gap_1,gap_2)
            if al[i,j]==match:
                tb[i,j]=[i-1,j-1]
            elif al[i,j]==gap_1:
                tb[i,j]=[i-1,j]
            elif al[i,j]==gap_2:
                tb[i,j]=[i,j-1]
    score=np.amax(al)
    index_score=np.unravel_index(np.argmax(al, axis=None), al.shape)
    #Traceback:
    x='' 
    y='' 
    m=index_score[0]
    n=index_score[1]
    while al[m][n]!=0:
        if tb[m,n]==[m-1,n-1]:
            x=seqA[m-1]+x
            y=seqB[n-1]+y
            m=m-1
            n=n-1
        elif tb[m,n]==[m,n-1]:
            x='-'+x
            y=seqB[n-1]+y
            n=n-1
        else:
            y='-'+y
            x=seqA[m-1]+x
            m=m-1
    alignment={'Sequence A: ': x, 'Sequence B: ':y}
    return score, alignment

def output_alignment(score, alignment,blosum):
    seq1=alignment['Sequence A: ']
    seq2=alignment['Sequence B: ']
    n = 80
    seq_1 = []
    seq_2 = []
    i = 0
    count_id=0
    count_si=0
    length=len(seq1)
    while i < len(seq1):
        if i+n < len(seq1):
            seq_1.append(seq1[i:i+n])
            seq_2.append(seq2[i:i+n])
        else:
            seq_1.append(seq1[i:len(seq1)])
            seq_2.append(seq2[i:len(seq2)])
        i += n
    print('Score = '+str(score))
    for i in range(len(seq_1)):
        print(seq_1[i])
        comparison=''
        for j in range(len(seq_1[i])):
            if seq_1[i][j]!='-' and seq_2[i][j]!='-':
                if seq_1[i][j]==seq_2[i][j]:
                    comparison+='|'
                    count_id+=1
                    count_si+=1
                elif blosum[seq_1[i][j]][seq_2[i][j]]>0:
                    comparison+=':'
                    count_si+=1
                else:
                    comparison+=' '
            else:
                comparison+=' '
        print(comparison)
        print(seq_2[i])
    print('Sequence identity= {} %, Sequence similarity= {} %'.format(round((count_id/length)*100,2),round((count_si/length)*100,2)))
    return 

if sys.argv[1].endswith('.fasta') and sys.argv[2].endswith('.fasta'):
    seqA=open_file(sys.argv[1])
    seqB=open_file(sys.argv[2])
else:
    seqA=sys.argv[1]
    seqB=sys.argv[2]
w=int(sys.argv[3])
blosum=get_scoring_matrix(sys.argv[4])
score,alignment=local_alignment(seqA,seqB,w,blosum)

if __name__=='__main__':
    output_alignment(score,alignment,blosum)

