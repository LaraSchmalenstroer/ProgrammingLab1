#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import itertools
get_ipython().run_line_magic('matplotlib','inline')

def open_file(alignment):
    with open(alignment,'r') as file:
        sequences=file.read().split('\n')
    matrix=[[char for char in seq] for seq in sequences]
    return matrix

def get_aminoacids(matrix):
    aa=set(''.join(matrix[0]))
    aminoacids=[]
    if len(aa)==20:
        seq="".join(map(str,aa))
        aminoacids=[char for char in seq]
    else:
        n=1
        while len(aminoacids)<20 or n==len(matrix):
            aa=set(''.join(matrix[n]))
            seq="".join(map(str,aa))
            for char in seq:
                if char not in aminoacids:
                    aminoacids.append(char)
    return aminoacids

def get_single_frequencies(matrix,aminoacids):
    number_aa=len(matrix[0])*len(matrix)
    fa={}
    pa={}
    for i in range(len(aminoacids)):
        aa=aminoacids[i]
        freq=0
        for j in range(len(matrix)): 
            freq+=matrix[j].count(aa)
        fa[aa]=freq
        pa[aa]=freq/number_aa
    return fa,pa

def get_pairwise_frequencies(matrix,aminoacids):
    fab={}
    amino=''.join(aminoacids)
    poss=itertools.combinations_with_replacement(amino,2)
    for po in poss:
        p=''.join(po)
        fab[p]=1 #add pseudocounts
    for i in range(len(matrix[0])):
        letters=''
        for n in range(len(matrix)):
            letters+=matrix[n][i]
        combi=itertools.combinations(letters,2)
        for com in combi:
#             if com[0]==com[1]:
#                 ind="".join(map(str,c))
#                 ind=ind*2
#                 fab[ind]+=1
            c=set(com)
            if len(c)==1:
                ind="".join(map(str,c))
                ind=ind*2
                fab[ind]+=1
            else:
                ind="".join(map(str,c))
                ind=ind[0]+ind[1]
                if ind not in fab.keys():
                    ind=ind[::-1]
                fab[ind]+=1
    number_pairwise=sum(fab.values())
    pab={key:(fab[key]/number_pairwise) for key in fab}
    return fab,pab

def get_eab(fab,pa):
    eab={}
    for key in fab:
        word=key
        l1=word[0]
        l2=word[0]
        if l1==l2:
            prob=pa[l1]**2
        else:
            prob=pa[l1]*pa[l2]*2
        eab[word]=prob
    return eab

def get_sab(pab,eab):
    sab={}
    from math import log2
    for key in pab:
        score=2*log2(pab[key]/eab[key])
        sab[key]=int(score)
    return sab
        
def create_score_matrix(sab,aminoacids):
    blosum=np.zeros((len(aminoacids),len(aminoacids)))
    #blosum=np.concatenate(names,blosum)
    pairwise=pd.DataFrame(blosum,index=aminoacids,columns=aminoacids)
    for i in range(len(aminoacids)):
        aa1=aminoacids[i]
        for j in range(len(aminoacids)):
            aa2=aminoacids[j]
            combi=aa1+aa2
            if combi not in sab.keys():
                combi=combi[::-1]
            pairwise[aa1][aa2]=sab[combi]
    print(pairwise)
    return pairwise

alignment=sys.argv[1]
matrix=open_file(alignment)
aminoacids=get_aminoacids(matrix)
fa,pa=get_single_frequencies(matrix,aminoacids)
fab,pab=get_pairwise_frequencies(matrix,aminoacids)
eab=get_eab(fab,pa)
sab=get_sab(pab,eab)

blosum62=create_score_matrix(sab,aminoacids)

if __name__=='__main__':
    blosum62


# In[ ]:




