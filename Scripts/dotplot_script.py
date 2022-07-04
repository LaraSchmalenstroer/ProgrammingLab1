#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np 
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib','inline')

def open_file(filename):
    with open(filename,'r') as seq:
        fasta=seq.read().split('\n')[1:]
    fasta_seq=''.join(fasta)
    return fasta_seq

def dotplot(seqA,seqB,w=1,s=1):
    assert w%2==1, 'Please make sure your window size is an odd integer.'
    window=int((w-1)/2)
    seqA_n=('|'*window)+seqA+('|'*window)
    seqB_n=('-'*window)+seqB+('-'*window)
    plot=np.zeros((len(seqA_n),len(seqB_n)),dtype=int)
    for i in range(len(seqA)):
        partA=seqA_n[i:i+w]
        for j in range(len(seqB)):
            partB=seqB_n[j:j+w]
            count=0
            sets=[{partA[m],partB[m]} for m in range(len(partA))]
            for el in sets:
                if len(el)<2:
                    count+=1
            if count>=s:
                plot[i,j]=1
    return plot[:-4,:-4]

def dotplot2ascii(dp, seqA, seqB, heading, filename):
    list_dp=[list(dp[i]) for i in range(len(dp))]
    with open(filename,'w') as output:
        output.write(heading.upper()+'\n\n')
        output.write(' |'+seqB+'\n')
        output.write('_+'+'_'*(len(seqB))+'\n')
        for i in range(len(list_dp)):
            string=''
            for j in range(len(list_dp[i])):
                if list_dp[i][j]==1:
                    string+='*'
                else:
                    string+=' '
            output.write(seqA[i]+'|'+string+'\n')
    return

def dotplot2Graphics_v1(dp,labelA,labelB,heading,filename):
    assert filename.endswith('.png') or filename.endswith('.ps') or filename.endswith('.pdf'),     'Please enter a valid file extension type, valid types: .pdf, .png or .ps'
    xaxis=[labelB[i] for i in range(len(labelB))]
    yaxis=[labelA[i] for i in range(len(labelA)-1,-1,-1)]
    dp=np.rot90(dp, k=3)
    figure=plt.figure(figsize=(10,10))
    for i in range(len(xaxis)):
        for j in range(len(yaxis)):
            if dp[i,j]==1:
                plt.scatter(i,j,c='black', marker='+')
    if len(labelA)<100 and len(labelB)<100:
        plt.xticks(range(len(xaxis)),xaxis)
        plt.yticks(range(len(yaxis)),yaxis)
    else: 
        plt.xticks(range(len(xaxis)))
        plt.yticks(range(len(yaxis)))
    plt.title(heading, fontsize=15)
    plt.xlabel(labelB, fontsize=12)
    plt.ylabel(labelA, fontsize=12)
    plt.show()
    figure.savefig(filename)
    return

def dotplot2Graphics_v2(dp,labelA,labelB,heading,filename):
    assert filename.endswith('.png') or filename.endswith('.ps') or filename.endswith('.pdf'),     'Please enter a valid file extension type, valid types: .pdf, .png or .ps'
    figure = plt.figure(figsize=(10,10))
    axes = figure.add_subplot(111)
    xaxis = [labelB[i] for i in range(len(labelB))]
    yaxis = [labelA[i] for i in range(len(labelA))]
    caxes = axes.matshow(dp, cmap ='binary')
    plt.xlabel(labelB, fontsize=12)
    plt.ylabel(labelA, fontsize=12)
    plt.xticks(range(len(xaxis)),xaxis)
    plt.yticks(range(len(yaxis)),yaxis)
    plt.title(heading, fontsize=15)
    plt.show()
    figure.savefig(filename)
    return 

w=int(sys.argv[1])
s=int(sys.argv[2])
seqA=open_file(sys.argv[3])
seqB=open_file(sys.argv[4])
labelA=sys.argv[3]
labelB=sys.argv[4]
heading=sys.argv[5]
filename=sys.argv[6]

dp=dotplot(seqA,seqB,w,s)

if __name__=='__main__':
    if filename.endswith('.txt'):
        dotplot2ascii(dp,seqA,seqB,heading,filename)
    else:
        dotplot2Graphics_v2(dp,labelA,labelB,heading,filename)

