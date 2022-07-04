#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
from fastatools import single_fasta_sequence, write_fasta

def get_cdna_seq(input_seq):
    comp_base={'A':'T','C':'G','G':'C','T':'A'}
    input_seq=input_seq.upper()
    comp_seq=[]
    for base in input_seq:
        comp_seq.append(comp_base[base])
    comp_seq.reverse()
    cdna=''.join(comp_seq)
    return cdna

def find_orfs(input_seq,cdna=False):
    start='ATG'
    stop={'TAA','TAG','TAA'}
    orfs=[]
    frame1=0
    frame2=1
    frame3=2
    for i in range(len(input_seq)):
        if input_seq[i:i+3] in stop:
            last=i+3
            if cdna==True:
                last=len(input_seq)-last+1
            if len(input_seq[frame1:i+3])%3==0:
                orf=''
                for j in range(len(input_seq[frame1:i+3])):
                    if input_seq[j:j+3]==start and len(orf)==0:
                        orf=input_seq[j:i+3]
                        first=j+1
                        if cdna==True:
                            first=len(input_seq)-first+1
                        orfs.append((f'frame: 3',orf, f' | :{first}-{last}'))
                        #first=last
                        frame1=i+3
                    elif len(orf)!=0:
                        break
            elif len(input_seq[frame2:i+3])%3==0:
                orf=''
                for j in range(len(input_seq[frame2:i+3])):
                    if input_seq[j:j+3]==start and len(orf)==0:
                        orf=input_seq[j:i+3]
                        first=j+1
                        if cdna==True:
                            first=len(input_seq)-first+1
                        orfs.append((f'frame: 3',orf, f' | :{first}-{last}'))
                        #first=last
                        frame2=i+3
                    elif len(orf)!=0:
                        break
            elif len(input_seq[frame3:i+3])%3==0:
                orf=''
                for j in range(len(input_seq[frame3:i+3])):
                    if input_seq[j:j+3]==start and len(orf)==0:
                        orf=input_seq[j:i+3]
                        first=j+1
                        if cdna==True:
                            first=len(input_seq)-first+1
                        orfs.append((f'frame: 3',orf, f' | :{first}-{last}'))
                        #first=last
                        frame3=i+1
                    elif len(orf)!=0:
                        break
    print(orfs)
    return orfs

header,input_seq=single_fasta_sequence(sys.argv[1])
cdna_seq=get_cdna_seq(input_seq)
output_file=sys.argv[2]
orfs=find_orfs(input_seq)
orfs_cdna=find_orfs(cdna_seq,cdna=True)

if __name__=='__main__':
    with open(output_file,'a') as file:
        for i in range(len(orfs)):
            header,input_seq=single_fasta_sequence(sys.argv[1])
            header=header+str(orfs[i][2])
            write_fasta(file,header,orfs[i][1])
        for i in range(len(orfs_cdna)):
            header,input_seq=single_fasta_sequence(sys.argv[1])
            header=header+str(orfs_cdna[i][2])
            write_fasta(file,header,orfs_cdna[i][1])
