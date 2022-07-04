#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
import re
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

def get_sequence_positions(fasta_file):
    list_seq=fasta_list(fasta_file)
    headers=[]
    positions={}
    for header,seq in list_seq:
        headers.append(header)
    for header in headers:
        x = re.search(r"[0-9]+[-][0-9]+", header)
        match=x.group(0)
        match=match.split('-')
        positions[int(match[0])]=int(match[1])
    return positions

def find_orfs(input_seq,threshold,cdna=False):
    start='ATG'
    stop={'TAA','TAG','TAA'}
    orfs=[]
    for i in range(3):
        orf=[]
        start_in=i
        end_in=i+3
        codon=input_seq[start_in:end_in]
        while codon:
            if codon==start and len(orf)==0:
                orf.append(codon)
                if cdna==False:
                    first=start_in+1
                elif cdna==True:
                    first=len(input_seq)-start_in
            elif len(orf)!=0:
                if codon not in stop:
                    orf.append(codon)
                elif codon in stop:
                    orf.append(codon)
                    if cdna==False:
                        last=end_in
                        if last-first>=threshold: #the approach counts the codons, so no multiplication * 3
                            orfs.append(('frame '+str(i+1),''.join(orf),f' | :{first}-{last}'))
                    elif cdna==True:
                        last=len(input_seq)-end_in+1
                        if first-last>=threshold: 
                            orfs.append(('frame '+str(i+1),''.join(orf),f' | :c{first}-{last}'))
                    orf=[]
            start_in=end_in
            end_in+=3
            codon=input_seq[start_in:end_in]
    print(orfs)
    return orfs

header,input_seq=single_fasta_sequence(sys.argv[1])
cdna_seq=get_cdna_seq(input_seq)
output_file=sys.argv[2]
threshold=int(sys.argv[3])
orfs=find_orfs(input_seq,threshold,cdna=False)
orfs_cdna=find_orfs(cdna_seq,threshold,cdna=True)

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
