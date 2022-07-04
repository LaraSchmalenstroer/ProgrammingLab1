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
                        orfs.append(('frame '+str(i+1),''.join(orf),f' | :{first}-{last}'))
                    elif cdna==True:
                        last=len(input_seq)-end_in+1
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
