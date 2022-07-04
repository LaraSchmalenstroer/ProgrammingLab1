#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
import orffinder
from fastatools import fasta_list
import re

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

def no_orfs(orfs):
    return len(orfs)

def no_genes(genes):
    return len(genes)

def corr_genes(orfs,genes):
    set_items=orfs.items()&genes.items()
    number=len(set_items)
    ratio=number/len(genes)
    return (number,ratio)

def corr_stop(orfs,genes):
    set_keys=orfs.keys()&genes.keys()
    number=len(set_keys)
    ratio=number/len(genes)
    missing_genes=len(genes)-number
    return (number,ratio,missing_genes)

def properties_orfs(number_orfs,number_genes,orfs_corr,ratio_orfs,stop_corr,ratio_stop,missing_genes):
    print(f'Number of orfs found: {number_orfs}')
    print(f'Number of genes: {number_genes}')
    print(f'Number of orfs correctly predicting genes: {orfs_corr}, ({ratio_orfs} % of genes predicted)')
    print(f'Number of orfs correctly predicting a stop codon: {stop_corr}, ({ratio_stop} % of stop \
    codons predicted)')
    print(f'Number of missing genes: {missing_genes}')
    return

file_orfs=sys.argv[1]
file_genes=sys.argv[2]
orfs=get_sequence_positions(sys.argv[1])
genes=get_sequence_positions(sys.argv[2])
number_orfs=no_orfs(orfs)
number_genes=no_genes(genes)
orfs_corr,ratio_orfs=corr_genes(orfs,genes)
stop_corr,ratio_stop,missing_genes=corr_stop(orfs,genes)

if __name__=='__main__':
    properties_orfs(number_orfs,number_genes,orfs_corr,ratio_orfs,stop_corr,ratio_stop,missing_genes)
    