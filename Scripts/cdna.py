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

header,input_seq=single_fasta_sequence(sys.argv[1])
cdna=get_cdna_seq(input_seq)
output_file=sys.argv[2]
header=header+', cDNA'

if __name__=='__main__':
    with open(output_file,'w') as file:
        write_fasta(file,header,cdna)


