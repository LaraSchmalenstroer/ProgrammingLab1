#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def single_fasta_sequence(filename):
    with open(filename,'r') as seq:
#         title=seq.read().split('\n')[0][1:]
        fasta=seq.read().split('\n')
    title=fasta[0][1:]
    fasta_seq=''.join(fasta[1:])
    return (title,fasta_seq)
    
def fasta_list(filename):
    with open(filename,'r') as seqs:
        fasta=seqs.read().split('>')[1:]
    list_fasta=[]
    for i in range(len(fasta)):
        s=fasta[i].split('\n')
        title=s[0]
        fasta_seq=''.join(s[1:])
        list_fasta.append((title,fasta_seq))
    return list_fasta

def fasta_sequences(filename):
    with open(filename,'r') as file:
        sequence=''
        while True:
            line = file.readline()
            if sequence=='':
                sequence+=line
            elif sequence!='' and line[0]!='>':
                sequence+=line
            else:
                sequence=sequence.split('\n')
                header=sequence[0][1:]
                seq=''.join(sequence[1:])
                sequence=line
                yield(header,seq)
            if not line:
                break
    return

def write_fasta(file,header,sequence):
    n=70
    seq=[]
    i=0
    while i < len(sequence):
        if i+n < len(sequence):
            seq.append(sequence[i:i+n])
        else:
            seq.append(sequence[i:len(sequence)])
        i += n
    if header[0]!='>':
        header='>'+header
    file.write(header+'\n')
    for line in seq:
        file.write(line+'\n')
    return
