if __name__=='__main__':
    import sys
    #print('hello')
    dp=new_d(seqA,seqB,w,s)
    if output.endswith('.txt'):
    	dotplot2ascii(seqA,seqB,title,output)
    elif output.endswith('.pdf') or output.endswith('.png') or output.endswith('.ps'):
        dotplot2Graphics_v2(dp,seqA,seqB,title,output)