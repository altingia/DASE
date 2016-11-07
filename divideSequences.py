'''
Created on 2016/05/31

@author: yonezawa
'''

def divideSequences (seqfile, seqdir = ''):
    import os
    
    if len(seqdir) > 0 and seqdir[-1] != '/':
        seqdir += '/'
        
    sequences = {}
    transcript = ''
    contig = ''
    result_dir = seqdir + 'sequences/'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
        
    for line in open(seqdir + seqfile):
        if line[0] == '>':
            elm = line.strip()[1:].split()
            parts = elm[0].split('_')
            transcript = '_'.join(parts[:-1]).replace('|', '_')
            contig = parts[-1]
            
            sequences.setdefault(transcript, {})
            sequences[transcript].setdefault(contig, '')
        else:
            sequences[transcript][contig] += line.strip()
            
    for transcript in sequences.keys():
        # print(transcript)
        contig_list = list(sequences[transcript].keys())
        
        if len(contig_list) > 1:
            writefile = open(result_dir + transcript + '.fasta', 'w')
            
            for contig in contig_list:
                writefile.write('>' + contig + '\n')
                writefile.write(sequences[transcript][contig] + '\n')
        
            writefile.close()
            
if __name__ == '__main__':
    seqdir = '/Users/yonezawa/research/data/S_cerevisiae/'
    seqfile = 'kissplice.fa'
    
    divideSequences(seqfile, seqdir)