'''
Created on 2015/02/05

@author: yonezawa
'''

def executeMafft (mafft_exe, directory = '', gap_penalty = 10.0):
    import os, sys
    from Bio.Align.Applications import MafftCommandline
        
    if len(directory) > 0 and directory[-1] != '/':
        directory += '/'
    
    if len(mafft_exe) == 0:        
        sys.stderr.write('Install mafft before execution.')
        sys.exit(-1)
    
    after = directory + 'aligned_contigs/'
    if not os.path.exists(after):
        os.mkdir(after)
    
    seq_dir = directory + 'sequences/'
    seqfiles = os.listdir(seq_dir)
    for seqfile in seqfiles:
        if seqfile[-6:] == '.fasta':
            sequences = {}
            seq_ids = []
            for line in open(seq_dir + seqfile, 'r'):
                if line[0] == '>':
                    seq_ids.append(line.strip()[1:])
                else:
                    sequences.setdefault(seq_ids[-1], '')
                    sequences[seq_ids[-1]] += line.strip()
            
            transcript = seqfile[:seqfile.find('.')]
            mafft_cline = MafftCommandline(mafft_exe, input = seq_dir + seqfile)
            mafft_cline.set_parameter('--op', gap_penalty)
            writefile = open(after + transcript + '_aligned.fasta', 'w')
            stdout, stderr = mafft_cline()
            writefile.write(stdout)
            writefile.close()
        
if __name__ == '__main__':
    species = 'S_cerevisiae'
    directory = '/Users/yonezawa/research/data/' + species + '/'
    executeMafft('/usr/local/bin/mafft', directory = directory, gap_penalty = 10.0)