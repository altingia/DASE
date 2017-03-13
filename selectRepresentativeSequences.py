'''
Created on 2017/03/09

@author: yonezawa
'''

def selectRepresentative (seqs):
    majorities = {}
    isoforms = list(seqs.keys())
    representatives = {}
    
    tmp_isoform = isoforms[0]
    for site in range(len(seqs[tmp_isoform])):
        count = {}
        for isoform in isoforms:
            if seqs[isoform][site] != '-':
                count.setdefault(seqs[isoform][site], 0)
                count[seqs[isoform][site]] += 1
                
        rep_nts = []
        max_count = 0
        for nt in count.keys():
            if max_count < count[nt]:
                max_count = count[nt]
                rep_nts = [nt]
            elif max_count == count[nt]:
                rep_nts.append(nt)
        
        for isoform in isoforms:
            majorities.setdefault(isoform, 0)
            if seqs[isoform][site] in rep_nts:
                majorities[isoform] += 1
                
    for isoform in sorted(majorities.keys(), key = lambda x: majorities[x], reverse = True):
        if len(list(representatives.keys())) == 0:
            representatives[isoform] = seqs[isoform]
            break
        
    return representatives

def summarizeRepresentative (directory):
    import glob
    aligned_files = glob.glob(directory + '*_aligned.fasta')
    representatives = {}
    
    for aligned_file in aligned_files:
        parts = aligned_file.split('/')[-1].split('_')
        transcript = '_'.join(parts[:-1])
        seqs = {}
        
        for line in open(aligned_file):
            if line[0] == '>':
                isoform = line.strip()[1:]
                seqs.setdefault(isoform, '')
            else:
                seqs[isoform] += line.strip()
        
        representatives[transcript] = selectRepresentative(seqs)        
    
    return representatives

if __name__ == '__main__':
    from time import time
    directory = '/Users/yonezawa/research/data/human_liver_brain/aligned_contigs/'
    result = directory + 'representative_sequences.fasta'
    
    start_time = time()
    representatives = summarizeRepresentative(directory)
    writefile = open(result, 'w')
    for transcript in sorted(representatives.keys()):
        isoform = list(representatives[transcript].keys())[0]
        writefile.write('>' + transcript + '_' + isoform + '\n')
        writefile.write(representatives[transcript][isoform] + '\n')
    writefile.close()
    end_time = time()
    print('Execution time: {0:.2f} seconds.'.format(end_time - start_time))