'''
Created on 2016/05/30

@author: yonezawa
'''

def calcNorm (p1, p2):
    result = 0
    
    for i in range(len(p1)):
        result += p1[i] * p2[i]
    
    return result

def calcInnerProduct (p1, p2):
    return calcNorm(p1, p2) / (calcNorm(p1, p1) * calcNorm(p2, p2)) ** 0.5

def calcAngles (expressions, coverages):
    from math import acos, pi
    angles = {}
    
    for transcript in expressions.keys():
        if transcript in coverages.keys():
            contigs = list(expressions[transcript].keys())
            
            for i in range(len(contigs) - 1):
                for j in range(i + 1, len(contigs)):
                    inner_product = calcInnerProduct(expressions[transcript][contigs[i]], expressions[transcript][contigs[j]])
                    if inner_product >= 1:
                        inner_product = 1
                    
                    if contigs[i] in coverages[transcript].keys() and contigs[j] in coverages[transcript][contigs[i]].keys():
                        angles.setdefault(transcript, {})
                        angles[transcript].setdefault(contigs[i], {})
                        angles[transcript][contigs[i]][contigs[j]] = acos(inner_product) / pi * 180
                        angles[transcript].setdefault(contigs[j], {})
                        angles[transcript][contigs[j]][contigs[i]] = acos(inner_product) / pi * 180
                    
    return angles

def calcCoverages (sequences, coverage_threshold):
    coverages = {}
    for transcript in sequences.keys():
        contig_list = list(sequences[transcript].keys())
        if len(contig_list) > 1:
            coverages.setdefault(transcript, {})
            for i in range(len(contig_list) - 1):
                seq1 = sequences[transcript][contig_list[i]]
                for j in range(i + 1, len(contig_list)):
                    seq2 = sequences[transcript][contig_list[j]]
                    total_length = 0
                    same_site = 0
                    
                    for k in range(len(seq1)):
                        if seq1[k] != '-' or seq2[k] != '-':
                            total_length += 1
                            if seq1[k] == seq2[k]:
                                same_site += 1
                    
                    if float(same_site) / total_length >= coverage_threshold:            
                        coverages[transcript].setdefault(contig_list[i], {})
                        coverages[transcript][contig_list[i]][contig_list[j]] = float(same_site) / total_length
                        coverages[transcript].setdefault(contig_list[j], {})
                        coverages[transcript][contig_list[j]][contig_list[i]] = float(same_site) / total_length
    return coverages

if __name__ == '__main__':
    from parseFiles import parseFiles  # @UnresolvedImport
    
    directory = '/Users/yonezawa/research/data/human_tissues/'
    expresionfiles = []
    for ERR in ['kallisto_out/brain', 'kallisto_out/heart']:
        expresionfiles.append(directory + ERR + '/abundance.tsv')
        
    log_ex, sequences = parseFiles(expresionfiles, file_format = 'kallisto', threshold = 2, seq_dir = '/Users/yonezawa/research/data/human_tissues/aligned_contigs')

    coverages = calcCoverages(sequences)
    angles = calcAngles(log_ex, coverages)
    
    writefile = open(directory + 'brain_heart_comparison.dat', 'w')
    
    for transcript in angles.keys():
        contig_list = list(sorted(angles[transcript].keys(), key = lambda x: int(x[1:])))
        
        if len(contig_list) > 2:
            print(transcript)
            for i in range(len(contig_list) - 1):
                c1 = contig_list[i]
                for j in range(i + 1, len(contig_list)):
                    c2 = contig_list[j]
                    
                    print('{0:s}\t{1:s}\t{2:.3f}\t{3:.3f}'.format(c1, c2, coverages[transcript][c1][c2], angles[transcript][c1][c2]))
                    writefile.write('{0:s}\t{1:s}\t{2:.3f}\t{3:.3f}\n'.format(c1, c2, coverages[transcript][c1][c2], angles[transcript][c1][c2]))
    writefile.close()