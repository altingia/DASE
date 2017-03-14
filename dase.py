'''
Created on 2016/05/31

@author: yonezawa
'''
def calcDistances (coverages, angles):
    from math import cos, radians
    
    distances = {}
    
    for transcript in angles.keys():
        contigs = sorted(list(angles[transcript].keys()), key = lambda x: int(x[1:]))
        distance = 0
        
        for i in range(len(contigs) - 1):
            for j in range(i + 1, len(contigs)):
                if contigs[i] in coverages[transcript].keys() and contigs[j] in coverages[transcript][contigs[i]].keys():
                    distance = 1 - (coverages[transcript][contigs[i]][contigs[j]] + cos(radians(angles[transcript][contigs[i]][contigs[j]]))) / 2
                    distances.setdefault(transcript, {})
                    distances[transcript].setdefault(contigs[i], {})
                    distances[transcript][contigs[i]][contigs[j]] = distance
                    distances[transcript].setdefault(contigs[j], {})
                    distances[transcript][contigs[j]][contigs[i]] = distance
    
    return distances

def rankASs (distances, min_transcripts = 3):
    from heap import heap, heappop  # @UnresolvedImport
    distance_list = []
    contig_list = []
    
    for transcript in distances.keys():
        contigs = list(distances[transcript].keys())
        
        if len(contigs) >= min_transcripts:        
            for i in range(len(contigs)):
                sum_dist = 0
                target_num = 0
                
                for j in range(len(contigs)):
                    if i != j and contigs[j] in distances[transcript][contigs[i]].keys():
                        sum_dist += distances[transcript][contigs[i]][contigs[j]]
                        target_num += 1
                
                contig_list.append(transcript + '_' + contigs[i])
                distance_list.append(sum_dist / target_num)
    
    heap_distances, heap_contigs  = heap(distance_list, contig_list)
        
    return heappop(heap_distances, heap_contigs)

def filterASs (distances, max_rank, min_transcripts = 3):
    ordered_list = rankASs(distances, min_transcripts)
    
    return ordered_list[:max_rank + 1]

if __name__ == '__main__':
    from parseFiles import parseFiles  # @UnresolvedImport
    from calculation import calcCoverages, calcAngles  # @UnresolvedImport
    from divideSequences import divideSequences  # @UnresolvedImport
    from executeMafft import executeMafft  # @UnresolvedImport
    from selectRepresentativeSequences import summarizeRepresentative  # @UnresolvedImport
    from determineBirths import determineBirth # @UnresolvedImport
    import argparse, os, sys
    
    species = 'human_liver_brain'
    directory = os.getcwd()
    # directory = '/Users/yonezawa/research/data/' + species
    usage = 'python {0:s}'.format(__file__)
    
    if directory[-1] != '/':
        directory += '/'
    
    parser = argparse.ArgumentParser(usage = usage)
    parser.add_argument('-ex', '--expression_files', type = str, dest = 'ex', help = 'expression files, comma separated', required = True)
    parser.add_argument('-mafft', '--mafft-dir', type = str, dest = 'mafft', help = 'where the executable MAFFT is', required = True)
    parser.add_argument('-trdir', '--trinity-dir', type = str, dest = 'trinity_dir', help = 'directory of the Trinity result', required = True)
    parser.add_argument('-seq', '--sequence-file', type = str, dest = 'seq', help = 'sequence file (FASTA format is necessary)', required = True)
    parser.add_argument('-ef', '--expression_file_format', type = str, dest = 'file_format', help = 'expression file format. (default: kallisto)')
    parser.add_argument('-th', '--expression_threshold', type = float, dest = 'expression_threshold', help = 'threshold of logarithm of the expression. (default: 2)')
    parser.add_argument('-cov', '--coverage_threshold', type = float, dest = 'coverage_threshold', help = 'threshold of the coverage. (dafault: 0.25)')
    parser.add_argument('-rep', '--representative_sequences', type = str, dest = 'representative_file', help = 'fasta-format output file for the representative sequence of each transcripts. (default: representative_sequences.fasta)')
    parser.add_argument('-o', '--output', type = str, dest = 'output_file', help = 'output file. (default: distance_list_{min_contigs}.dat)')
    parser.add_argument('-nc', '--number_of_contigs', type = int, dest = 'min_contigs', help = 'min. contigs in one transcript considered. (default: 3)')
    parser.add_argument('-gap', '--gap_penalty', type = float, dest = 'gap_penalty', help = 'gap penalty for MAFFT. (default: 10.0)')
    args = parser.parse_args()
    
    seqfile = args.seq
    
    trinity_dir = args.trinity_dir
    if trinity_dir[-1] != '/':
        trinity_dir += '/'
    
    expressionfiles = args.ex.split(',')
    if len(expressionfiles) < 2:
        sys.stderr.write('More than 1 expression files are necessary.')
        sys.exit(-1)
        
    for i in range(len(expressionfiles)):
        expressionfiles[i] = directory + expressionfiles[i]
    
    min_transcripts = 3
    if args.min_contigs:
        min_transcripts = args.min_contigs
        
    coverage_threshold = 0.25
    if args.coverage_threshold:
        coverage_threshold = args.coverage_threshold

    representative_file = directory + 'representative_sequences.fasta'
    if args.representative_file:
        representative_file = args.representative_file
    
    resultfile = directory + 'transcript_list_' + str(min_transcripts) + '_with_coverage_' + str(coverage_threshold) + '.dat'
    if args.output_file:
        resultfile = directory + args.output_file
    
    file_format = 'kallisto'
    if args.file_format:
        file_format = args.file_format
        
    gap_penalty = 10.0
    if args.gap_penalty:
        gap_panalty = args.gap_penalty
        
    threshold = 2
    if args.expression_threshold:
        threshold = args.threshold

    divideSequences(seqfile, seqdir = directory)
    print("Loaded the sequence file {0:s}.".format(seqfile))
    executeMafft(args.mafft, directory = directory, gap_penalty = gap_penalty)
    print("Executed MAFFT for each transcript containing more than 1 contigs.")
    representatives = summarizeRepresentative(directory + 'aligned_contigs/')
    rep_file = open(representative_file, 'w')
    for transcript in sorted(representatives.keys()):
        isoform = list(representatives[transcript].keys())[0]
        rep_file.write('>' + transcript + '_' + isoform + '\n')
        rep_file.write(representatives[transcript][isoform] + '\n')
    rep_file.close()
    print("Selected representative sequences.")
        
    log_ex, sequences = parseFiles(expressionfiles, file_format = 'kallisto', threshold = threshold, seq_dir = directory + 'aligned_contigs')
    print("Loaded the expression files.")

    coverages = calcCoverages(sequences, coverage_threshold)
    angles = calcAngles(log_ex, coverages)
    
    distances = calcDistances(coverages, angles)
    print("Calculated the distances.")
    
    ordered_list = filterASs(distances, len(distances) - 1, min_transcripts)
    print("Ranked contigs according to the above distances.")
    
    trinity_ids = ['' for _ in range(len(ordered_list))]
    for i in range(len(ordered_list)):
        trinity_ids[i] = ordered_list[i][0]
    read_counts = determineBirth(trinity_ids, trinity_dir)
    print("Summarized the evidence reads for each contig.")
    
    writefile = open(resultfile, 'w')
    for pair in ordered_list:
        writefile.write('{0:s}\t{1:.4f}'.format(pair[0], pair[1]))
        for sample in read_counts[pair[0]].keys():
            writefile.write('\t' + sample + '\t' + str(read_counts[pair[0]][sample]))
        writefile.write('\n')
    writefile.close()
    
    print('All tasks finished.')