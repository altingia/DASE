'''
Created on 2017/03/09

@author: yonezawa
'''
def determineBirth (order_list, directory):
    parent_dir = directory + 'read_partitions/Fb_0/'    
    counts = {}
    
    for trinity_id in order_list:
        parts = trinity_id.split('_')
        base_num = int(int(parts[1][2:]) / 100)
        transcript_id = 'c' + parts[1][2:]
        contig_id = '_'.join(parts[2:])
        base_dir = parent_dir + 'CBin_' + str(base_num) + '/'
        read_seqs = {}
        flag = False
        seq = ''
        
        for line in open(base_dir + transcript_id + '.trinity.reads.fa.out.Trinity.fasta'):
            if line[0] == '>':
                elm = line[1:].strip().split()
                if elm[0] == contig_id:
                    flag = True
                else:
                    flag = False
            elif flag:
                seq += line.strip()
                
        read_id = ''
        for line in open(base_dir + transcript_id + '.trinity.reads.fa'):
            if line[0] == '>':
                read_id = line.strip()[1:line.find('/')]
                read_seqs.setdefault(read_id, '')
            else:
                read_seqs[read_id] += line.strip()
                
        counts.setdefault(trinity_id, {})
        for read_id in read_seqs.keys():
            sample_id = read_id[:read_id.find('.')]
            counts[trinity_id].setdefault(sample_id, 0)
            counts[trinity_id][sample_id] += 1
                
    return counts

if __name__ == '__main__':
    from time import time
    start_time = time()
    
    directory = '/Users/yonezawa/research/data/human_liver_brain/'
    distance_file = directory + 'distance_list_3.dat'
    result = directory + 'distance_list_3_with_birth.dat'
    order_list = []
    distances = {}
    for line in open(distance_file):
        elm = line.strip().split('\t')
        order_list.append(elm[0])
        distances[elm[0]] = elm[1]
    
    counts = determineBirth(order_list, directory + 'trinity_result/')
    
    writefile = open(result, 'w')
    for contig in counts.keys():
        writefile.write(contig + '\t' + distances[contig])
        for sample in counts[contig].keys():
            writefile.write('\t' + sample + '\t' + str(counts[contig][sample]))
        writefile.write('\n')
    writefile.close()
    end_time = time()
    print('Execution time: {0:.2f} seconds.'.format(end_time - start_time))