'''
Created on 2016/06/02

@author: yonezawa
'''
def heap (value_list, label_list):
    from math import floor
    
    value_heap = []
    label_heap = []
    
    for i in range(len(value_list)):
        value_heap.append(value_list[i])
        label_heap.append(label_list[i])
        
        index = i
        parent = int(floor((index - 1) / 2))
        
        while index > 0 and value_heap[parent] < value_heap[index]:
            tmp = value_heap[index]
            value_heap[index] = value_heap[parent]
            value_heap[parent] = tmp
            tmp = label_heap[index]
            label_heap[index] = label_heap[parent]
            label_heap[parent] = tmp
            
            index = parent
            parent = int(floor((index - 1) / 2))
    
    return value_heap, label_heap

def heappop (value_heap, label_heap):
    result = []
    
    while len(value_heap) > 0:
        result.append([label_heap[0], value_heap[0]])
        
        label_heap[0] = label_heap[-1]
        value_heap[0] = value_heap[-1]
        label_heap = label_heap[:-1]
        value_heap = value_heap[:-1]
        
        index = 0
        while ((index + 1) * 2 == len(value_heap) and value_heap[(index + 1) * 2 - 1] > value_heap[index]) or ((index + 1) * 2 < len(value_heap) and (value_heap[(index + 1) * 2 - 1] > value_heap[index] or value_heap[(index + 1) * 2] > value_heap[index])):
            swap_index = (index + 1) * 2
            if (index + 1) * 2 == len(value_heap) or value_heap[(index + 1) * 2 - 1] > value_heap[(index + 1) * 2]:
                swap_index -= 1

            tmp = value_heap[index]
            value_heap[index] = value_heap[swap_index]
            value_heap[swap_index] = tmp
            tmp = label_heap[index]
            label_heap[index] = label_heap[swap_index]
            label_heap[swap_index] = tmp
            
            index = swap_index
            
    return result

def findBoundary (ordered_list, threshold):
    boundary = int(len(ordered_list) / 2)
    max_bound = len(ordered_list) - 1
    min_bound = 0
    
    if threshold > ordered_list[0][1]:
        return 0
    if threshold < ordered_list[-1][1]:
        return len(ordered_list) - 1
    
    while ordered_list[boundary][1] < threshold or ordered_list[boundary + 1][1] > threshold:
        if ordered_list[boundary][1] < threshold:
            max_bound = boundary
        else:
            min_bound = boundary
        boundary = int((min_bound + max_bound) / 2)
    
    return boundary

if __name__ == '__main__':
    value_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    label_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    
    value_heap, label_heap = heap(value_list, label_list)
    
    result = heappop(value_heap, label_heap)
    
    for pair in result:
        print('{0:s}\t{1:d}'.format(pair[0], pair[1]))
        
    print(findBoundary(result, 8.43))