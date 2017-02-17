'''
Created on 2017/02/13z

@author: yonezawa
'''

def distanceFromEven (point):
    first_term = 0
    denominator = len(point) ** 0.5
    numerator = 0
    
    for i in range(len(point)):
        first_term += point[i] ** 2
        numerator += point[i]

    second_term = (numerator / denominator) ** 2
    
    distance = (first_term - second_term) ** 0.5
        
    return distance

def oneClassSVM (data, nu = 0.1, kernel = 'rbf'):
    import numpy as np
    from sklearn.svm import OneClassSVM
    # import matplotlib.pyplot as plt
    
    original_data = []    
    tmp_data = []
    labels = []
    max_expression = 0
    
    for transcript in data.keys():
        for variant in sorted(data[transcript].keys()):
            labels.append(transcript + '_' + variant)
            
            this_data = data[transcript][variant]
            for i in range(len(this_data)):
                if max_expression < this_data[i]:
                    max_expression = this_data[i]
                    
            tmp_data.append([distanceFromEven(this_data)])
            original_data.append(this_data)

    data_input = np.array(tmp_data)
    
    clf = OneClassSVM(nu = nu, kernel = kernel)
    _ = clf.fit(data_input)
    prediction = clf.predict(data_input)
    
    x, y = np.meshgrid(np.linspace(0, max_expression, 100), np.linspace(0, max_expression, 100))
    Z = clf.decision_function(np.c_[abs(x.ravel() - y.ravel())]).reshape(x.shape)
    """
    max_expression = 15
    s = 40
    normal = data_input[prediction == 1]
    anomaly = data_input[prediction == -1]
    plt.scatter(normal[:, 0], np.zeros(normal.shape), c = 'blueviolet', s = s)
    plt.scatter(anomaly[:, 0], np.zeros(anomaly.shape), c = 'gold', s = s)
    plt.axis('tight')
    min_expression = 0
    plt.xlim((min_expression, max_expression))
    x, y = np.meshgrid(np.linspace(0, max_expression, 100), np.linspace(0, max_expression, 100))
    Z = clf.decision_function(np.c_[abs(x.ravel() - y.ravel())]).reshape(x.shape)
    plt.contourf(x, y, Z, level = [0, Z.max()])
    normal = []
    anomaly = []
    for i in range(prediction.shape[0]):
        if prediction[i] == 1:
            normal.append(original_data[i])
        else:
            anomaly.append(original_data[i])
    normal = np.array(normal)
    anomaly = np.array(anomaly)
    # normal = data_input[prediction == 1]
    # anomaly = data_input[prediction == -1]
    plt.scatter(normal[:, 0], normal[:, 1], c = 'blueviolet', s = s)
    plt.scatter(anomaly[:, 0], anomaly[:, 1], c = 'gold', s = s)
    plt.axis('tight')
    plt.xlim((0, max_expression))
    plt.ylim((0, max_expression))
    plt.show()
    """
    
    return prediction, labels, Z


if __name__ == '__main__':
    from math import log2
    from sys import maxsize
    
    directory = '/Users/yonezawa/research/data/paralvinella/'
    expressions = {}
    tissues = ['hot', 'cold']
    threshold = 2
    nus = [1e-4 * (i + 1) for i in range(9)]
    nus.extend([1e-3 * (i + 1) for i in range(10)])
    nu_value = {}
    
    datafile = directory + 'Tri_tra_Phe.TMM.fpkm.matrix'
        
    for line in open(datafile):
        elm = line.strip().split('\t')
        
        if elm[0] != 'result_hot':
            if log2(float(elm[1]) + 1) + log2(float(elm[2]) + 1) >= threshold:
                parts = elm[0].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                expressions.setdefault(transcript, {})
                expressions[transcript].setdefault(variant, {})            
                for i in range(1, 3):
                    tissue = tissues[2 - i]
                    expressions[transcript][variant][tissue] = log2(float(elm[i]) + 1)
    
    for nu in nus:
        prediction, labels, _ = oneClassSVM(expressions, nu = nu, kernel = 'rbf')
        upper_bound = -maxsize + 1
        
        for i in range(len(labels)):
            if prediction[i] > 0:
                parts = labels[i].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                
                this_data = []
                for tissue in tissues:
                    this_data.append(expressions[transcript][variant][tissue])
                distance = distanceFromEven(this_data)
                if upper_bound < distance:
                    upper_bound = distance
        
        for i in range(len(labels)):
            if prediction[i] < 0:
                parts = labels[i].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                
                this_data = []
                for tissue in tissues:
                    this_data.append(expressions[transcript][variant][tissue])
                distance = distanceFromEven(this_data)
                
                if distance > upper_bound:
                    nu_value.setdefault(labels[i], nu)
                    if nu_value[labels[i]] > nu:
                        nu_value[labels[i]] = nu
        
        print('nu: {0:.4f}'.format(nu))
        
    for label in sorted(nu_value.keys(), key = lambda x: nu_value[x]):
        parts = label.split('_')
        transcript = '_'.join(parts[:-1])
        variant = parts[-1]
        print('{0:s}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(label, nu_value[label], expressions[transcript][variant]['hot'], expressions[transcript][variant]['cold']))