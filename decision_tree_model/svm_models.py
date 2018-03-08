# Title: svm_models
# Author: Zhaozhen Xu, xuuuzzhen520@gmail.com
# Affiliation: University of Bristol
# Run Using: Python 3.5
# Last Updated: 06/12/2017

from sklearn import svm
import pandas as pd
import os
import glob
import numpy as np
from scipy import stats


# train the svm models to classified the ploidy, protein, and rna
def return_models():
    # label of the training data
    ploidy_label = [0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0]
    protein_label = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    rna_label = [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0]

    ploidy = pd.read_pickle('data_ploidy.pickle')
    data_x = ploidy.as_matrix()
    ploidy_svm = svm.SVC()
    ploidy_svm.fit(data_x, ploidy_label)

    protein = pd.read_pickle('data_protein.pickle')
    data_x = protein.as_matrix()
    protein_svm = svm.SVC()
    protein_svm.fit(data_x, protein_label)

    rna = pd.read_pickle('data_rna.pickle')
    data_x = rna.as_matrix()
    rna_svm = svm.SVC()
    rna_svm.fit(data_x, rna_label)

    return ploidy_svm, protein_svm, rna_svm



# get the data by column
def data_byColumn(example_data, column_name):
    data = example_data[column_name]
    d = list(data)

    return d


# get the length of the simulation
def get_length(example_data):
    data = example_data['chromosome_ploidy']
    d = list(data)
    end_point = len(d)

    return end_point


# find the first and last value of the list(matrix)
def get_value(data):
    first_value = data[0]
    last_value = data[-1]

    return first_value, last_value


# calculate the standard deviation of the linear regression
def cal_std(data):
    std = np.std(data)

    return std


# calculate the gradient of the linear regression
def cal_gradient(data):
    x = range(1, len(data)+1)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,data)

    return gradient


# normalise the data with scale
def normalise_data(data, Max, Min):
    value = []
    for index in range(0, 4):
        if Max[index] != Min[index]:
            value.append((data[index]-Min[index])/(Max[index]-Min[index]))
        else:
            value.append(data)

    return value


# process the raw data
def process_data(data, Max, Min):
    processed_data = []
    for d in data:
        p = []
        first_value, last_value = get_value(d)
        p.append(first_value)
        p.append(last_value)
        p.append(cal_std(d))
        p.append(cal_gradient(d))

        normalised = normalise_data(p, Max, Min)
        processed_data.append(normalised)

    return processed_data


# get the raw data of the simulation
def get_all_data(path):
    data_pinched = []
    data_protein = []
    data_rna = []
    data_ploidy = []
    length = []
    filename = []

    column = ['geometry_pinchedDiameter', 'mass_proteinWt', 'mass_rnaWt', 'chromosome_ploidy']
    array_name = ['data_pinched', 'data_protein', 'data_rna', 'data_ploidy']

    for f in glob.glob(os.path.join(path, '*.pickle')):
        filename.append(f)
        example_data = pd.read_pickle(f)
        count = 0
        l = get_length(example_data)
        length.append(l)
        for column_name in column:
            data = data_byColumn(example_data, column_name)
            vars()[array_name[count]].append(data)
            count += 1
    scale = pd.read_pickle('scale.pickle')
    processed_protein = process_data(data_protein, list(scale['Max_protein']), list(scale['Min_protein']))
    processed_rna = process_data(data_rna, list(scale['Max_rna']), list(scale['Min_rna']))
    processed_ploidy = []

    # the length of training data(i.e. run1) is 50,000
    # the length of testing data(i.e. other runs) is 50,001
    # remove the first data point of the testing data
    for d in data_ploidy:
        d = d[1:]
        processed_ploidy.append(d)

    return data_pinched, processed_protein, processed_rna, processed_ploidy, length, filename
    # return data_pinched, processed_protein, processed_rna, data_ploidy, length, filename
    

# return the label of data
def return_labels(data, model):
    label = model.predict(data)

    return label
