# Title: tree_dataprocessing
# Author: Zhaozhen Xu, xuuuzzhen520@gmail.com
# Affiliation: University of Bristol
# Run Using: Python 3.5
# Last Updated: 06/12/2017

import pandas as pd
import glob
import os
from scipy import stats
import numpy as np


# find the first and last value of the list(matrix)
def get_value(data):
    # print('first value')
    # print('last value')
    first_value = data[0]
    last_value = data[-1]
    return first_value, last_value


# calculate the standard deviation of the linear regression
def cal_std(data):
    # print('std')
    std = np.std(data)
    return std


# calculate the gradient of the linear regression
def cal_gradient(data):
    # print('gradient')
    # gradient = np.gradient(data)
    # print(gradient)
    x = range(1, len(data)+1)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,data)
    return gradient


# get the data by column
def data_byColumn(example_data, column_name):
    # example_data = pd.read_pickle(file_name)
    data = example_data[column_name]
    d = list(data)
    # print(d)
    return d


# data processing for protein and RNA
def process_data(data, Max, Min):
    procdssed_data = []
    for d in data:
        p = []
        first_value, last_value = get_value(d)
        p.append(first_value)
        p.append(last_value)
        p.append(cal_std(d))
        p.append(cal_gradient(d))

        normalised = normalise_data(p, Max, Min)
        procdssed_data.append(normalised)
    return procdssed_data


# normalise the data using wild-type simulation
def normalise_data(data, Max, Min):
    value = []
    for index in range(0, 4):
        if Max[index] != Min[index]:
            value.append((data[index]-Min[index])/(Max[index]-Min[index]))
        else:
            value.append(data)
    return value


# find the scale of the normalisation
def normalised_scale(path):
    # get the protein and rna data from the wild-type simulations
    processed_protein, processed_rna = get_wildtype_data(path)
    data_protein = np.array(processed_protein)

    Max_protein = []
    Min_protein = []
    data_rna = np.array(processed_rna)
    Max_rna = []
    Min_rna = []

    # find the maximum and minimum value for each features in all wild type simulations
    for index in range(0, 4):
        max_value = data_protein[:, index].max()
        # print(max_value)
        Max_protein.append(max_value)
        max_value = data_rna[:, index].max()
        Max_rna.append(max_value)
        min_value = data_protein[:, index].min()
        Min_protein.append(min_value)
        min_value = data_rna[:, index].min()
        Min_rna.append(min_value)
    # print(Max)
    # print(Min)
    return Max_protein, Min_protein, Max_rna, Min_rna


def transpose_scale(Max_protein, Min_protein, Max_rna, Min_rna):
    scale = []
    scale.append(Max_protein)
    scale.append(Min_protein)
    scale.append(Max_rna)
    scale.append(Min_rna)
    s = np.array(scale)
    return s.transpose()


# find the length of the simulation
def get_length(example_data):
    # print('length of simulation')
    data = example_data['chromosome_ploidy']
    d = list(data)
    end_point = len(d)
    return end_point


# get the wild type data
def get_wildtype_data(path):
    data_protein = []
    data_rna = []

    column = ['mass_proteinWt', 'mass_rnaWt']
    array_name = ['data_protein', 'data_rna']

    for filename in glob.glob(os.path.join(path, '*.pickle')):
        example_data = pd.read_pickle(filename)
        count = 0
        for column_name in column:
            data = data_byColumn(example_data, column_name)
            vars()[array_name[count]].append(data)
            count += 1
    process_protein = []
    process_rna = []
    for d in data_protein:
        p = []
        first_value, last_value = get_value(d)
        p.append(first_value)
        p.append(last_value)
        p.append(cal_std(d))
        p.append(cal_gradient(d))
        process_protein.append(p)
    for d in data_rna:
        p = []
        first_value, last_value = get_value(d)
        p.append(first_value)
        p.append(last_value)
        p.append(cal_std(d))
        p.append(cal_gradient(d))
        process_rna.append(p)
    return process_protein, process_rna


# get the data from run1
def get_all_data(path):
    # data_growth = []
    MAX_LENGTH = 50000
    # MAX_LENGTH = 50001
    data_protein = []
    data_rna = []
    data_ploidy = []
    data_pinched = []
    length = []
    essential = []

    column = ['mass_proteinWt', 'mass_rnaWt', 'chromosome_ploidy',
              'geometry_pinchedDiameter']
    array_name = ['data_protein', 'data_rna', 'data_ploidy', 'data_pinched']
    pickle = 0
    for filename in glob.glob(os.path.join(path, '*.pickle')):
        example_data = pd.read_pickle(filename)
        count = 0
        l = get_length(example_data)
        length.append(l)
        if l == MAX_LENGTH:
            for column_name in column:
                data = data_byColumn(example_data, column_name)
                vars()[array_name[count]].append(data)
                count += 1
            essential.append(pickle)
        pickle += 1
    print('get all data')
    return data_protein, data_rna, data_ploidy, data_pinched, length, essential


# save the processed data from the run1 and the scale for normalising the data
# get the raw data
run1_path = '/Users/superxzz/Documents/project/panda_dataframe/1kos/run1'
data_protein, data_rna, data_ploidy, data_pinched, length, essential = get_all_data(run1_path)

# normalisation
wildtype_path = '/Users/superxzz/Documents/project/panda_dataframe/wildtype'
Max_protein, Min_protein, Max_rna, Min_rna = normalised_scale(wildtype_path)
scale = transpose_scale(Max_protein, Min_protein, Max_rna, Min_rna)

# processing the data of protein and rna
processed_rna = process_data(data_rna, Max_rna, Min_rna)
processed_protein = process_data(data_protein, Max_protein, Min_protein)

# save to pickle file
data = pd.DataFrame(data=processed_protein)
data.to_pickle('data_protein.pickle')
data = pd.DataFrame(data=processed_rna)
data.to_pickle('data_rna.pickle')
data = pd.DataFrame(data=data_ploidy)
data.to_pickle('data_ploidy.pickle')
data = pd.DataFrame(data=data_pinched)
data.to_pickle('data_pinched.pickle')
data = pd.DataFrame(data=scale, columns=['Max_protein', 'Min_protein', 'Max_rna', 'Min_rna'])
data.to_pickle('scale.pickle')
