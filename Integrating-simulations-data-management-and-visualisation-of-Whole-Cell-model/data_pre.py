# Carina Xu
# Python 3.5
# Last updated: 15-09-2017

import pandas as pd
import numpy as np
import os
import glob
from scipy import stats

# This is the data preprocessing methods without normalisation

# read the data by time series
def get_data(example_data, column_name):
    # example_data = pd.read_pickle(file_name)
    data = example_data[column_name]
    d = list(data)
    # print(d)
    return d

# calculate the standard deviation
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

# get the first and last value of the input data
def get_value(data):
    # print('first value')
    # print('last value')
    first_value = data[0]
    last_value = data[-1]
    return first_value, last_value

# get the time step the data start to increase
def get_changepoint_increase(data):
    # print('data changing point')
    d = pd.Series(data)
    change = d[d.pct_change() > 0].index.tolist()
    if not change:
        change_point = -1
    else:
        change_point = change[0]+1
    return change_point

# get the time step the data start to decrease
def get_changepoint_decrease(data):
    # print('initial data')
    # print('data changing point')
    d = pd.Series(data)
    initial_value = d[0]
    change = d[d.pct_change() < 0].index.tolist()
    if not change:
        change_point = -1
    else:
        change_point = change[0]+1
    return initial_value, change_point

# get the time step ploidy reach to 2, and the finaly value of ploidy
def get_changepoint_ploidy(data):
    # print('data changing point')
    if data[-1] == 2:
        change_point = data.index(2)+1
    else:
        change_point = -1
    end_value = data[-1]
    return change_point, end_value

# get the length of the simulation
def get_end(data):
    # print('length of simulation')
    end_point = len(data)+1
    if end_point >= 50000:
        end_point = -1
    return end_point

# count the number of non-zero data
def count_nonzero(data):
    # print('number of non zero value')
    count = np.count_nonzero(data)
    return count

# count the number of values equal to 1,2,3,4 respectively
def count_values(data):
    # print('number of 0')
    # print('number of 1')
    # print('number of 2')
    # print('number of 3')
    count = []
    # unique = np.unique(data)
    for num in [0,1,2,3]:
        count.append(data.count(num))
    return count

# get the pinched value (1 or 0)
def pinched_value(data):
    # print('pinched value')
    value = data[-1]
    return value

# data processing
def data_preprocess(path):
    # path = '/Users/superxzz/Documents/project/panda_dataframe/wildtype'
    # column = ['mass_total', 'mass_cell', 'mass_cellDry', 'mass_media', 'mass_waterWt', 'mass_metaboliteWt', 'mass_dnaWt', 'mass_rnaWt', 'mass_proteinWt', 'chromosome_ploidy', 'chromosome_segregated', 'geometry_pinchedDiameter', 'geometry_pinched', 'ftszRing_numEdgesOneStraight', 'ftszRing_numEdgesTwoStraight', 'ftszRing_numEdgesTwoBent', 'ftszRing_numResidualBent', 'ftszRing_numEdges', 'ribosome_nActive', 'ribosome_nNotExist', 'ribosome_nStalled', 'rnaPolymerase_nActive', 'rnaPolymerase_nSpecificallyBound', 'rnaPolymerase_nNonSpecificallyBound', 'rnaPolymerase_nFree', 'metabolicReaction_growth']
    # time series used in clustering
    column = ['metabolicReaction_growth', 'mass_proteinWt', 'mass_rnaWt', 'chromosome_ploidy', 'geometry_pinchedDiameter']
    x = []
    # y = []

    for filename in glob.glob(os.path.join(path, '*.pickle')):
        # print(filename)
        data_x = []

        # file_path = path+'/'+filename
        # print(file_path)
        example_data = pd.read_pickle(filename)
        for column_name in column:
            # print(column_name+':')
            data = get_data(example_data, column_name)
            # data = pd.Series(d)
            # print(data[0])
            if column_name =='geometry_pinched':
                value = pinched_value(data)
                data_x.append(value)
            elif column_name == 'ribosome_nStalled':
                value = count_nonzero(data)
                # print(value)
                data_x.append(value)
            elif column_name == 'rnaPolymerase_nSpecificallyBound':
                value = count_values(data)
                # print(value)
                for v in value:
                    data_x.append(v)
            elif column_name == 'chromosome_ploidy':
                value, end = get_changepoint_ploidy(data)
                data_x.append(value)
                data_x.append(end)
            elif column_name =='geometry_pinchedDiameter' or column_name == 'ftszRing_numEdges':
                # print(data)
                initial_value, changepoint = get_changepoint_decrease(data)
                data_x.append(initial_value)
                data_x.append(changepoint)
            elif column_name == 'chromosome_segregated' or column_name == 'ftszRing_numEdgesOneStraight' or column_name == 'ftszRing_numEdgesTwoStraight' or column_name == 'ftszRing_numEdgesTwoBent' or column_name == 'ftszRing_numResidualBent':
                # print(data)
                value = get_changepoint_increase(data)
                data_x.append(value)
            elif column_name =='rnaPolymerase_nFree':
                start, end = get_value(data)
                data_x.append(end)
                gradient = cal_gradient(data)
                data_x.append(gradient)
                # data_x.append(intercept)
                std = cal_std(data)
                data_x.append(std)
            elif column_name != 'mass_total' and column_name != 'mass_cellDry' and column_name != 'mass_dnaWt':
                start, end = get_value(data)
                data_x.append(start)
                data_x.append(end)
                gradient = cal_gradient(data)
                data_x.append(gradient)
                # data_x.append(intercept)
                std = cal_std(data)
                data_x.append(std)
            # print(len(data_x))
        timelength = get_end(data)
        data_x.append(timelength)
        # print(data_x)
        # print(len(data_x))
        x.append(data_x)
    # print(x)
    return x


