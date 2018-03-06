# Carina Xu
# Python 3.5
# Last updated: 15-09-2017

import pandas as pd
import numpy as np
import os
import glob
from scipy import stats
import normalise
import math

# This is the data preprocessing methods with normalisation and feature weighted

# normalised the data using (data-min)/(max-min)
def normalise_data(count, data, Max, Min):
    index = count
    if Max[index] != Min[index]:
        value = (data-Min[index])/(Max[index]-Min[index])
    else:
        value = data
    return value

# get the scale used to normalise the data
def get_normalisation_scale():
    Max, Min = normalise.normalisation()
    return Max, Min

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
    column = ['metabolicReaction_growth', 'mass_proteinWt', 'mass_rnaWt', 'chromosome_ploidy', 'geometry_pinchedDiameter']
    x = []
    # y = []
    Max, Min = get_normalisation_scale()
    print(Max)
    print(Min)
    for filename in glob.glob(os.path.join(path, '*.pickle')):
        # print(filename)
        data_x = []
        count = 0
        # scales = []
        # file_path = path+'/'+filename
        # print(file_path)
        example_data = pd.read_pickle(filename)
        for column_name in column:
            # print(count)
            # print(column_name+':')
            data = get_data(example_data, column_name)
            # data = pd.Series(d)
            # print(data[0])
            if column_name =='geometry_pinched':
                v = pinched_value(data)
                value = normalise_data(count,v,Max, Min)
                data_x.append(value)
                count +=1
            elif column_name == 'ribosome_nStalled':
                v = count_nonzero(data)
                # print(value)
                value = normalise_data(count, v, Max, Min)
                data_x.append(value)
                count +=1
            elif column_name == 'rnaPolymerase_nSpecificallyBound':
                value = count_values(data)
                # print(value)
                for v in value:
                    val = normalise_data(count, v, Max, Min)
                    data_x.append(val)
                    count += 1
            elif column_name == 'chromosome_ploidy':
                v,e= get_changepoint_ploidy(data)
                value = normalise_data(count, v, Max, Min)
                data_x.append(value)
                count += 1
                end = 6*normalise_data(count,e, Max, Min)
                # end = normalise_data(count, e, Max, Min)
                data_x.append(end)
                # scales.append(value)
                count += 1
            elif column_name =='geometry_pinchedDiameter' or column_name == 'ftszRing_numEdges':
                # print(data)
                initial_value, changepoint = get_changepoint_decrease(data)
                initial = normalise_data(count, initial_value, Max, Min)
                # value = scale_data(initial, scales[3], scales[0], scales[1], scales[2])
                data_x.append(initial)
                count += 1
                change = normalise_data(count, changepoint, Max, Min)
                # value = scale_data(change, scales[3], scales[0], scales[1], scales[2])
                data_x.append(change)
                count += 1
            elif column_name == 'chromosome_segregated' or column_name == 'ftszRing_numEdgesOneStraight' or column_name == 'ftszRing_numEdgesTwoStraight' or column_name == 'ftszRing_numEdgesTwoBent' or column_name == 'ftszRing_numResidualBent':
                # print(data)
                v = get_changepoint_increase(data)
                value = normalise_data(count, v, Max, Min)
                data_x.append(value)
                count += 1
            elif column_name == 'mass_rnaWt':
                start, end = get_value(data)
                start_value = normalise_data(count, start, Max, Min)
                data_x.append(start_value)
                count += 1
                end_value = 3*normalise_data(count, end, Max, Min)
                # end_value = normalise_data(count, end, Max, Min)
                data_x.append(end_value)
                count += 1
                gradient = cal_gradient(data)
                value = 3*normalise_data(count, gradient, Max, Min)
                # value = normalise_data(count, gradient, Max, Min)
                data_x.append(value)
                count += 1
                # data_x.append(intercept)
                std = cal_std(data)
                value = normalise_data(count, std, Max, Min)
                data_x.append(value)
                # scales.append(value)
                count += 1
            elif column_name != 'mass_total' and column_name != 'mass_cellDry' and column_name != 'mass_dnaWt':
                start, end = get_value(data)
                start_value = normalise_data(count, start, Max, Min)
                data_x.append(start_value)
                count += 1
                end_value = normalise_data(count, end, Max, Min)
                data_x.append(end_value)
                count += 1
                gradient = cal_gradient(data)
                value = normalise_data(count, gradient, Max, Min)
                data_x.append(value)
                count += 1
                # data_x.append(intercept)
                std = cal_std(data)
                value = normalise_data(count, std, Max, Min)
                data_x.append(value)
                count += 1
            # print(len(data_x))
        timelength = get_end(data)
        value = normalise_data(count, timelength, Max, Min)
        data_x.append(value)
        # print(data_x)
        # print(count)
        # print(len(data_x))
        x.append(data_x)
    # print(x)
    return x


