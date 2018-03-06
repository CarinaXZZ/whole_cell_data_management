# Carina Xu
# Python 3.5
# Last updated: 15-09-2017

import data_pre
import numpy as np

# the function using wild type simulation to find the scale for the normalisation
def normalisation():
    cluster_data = data_pre.data_preprocess('/Users/superxzz/Documents/project/panda_dataframe/wildtype')
    data = np.array(cluster_data)
    Max = []
    Min = []

    # find the maximum and minimum value for each features in all wild type simulations
    for index in range(0, 17):
        max_value = data[:, index].max()
        # print(max_value)
        Max.append(max_value)
        min_value = data[:, index].min()
        if index == 13:
            min_value = 1
        # print(min_value)
        Min.append(min_value)
    # print(Max)
    # print(Min)
    return Max, Min
