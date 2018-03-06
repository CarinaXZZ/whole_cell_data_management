# Carina Xu
# Python 3.5
# Last updated: 15-09-2017

import pandas as pd
import matplotlib.pyplot as plt


# visualisation of the raw data

column = ['mass_total', 'mass_cell', 'mass_cellDry', 'mass_media', 'mass_waterWt', 'mass_metaboliteWt', 'mass_dnaWt', 'mass_rnaWt', 'mass_proteinWt', 'chromosome_ploidy', 'chromosome_segregated', 'geometry_pinchedDiameter', 'geometry_pinched', 'ftszRing_numEdgesOneStraight', 'ftszRing_numEdgesTwoStraight', 'ftszRing_numEdgesTwoBent', 'ftszRing_numResidualBent', 'ftszRing_numEdges', 'ribosome_nActive', 'ribosome_nNotExist', 'ribosome_nStalled', 'rnaPolymerase_nActive', 'rnaPolymerase_nSpecificallyBound', 'rnaPolymerase_nNonSpecificallyBound', 'rnaPolymerase_nFree', 'metabolicReaction_growth']

example_data = pd.read_pickle('/Users/superxzz/Documents/project/panda_dataframe/wildtype/basic_summary_1.pickle')


for column_name in column:
    # print(column_name)
    title_name = 'Raw data of ' + column_name
    data = example_data[column_name]
    x = range(1, len(data) + 1)
    plt.figure()
    plt.plot(x, data)
    plt.xlabel('Time(s)')
    plt.ylabel(column_name)
    plt.title(column_name)
    plt.savefig(column_name)
