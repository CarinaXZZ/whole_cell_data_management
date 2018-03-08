# Title: tree_classifier
# Author: Zhaozhen Xu, xuuuzzhen520@gmail.com
# Affiliation: University of Bristol
# Run Using: Python 3.5
# Last Updated: 06/12/2017

import csv
import svm_models

# recive the path of the pickle data file
path = input('Please input the path of the data:')
print(path)

label = []
result = []

# get data and svm models
data_pinched, data_protein, data_rna, data_ploidy, length, filename = svm_models.get_all_data(path)
print('get data')
ploidy_svm, protein_svm, rna_svm = svm_models.return_models()
print('get models')

# tree model
# mutant code: 0:DNA, 1:Septum, 2:Protein, 3:Metabolic, 4:Non-essential, 5:Slow-growing, 6:RNA
# label: 0:mutant, 1:normal
for i in range(0, len(length)):
    if length[i] >= 50000:
        if min(data_pinched[i]) == max(data_pinched[i]):
            if svm_models.return_labels(data_ploidy[i], ploidy_svm) == 1:
                if svm_models.return_labels(data_rna[i], rna_svm) == 1:
                    if svm_models.return_labels(data_protein[i], protein_svm) == 1:
                        label.append(1)
                        print('The simulation' + filename[i] + 'is Septum mutant.')
                        result.append('septum')
                    else:
                        label.append(2)
                        print('The simulation' + filename[i] + 'is Protein mutant.')
                        result.append('protein')
                else:
                    label.append(6)
                    print('The simulation' + filename[i] + 'is RNA mutant.')
                    result.append('RNA')
            else:
                if svm_models.return_labels(data_protein[i], protein_svm) == 1:
                    label.append(0)
                    print('The simulation' + filename[i] + 'is DNA mutant.')
                    result.append('DNA')
                else:
                    if svm_models.return_labels(data_rna[i], rna_svm) == 1:
                        label.append(2)
                        print('The simulation' + filename[i] + 'is Protein mutant.')
                        result.append('protein')
                    else:
                        label.append(3)
                        print('The simulation' + filename[i] + 'is Metabolic mutant.')
                        result.append('metabolic')
        else:
            label.append(5)
            print('The simulation' + filename[i] + 'is Slow growing.')
            result.append('slow-growing')
    else:
        label.append(4)
        print('The simulation' + filename[i] + 'is non-essential.')
        result.append('non-essential')

print('The mutant is:')
print(label)

# save the result to csv file
name = path.split('/')
csvfilename = name[-1] + '.csv'
with open(csvfilename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='|')
    for i in range(0, len(length)):
        writer.writerow([filename[i], result[i]])
