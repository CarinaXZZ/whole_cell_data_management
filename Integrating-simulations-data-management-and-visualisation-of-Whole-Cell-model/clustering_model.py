# Carina Xu
# Python 3.5
# Last updated: 15-09-2017

import data_pre_normal
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.decomposition import PCA
import pylab as pl
from sklearn import cluster


# ground truth of the clustering models
cluster_true = [0,3,4,3,2,4,3,4,4,2,3,4,4,7,2,4,3,4,2,2,2,2,2,4,2,2,2,2,2,2,2,2,6,2,1,2,2,2,2,2,4,3,4,6,2,4,2,2,7,2,3,3,3,2,4,4,2,4,4,4,4,4,4,4,2,2,2,3,2,2,4,1,1,1,4,4,4,4,3,4,2,4,3,4,4,3,3,4,4,2,1,1,4,4,4,3,3,3,3,6,4,2,4,4,5,4,4,4,3,3,0,1,4,2,4,6,2,2,2,4,0,4,4,4,4,4,6,3,3,3,3,3,3,3,0,4,3,7,2,3,4,4,4,4,2,4,4,2,4,4,3,3,3,3,3,3,3,4,4,4,2,4,0,4,4,4,3,2,3,3,4,2,4,4,3,4,2,4,2,3,4,4,7,7,3,4,2,4,2,3,4,3,4,4,4,4,3,4,4,2,3,2,2,6,2,2,3,4,4,2,2,4,4,2,2,4,3,3,4,4,4,1,1,4,4,4,4,4,3,3,4,4,5,4,3,4,4,4,4,3,4,4,4,4,4,3,2,6,0,1,6,4,2,4,3,3,3,3,2,3,3,4,4,2,2,3,2,4,4,2,3,4,2,3,4,3,2,4,2,2,4,4,4,0,1,4,3,4,4,2,4,5,4,3,4,2,3,4,4,4,3,4,4,4,3,3,2,3,4,4,4,0,3,3,3,3,6,2,4,2,4,4,2,2,2,2,6,0,2,2,4,2,4,2,2,3,4,4,4,2,3,4,4,4,3,3,2,3,3,4,4,4,4]
# ground truth of the prediciton
ground_true = [0, 3, 4, 3, 2, 4, 3, 4, 4, 2, 3, 4, 4, 7, 2, 4, 3, 4, 2, 2, 2, 2, 2, 4, 2, 6, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 6, 4, 3, 4, 2, 2, 4, 2, 2, 7, 2, 3, 3, 3, 2, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 2, 6, 3, 2, 2, 4, 1, 1, 1, 4, 0, 4, 4, 3, 4, 2, 4, 3, 4, 4, 3, 3, 4, 4, 2, 1, 1, 4, 4, 4, 3, 3, 3, 3, 2, 4, 2, 4, 4, 1, 5, 4, 4, 3, 7, 0, 1, 4, 2, 4, 2, 1, 2, 2, 4, 0, 4, 4, 4, 4, 4, 2, 3, 3, 3, 3, 3, 3, 3, 0, 4, 3, 7, 2, 3, 4, 4, 4, 4, 2, 4, 4, 2, 4, 4, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 2, 4, 0, 4, 4, 4, 3, 2, 3, 3, 4, 2, 4, 4, 3, 2, 3, 3, 0, 1, 4, 3, 3, 3, 3, 3, 4, 3, 4, 2, 2, 4, 4, 4, 3, 4, 2, 3, 4, 4, 4, 3, 4, 4, 4, 3, 3, 2, 3, 4, 4, 4, 0, 3, 3, 3, 3, 2, 2, 4, 4, 4, 4, 2, 2, 6, 2, 2, 0, 6, 2, 0, 4, 2, 4, 2, 2, 4, 3, 4, 4, 4, 6, 3, 4, 4, 4, 3, 3, 2, 3, 3, 4, 4, 4, 4]

cluster_pred = []

# get the data for clustering (run 1 of single gene knockouts)
cluster_data = data_pre_normal.data_preprocess('/Users/superxzz/Documents/project/panda_dataframe/1kos/run1')
# get the data for prediction (run 2 of single gene knockouts)
predict_data = data_pre_normal.data_preprocess('/Users/superxzz/Documents/project/panda_dataframe/1kos/run2')

# K-Means
print('K-Means')
kmeans = KMeans(n_clusters=8).fit(cluster_data)
cluster_pred = kmeans.labels_
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)

# the visualisation of the result of clustering using PCA
pca = PCA(n_components=2).fit(cluster_data)
pca_2d = pca.transform(cluster_data)
for i in range(0, pca_2d.shape[0]):
    if cluster_pred[i] == 0:
        c1 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='c')
    elif cluster_pred[i] == 1:
        c2 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='m')
    elif cluster_pred[i] == 2:
        c3 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='g')
    elif cluster_pred[i] == 3:
        c4 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='r')
    elif cluster_pred[i] == 4:
        c5 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='y')
    elif cluster_pred[i] == 5:
        c6 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='k')
    elif cluster_pred[i] == 6:
        c7 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='b')
    elif cluster_pred[i] == 7:
        c8 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='r', marker='*')

pl.legend([c1, c2, c3, c4, c5, c6, c7, c8], ['Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 8'])
pl.title('K-means clusters the 1kos dataset into 8 clusters')
# pl.ylim((-10000,5000))
pl.savefig('plot_kmeans')

# plot the original label of dataset as the ground truth
for i in range(0, pca_2d.shape[0]):
    if cluster_true[i] == 0:
        c1 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='r')
    elif cluster_true[i] == 1:
        c2 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='b')
    elif cluster_true[i] == 2:
        c3 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='g')
    elif cluster_true[i] == 3:
        c4 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='c')
    elif cluster_true[i] == 4:
        c5 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='m')
    elif cluster_true[i] == 5:
        c6 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='y')
    elif cluster_true[i] == 6:
        c7 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='k')
    elif cluster_true[i] == 7:
        c8 = pl.scatter(pca_2d[i, 0], pca_2d[i, 1], c='r', marker='*')

pl.legend([c1, c2, c3, c4, c5, c6, c7, c8], ['DNA', 'Septum', 'Protein', 'Metabolic', 'Non-essential', 'Slow Growing', 'Protein and DNA', 'RNA'])
pl.title('Original label of mutants for 1kos dataset')
# pl.ylim((-10000,5000))
pl.savefig('plot_original')


# spectral clustering
print('spectral clustering')
spectral = cluster.SpectralClustering(n_clusters= 8 , eigen_solver='arpack', affinity="nearest_neighbors").fit(cluster_data)
cluster_pred = spectral.labels_
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)

# agglomerative clustering with complete linkage
print('agglomerative clustering with complete linkage')
hierarchical = cluster.AgglomerativeClustering(n_clusters=8, linkage="complete").fit(cluster_data)
cluster_pred = hierarchical.labels_
# print(cluster_pred)
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)


# agglomerative clustering with ward linkage
print('agglomerative clustering with ward linkage')
hierarchical = cluster.AgglomerativeClustering(n_clusters=8, linkage="ward").fit(cluster_data)
cluster_pred = hierarchical.labels_
# print(cluster_pred)
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)


# agglomerative clustering with average linkage
print('agglomerative clustering with average linkage')
hierarchical = cluster.AgglomerativeClustering(n_clusters=8, linkage="average")
cluster_h = hierarchical.fit(cluster_data)
cluster_pred = cluster_h.labels_
# print(cluster_pred)
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)

# prediction by agglomerative clustering with average linkage
hierarchical_label=[]
hierarchical_test = hierarchical.fit_predict(predict_data)
# cluster_result = hierarchical_test.labels_
print('result of prediction')
accuracy2 = metrics.adjusted_rand_score(ground_true, hierarchical_test)
print(accuracy2)

# DBSCAN with eps =0.9 MinPt =10
print('DBSCAN')
dbscan = cluster.DBSCAN(eps=0.9, min_samples=10).fit(cluster_data)
cluster_pred = dbscan.labels_
accuracy1 = metrics.adjusted_rand_score(cluster_true, cluster_pred)
print(accuracy1)
