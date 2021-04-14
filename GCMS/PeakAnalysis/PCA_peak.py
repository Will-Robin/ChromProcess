import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from ChromProcess import Classes
from ChromProcess import mass_spectra
from ChromProcess import file_import as f_i
from ChromProcess import chromatogram_operations as chrom_ops
from ChromProcess import deconvolution
from ChromProcess import simple_functions as simp_func

directory = Path(r'C:\Users\willi\Documents\Data\GCMS\FRN\FRN107\A')
region = [10.52,10.8]
ic_thresh = 0.2
features = 10

chrom = Classes.Chromatogram(str(directory/'FRN107_002.cdf'), mass_spec = True)
idx = np.where((chrom.time > region[0])&(chrom.time < region[1]))[0]

ic_dict = chrom_ops.getIonChromatogramsFromRegion(chrom,
                                                  region[0], region[1],
                                                  threshold = ic_thresh)
mass_ax = np.array([*ic_dict])
ic_stack = np.array(list(ic_dict.values()))
ic_stack = np.vstack((ic_stack, chrom.time[idx]))

corr_mat = np.corrcoef(ic_stack)

results = simp_func.runPCA(corr_mat, n_components = features)

from sklearn.cluster import KMeans

# perform k_means clustering on dimension-reduced data
k_means = KMeans(n_clusters=features, init='k-means++', n_init=10, max_iter=300,
                 tol=0.00001, precompute_distances='deprecated', verbose=0,
                 random_state=42, copy_x=True, n_jobs='deprecated',
                 algorithm='elkan')

kmc = k_means.fit(corr_mat)
labels = k_means.labels_
cluster_colour_map = ['#DDCC77', '#CC6677', '#E69F00', '#AA4499', '#117733',
                      '#332288', '#88CCEE', '#000000', '#882255', '#44AA99']

for x in range(0,len(ic_stack)):
    ind = np.where(labels == x)[0]
    for i in ind:
        plt.plot(chrom.time[idx], ic_stack[i],c=cluster_colour_map[x])
plt.show()
