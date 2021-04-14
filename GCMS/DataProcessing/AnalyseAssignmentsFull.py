from pathlib import Path
from ChromProcess import info_params
from ChromProcess import plotting
import matplotlib.pyplot as plt
from NorthNet import data_processing
import numpy as np
import pickle
import math

import umap
from sklearn.preprocessing import minmax_scale
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


# database root folder
database_root = Path(r'C:\Users\willi\Documents\PrebioticDatabase')
collated_folder = database_root/'AggregatedData'

# Import data
with open(collated_folder/'All_mass_spectra', 'rb') as f:
    data_mat = pickle.load(f)

with open(collated_folder/'mass_spectra_labels', 'rb') as f:
    compound_labels = pickle.load(f)

with open(collated_folder/'all_masses', 'rb') as f:
    all_masses = pickle.load(f)

with open(collated_folder/'integral_mat', 'rb') as f:
    integral_mat = pickle.load(f)

# Scale the data ( (x - mean) / stdev)
scaled_data_mat = minmax_scale(data_mat)

clrs = np.array([info_params.colour_assignments[x] for x in compound_labels])

accuracies = []
retention_times = {}
mass_spectra = {}

# get data according to assigned carbon chain length
local_rts = data_mat[:,0]

# local colour scheme and compound labels
indiv_comps = list(set(comps))

# Number of features to search for in the data is expected to be
# the number of assigned compounds
n_clust = len(indiv_comps)

# reduce dimentions of data (hopefully removing noise)
pca = PCA(n_components = int(len(compound_labels)/3))
dim_red = pca.fit_transform(scaled_data_mat)

exp_var_tot = np.cumsum(pca.explained_variance_ratio_)
plt.scatter(np.linspace(1,len(exp_var_tot), len(exp_var_tot)), exp_var_tot, c = 'k')
plt.savefig('Explained_variance_full_assigns.png', dpi = 600)
plt.show()
plt.close()

# try to spread out the data a bit more using LDA
lda = LinearDiscriminantAnalysis(n_components=n_clust-1,store_covariance=True)
X_r2 = lda.fit(dim_red, compound_labels).transform(dim_red)

# build correlation matrix as an estimate of data similarity
corr_mat = np.corrcoef(X_r2)

# Choose data to cluster
clustering_basis = X_r2

# Perform clustering
from sklearn.cluster import KMeans
cluster_tag = 'KMeans'
clustering = KMeans(n_clusters=n_clust).fit(clustering_basis,
                                            sample_weight = integral_mat)

n_lab = clustering.labels_

# Create UMAP layout of correlation matrix section
reducer = umap.UMAP(n_components = n_clust)
X_embedded = reducer.fit_transform(clustering_basis)

# Plot the results
sq_rt = np.sqrt(n_clust)
square = math.ceil(sq_rt)

fig,ax = plt.subplots(figsize = (2*square, 2*square),
                        sharey = True, sharex = True)
cluster_assigns = {}
for z in range(n_clust):
    inds = np.where(n_lab == z)[0]
    if len(inds) == 0:
        continue

    loc_comps = [compound_labels[i] for i in inds]
    lib = {l:loc_comps.count(l) for l in loc_comps}
    lib_counts = [lib[l] for l in lib]
    lib_clrs = [info_params.colour_assignments[cmp] for cmp in lib]

    com_lab = max(set(loc_comps), key = loc_comps.count)

    max_c = info_params.smiles_to_names[com_lab]

    cluster_assigns[z] = com_lab

    if com_lab in retention_times:
        retention_times[com_lab] = np.hstack((retention_times[com_lab], local_rts[inds]))
    else:
        retention_times[com_lab] = local_rts[inds]

    accuracy = 100*(loc_comps.count(com_lab)/len(loc_comps))
    accuracies.append(accuracy)

    X = X_embedded[inds,0]
    Y = X_embedded[inds,1]

    ax.scatter(X,Y, c = clrs[inds],
                  zorder = 1, s= 2, alpha = 0.3)

    ax.annotate('{}\n{}%'.format(max_c, np.round(accuracy,2)), xy = (X.mean(), Y.mean()),
                horizontalalignment = 'center', fontsize = 6,
                fontweight = 'bold')

    plotting.confidence_ellipse(X_embedded[inds,0], X_embedded[inds,1], ax,
                        n_std=3.0, weights = integral_mat[inds],
                        facecolor=info_params.colour_assignments[com_lab],
                        alpha = 0.5, zorder = 0)

ax.set_xlabel('UMAP 1', fontsize = 6)
ax.set_ylabel('UMAP 2', fontsize = 6)
ax.set_xticklabels([])
ax.set_yticklabels([])

plt.savefig(collated_folder/'Full_assignment_clustering.png', dpi = 600)
plt.close()

print('Overall accuracy: ({}+-{})%'.format(np.mean(accuracies),np.std(accuracies, ddof = 1)))

for r in retention_times:
    print(info_params.smiles_to_names[r], np.mean(retention_times[r]))
    plt.hist(retention_times[r], bins = len(retention_times[r]), color = info_params.colour_assignments[r])

plt.savefig(collated_folder/'retention_time_histogram.png', dpi = 600)
