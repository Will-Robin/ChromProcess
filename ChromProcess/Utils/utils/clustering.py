import numpy as np

def cluster(values, bound = 0.1):

    values = np.sort(values)

    cluster = []
    for m in range(0,len(values)):
        if len(cluster) > 0:
            clust_av = np.average(cluster)

            if abs(values[m]-clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(values[m])

    yield cluster

def cluster_indices(values, bound = 0.1):

    sortedvalues = np.sort(values)

    cluster = []
    for m in range(0,len(sortedvalues)):
        if len(cluster) == 0:
            pass
        else:
            clust_av = np.average(sortedvalues[cluster])

            if  abs(sortedvalues[m]-clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(m)

    yield cluster
