import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt
from ChromProcess import mass_spectra, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations
import pickle

def _1gaussian(x, amp1, cen1, sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2)))

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN088\B"
chromatogram_number = 4
region = [13.14,13.36]

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv" # 3/16 vs 2/18
IS_pos = 6.74 # 6.925/ 6.74

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 1000 # to cut off noise if needed.
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.
os.chdir(directory)
f_list = os.listdir()
c_names = []
for f in f_list:
    if ".cdf" in f:
        c_names.append(f)
    if "conditions" in f:
        cond_file = Classes.Information(f)

c_names.sort()
chrom = Classes.Chromatogram(c_names[chromatogram_number-1], mass_spec = use_mass_spectra)

import numpy as np
processing_functions.MS_intensity_threshold_chromatogram(chrom, threshold = 500)

inds = np.where((chrom.time > cond_file.internal_ref_region[0]) & (chrom.time < cond_file.internal_ref_region[1]))[0]
max = np.amax(chrom.signal[inds])
inds = chrom.signal == max
time_base = chrom.time[inds]
chrom.time = chrom.time - time_base
chrom.time = chrom.time + IS_pos

os.makedirs("inspect_regions", exist_ok = True)
os.chdir("inspect_regions")

inds = np.where((chrom.time > region[0])&(chrom.time< region[1]))
time, ion_dict = mass_spectra.ion_chromatogram_region(chrom,region[0], region[1])
ion_dict = series_operations.bin_dictionary(ion_dict, bound = 0.1, rounding = 2)

tensor = np.zeros(len(time))
masses = np.array([])
for i in ion_dict:
    if np.amax(ion_dict[i]) < 1000:
        pass
    elif int(i) == 73:
        pass
    else:
        tensor = np.vstack((tensor, ion_dict[i]))
        masses = np.hstack((masses,i))

tensor = tensor[1:]

from sklearn.decomposition import FastICA
import pylab as pl
S = tensor.T

S /= S.std(axis=0)
ica = FastICA(n_components = 5)
S_ = ica.fit(tensor).transform(tensor)  # Get the estimated sources
A_ = ica.mixing_  # Get estimated mixing matrix

pl.figure()
pl.subplot(3, 1, 1)
pl.plot(S)
pl.title('True Sources')
pl.subplot(3, 1, 2)
pl.plot(tensor.T)
pl.title('Observations (mixed signal)')
pl.subplot(3, 1, 3)
pl.plot(S_)
pl.title('ICA estimated sources')
pl.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.36)
pl.show()

quit()
roll_vals = np.arange(-(len(time))+1,len(time))
t_step = time[1]-time[0]
t_lags = np.zeros(len(roll_vals))
for r in range(0,len(roll_vals)):
    t_lags[r] = roll_vals[r]*t_step
t_lag_corr_mat = np.zeros((len(tensor), len(tensor), len(t_lags)))

for x in range(0,len(tensor)):
    x_i = tensor[x] - tensor[x].mean()
    for y in range(0,len(tensor)):
        x_j = tensor[y] - tensor[y].mean()
        con = np.correlate(x_i,x_j, mode = "full")
        ccor = con/(len(time)*x_i.std()*x_j.std())
        t_lag_corr_mat[x,y] = np.nan_to_num(ccor)

ind = np.where(t_lags == 0.0)[0]

corr_mat  = np.zeros(len(masses))
for x in range(0,len(t_lag_corr_mat)):
    corr_mat = np.vstack((corr_mat,  t_lag_corr_mat[:,x,ind].T))
corr_mat = corr_mat[1:]

from NetFit import data_processing
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

linkage_matrix = linkage(corr_mat, "single")

fl = fcluster(linkage_matrix,3,criterion='maxclust')

clusters = {}
for c,f in enumerate(fl):
    if f in clusters:
        clusters[f].append(tensor[c])
    else:
        clusters[f] = [tensor[c]]

c_seq = ["r","b","g"]
TIC = np.zeros(len(time))
for c,f in enumerate(clusters,1):
    sum = np.zeros(len(time))
    for Y in clusters[f]:
        sum += Y
        #plt.plot(time, Y)
    TIC = np.vstack((TIC,sum))

TIC = TIC[1:]
TIC = TIC.T

roll_vals = np.arange(-(len(TIC[0]))+1,len(TIC[0]))
t_step = TIC[0][1]-TIC[0][0]
m_lags = np.zeros(len(roll_vals))
for r in range(0,len(roll_vals)):
    m_lags[r] = roll_vals[r]*t_step
tlcm = np.zeros((len(TIC), len(TIC), len(m_lags)))

for x in range(0,len(TIC)):
    x_i = TIC[x] - TIC[x].mean()
    for y in range(0,len(TIC)):
        x_j = TIC[y] - TIC[y].mean()
        con = np.correlate(x_i,x_j, mode = "full")
        ccor = con/(len(time)*x_i.std()*x_j.std())
        tlcm[x,y] = np.nan_to_num(ccor)

X,Y = np.meshgrid(m_lags,time)
fig, ax = plt.subplots(ncols=2)
t_ic = TIC.T

plt.ion()
for c in range(0,len(t_ic)):
    ax[1].plot(time, t_ic[c])

for x in range(0,len(tlcm)):
    ax1 = ax[0].contourf(X,Y,tlcm[x], cmap = "seismic")
    ax[1].scatter([time[x] for x in range(0,len(TIC[x]))], TIC[x])
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.show()

quit()

tensor = tensor.T

roll_vals = np.arange(-(len(masses))+1,len(masses))
t_step = masses[1]-masses[0]
m_lags = np.zeros(len(roll_vals))
for r in range(0,len(roll_vals)):
    m_lags[r] = roll_vals[r]*t_step
t_lag_corr_mat = np.zeros((len(tensor), len(tensor), len(m_lags)))

for x in range(0,len(tensor)):
    x_i = tensor[x] - tensor[x].mean()
    for y in range(0,len(tensor)):
        x_j = tensor[y] - tensor[y].mean()
        con = np.correlate(x_i,x_j, mode = "full")
        ccor = con/(len(time)*x_i.std()*x_j.std())
        t_lag_corr_mat[x,y] = np.nan_to_num(ccor)


X,Y = np.meshgrid(m_lags,time)
X2,Y2 = np.meshgrid(masses,time )
for x in range(0,len(t_lag_corr_mat)):
    fig, ax = plt.subplots(nrows = 2)
    ax[1].contourf(X,Y,t_lag_corr_mat[x], zorder = 0 , cmap = "seismic")
    ax[0].contourf(X2,Y2,tensor, zorder = 0 , cmap = "seismic")
    ax[0].plot(masses,np.full(len(masses),time[x]),c = "k")
    plt.title
    plt.show()
    plt.close()
quit()
for x in range(0,len(t_lag_corr_mat)):
    fig = plt.figure()
    ax = fig.add_subplot(projection = "3d")
    for y in range(0,len(t_lag_corr_mat)):
        ax.plot(m_lags, np.full(len(m_lags), time[y]), t_lag_corr_mat[x,y], c = "k", alpha = 0.5)
    ax.set_xlabel("m/z")
    ax.set_label("time/ min")
    ax.set_label("corr")
    plt.show()
    plt.close()

quit()


pass_mat = np.zeros((len(masses), len(time)))
x_i = chrom.signal[inds] - chrom.signal[inds].mean()
for x in range(0,len(tensor)):
    x_j = tensor[x] - tensor[x].mean()
    con = np.correlate(x_i,x_j, mode = "same")
    ccor = con/(len(time)*x_i.std()*x_j.std())
    pass_mat[x] = np.nan_to_num(ccor)

X,Y = np.meshgrid(time, masses)
#plt.contourf(X,Y, pass_mat,cmap = "seismic", levels = 1000)
fig = plt.figure()
fig.add_subplot(projection = "3d")
for x in range(0,len(pass_mat)):
    plt.plot(time, np.full(len(time),masses[x]), pass_mat[x])
plt.show()
quit()

tensor = tensor.T

diff =np.zeros(len(masses))
for x in range(0,len(tensor)-1):
    diff = np.vstack((diff,tensor[x]-tensor[x+1]))
diff = diff[1:]
diff = diff.T


diff_stack = np.zeros(( len(masses),len(diff[0]) ))
for x in range(0,len(diff)):
    diff_stack[x] =  processing_functions.savitzky_golay(diff[x]/np.amax(diff[x]), 21, 3, deriv=0, rate=1)

X,Y = np.meshgrid( time[:-1], masses )

plt.contourf(X,Y, diff_stack, cmap = "seismic")
plt.show()

quit()

tensor = tensor.T

score_mat = np.zeros((len(tensor), len(tensor)))
for x in range(0,len(tensor)):
    x_i = tensor[x] - tensor[x].mean()
    for y in range(0,len(tensor)):
        x_j = tensor[y] - tensor[y].mean()
        con = np.correlate(x_i,x_j, mode = "valid")
        ccor = con/(len(time)*x_i.std()*x_j.std())
        score_mat[x,y] = np.nan_to_num(ccor)


fig,ax = plt.subplots(nrows = 1)
X,Y = np.meshgrid(time, time)

ax1 = ax.contourf(X,Y,score_mat,cmap = "seismic")
fig.colorbar(ax1, ax=ax)
ax.set_xlabel("time/ min")
ax.set_ylabel("scan time/ min")
plt.show()

quit()


fig = plt.figure()
fig.add_subplot(projection = "3d")
for x in range(0,len(score_mat),2):
    plt.plot(time, np.full(len(time),x),score_mat[x])
    #plt.plot(time, np.full(len(time),x),tensor[:,x]/np.amax(tensor[:,x]), c = "k")

#plt.plot(time,np.full(len(time),time[0]), chrom.signal[inds]/np.amax(chrom.signal[inds]), c = "k")
plt.show()
quit()


from sklearn.decomposition import PCA

pca = PCA(n_components=3)

pca.fit(tensor)

for x in range(0,len(tensor)):
    if x == 0:
        plt.plot(time, tensor[x], label = "signal", c = "k", alpha = 0.5)
    else:
        plt.plot(time, tensor[x], c = "k", alpha = 0.5)

count = 1
for length, vector in zip(pca.explained_variance_, pca.components_):
    v = vector * 3 * np.sqrt(length)
    plt.plot(time,pca.mean_ + v, label = "pca {} {}%".format(count, round(100*pca.explained_variance_ratio_[count-1],2)), linewidth = 3)
    count +=1

v1 = np.dot( pca.components_[0] * 3 * np.sqrt(pca.explained_variance_[0]),  pca.components_[2] * 3 * np.sqrt(pca.explained_variance_[2]) )
plt.plot(time,  pca.mean_ + v1, label = "{}".format("LC 1, 3"), linewidth = 3, c = "r")

plt.legend()
plt.show()

quit()

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
fig.add_subplot(projection = "3d")
for x in range(0,len(tensor)):
    plt.plot(time,  np.full(len(time),masses[x]),tensor[x])
plt.show()
quit()


import tensorly as tl
from tensorly.decomposition import parafac
from tensorly.decomposition import tucker

ts = tl.tensor(tensor)

weights, factors = tucker(ts, 10)

[print(f.shape)for f in factors]
print(len(time))

plt.plot(factors[0][:,0])
plt.plot(factors[0][:,1])

plt.show()
plt.close()
#plt.plot(chrom.time[inds], chrom.signal[inds], c = "k")
for n in range(0,len(factors[1][0])):
    plt.plot(time, factors[1][:,n])

plt.show()
