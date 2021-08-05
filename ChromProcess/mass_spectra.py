import os
import numpy as np
import matplotlib.pyplot as plt

from ChromProcess import Classes
from ChromProcess import info_params
from ChromProcess import mass_spectra as ms
from ChromProcess import simple_functions as s_f
from ChromProcess import series_operations as s_o

'''
Functions for dealing with mass spectra.
'''

def ion_chromatogram_region(chromatogram, lower, upper, threshold = 0.1):

    if len(chromatogram.scan_indices) == 0:
        pass
    else:
        ion_dict = {}
        inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

        time = chromatogram.time[inds]
        signal = chromatogram.signal[inds]
        scan_inds = chromatogram.scan_indices[inds]
        p_counts = chromatogram.point_counts[inds]

        for s in range(0,len(time)):

            inten = chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
            masses = chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

            if len(inten) > 0:

                filt_inds = np.where(inten > threshold*np.amax(inten))[0]
                inten = inten[filt_inds]
                masses = masses[filt_inds]

                round = np.round(masses, 2)

                for m in range(0,len(round)):
                    if round[m] in ion_dict:
                        ion_dict[round[m]][s] = inten[m]
                    else:
                        ion_dict[round[m]] = np.zeros(len(time))
                        ion_dict[round[m]][s] = inten[m]
            else:
                pass

    return time, ion_dict

def peak_ion_chromatograms(chromatogram, spectrum_filter = 0.1):

    if len(chromatogram.scan_indices) == 0:
        pass
    else:
        for p in chromatogram.peaks:

            inds = chromatogram.peaks[p].indices

            time = chromatogram.time[inds]
            signal = chromatogram.signal[inds]
            scan_inds = chromatogram.scan_indices[inds]
            p_counts = chromatogram.point_counts[inds]

            for s in range(0,len(time)):

                inten = chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
                masses = chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

                if len(inten) > 0:

                    filt_inds = np.where(inten > spectrum_filter*np.amax(inten))[0]
                    inten = inten[filt_inds]
                    masses = masses[filt_inds]

                    round = np.round(masses, 2)

                    for m in range(0,len(round)):
                        if round[m] in chromatogram.peaks[p].ion_chromatograms:
                            chromatogram.peaks[p].ion_chromatograms[round[m]][s] = inten[m]
                        else:
                            chromatogram.peaks[p].ion_chromatograms[round[m]] = np.zeros(len(time))
                            chromatogram.peaks[p].ion_chromatograms[round[m]][s] = inten[m]
                else:
                    pass

            ms.bin_ion_chromatograms(chromatogram.peaks[p], stdev = 0.1)

def bin_ion_chromatograms(peak, stdev = 0.001):

    sorted_masses = sorted([*peak.ion_chromatograms])

    clusters = []
    for c in s_o.cluster(sorted_masses, bound = stdev):
        clusters.append(c)

    out_log = {}
    for c in clusters:

        position = round(np.average(c),2)

        out_log[position] = []

        for o in range(0,len(sorted_masses)):
            if sorted_masses[o] in c:
                if len(out_log[position]) == 0:
                    out_log[position] = peak.ion_chromatograms[ sorted_masses[o] ]
                else:
                    for count,p in enumerate(out_log[position]):
                        if p == 0:
                            out_log[position][count] = peak.ion_chromatograms[ sorted_masses[o] ][count]
                        else:
                            pass

    peak.ion_chromatograms = out_log

def integrate_ion_chromatograms(chromatogram, threshold = 0.1):

    for p in chromatogram.peaks:
        inds = chromatogram.peaks[p].indices

        time = chromatogram.time[inds]

        for i in chromatogram.peaks[p].ion_chromatograms:
            if np.amax(chromatogram.peaks[p].ion_chromatograms[i]) > threshold*np.amax(chromatogram.signal[chromatogram.peaks[p].indices]):
                if chromatogram.internal_reference == False:
                    integral = np.trapz(chromatogram.peaks[p].ion_chromatograms[i], x = time)
                else:
                    integral = np.trapz(chromatogram.peaks[p].ion_chromatograms[i], x = time)/chromatogram.internal_reference.integral
                chromatogram.peaks[p].ion_integrals[i] = integral

def ion_chromatograms_relative_to_peak(series):

    best_peaks = []
    for p in series.peak_series:
        max_integral = 0
        best_peak = False
        for c,pk in enumerate(series.peak_series[p]):
            if pk.integral > max_integral:
                best_peak = (c,pk)
                max_integral = pk.integral
        if best_peak:
            best_peaks.append(best_peak)

    with open("ion_fractions.csv", "w") as f:
        f.write("peak, ion, fraction\n")
        for p in best_peaks:
            for i in p[1].ion_integrals:
                f.write("{},".format(p[1].retention_time))
                f.write("{},".format(i))
                f.write("{}".format(p[1].ion_integrals[i]/p[1].integral))
                f.write("\n")

def peak_character(series):
    '''
    Express peak identities as a function of residence times and ion chromatograms.
    '''

    for c in series.chromatograms:
        for p in c.peaks:

            characteristics = [str(c.peaks[p].retention_time), str(c.peaks[p].integral)]

            [characteristics.append(str(x)) for x in [*c.peaks[p].ion_integrals]]

            c.peaks[p].character = ";".join(characteristics)

def ion_chromatogram_integral_series(series):

    for p in series.peak_series:
        track_masses = {}
        for count, tp in enumerate(series.peak_series[p]):
            for m in tp.ion_integrals:
                if m in track_masses:
                    track_masses[m][count] = tp.ion_integrals[m]
                else:
                    track_masses[m] = [0 for x in series.x_series]
                    track_masses[m][count] = tp.ion_integrals[m]

        series.ion_series[p] = track_masses

def get_peak_mass_spectra(series):
    for c in series.chromatograms:
        if len(c.scan_indices) == 0:
            pass
        else:
            for p in c.peaks:
                ind = np.where(c.time == c.peaks[p].retention_time)[0]

                start = c.scan_indices[ind][0]
                end = start + c.point_counts[ind][0]

                c.peaks[p].mass_spectrum = [np.round(c.mass_values[start:end],2), c.mass_intensity[start:end]]

def write_mass_spectra(series):
    for count,c in enumerate(series.chromatograms):
        os.makedirs("{}".format(count), exist_ok = True)
        os.chdir("{}".format(count))

        for p in c.peaks:
            if c.peaks[p].mass_spectrum:

                masses = c.peaks[p].mass_spectrum[0]
                intensities = c.peaks[p].mass_spectrum[1]

                intensities = 100*intensities/np.amax(intensities)

                fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
                ax = plt.subplot(111)

                (markers, stemlines, baseline) = ax.stem(masses, intensities, use_line_collection = True)
                plt.setp(markers, marker=' ')
                plt.setp(baseline, visible=False)
                plt.setp(stemlines, color = "k")
                positions = []
                for x,y in zip(masses, intensities):
                    pos = x
                    if y > 8 and round(pos,1) not in positions:
                        ax.annotate('{}'.format(pos),
                                    xy=(x,y),
                                    xytext=(0*3, 3),  # use 3 points offset
                                    textcoords="offset points",  # in both directions
                                    ha='center', va='bottom',
                                    color = "b")
                        positions.append(round(pos,1))

                ax.set_xlabel("m/z")
                ax.set_ylabel("% abundance")
                ax.set_ylim(0,110)
                plt.savefig("peak_{}_{}.png".format(p, c.filename))
                plt.clf()
                plt.close()

                with open("peak_{}_{}.csv".format(p, c.filename), "w") as f:
                    f.write("m/z, intensity\n")
                    for m,i in zip(masses,intensities):
                        f.write("{},{}\n".format(m,i))

        os.chdir("..")

def mass_spectrum_peak_picker(series, mass_split = 5):
    for count,c in enumerate(series.chromatograms):
        os.makedirs("{}".format(count), exist_ok = True)
        os.chdir("{}".format(count))

        for p in c.peaks:
            masses = c.peaks[p].mass_spectrum[0]
            intensities = c.peaks[p].mass_spectrum[1]

            i = np.argsort(masses)
            masses = masses[i]
            intensities = intensities[i]

            clusters = []
            for clust in ms.cluster(masses, n = mass_split):
                clusters.append((clust).astype(int))

            mass_out = []
            for counter, clust in enumerate(clusters):
                om = masses[clust]
                max_ind = np.argmax(intensities[clust])
                mass_out.append(om[max_ind])

            with open("peak_maxima_{}_{}.csv".format(p, c.filename), "w") as f:
                f.write("m/z\n")
                for m in mass_out:
                    f.write("{}\n".format(m))

            mass_out = []

        os.chdir("..")

def bin_masses(all_masses, data_mat):
    clusters = []
    for c in s_f.cluster_indices(all_masses, bound = 0.1):
        clusters.append(c)

    binned_masses = np.zeros(len(clusters))
    binned_data = np.zeros((len(data_mat),len(clusters)+1))
    binned_data[:,0] = data_mat[:,0]

    for c,cl in enumerate(clusters):
        binned_masses[c] = np.average(all_masses[cl])
        shft = [x+1 for x in cl]
        sum_inten = np.sum(data_mat[:,shft], axis = 1)
        binned_data[:,c+1] = sum_inten

    # make sure mass spectra a reported as relative abundance
    for x in range(0,len(binned_data)):
        binned_data[x,1:] = binned_data[x,1:]/np.max(binned_data[x,1:])

    return binned_masses, binned_data
