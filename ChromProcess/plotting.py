import os
import numpy as np
import matplotlib.pyplot as plt # plotting

from ChromProcess import info_params
from ChromProcess import processing_functions
from ChromProcess import deconvolution as d_c


def plot_conc_series(series, bound_dict):
    '''
    Plot a time series of integrals for a Chromatogram_Series.

    Parameters
    ----------

    Returns
    -------

    '''

    out_unit = " M"

    y_data = series.conc_series
    x_data = np.array([x.timepoint for x in series.chromatograms])

    name = '{}_conc_series'.format(series.set_name)

    fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
    ax = plt.subplot(111)

    for p in series.conc_series:

        label_name = processing_functions.name_peak(float(p),bound_dict)
        if label_name in info_params.colour_assignments:
            plot_col = info_params.colour_assignments[label_name]
        else:
            plot_col = "0.75"

        ax.plot(x_data,series.conc_series[p], "-o", label = "{} ({})".format(label_name,p), c = plot_col)

    ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
    ax.set_xlabel("{}".format(series.x_name), fontsize = info_params.font)
    ax.set_ylabel("{}/ {}".format("concentration", out_unit), fontsize = info_params.font)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              fancybox=True, shadow=True, ncol=3, fontsize = info_params.legend)
    plt.savefig('{}.png'.format(name))
    plt.close()
    plt.clf()

def plot_integral_series(series, bound_dict):
    '''
    Plot a time series of integrals for a Chromatogram_Series.

    Parameters
    ----------

    Returns
    -------

    '''

    out_unit = "I$_p$/I$_I$$_S$"

    y_data = series.integral_series
    x_data = np.array([x.timepoint for x in series.chromatograms])

    name = '{}_integral_series'.format(series.set_name)

    fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
    ax = plt.subplot(111)

    for p in series.integral_series:

        label_name = processing_functions.name_peak(float(p),bound_dict)

        if label_name in info_params.colour_assignments:
            plot_col = info_params.colour_assignments[label_name]
        else:
            plot_col = "0.75"

        ax.plot(x_data,series.integral_series[p], "-o", label = "{} ({})".format(label_name,p), c = plot_col)

    ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
    ax.set_xlabel("{}".format(series.x_name), fontsize = info_params.font)
    ax.set_ylabel("{}".format("I$_S$/I$_I$$_S$"), fontsize = info_params.font)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              fancybox=True, shadow=True, ncol=4, fontsize = info_params.legend)
    plt.savefig('{}.png'.format(name))
    plt.close()
    plt.clf()

def plot_chromatograms(series, information):
    '''
    Plot the chromatograms for a Chromatogram_Series
    Parameters
    ----------

    Returns
    -------

    '''
    for v,b in enumerate(series.regions):
        plt.figure(figsize = (info_params.across,info_params.up))
        ax = plt.subplot(111)

        for z,c in enumerate(series.chromatograms):
            region_inds = np.where((c.time > b[0]) & (c.time < b[1]))[0]
            time = c.time[region_inds]
            signal = c.signal[region_inds]

            for p in c.peaks:
                if len(c.peaks[p].indices) == 0:
                    pass

                elif b[0] < p < b[1]:
                    ind = np.where(time == p)[0]

                    ax.scatter(time[ind],signal[ind], marker =  "x", s= 20, c = "k", zorder = 2)
                    ax.scatter(c.time[c.peaks[p].indices[0]],c.signal[c.peaks[p].indices[0]], marker =  "|", s= 100, c = "k", zorder = 2)
                    ax.scatter(c.time[c.peaks[p].indices[-1]],c.signal[c.peaks[p].indices[-1]], marker =  "|", s= 200, c = "k", zorder = 2)

            ax.plot(time,signal, label ='{} {}'.format(series.x_series[z], series.x_name.split("/")[1]), linewidth  = info_params.lines, zorder = 1)


        for g in information.boundaries:
            region_inds = np.where((c.time > b[0]) & (c.time < b[1]))[0]
            time = c.time[region_inds]
            signal = c.signal[region_inds]
            if information.boundaries[g][1] < b[1] and information.boundaries[g][0] > b[0]:

                lower_bound_y = np.linspace(np.amin(signal),np.amax(signal), num = 100)
                lower_bound_x = np.full(len(lower_bound_y),information.boundaries[g][0])

                upper_bound_y = np.linspace(np.amin(signal),np.amax(signal), num = 100)
                upper_bound_x = np.full(len(upper_bound_y),information.boundaries[g][1])

                ax.scatter(lower_bound_x, lower_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                ax.scatter(upper_bound_x, upper_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                ax.text(np.mean(information.boundaries[g]),np.amax(signal),g, fontsize = 10, rotation = 90, rotation_mode = "anchor")

        ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
        ax.set_xlabel("time/ min", fontsize = info_params.font)
        ax.set_ylabel("Intensity", fontsize = info_params.font)
        ax.set_xlim(time[0],time[-1])
        box = ax.get_position()
        ax.set_position([box.x0 + box.width*0.2, box.y0 + box.height * 0.2, box.width*0.8, box.height * 0.8])
        # Put a legend below current axis
        '''
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                  fancybox=True, shadow=True, ncol=3, fontsize = info_params.labels-6)
        '''
        plt.savefig('{}_{}_to_{}_{}_chromatogram_inferred_integrals.png'.format(c.filename,b[0],b[1],z))
        plt.close()
        plt.clf()

def region_heatmap(series, information):
    '''
    Parameters
    -----------


    Returns
    -------

    '''

    min_length = 1e100
    for c in series.chromatograms:
        if len(c.time) < min_length:
            min_length = len(c.time)

    heat_stack = np.empty((len(series.chromatograms),len(series.chromatograms[0].signal[:min_length])))


    for c in range(0,len(series.chromatograms)):
        heat_stack[c] = series.chromatograms[c].signal[:min_length]

    for r,b in enumerate(series.regions):

        region_inds = np.where((series.chromatograms[0].time > b[0]) & (series.chromatograms[0].time < b[1]))[0]
        fig, ax = plt.subplots(figsize=(15,10))

        X = series.chromatograms[0].time[region_inds]
        Y = series.x_series
        zi = heat_stack[:,region_inds]

        # create x-y points to be used in heatmap
        xi, yi = np.meshgrid(X, Y)
        CS = ax.contourf(xi, yi, zi, 200, cmap=plt.cm.rainbow)

        for g in information.boundaries:
            if information.boundaries[g][1] < b[1] and information.boundaries[g][0] > b[0]:

                lower_bound_y = np.linspace(Y[0],Y[-1], num = 100)
                lower_bound_x = np.full(len(lower_bound_y),information.boundaries[g][0])

                upper_bound_y = np.linspace(Y[0],Y[-1], num = 100)
                upper_bound_x = np.full(len(upper_bound_y),information.boundaries[g][1])

                ax.scatter(lower_bound_x, lower_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                ax.scatter(upper_bound_x, upper_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                ax.text(np.mean(information.boundaries[g]),np.amin(Y),g, fontsize = 10, rotation = 90, rotation_mode = "anchor")

        for p in series.chromatograms[-1].peaks:
            if b[0] < p < b[1]:
                ax.annotate(round(p,3), xy = (p,np.amax(Y)), horizontalalignment = 'center')

        cbar = plt.colorbar(mappable = CS)
        cbar.set_label('signal', rotation=270, fontsize = 8, labelpad = 20)
        ax.tick_params(labelsize = info_params.font-2, axis = "both")
        ax.set_xlabel("retention time/ min", fontsize = info_params.font)
        ax.set_ylabel("{}".format(series.x_name), fontsize = info_params.font)
        #ax.set_title("{} concentrations".format(n))
        #plt.legend()
        plt.savefig("{}_{}_contour_plot.png".format(b[0],b[1]))
        plt.clf()
        plt.close()

    return series.chromatograms[0].time[:min_length], heat_stack

def plot_peak_ion_chromatogram(peak, chromatogram, name):
    '''
    '''

    fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
    ax = plt.subplot(111)
    if len(peak.ion_chromatograms) == 0:
        pass
    else:
        #ax.plot(series.chromatograms[c].time[peak.indices],series.chromatograms[c].signal[peak.indices], linewidth = info_params.lines, c = "k", alpha = 0.5)

        for i in peak.ion_chromatograms:

            if int(i) in info_params.frag_colours:
                ax.plot(chromatogram.time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines, c = info_params.frag_colours[int(i)])
            else:
                ax.plot(chromatogram.time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines)

            if np.amax(peak.ion_chromatograms[i]) > 0.1*peak.integral:
                ax.annotate(i, xy = (chromatogram.time[peak.indices][np.argmax(peak.ion_chromatograms[i])], np.amax(peak.ion_chromatograms[i])), fontweight = "bold", horizontalalignment = "center", fontsize = 12)

        ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
        ax.set_xlabel("time/ min.", fontsize = info_params.font)
        ax.set_ylabel("Ion count", fontsize = info_params.font)

        plt.savefig('{}_{}_ion_chromatogram.png'.format(name,round(peak.retention_time,3)))
        plt.close()
        plt.clf()


def plot_peak_ion_chromatogram_from_series(series, information, cutoff = 0.1):

    best_peaks = []
    for p in series.peak_series:
        max_integral = 0
        best_peak = False
        for c,pk in enumerate(series.peak_series[p]):
            if pk.integral > max_integral and pk.integral > cutoff:
                best_peak = (c,pk)
                max_integral = pk.integral
        if best_peak:
            best_peaks.append(best_peak)

    for c,peak in best_peaks:

        fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
        ax = plt.subplot(111)
        if len(peak.ion_chromatograms) == 0:
            pass
        else:
            for i in peak.ion_chromatograms:

                if int(i) in info_params.frag_colours:
                    ax.plot(series.chromatograms[c].time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines, c = info_params.frag_colours[int(i)])
                else:
                    ax.plot(series.chromatograms[c].time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines)

                if np.amax(peak.ion_chromatograms[i]) > 0.1*peak.integral:
                    ax.annotate(i, xy = (series.chromatograms[c].time[peak.indices][np.argmax(peak.ion_chromatograms[i])], np.amax(peak.ion_chromatograms[i])), fontweight = "bold", horizontalalignment = "center", fontsize = 12)

            ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
            ax.set_xlabel("time/ min.", fontsize = info_params.font)
            ax.set_ylabel("Ion count", fontsize = info_params.font)

            peak_name = processing_functions.name_peak(peak.retention_time, information.boundaries)
            if peak_name in info_params.class_assignments:
                class_assign = info_params.class_assignments[peak_name]
            else:
                class_assign = peak_name

            plt.savefig('{}_{}_ion_chromatograms.png'.format(class_assign, round(peak.retention_time,3)))
            plt.close()
            plt.clf()

def write_peak_ion_chromatogram_from_series(series):

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

    for c,peak in best_peaks:

        output = series.chromatograms[c].time[peak.indices]

        for i in peak.ion_chromatograms:
            output = np.vstack((output,peak.ion_chromatograms[i]))

        output = output.T

        with open('{}_ion_chromatograms.csv'.format(round(peak.retention_time,3)), "w") as f:
            f.write("time/ min,")
            for i in peak.ion_chromatograms:
                f.write("{},".format(i))

            f.write("\n")

            for x in range(0,len(output)):
                for y in range(0,len(output[x])):
                    f.write("{},".format(output[x,y]))

                f.write("\n")



def plot_ion_channels(series):

    mass_track = {}
    for p in series.peak_series:
        trace = series.peak_series[p]

        for pk in trace:
            for x in [*pk.ion_integrals]:
                mass_track[x] = 1

    peaks = sorted([*series.peak_series])

    for m in [*mass_track]:
        fig, ax = plt.subplots(figsize=(15,10))

        if len(peaks) < 2:
            grid = np.zeros((len(series.x_series), 2))
            xaxis = [peaks[0],peaks[0]+1]

        else:
            grid = np.zeros((len(series.x_series), len(peaks)))
            xaxis = peaks

        for p_count,p in enumerate(peaks):
            trace = series.peak_series[p]

            for t_count, pk in enumerate(trace):
                for mass in pk.ion_integrals:
                    if round(mass,1) == round(m,1):
                        grid[t_count, p_count] = 100*pk.ion_integrals[mass]/pk.integral
                    else:
                        pass

        X = xaxis
        Y = series.x_series
        zi = grid

        # create x-y points to be used in heatmap
        xi, yi = np.meshgrid(X, Y)
        CS = ax.contourf(xi, yi, zi, 200, cmap=plt.cm.rainbow, zorder = 0)
        cbar = plt.colorbar(mappable = CS)
        cbar.set_label('percentage of peak integral', rotation=270, fontsize = 8, labelpad = 20)
        ax.tick_params(labelsize = info_params.font-2, axis = "both")
        ax.set_xlabel("retention time/ min", fontsize = info_params.font)
        ax.set_ylabel("{}".format(series.x_name), fontsize = info_params.font)
        ax.set_title("mz_{}_tracking".format(round(m,1)))
        plt.savefig("mz_{}_tracking.png".format(round(m,1)))
        plt.clf()
        plt.close()

def plot_2D_grid(mass_spectra_dict,x2,y2,filename):
    '''
    Plots the total intensity chromatogram and an 2D scatter of time vs m/z
    with intensities represented by a colormap.

    Parameters
    ----------

    Returns
    -------

    '''
    bounds = [0,float([*mass_spectra_dict][-1])]
    fig   = plt.figure(figsize = (10,10))
    gs    = gridspec.GridSpec(10,10)
    ax1   = fig.add_subplot(gs[3:-1,:]) # mass spectra
    ax2   = fig.add_subplot(gs[0:3,:-2]) # chromatogram
    for k in mass_spectra_dict:
        if bounds[0] < float(k) and float(k) < bounds[1]:

            indices = np.where(mass_spectra_dict[k]['Mass'] > 100)[0]
            x = [float(k) for m in mass_spectra_dict[k]['Mass'][indices]]

            if len(x) == 0:
                pass
            else:
                y = mass_spectra_dict[k]['Mass'][indices]
                z = mass_spectra_dict[k]['Intensity'][indices]
                sc = ax1.scatter(x, y, c = z/np.amax(z), alpha = 1, cmap = "viridis", s=20)
        else:
            pass

    ax2.plot(x2,y2)
    ax2.set_xlim(bounds[0],bounds[1])
    ax1.set_xlim(bounds[0],bounds[1])
    plt.tight_layout()
    ax1.set_ylabel("mass")
    #ax2.set_ylabel("total intensity")
    ax1.set_xlabel("time")
    plt.colorbar(sc,ax=ax1,pad = 0.054)
    plt.savefig("out_{}.png".format(filename),dpi = 900)

    return 'This is for making a graph!'

def plot_mass_spectra(series):
    '''
    Plotting mass spectra
    Parameters
    ----------
    series

    Returns
    -------

    '''

    os.makedirs("Mass_spectra", exist_ok=True)
    os.chdir("Mass_spectra")

    for c in series.chromatograms:
        for x,v in enumerate(series.regions):
            region_inds = np.where((c.time > v[0]) & (c.time < v[1]))[0]
            inds = c.peaks[x]['Peak_indices']
            newt = c.time[region_inds]
            peak_list = newt[inds]
            for p in peak_list:
                fig = plt.figure(figsize = (info_params.across,info_params.up) )
                ax = plt.subplot(111)

                ax.bar(c.mass_spectra[str(round(p,3))]["Mass"],100*c.mass_spectra[str(round(p,3))]["Intensity"]/max(c.mass_spectra[str(round(p,3))]["Intensity"]), color = 'k')
                ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
                ax.set_xlabel("m/z", fontsize = info_params.font)
                ax.set_ylabel("% abundance", fontsize = info_params.font)
                plt.savefig('mass_spectrum_{}.png'.format(str(round(p,3))))
                plt.clf()
                plt.close()

    os.chdir("..")

def plot_mass_spectra_in_time(mass_spectra_dict):
    '''
    Generates a plot of how the mass spectra change over time.

    Parameters
    ----------

    Returns
    -------

    '''

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111, projection='3d')
    for d in mass_spectra_dict:
        for k in mass_spectra_dict[d]:
            if 14 < float(k) and float(k) < 20:
                indices = np.where(mass_spectra_dict[d][k]['Mass'] > 300)[0]
                time = [float(d) for m in mass_spectra_dict[d][k]['Mass'][indices]]
                retention_time = [float(k) for m in mass_spectra_dict[d][k]['Mass'][indices]]
                if len(retention_time) == 0:
                    pass
                else:
                    masses = mass_spectra_dict[d][k]['Mass'][indices]
                    intensities = mass_spectra_dict[d][k]['Intensity'][indices]
                    ax.scatter(time, retention_time, masses, s= 2, c = intensities/np.amax(intensities), alpha = 0.2)
                    #ax.scatter(y, x, time, s= 2, c = z/np.amax(z), alpha = 0.1, cmap = "viridis")
            else:
                pass

    ax.set_xlabel('time/ min')
    ax.set_ylabel('retention time/ min')
    ax.set_zlabel('m/z')

    print("making graph")
    plt.savefig('test_cloud_plot.png')

def plot_mass_spectra_in_time_2(mass_spectra_dict1,mass_spectra_dict2):
    '''
    Alternate function for plotting mass spectra over time.
    Parameters
    ----------

    Returns
    -------

    '''

    fig = plt.figure(figsize = (10,10))
    ax1 = fig.add_subplot(111, projection='3d')
    for d in mass_spectra_dict1:
        for k in mass_spectra_dict1[d]:
            if 14 < float(k) and float(k) < 20:
                indices = np.where(mass_spectra_dict1[d][k]['Mass'] > 300)[0]
                time = [float(d) for m in mass_spectra_dict1[d][k]['Mass'][indices]]
                retention_time = [float(k) for m in mass_spectra_dict1[d][k]['Mass'][indices]]
                if len(retention_time) == 0:
                    pass
                else:
                    masses = mass_spectra_dict1[d][k]['Mass'][indices]
                    intensities = mass_spectra_dict1[d][k]['Intensity'][indices]
                    ax1.scatter(time, retention_time, masses, s= 2, c = intensities/np.amax(intensities), alpha = 0.2)
                    #ax.scatter(y, x, time, s= 2, c = z/np.amax(z), alpha = 0.1, cmap = "viridis")
            else:
                pass

    for d in mass_spectra_dict2:
        for k in mass_spectra_dict2[d]:
            if 14 < float(k) and float(k) < 20:
                indices = np.where(mass_spectra_dict2[d][k]['Mass'] > 300)[0]
                time = [float(d) for m in mass_spectra_dict2[d][k]['Mass'][indices]]
                retention_time = [float(k) for m in mass_spectra_dict2[d][k]['Mass'][indices]]
                if len(retention_time) == 0:
                    pass
                else:
                    masses = mass_spectra_dict2[d][k]['Mass'][indices]
                    intensities = mass_spectra_dict2[d][k]['Intensity'][indices]
                    ax1.scatter(time, retention_time, masses, s= 2, c = intensities/np.amax(intensities), alpha = 0.2, cmap = 'magma')
                    #ax.scatter(y, x, time, s= 2, c = z/np.amax(z), alpha = 0.1, cmap = "viridis")
            else:
                pass

    ax1.set_ylabel('retention time/ min')
    ax1.set_zlabel('m/z')
    ax1.set_xlabel('time/ min')

    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False

    ax1.xaxis.pane.set_edgecolor('w')
    ax1.yaxis.pane.set_edgecolor('w')
    ax1.zaxis.pane.set_edgecolor('w')

    ax.grid(False)

    print("making graph")
    plt.savefig('test_cloud_plot.png')

def plot_3D_scatter(mass_spectra_dict,x2,y2,filename):
    '''
    A 3D scatter plot

    Parameters
    ----------

    Returns
    -------

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for k in mass_spectra_dict:
        if 7 < float(k) and float(k) < 20:
            indices = np.where(mass_spectra_dict[k]['Mass'] > 100)[0]
            x = [float(k) for m in mass_spectra_dict[k]['Mass'][indices]]

            if len(x) == 0:
                pass
            else:
                y = mass_spectra_dict[k]['Mass'][indices]
                z = mass_spectra_dict[k]['Intensity'][indices]

                ax.scatter(x, y, z, s= 5, color = 'b', alpha = 0.3)
        else:
            pass

    ax.set_xlabel('time/ min')
    ax.set_ylabel('m/z')
    ax.set_zlabel('intensity')

    print("making graph")
    plt.savefig('test_3D_plot.png')

def plot_series_internal_standard_log(series, attribute, bound_dict):
    '''
    Plot a time series of integrals for a Chromatogram_Series.

    Parameters
    ----------

    Returns
    -------

    '''

    out_unit = " unknown unit"

    if "integral" in attribute:
        out_unit = " integral"

    if "conc" in attribute:
        out_unit = " M"

    y_data = getattr(series,attribute)
    x_data = np.array([x.timepoint for x in series.chromatograms])

    for c,r in enumerate(y_data):
        fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
        ax = plt.subplot(111)

        for p in y_data[r]:

            label_name = processing_functions.name_peak(float(p),bound_dict)
            name = '{}'.format(label_name)
            if label_name in info_params.colour_assignments:
                plot_col = info_params.colour_assignments[label_name]
            else:
                plot_col = "0.75"

            ax.plot(x_data,y_data[r][p], "-o", label = "{} ({})".format(label_name,p), c = plot_col)

        ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
        ax.set_xlabel("{}".format(series.x_name), fontsize = info_params.font)
        ax.set_ylabel("{}/ {}".format(attribute, out_unit), fontsize = info_params.font)
        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                  fancybox=True, shadow=True, ncol=3, fontsize = info_params.labels-4)
        plt.savefig('peak_timecourses_{}.png'.format(name))
        plt.close()
        plt.clf()

def plot_chromatogram_with_assignments(chromatogram,series,bound_dict, labelling = False):
    '''
    Parameters
    -----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram to be plotted.
    series:  ChromProcess Chromatogram_Series object
        just for the regions

    Returns
    -------
    '''
    if chromatogram.c_type == "GCMS":
        C_chain_regions = {"tetradecane" : [6.6,7.05],
                           "C3": [5.2,6.6],
                           "C4": [7.5,10],
                           "C5": [10,12],
                           "C6": [12,18],
                           "C7": [18,25]}

        C_chain_colors = { "tetradecane" : "k",
                           "C3": "#f09c08",
                           "C4": "#2738e7",
                           "C5": "#cb340b",
                           "C6": "#30bd37",
                           "C7": "#592387"}
    else:
        C_chain_regions = {"sugars" : [0.5,0.9],
                           "DHA": [0.95,1.1],
                           "DNPH": [1.1,1.4],
                           "HCHO": [2,2.5]}

        C_chain_colors = {"sugars" : "#2738e7",
                           "DHA": "#f09c08",
                           "DNPH": "#ff7800",
                           "HCHO": "b"}

    if series.internal_ref_region:
        internal_ref_height = np.amax(chromatogram.signal[np.where((chromatogram.time > series.internal_ref_region[0])&(chromatogram.time < series.internal_ref_region[1]))[0]])
    else:
        internal_ref_height = 1

    fig = plt.figure(figsize = (10,5))
    ax = plt.subplot(111)

    for c in C_chain_regions:
        inds = np.where((chromatogram.time > C_chain_regions[c][0])&(chromatogram.time < C_chain_regions[c][1]))[0]
        ax.plot(chromatogram.time[inds],chromatogram.signal[inds]/internal_ref_height, c = C_chain_colors[c])

    if labelling:
        for p in chromatogram.peaks:
            label_name = processing_functions.name_peak(float(p),bound_dict)
            ind = np.where(chromatogram.time == float(p))[0]
            if label_name.split(" ")[0] in info_params.numbering_scheme:
                plt.annotate(info_params.numbering_scheme[label_name.split(" ")[0]], xy = (float(p), 0.01+chromatogram.signal[ind]/internal_ref_height), fontweight = "bold", horizontalalignment = "center", fontsize = 12)
            else:
                plt.annotate(label_name, xy = (float(p), 0.01+chromatogram.signal[ind]/internal_ref_height),fontsize = 8, fontweight = "bold", horizontalalignment = "left", verticalalignment = "bottom", rotation = 45)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xlabel("time/ min", fontsize = 14)
    ax.set_ylabel("signal/ tetradecane peak height", fontsize = 14)
    plt.savefig("{}_chromatogram.png".format(chromatogram.filename))
    plt.close()
    plt.clf()

def region_heatmap_with_ion_chroms(series, information, cutoff = 0.1):
    '''
    Parameters
    -----------


    Returns
    -------

    '''

    # find best peaks
    best_peaks = []
    for p in series.peak_series:
        max_integral = 0
        best_peak = False
        for c,pk in enumerate(series.peak_series[p]):
            if pk.integral > max_integral and pk.integral > cutoff:
                best_peak = (c,pk)
                max_integral = pk.integral
        if best_peak:
            best_peaks.append(best_peak)




    min_length = 1e100
    for c in series.chromatograms:
        if len(c.time) < min_length:
            min_length = len(c.time)

    heat_stack = np.empty((len(series.chromatograms),len(series.chromatograms[0].signal[:min_length])))


    for c in range(0,len(series.chromatograms)):
        heat_stack[c] = series.chromatograms[c].signal[:min_length]

    for r,b in enumerate(series.regions):

        region_inds = np.where((series.chromatograms[0].time > b[0]) & (series.chromatograms[0].time < b[1]))[0]
        fig, ax = plt.subplots(figsize=(15,10))

        X = series.chromatograms[0].time[region_inds]
        Y = series.x_series
        zi = heat_stack[:,region_inds]

        # create x-y points to be used in heatmap
        xi, yi = np.meshgrid(X, Y)
        CS = ax.contourf(xi, yi, zi, 200, cmap=plt.cm.rainbow)

        for g in information.boundaries:
            if information.boundaries[g][1] < b[1] and information.boundaries[g][0] > b[0]:

                lower_bound_y = np.linspace(Y[0],Y[-1], num = 100)
                lower_bound_x = np.full(len(lower_bound_y),information.boundaries[g][0])

                upper_bound_y = np.linspace(Y[0],Y[-1], num = 100)
                upper_bound_x = np.full(len(upper_bound_y),information.boundaries[g][1])

                ax.scatter(lower_bound_x, lower_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                ax.scatter(upper_bound_x, upper_bound_y, marker =".", c= "k", alpha = 0.5, edgecolors = None, s = 2)
                c_a = info_params.class_assignments[g]
                ax.text(np.mean(information.boundaries[g]),np.amin(Y),'{} ({})'.format(g,c_a), fontsize = 10, rotation = 90, rotation_mode = "anchor")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height * 0.9])

        for c,peak in best_peaks:
            if b[0] < peak.retention_time < b[1]:
                left = series.chromatograms[c].time[peak.indices][0]
                right = series.chromatograms[c].time[peak.indices][-1]
                height = (np.amax(Y)-np.amin(Y))/10
                bottom = np.amax(Y)
                ion_ax = ax.inset_axes([left, bottom, right-left, height], transform=ax.transData)

                if len(peak.ion_chromatograms) == 0:
                    pass
                else:
                    for i in peak.ion_chromatograms:
                        if int(i) in info_params.frag_colours:
                            ion_ax.plot(series.chromatograms[c].time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines, c = info_params.frag_colours[int(i)])
                        else:
                            ion_ax.plot(series.chromatograms[c].time[peak.indices], peak.ion_chromatograms[i], linewidth = info_params.lines)

                        if np.amax(peak.ion_chromatograms[i]) > 0.1*peak.integral:
                            ion_ax.annotate(i, xy = (series.chromatograms[c].time[peak.indices][np.argmax(peak.ion_chromatograms[i])], np.amax(peak.ion_chromatograms[i])),
                                            fontweight = "bold", horizontalalignment = "center", fontsize = 6)

                    ion_ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
                    #ion_ax.set_ylabel("Ion count", fontsize = info_params.font/2)
                    ion_ax.set_xticklabels([])
                    ion_ax.set_yticklabels([])
                    ion_ax.set_axis_off()

        cbar = plt.colorbar(mappable = CS)
        cbar.set_label('signal', rotation=270, fontsize = 8, labelpad = 20)
        ax.tick_params(labelsize = info_params.font-2, axis = "both")
        ax.set_xlabel("retention time/ min", fontsize = info_params.font)
        ax.set_ylabel("{}".format(series.x_name), fontsize = info_params.font)

        plt.savefig("{}_{}_contour_plot.png".format(b[0],b[1]))
        plt.clf()
        plt.close()

def confidence_ellipse(x, y, ax = None, weights = None, n_std=3.0,
                        facecolor='none', edgecolor = 'k',**kwargs):
    """
    Adapted from matplotlib documentation.

    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    import numpy as np
    from ChromProcess import plotting

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y, aweights = weights)

    ellipse = plotting.create_ellipse((0, 0), cov,
                              facecolor = facecolor,
                              edgecolor = edgecolor,
                              **kwargs)


    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def confidence_ellipse_from_mean_cov(mean, cov, ax = None, weights = None, n_std=3.0,
                        facecolor='none', edgecolor = 'k',**kwargs):
    """
    Adapted from matplotlib documentation.

    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    import numpy as np
    from ChromProcess import plotting

    ellipse = plotting.create_ellipse((0, 0), cov,
                              facecolor = facecolor,
                              edgecolor = edgecolor,
                              **kwargs)


    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mean[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def create_ellipse(mean, cov,
                  facecolor = 'k',edgecolor = 'k', **kwargs):

    from matplotlib.patches import Ellipse

    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = Ellipse(mean, width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, edgecolor = edgecolor, **kwargs)

    return ellipse
