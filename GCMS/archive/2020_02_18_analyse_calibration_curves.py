import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, np
from ChromProcess import curve_fit
from NetFit import Classes as NClass
from NetFit import file_in_out as N_io

def quadratic(x, A,B):
    return A*(x**2) + B*x

os.chdir(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations")

paths = {}
with open("Calibration_paths.csv", "r") as f:
    next(f)
    for line in f:
        l = line.strip("\n")
        l = l.split(",")

        paths[l[0]] = [x for x in l[1:] if x != ""]

calibration_parameters = {"Sugar":["A","B","C"]}

for p in paths:

    datasets = []
    print(p)
    for address in paths[p]:
        os.chdir(r"{}/Data_reports".format(address))
        f_lst = os.listdir()
        for f in f_lst:
            if "all_integrals" in f and not "._" in f:
                #create a dataset
                print(address)
                datasets.append( N_io.get_data(f, x_axis_key = "[sample]/[internal standard]",  flow_data = False) )

    conc_extract = np.array([])
    peaks = {}
    coords = []
    for d in datasets:
        del_list = []
        A = 0
        for dep in d.dependents:
            if float(dep.split("(")[0]) > A:
                A = float(dep.split("(")[0])

        for dep in [*d.dependents]:
            if float(dep.split("(")[0]) == A:
                d.dependents["A"] = d.dependents[dep][:]
            else:
                d.dependents["B"] = d.dependents[dep][:]

            del_list.append(dep)

        for r in del_list:
            del d.dependents[r]


    for d in datasets:

        conc_extract = np.hstack((conc_extract,d.time))

        for dep in d.dependents:
            if dep.split("(")[0] in peaks:
                for x,y in zip(d.time,d.dependents[dep]):
                    peaks[dep.split("(")[0]].append((x,y))
            else:
                peaks[dep.split("(")[0]] = []
                for x,y in zip(d.time,d.dependents[dep]):
                    peaks[dep.split("(")[0]].append((x,y))

    for p in peaks:
        peaks[p] = sorted(peaks[p])

    os.chdir(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\Calibration_results")
    for c,d in enumerate(peaks,1):
        conc_ax = np.array([x for x,_ in peaks[d]])
        integral_ax  = np.array([x for _,x in peaks[d]])
        popt,pcov = curve_fit(quadratic, conc_ax, integral_ax, p0 = [0.5,0.5])
        perr = np.sqrt(np.diag(pcov))
        p_high = popt + perr
        p_low = popt - perr
        upper = quadratic(conc_ax,*p_high)# upper bound on fit
        lower = quadratic(conc_ax,*p_low) # lower bound on fit

        #plt.fill_between(conc_ax,upper, lower, facecolor = "k", alpha = 0.3)
        plt.scatter(conc_ax, integral_ax)
        plt.plot(conc_ax, quadratic(conc_ax,*popt), c = "k") # best fit
        print(popt)
        calibration_parameters["{} {}".format(datasets[0].name.replace("calibration", "").replace("_",""),c)] = popt[:2]

    plt.title(datasets[0].name)
    plt.savefig("{}.png".format(datasets[0].name))
    plt.clf()


print(calibration_parameters)
os.chdir(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations")
os.makedirs("Calibration_results", exist_ok = True)
os.chdir("Calibration_results")

with open("calibration_results.csv", "w") as f:
    for c in calibration_parameters:
        f.write("{},{},{},{}".format(c, calibration_parameters[c][0],calibration_parameters[c][1],0))
        f.write("\n")
