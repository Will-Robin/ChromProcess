import sys
sys.path.append(r'C:\Users\willi\Documents/Packages')
from ChromProcess import Classes, os, file_import, plt, series_builder, plotting, file_output

folder = r"F:\2020_05_18" # You  can paste the folder path here, but it must be preceded by an 'r' (in Windows)
os.chdir(folder) # changes the working directory
drlist = os.listdir() # get list of files in directory
drlist.sort() # sorts the files

'''Please make sure that your files are named so that they are read in the right
order to be aligned with the time or concentration series array. Otherwise the
automatic plotting/ series generation will not work properly!'''

chromatograms = []
for f in drlist:
    if f.endswith(".txt"): # this is how the programs finds the data files: .txt is assumed to be a HPLC file, .cdf is assumed to be GCMS
        print(f)
        chromatograms.append(Classes.Chromatogram(f))
    if "conditions" in f:
        cond_file = f

n= 1
fig,ax= plt.subplots(nrows = len(chromatograms), sharex=True)
for count,c in enumerate(chromatograms):
     ax[count].plot(c.time, c.signal, label = n)
     ax[count].set_title("{}".format(c.filename))
     n+=1
plt.legend()
plt.show()
plt.close()
