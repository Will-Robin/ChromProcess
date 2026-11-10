import numpy as np
from ChromProcess.Classes import Chromatogram


def chrom_from_rainbow(datafile, load_ms: bool = True) -> Chromatogram:
    fid = next(df for df in datafile.datafiles if df.detector == "FID")

    chromatogram = Chromatogram()

    chromatogram.time = np.asarray(fid.xlabels)

    chromatogram.signal = np.asarray(fid.data[:, 0])

    chromatogram.filename = datafile.path

    chromatogram.x_unit = "minutes"
    chromatogram.y_unit = "pA"

    if load_ms:
        ms = next(df for df in datafile.datafiles if df.detector == "MS")

        mz_values = []
        mz_intensity = []
        scan_indices = []
        point_counts = []

        offset = 0

        for spectrum in ms.data:
            scan_indices.append(offset)

            mz_values.extend(ms.ylabels)
            mz_intensity.extend(spectrum)

            point_counts.append(len(ms.ylabels))

            offset += len(ms.ylabels)

        chromatogram.ms_time = np.asarray(ms.xlabels, dtype=float)

        chromatogram.mz_values = np.asarray(mz_values, dtype=float)

        chromatogram.mz_intensity = np.asarray(mz_intensity, dtype=float)

        chromatogram.scan_indices = np.asarray(scan_indices, dtype=int)

        chromatogram.point_counts = np.asarray(point_counts, dtype=int)

    return chromatogram
