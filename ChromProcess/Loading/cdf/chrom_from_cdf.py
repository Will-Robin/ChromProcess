from pathlib import Path
from ChromProcess import Classes
from cdf_loading import load_from_cdf

JEOL_keys = [
   'error_log',
   'a_d_sampling_rate',
   'a_d_coaddition_factor',
   'scan_acquisition_time',
   'scan_duration',
   'inter_scan_time',
   'resolution',
   'actual_scan_number',
   'total_intensity',
   'mass_range_min',
   'mass_range_max',
   'time_range_min',
   'time_range_max',
   'scan_index',
   'point_count',
   'flag_count',
   'mass_values',
   'intensity_values',
   'instrument_name',
   'instrument_id',
   'instrument_mfr',
   'instrument_model',
   'instrument_serial_no',
   'instrument_sw_versi on',
   'instrument_fw_version',
   'instrument_os_version',
   'instrument_app_version',
   'instrument_comments'
]

class JEOL:
    TIME_CONVERSION = 60
    TIME_KEY = 'scan_acquisition_time'
    TIC_KEY = 'total_intensity'
    MS_INTENSITY_KEY = 'intensity_values'
    MZ_KEY = 'mass_values' 
    SCAN_INDEX_KEY = 'scan_index'
    POINT_COUNT_KEY = 'point_count'

def chrom_from_cdf(filename, load_ms = False):
    '''
    Parameters
    ----------
    filename: str or pathlib Path

    Returns
    -------
    chrom: ChromProcess Chromatogram
    '''

    if isinstance(filename, str):
        fname = Path(filename)
    else:
        fname = filename

    data_keys = [JEOL.TIME_KEY, JEOL.TIC_KEY]

    if load_ms is True:
        ms_extensions = [
                        JEOL.MS_INTENSITY_KEY, 
                        JEOL.MZ_KEY,
                        JEOL.MS_INTENSITY_KEY,
                        JEOL.SCAN_INDEX_KEY,
                        JEOL.POINT_COUNT_KEY
                        ]

        data_keys.extend(ms_extensions)

    data_container = load_from_cdf(fname, data_keys)

    chrom = Classes.Chromatogram()

    chrom.time = data_container[JEOL.TIME_KEY]/JEOL.TIME_CONVERSION
    chrom.signal = data_container[JEOL.TIC_KEY]

    chrom.filename = fname.name

    chrom.x_unit = 'total ion counts'
    chrom.y_unit = 'min.'

    if load_ms is True:
        chrom.mz_values = data_container[JEOL.MZ_KEY]
        chrom.mz_intensity = data_container[JEOL.MS_INTENSITY_KEY]
        chrom.scan_indices = data_container[JEOL.SCAN_INDEX_KEY]
        chrom.point_counts = data_container[JEOL.POINT_COUNT_KEY]

    return chrom

