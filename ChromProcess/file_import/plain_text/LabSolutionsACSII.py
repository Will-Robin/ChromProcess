import re
import numpy as np

start_time_fieldname = 'Start Time'
end_time_fieldname = 'End Time'
data_start_fieldname = 'R.Time (min)'
wavelength_fieldname = 'Wavelength(nm)'
units_fieldname = 'Intensity Units'

block_regex = f'(?<={start_tag}).*?(?={end_tag})'
field_regex = f'^{field_name}'

def load_from_text_blocks(filename):
    '''
    A specific parser for .txt files exported from Shimadzu LabSolutions
    software.

    Extracts data from the chromatogram file into a dictionary.

    Parameters
    ----------
    file: str or pathlib Path
        name of chromatogram file
        (ASCII exported from Shimadzu LC software)
    dat_dict: dict
        dictionary of information from chromatogram file.

    Returns
    -------
    data_container: dict
        Dictionary containing data scraper from the file.
    '''

    with open(filename, 'r') as file:
        text = file.read()

    data_container = {}
    regexes = []
    for reg in regexes:
        matches = re.findall(reg, text)
        data_container[reg] = matches
        
    return data_container

