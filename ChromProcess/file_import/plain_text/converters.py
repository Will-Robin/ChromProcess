import numpy as np
from ChromProcess import Classes

def parse_text_columns(text, point_delimiter, ordinal_delimiter):
    '''
    Converts data organised in columns in a text file and converts it to a list
    of lists [[x values],[y values]].

    Parameters
    ----------
    text: str
        Text to be parsed.
    point_delimiter: str
        Delimiter for separating point coordintes. 
    ordinal_delimiter: str
        Delimiter for separating x,y pairs
    Returns
    -------
    data_container: dict
        Dictionary containing data scraper from the file.
    '''

    lines = text.split(point_delimiter)
    scrub_lines = [x for x in lines if x != '']
    points = [list(x.split(ordinal_delimiter)) for x in scrub_lines]
    data = list(map(list, zip(*points)))

    return data

def chrom_from_text(x_values, y_values, x_unit, y_unit, filename):
    '''
    Create a chromatogram object and insert time and signal
    information into it.

    Parameters
    ----------
    x_values: list
    y_values: list
    x_unit: str
    y_unit: str

    Returns
    -------
    chrom: ChromProcess Chromatogram object
    '''

    chrom = Classes.Chromatogram()
    chrom.filename = filename

    chrom.x_unit = x_unit
    chrom.y_unit = y_unit

    x_values = [float(x) for x in x_values]
    y_values = [float(y) for y in y_values]

    chrom.time = np.array(x_values)
    chrom.signal = np.array(y_values)

    return chrom
