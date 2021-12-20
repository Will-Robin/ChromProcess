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

    assert isinstance(text, str), "text arg should be a string"
    assert isinstance(point_delimiter, str), "point_delimiter arg should be str"
    assert isinstance(ordinal_delimiter, str), "ordinal_delimiter arg should be str"

    lines = text.split(point_delimiter)
    scrub_lines = [x for x in lines if x != '']
    points = [list(x.split(ordinal_delimiter)) for x in scrub_lines]
    data = list(map(list, zip(*points)))

    return data
