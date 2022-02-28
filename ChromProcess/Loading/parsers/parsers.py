def parse_text_columns(text, point_delimiter, ordinal_delimiter):
    """
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
    data: list
    """

    assert isinstance(text, str), "text arg should be a string"
    assert isinstance(point_delimiter, str), "point_delimiter arg should be str"
    assert isinstance(ordinal_delimiter, str), "ordinal_delimiter arg should be str"

    lines = text.split(point_delimiter)
    scrub_lines = [x for x in lines if x != ""]
    points = [list(x.split(ordinal_delimiter)) for x in scrub_lines]
    data = list(map(list, zip(*points)))

    return data


def import_file_section(file, start_token, end_token):
    """
    Load a section of lines between two tokens.

    Parameters
    ----------
    file: str or pathlib Path
        Path to the file.
    start_token: str
        String in line to start reading file from.
    end_token:
        String in line to end reading file from.

    Returns
    -------
    c_set: list
    """

    spl_lin = lambda x: [e for e in x.strip("\n").split(",") if e != ""]
    readstate = False
    c_set = []
    with open(file, "r", encoding="latin-1") as f:
        for _, line in enumerate(f):
            if start_token in line:
                readstate = True
                line = next(f)
            if end_token in line:
                readstate = False
            if readstate:
                newline = spl_lin(line)
                c_set.append(newline)

    return c_set
