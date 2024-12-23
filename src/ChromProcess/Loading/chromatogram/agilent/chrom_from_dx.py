"""
File contents.

A .dx file is a zipped archive containing a set of files similar to the list
below:

```
263d4215-846a-4338-8581-9284bb43a0a5.CH
3b5c4044-e34c-43ee-bd51-1b0df7b2063a.CH
[Content_Types].xml
a0eb515a-2b3e-45fd-8a20-625671b8912a.CH
bc602eb9-65a5-4129-af95-7db006b4ebf4.CH
injection.acmd
```

`[Content_Types].xml` gives an overview of the contents of each filetype.

`injection.acmd` is an xml formatted file containing information about the data
contained in the .CH files.

Files with a .CH suffix contain data traces (those which are recorded as
specified in the GC Method, such as the detector trace and oven temperature).
They are 'flat' binaries (a long list of ones and zeros) which contain several
data types side-by-side in the sequences (segment for string-based information,
segment for detector traces, and so on). Each of these components can be
snipped out of the file with knowledge of their offset (where they begin in the
binary). The offsets can be determined using a binary editor:

- utf-16 encoded values can be sliced out and re-encoded into strings. These
  strings are length-prefixed, meaning that their fist byte indidates the
  length of the content string. This info can be used to cut out the string
  from the bytes. There may be some null characters left over in the strings,
  but they can be removed.
- >i4, >f4, >f8 (> means 'big endian') encoded values can be decoded using
  numpy.frombuffer. Numpy takes care of the details here in terms of clipping
  the data out of the container bytes.

|  Offset  | Encoding |
|----------|----------|
| 0x0      | "utf-16" |
| 0xfc     | ">H"     |
| 0x100    | ">H"     |
| 0x108    | ">i4"    |
| 0x11a    | ">f4"    |
| 0x11e    | ">f4"    |
| 0x146    | "utf-16" |
| 0x15b    | "utf-16" |
| 0x35a    | "utf-16" |
| 0x559    | "utf-16" |
| 0x758    | "utf-16" |
| 0x957    | "utf-16" |
| 0x9bc    | "utf-16" |
| 0x9e5    | "utf-16" |
| 0xa0e    | "utf-16" |
| 0xc11    | "utf-16" |
| 0xe11    | "utf-16" |
| 0xeda    | "utf-16" |
| 0xfd7    | "utf-16" |
| 0x104c   | "utf-16" |
| 0x1075   | "utf-16" |
| 0x1274   | ">f8"    |
| 0x127c   | ">d"     |
| 0x1800   | ">f"     |

Credit goes to the authors of yadg(1)

References
----------

(1) Kraus et al., (2022). yadg: yet another datagram. Journal of Open Source
    Software, 7(72), 4166, https://doi.org/10.21105/joss.04166
    https://github.com/dgbowl/yadg

(2) https://github.com/chemplexity/chromatography

(3) https://github.com/evanyeyeye/rainbow
"""

import zipfile
import numpy as np
from typing import Any, Union
from pathlib import Path
from ChromProcess.Classes import Chromatogram

sections = {
    "encoding": (0x0, "utf-16"),
    "data_offset": (0x108, ">i4"),
    "time_start": (0x11A, ">f4"),
    "time_end": (0x11E, ">f4"),
    "version_string": (0x146, "utf-16"),
    "description": (0x15B, "utf-16"),
    "sample": (0x35A, "utf-16"),
    "username": (0x758, "utf-16"),
    "date": (0x957, "utf-16"),
    "method_file": (0xA0E, "utf-16"),
    "units": (0x104C, "utf-16"),
    "channel_name": (0x1075, "utf-16"),
    "intercept": (0x1274, ">f8"),
    "slope": (0x127C, ">d"),
    "data_start": (0x1800, ">f"),
}


def decode_utf(binary_data: bytes, dtype: str) -> str:
    """
    Decode utf-8 or utf-16 encoded bytes to a string.

    Parameters
    ----------
    binary_data: bytes
        Data to be decoded.
    dtype: str
        Type of data.

    Returns
    -------
    str
    """
    # Find the length in bytes from the beginning of the string
    length_byte = binary_data[0:1]
    bytes_length = int.from_bytes(length_byte, byteorder="big")
    if dtype == "utf-16":
        # UTF-16 uses 16 bytes per character, not 8
        bytes_length *= 2
    # Cut out the string bytes and return as a re-encoded string
    string_bytes = binary_data[1 : bytes_length + 1]
    new_string = string_bytes.decode("utf-8")
    return new_string.replace("\0", "")


def read_header(
    data: bytes,
    offset: int,
    dtype: str,
) -> Union[str, np.ndarray]:
    """
    Decode a section of a string of bytes using an offset and expected data
    type.

    Parameters
    ----------
    data: bytes
        Data containing section to be decoded.
    offset: int
        Where the data starts in the bytes.
    dtype: str
        The expected type of data.

    Returns
    -------
    Union[str , np.ndarray]
    """

    sliced = data[offset:]

    if "utf-" in dtype:  # If the data are a string
        return decode_utf(sliced, dtype)

    else:
        return np.frombuffer(sliced, dtype=dtype, count=1)


def chrom_from_ch(file_bytes: bytes) -> Chromatogram:
    """
    Create a chromatogram object from the byte contents of an Agilent .ch file.

    """
    info: dict[str, Any] = {}
    for f in sections:
        info[f] = read_header(file_bytes, sections[f][0], sections[f][1])

    start = sections["data_start"][0]
    end = len(file_bytes)
    npoints = (end - start) // 8

    time = (
        np.linspace(info["time_start"][0], info["time_end"][0], num=npoints) / 60000
    )  # converted to minutes

    signal = (
        np.frombuffer(
            file_bytes,
            offset=start,
            dtype="<f8",
            count=npoints,
        )
        * info["slope"][0]
    )

    chrom = Chromatogram()

    chrom.time = time
    chrom.signal = signal

    chrom.x_unit = "minutes"
    chrom.y_unit = info["units"]

    return chrom


def chrom_from_dx(
    filename: Union[str, Path], trace_name: str = "FID1A"
) -> Chromatogram:
    """
    Load a chromatogram from an Agilent .dx file.

    Parameters
    ----------
    filename: Union[str, Path]

    Returns
    -------
    chrom: Chromatogram
    """
    # Unzip file to get hold of contents
    # Each file is read to find out what its contents are. Alternatively, the
    # information can be gleaned from injection.acmd, but parsing a load of XML
    # might be slower than just peeking at bytes.
    trace_name = trace_name.encode("utf-16").decode("utf-16")
    with zipfile.ZipFile(filename) as file:
        filenames = [x for x in file.namelist() if x.endswith("CH")]
        for name in filenames:
            content_bytes = file.read(name)
            trace = read_header(
                content_bytes, sections["channel_name"][0], sections["channel_name"][1]
            )
            if trace == trace_name:
                chrom = chrom_from_ch(content_bytes)
                if isinstance(filename, Path):
                    filename = filename.name
                else:
                    filename = Path(filename).name

                chrom.filename = filename

    return chrom
