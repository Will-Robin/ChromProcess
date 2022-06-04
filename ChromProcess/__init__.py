"""
# Welcome to the ChromProcess Documentation

For example usage of ChromProcess, see the `Tutorials` folder in the
repository. The examples given are one possible analysis workflow, and you can
build your own using ChromProcess as a library, which can be explored in these
docs.

Structures for data are found in ChromProcess.Classes. Using ChromProcess is
centred around the data contained in these objects. The general idea behind the
directory structure is that when one wants to *Load* a *chromatogram* from a
*text* file, the code can be found in `Loading/chromatogram/text`, and so on.

The data structures are all mutable, meaning the 'primary' data loaded into
them can be changed. It is best to consider source data files as 'read only'.
Once they are loaded into ChromProcess, they should be considered to have been
modified, even if no operations have been performed in the software. Data from
inside ChromProcess must be saved to new files, rather than regarded as primary
data. Keep primary data well organised and unmodified!

This documentation was generated from the ChromProcess source files using
[pdoc](https://pdoc.dev).
"""
