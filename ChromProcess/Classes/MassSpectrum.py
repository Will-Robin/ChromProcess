import sys

class MassSpectrum:
    def __init__(self, fname, mz, inten, pos = None):
        '''
        Parameters
        ----------
        fname: str or pathlib Path
        mz: array
        inten: array
        pos: None or float
        '''

        from pathlib import Path

        filename = fname
        if isinstance(fname, str):
            filename = Path(fname)
        elif isinstance(fname, Path):
            pass
        else:
            sys.exit('''''')

        self.filename = filename
        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

    def write_to_string(self):
        '''
        Write a row-wise csv for the mass spectrum.

        Returns
        -------
        ms_string: str
        '''
        ms_string = 'm/z,'
        for x in range(0,len(self.mz)):
            ms_string += f'{self.mz[x]},'

        ms_string += "\nrelative abundance,"

        for x in range(0,len(self.relative_abundances)):
            ms_string += f'{self.mz[x]}'
            
        ms_string += "\n"

        return ms_string
