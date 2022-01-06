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
        import ChromProcess.Writers as write

        ms_string = write.mass_spectrum_to_string(self)

        return ms_string
