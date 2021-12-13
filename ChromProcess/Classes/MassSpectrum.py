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
