
class DataPath:
    '''
    Contains information paths to data files.
    '''
    def __init__(self, exp_code, data_type, path):
        '''
        Parameters
        ----------
        exp_code: str
        data_type: str
        path: str or pathlib Path
        '''

        from pathlib import Path

        self.experiment_code = exp_code
        self.data_type = data_type

        if isinstance(path, str):
            self.path = Path(path)
        if isinstance(path, Path):
            self.path = path

