
class CalibrationAllocations:
    '''
    Contains information about which calibration file is assigned to which data
    set.
    '''
    def __init__(self,filename):

        '''
        Parameters
        ----------
        filename: str or pathlib Path
        '''

        self.experiments = []
        self.GCMS_allocations = {}
        self.HPLC_allocations = {}
        self.GC_allocations = {}

        if filename == '':
            pass
        else:
            self.import_file(filename)

    def import_file(self,filename):
        '''
        Parameters
        ----------
        filename: str or pathlib Path
        '''

        rdln = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        GCMS_idx = 0
        HPLC_idx = 0
        GC_idx = 0
        with open(filename, 'r') as f:
            for c, line in enumerate(f):
                if c == 0:
                    header = rdln(line)
                    GCMS_idx = header.index('GCMS_Calibration_file')
                    HPLC_idx = header.index('HPLC_Calibration_file')
                    GC_idx = header.index('GC_Calibration_file')
                else:
                    ins = rdln(line)
                    self.experiments.append(ins[0])
                    self.GCMS_allocations[ins[0]] = ins[GCMS_idx]
                    self.HPLC_allocations[ins[0]] = ins[HPLC_idx]
                    self.GC_allocations[ins[0]] = ins[GC_idx]

