
class ExperimentConditions:
    def __init__(self):
        '''
        Parameters
        ----------
        information_file: str or pathlib Path
        '''

        self.experiment_code = 'not specified'
        self.series_values = []
        self.series_unit = 'not specified'
        self.conditions = {}

    def write_to_file(self, directory = ''):
        '''
        Parameters
        ----------
        directory: str or pathlib Path
        '''

        from pathlib import Path

        fname = self.experiment_code
        exp_code = self.experiment_code
        if isinstance(directory, str):
            if directory == '':
                fname = f'{exp_code}_conditions.csv'
            else:
                fname = f'{directory}/{exp_code}_conditions.csv'

        elif isinstance(directory, Path):
            fname = directory/f'{exp_code}_conditions.csv'
        else:
            pass

        with open(fname, 'w') as f:
            f.write(f'Dataset,{exp_code}\n')
            f.write('start_experiment_information\n')

            f.write('series_values,')
            for s in self.series_values:
                f.write(f'{s},')

            f.write('\n')
            f.write(f'series_unit,{self.series_unit}\n')

            f.write('end_experiment_information\n')
            f.write("start_conditions\n")
            for c in self.conditions:
                f.write(f"{c},")
                [f.write(f"{x},") for x in self.conditions[c]]
                f.write("\n")
            f.write("end_conditions\n")

