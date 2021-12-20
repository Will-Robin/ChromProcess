from ChromProcess import Classes

def load_conditions_from_csv(filename):
    '''
    Parameters
    ----------
    filename: str or pathlib Path
    '''

    from ChromProcess.Utils import simple_functions as s_f

    conditions = Classes.ExperimentConditions()

    with open(filename, 'r', encoding = 'latin-1') as f:
        for line in f:
            if "Dataset" in line:
                ins = line.strip("\n")
                conditions.experiment_code = ins.split(",")[1]
            if "series_values" in line:
                ins = line.strip("\n")
                spl = ins.split(",")
                conditions.series_values =  [float(x) for x in spl[1:] if x != ""]
            if "series_unit" in line:
                ins = line.strip("\n")
                conditions.series_unit = ins.split(",")[1]

    # Read conditions
    condset = []
    readstate = False
    with open(filename, "r", encoding = 'latin-1') as f:
        for c,line in enumerate(f):
            if "start_conditions" in line:
                readstate = True
                line = next(f)
            if "end_conditions" in line:
                readstate = False
            if readstate:
                newline = line.strip("\n")
                condset.append([x for x in newline.split(",") if x != ""])
    c_out = {}
    for c in condset:
        c_out[c[0]] = []
        for x in c[1:]:
            if s_f.isfloat(x):
                c_out[c[0]].append(float(x))
            else:
                c_out[c[0]].append(x)

    conditions.conditions = c_out

    return conditions

