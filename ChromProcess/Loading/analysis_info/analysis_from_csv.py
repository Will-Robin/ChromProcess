from ChromProcess import Classes

def analysis_from_csv(fname):
    '''
    Parameters
    ----------
    fname: str or pathlib Path
    '''

    analysis = Classes.AnalysisInformation()

    rdlin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
    with open(fname, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if 'Dataset' in line:
            ins = rdlin(line)
            analysis.experiment_code = ins[1]

        if 'Method' in line:
            ins = rdlin(line)
            analysis.analysis_type = ins[1]

        if 'regions' in line:
            ins = rdlin(line)
            reg = [float(x) for x in ins[1:]]
            analysis.regions = [reg[x:x+2] for x in range(0, len(reg), 2)]

        if 'internal_standard_region' in line:
            ins = rdlin(line)
            reg = [float(x) for x in ins[1:]]
            analysis.internal_ref_region = reg

        if 'extract_mass_spectra' in line:
            ins = rdlin(line)
            use_ms = ins[1].lower()
            if ins[1] == 'true':
                analysis.use_MS = True
            if ins[1] == 'false':
                analysis.use_MS = False

        if 'mass_spectra_filter' in line:
            ins = rdlin(line)
            analysis.MS_cutoff = float(ins[1])

        if 'peak_pick_threshold' in line:
            ins = rdlin(line)
            analysis.peak_pick_threshold = float(ins[1])

        if 'dilution_factor,' in line:
            ins = rdlin(line)
            analysis.dilution_factor = float(ins[1])

        if 'dilution_factor_error' in line:
            ins = rdlin(line)
            analysis.dilution_factor_error = float(ins[1])

        if "internal_ref_concentration," in line:
            ins = rdlin(line)
            analysis.internal_ref_concentration = float(ins[1])

        if 'internal_ref_concentration_error' in line:
            ins = rdlin(line)
            analysis.internal_ref_concentration_error = float(ins[1])

    return analysis

