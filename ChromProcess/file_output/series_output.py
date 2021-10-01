from ChromProcess import info_params as i_p

def data_report_template_convert(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''
    name = "{}_{}_{}_data_report.csv".format(series.set_name,"conc_series",series.chromatograms[0].c_type)
    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = [series.x_series]
        p_header = [series.x_name]

        for p in series.conc_series:
            label_name = processing_functions.name_peak(float(p),information.boundaries)
            if processing_functions.int_test(label_name[-1]):
                label_name = label_name[:-2]

            species_smiles = i_p.canonical_SMILES[label_name]
            p_header.append(species_smiles+" M")
            out.append(series.conc_series[p])

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")

def peak_ion_chromatograms_report(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''

    for s in series.ion_series:

        name = "{}_{}_peak_report.csv".format(series.set_name,s)

        with open(name, "w") as outfile:
            outfile.write("Dataset,{}\n".format(series.set_name))
            outfile.write("Peak,{}\n".format(s))
            outfile.write("start_conditions\n")
            for c in series.conditions:
                outfile.write("{},".format(c))
                [outfile.write("{},".format(x)) for x in series.conditions[c]]
                outfile.write("\n")
            outfile.write("end_conditions\n")

            outfile.write("start_analysis_details\n")

            outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
            outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
            outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

            outfile.write("end_analysis_details\n")

            outfile.write("start_data\n")

            # create ion integral traces
            p_header = [series.x_name]
            out = [series.x_series]

            for m in series.ion_series[s]:
                p_header.append(m)
                out.append(series.ion_series[s][m])

            out = [list(i) for i in zip(*out)]

            [outfile.write("{},".format(x)) for x in p_header]
            outfile.write("\n")
            for o in out:
                [outfile.write("{},".format(x)) for x in o]
                outfile.write("\n")

            outfile.write("end_data\n")

def report_all_peak_integrals(series, information):

    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''

    y_data = series.integral_series
    x_data = np.array([x.timepoint for x in series.chromatograms])
    out_unit = " integral"
    name = '{}'.format("{}_{}_all_integrals_report.csv".format(series.set_name,series.chromatograms[0].c_type))

    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        outfile.write("dilution_factor, {}\n".format(series.dilution))
        outfile.write("internal_ref_concentration/ M, {}\n".format(series.internal_ref_concentration))
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = x_data
        p_header = [series.x_name]

        for p in series.integral_series:
            label_name = processing_functions.name_peak(float(p),information.boundaries)
            out = np.vstack((out,series.integral_series[p]))
            p_header.append("{}({})  {}".format(label_name,p,out_unit))

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")

def report_integral_log(series, information):
    '''
    Output integral log report.
    Parameters
    ----------

    Returns
    -------

    '''

    y_data = series.integral_series
    x_data = np.array([x.timepoint for x in series.chromatograms])
    out_unit = " integral"
    name = '{}'.format("{}_internal_reference_report.csv".format(series.set_name))

    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))
        outfile.write("Start_Date,{}\n".format(i_p.dates[series.set_name]))
        out = x_data
        p_header = [series.x_name]

        for c,r in enumerate(y_data):
            for p in y_data[r]:
                label_name = processing_functions.name_peak(float(p),information.boundaries)
                out = np.vstack((out,y_data[r][p]))
                p_header.append("{}({})  {}".format(label_name,p,out_unit))
        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for c,o in enumerate(out):
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

    return 0

def write_deconvoluted_output(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''
    name = "{}_deconvolution_report.csv".format(series.set_name)
    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = [series.x_series]
        p_header = [series.x_name]

        for p in series.deconvoluted_series:
            print(p)
            p_header.append("{} integral".format(p))
            out.append(series.deconvoluted_series[p])

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")
