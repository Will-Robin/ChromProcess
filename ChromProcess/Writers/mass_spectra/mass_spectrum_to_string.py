
def mass_spectrum_to_string(mass_spectrum):
    '''
    Write a row-wise csv for the mass spectrum.

    Parameters
    ----------
    mass_spectrum: Classes.MassSpectrum

    Returns
    -------
    ms_string: str
    '''

    ms_string = 'm/z,'
    for x in range(0,len(mass_spectrum.mz)):
        ms_string += f'{mass_spectrum.mz[x]},'

    ms_string += "\nrelative abundance,"

    for x in range(0,len(mass_spectrum.relative_abundances)):
        ms_string += f'{mass_spectrum.mz[x]}'

    ms_string += "\n"

    return ms_string
