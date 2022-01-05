
def write_to_csv(chromatogram, filename = 'chromatogram.csv'):
    '''
    Write chromatogram to a .csv file.
    '''

    time   = chromatogram.time
    signal = chromatogram.signal

    with open(filename, 'w') as f:
        f.write(f'{chromatogram.x_unit},{chromatogram.y_unit}\n')
    for x in range(0,len(time)):
        f.write(f'{time[x]},{signal[x]}\n')
