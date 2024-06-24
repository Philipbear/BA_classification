"""
this is to control the spectral quality of bile acid spectra
"""
import pandas as pd
from matplotlib import pyplot as plt


def generate_library_df(library_mgf):
    """
    Generate metadata dataframe for the mgf file
    """
    with open(library_mgf, 'r') as file:
        spectrum_list = []
        for line in file:
            # empty line
            _line = line.strip()  # remove leading and trailing whitespace
            if not _line:
                continue
            elif line.startswith('BEGIN IONS'):
                spectrum = {}
                # initialize spectrum
                mz_list = []
                intensity_list = []
            elif line.startswith('END IONS'):
                if len(mz_list) == 0:
                    continue
                spectrum['mz_ls'] = mz_list
                spectrum['intensity_ls'] = intensity_list
                spectrum_list.append(spectrum)
                continue
            else:
                # if line contains '=', it is a key-value pair
                if '=' in _line:
                    # split by first '='
                    key, value = _line.split('=', 1)
                    spectrum[key] = value
                else:
                    # if no '=', it is a spectrum pair
                    this_mz, this_int = _line.split()
                    try:
                        mz_list.append(float(this_mz))
                        intensity_list.append(float(this_int))
                    except:
                        continue

    df = pd.DataFrame(spectrum_list)

    return df


def calc_spec_quality(mz_ls, int_ls, mz_range=(50, 300)):
    """
    Calculate the spectral quality of the library
    specifically, the sum intensity % of the peaks in the mz_range
    """
    total_int = sum(int_ls)
    mz_int = 0
    for mz, intensity in zip(mz_ls, int_ls):
        if mz_range[0] <= mz <= mz_range[1]:
            mz_int += intensity
    return mz_int / total_int * 100


def analyze_spec_quality(library_mgf):

    df = generate_library_df(library_mgf)

    # remove column: SCANS
    df = df.drop(columns=['SCANS'])

    # add col: Scans, from 1 to len(df)
    df['Scans'] = range(1, len(df) + 1)

    # assess the quality of the spectra
    df['50_300_int_pct'] = df.apply(lambda x: calc_spec_quality(x['mz_ls'], x['intensity_ls']), axis=1)

    # write to pkl and tsv
    df.to_pickle('data/all_msms.pkl')
    df.to_csv('data/all_msms.tsv', sep='\t', index=False)

    # plot the distribution of the quality
    df['50_300_int_pct'].plot.hist(bins=20)
    plt.xlabel('50-300 m/z intensity %')
    plt.ylabel('Count')
    plt.title('Distribution of 50-300 m/z total intensity %')

    # save the plot
    plt.savefig('data/msms_quality.svg', format='svg', bbox_inches='tight', transparent=True)

    return df


if __name__ == '__main__':

    analyze_spec_quality('data/BA_Spectra_for_FDR_2.mgf')

