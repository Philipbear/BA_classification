import pandas as pd

from massql import msql_engine


NEW_QUERIES = {

    # monohydroxy
    'Monohydroxy': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=341.28:TOLERANCEMZ=0.01:INTENSITYPERCENT=2 AND MS2PROD=323.27:TOLERANCEMZ=0.01:INTENSITYPERCENT=2 AND MS2PREC=X AND MS2PROD=X-358.2871:TOLERANCEMZ=0.01:INTENSITYPERCENT=2",
    '3a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=107.086:TOLERANCEPPM=20:INTENSITYPERCENT=10:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=109.101:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.8:INTENSITYMATCHPERCENT=20",
    '7b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=177.127:TOLERANCEPPM=20:INTENSITYPERCENT=10:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=241.194:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.2:INTENSITYMATCHPERCENT=30",

    # dihydroxy
    'Dihydroxy': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=339.27:TOLERANCEMZ=0.01:INTENSITYPERCENT=5 AND MS2PROD=321.26:TOLERANCEMZ=0.01:INTENSITYPERCENT=5 AND MS2PREC=X AND MS2PROD=X-374.2894:TOLERANCEMZ=0.01:INTENSITYPERCENT=5",
    '1-OH-Sidechain; 1-OH-core_1': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=257.226:TOLERANCEPPM=10:INTENSITYPERCENT=40",
    '1-OH-Sidechain; 1-OH-core_2': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=257.226:TOLERANCEPPM=10:INTENSITYPERCENT=2 AND MS2MZ=147.117:TOLERANCEPPM=10:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=149.131:TOLERANCEPPM=10:INTENSITYMATCH=Y*1:INTENSITYMATCHPERCENT=30 AND MS2MZ=257.226:TOLERANCEPPM=10:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=161.132:TOLERANCEPPM=10:INTENSITYMATCH=Y*1:INTENSITYMATCHPERCENT=30",
    '1-OH-Sidechain; 1-OH-core_3': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=257.226:TOLERANCEPPM=10:INTENSITYPERCENT=2 AND MS2MZ=147.117:TOLERANCEPPM=10:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=149.131:TOLERANCEPPM=10:INTENSITYMATCH=Y*3:INTENSITYMATCHPERCENT=35",
    '1-ketone': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=161.132:TOLERANCEPPM=20:INTENSITYPERCENT=95",
    '3,12a-OH; 7,12a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=211.147:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=201.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.6:INTENSITYMATCHPERCENT=80",
    '3,12a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=243.174:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=239.179:TOLERANCEPPM=20:INTENSITYMATCH=Y*1:INTENSITYMATCHPERCENT=30",
    '7,12a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=243.174:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=239.179:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.25:INTENSITYMATCHPERCENT=30",
    '3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=211.147:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=201.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*4.2:INTENSITYMATCHPERCENT=40",
    '3,7-OH; 3,12b-OH; 7,12b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=161.132:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=201.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*2:INTENSITYMATCHPERCENT=20",
    '3,6-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=161.132:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=201.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.7:INTENSITYMATCHPERCENT=60",
    '3,7-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=175.148:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=173.132:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.5:INTENSITYMATCHPERCENT=30",
    '3,12b-OH; 7,12b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=175.148:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=173.132:TOLERANCEPPM=20:INTENSITYMATCH=Y*1:INTENSITYMATCHPERCENT=20",

    # trihydroxy
    'Trihydroxy': "QUERY scaninfo(MS2DATA) WHERE MS2PROD=337.25:TOLERANCEMZ=0.01:INTENSITYPERCENT=2 AND MS2PROD=319.24:TOLERANCEMZ=0.01:INTENSITYPERCENT=2 AND MS2PREC=X AND MS2PROD=X-390.277:TOLERANCEMZ=0.01:INTENSITYPERCENT=2",
    'All_ketone': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=107.086:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=105.07:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.55:INTENSITYMATCHPERCENT=40 AND MS2MZ=199.148:TOLERANCEPPM=20:INTENSITYMATCH=X:INTENSITYMATCHREFERENCE AND MS2MZ=209.132:TOLERANCEPPM=20:INTENSITYMATCH=X*0.45:INTENSITYMATCHPERCENT=40",
    '3k_7or12_OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=199.148:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=213.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*2:INTENSITYMATCHPERCENT=20",
    '3a_12k': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=319.242:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=309.257:TOLERANCEPPM=20:INTENSITYMATCH=Y*5.5:INTENSITYMATCHPERCENT=20",

    '3,7a,12a': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=209.132:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=199.148:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.7:INTENSITYMATCHPERCENT=30",
    '3,6,7-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=191.143:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=225.163:TOLERANCEPPM=20:INTENSITYMATCH=Y*1.2:INTENSITYMATCHPERCENT=60 AND MS2MZ=199.148:TOLERANCEPPM=20:INTENSITYMATCH=X:INTENSITYMATCHREFERENCE AND MS2MZ=213.163:TOLERANCEPPM=20:INTENSITYMATCH=X*0.7:INTENSITYMATCHPERCENT=10",
    '3,6b,7b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=173.132:TOLERANCEPPM=20:INTENSITYMATCH=X:INTENSITYMATCHREFERENCE AND MS2MZ=185.132:TOLERANCEPPM=20:INTENSITYMATCH=X*0.8:INTENSITYMATCHPERCENT=10",
    '3,6a,7a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=173.132:TOLERANCEPPM=20:INTENSITYMATCH=X:INTENSITYMATCHREFERENCE AND MS2MZ=185.132:TOLERANCEPPM=20:INTENSITYMATCH=X*1:INTENSITYMATCHPERCENT=5",
    '3,6b,7a-OH; 3,6a,7b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=159.117:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=161.132:TOLERANCEPPM=20:INTENSITYMATCH=Y*0.9:INTENSITYMATCHPERCENT=10",
    '3,6b,7a-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=125.096:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=309.257:TOLERANCEPPM=20:INTENSITYMATCH=Y*1:INTENSITYMATCHPERCENT=20",
    '3,6a,7b-OH': "QUERY scaninfo(MS2DATA) WHERE MS2MZ=125.096:TOLERANCEPPM=20:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2MZ=309.257:TOLERANCEPPM=20:INTENSITYMATCH=Y*2.5:INTENSITYMATCHPERCENT=20"
}


def correct_scans(input_mgf='data/BILELIB19.mgf'):
    """
    Correct the scans in the library
    """

    scan_cnt = 0
    new_lines = []
    # read mgf line by line
    with open(input_mgf, 'r') as file:
        for line in file:
            # empty line
            _line = line.strip()
            if not _line:
                new_lines.append(line)
                continue
            elif line.startswith('SCANS'):
                new_lines.append(f'SCANS={scan_cnt}\n')
                scan_cnt += 1
            else:
                new_lines.append(line)

    # write the new mgf
    out_name = input_mgf.replace('.mgf', '_corrected.mgf')
    with open(out_name, 'w') as file:
        for line in new_lines:
            file.write(line)


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

    # drop cols of mz_ls and intensity_ls
    df.drop(columns=['mz_ls', 'intensity_ls'], inplace=True)

    # save the result
    out_name = library_mgf.replace('.mgf', '.tsv')
    df.to_csv(out_name, sep='\t', index=False)


def massql_filter(input_mgf, massql_queries):
    """
    Filter the library for BA
    """
    # read the library
    df = pd.read_csv(input_mgf.replace('.mgf', '.tsv'), sep='\t')

    # process the queries
    for query_name, input_query in massql_queries.items():
        results_df = msql_engine.process_query(input_query, input_mgf)
        if len(results_df) == 0:
            df[query_name] = 0
            continue
        passed_scan_ls = results_df['scan'].values.tolist()
        passed_scan_ls = [str(x) for x in passed_scan_ls]

        df[query_name] = df['SCANS'].apply(lambda x: 1 if str(x) in passed_scan_ls else 0)

    # merge 1-OH-Sidechain; 1-OH-core_1 and 1-OH-Sidechain; 1-OH-core_2 and 1-OH-Sidechain; 1-OH-core_3
    df['1-OH-Sidechain; 1-OH-core'] = df['1-OH-Sidechain; 1-OH-core_1'] | df['1-OH-Sidechain; 1-OH-core_2'] | df[
        '1-OH-Sidechain; 1-OH-core_3']
    df.drop(['1-OH-Sidechain; 1-OH-core_1', '1-OH-Sidechain; 1-OH-core_2', '1-OH-Sidechain; 1-OH-core_3'], axis=1,
            inplace=True)

    # save the result
    out_name = input_mgf.replace('.mgf', '_massql.tsv')
    df.to_csv(out_name, sep='\t', index=False)


if __name__ == '__main__':
    # correct SCANS in the mgf file
    # correct_scans('data/BILELIB19.mgf')

    # generate_library_df('data/BILELIB19_corrected.mgf')

    massql_filter('data/BILELIB19_corrected.mgf', NEW_QUERIES)

