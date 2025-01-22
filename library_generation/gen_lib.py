"""
Prepare library mgf and tsv files that are ready to be uploaded to GNPS.
"""
import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from requests import get
from json import loads
from matchms.importing import load_from_mgf


def load_csv(dir_path):
    """
    Libhit files
    """

    dfs = os.listdir(dir_path)
    dfs = [x for x in dfs if x.endswith('.csv') and not x.startswith('.')]
    dfs = [os.path.join(dir_path, x) for x in dfs]

    all_df = pd.DataFrame()
    for df in dfs:
        print(f'Loading {df}')
        this_df = pd.read_csv(df)
        print('rows:', this_df.shape[0])

        if '_mono_' in df:
            this_df['id'] = 'mono_' + this_df['original_Scan'].astype(str)
        elif '_di_' in df:
            this_df['id'] = 'di_' + this_df['original_Scan'].astype(str)
        else:
            this_df['id'] = 'tri_' + this_df['original_Scan'].astype(str)

        all_df = pd.concat([all_df, this_df])

    # remove rows whose id is duplicated
    all_df = all_df[~all_df.duplicated(subset='id')]

    return all_df


def load_db_df():

    mono = pd.read_csv('data/mono_nowaterloss_scan_index.tsv', sep='\t')
    mono['id'] = 'mono_' + mono['stage2_merged_scan'].astype(str)

    di = pd.read_csv('data/di_nowaterloss_scan_index.tsv', sep='\t')
    di['id'] = 'di_' + di['stage2_merged_scan'].astype(str)

    tri = pd.read_csv('data/tri_nowaterloss_scan_index.tsv', sep='\t')
    tri['id'] = 'tri_' + tri['stage2_merged_scan'].astype(str)

    db_df = pd.concat([mono, di, tri])

    return db_df


def merge_df():
    """
    Merge libhit and db_df
    """

    libhit_df = load_csv('data')
    db_df = load_db_df()

    merged_df = pd.merge(libhit_df, db_df, on='id', how='left')

    merged_df['isomer_label'] = merged_df['isomer_label'].apply(lambda x: x.replace('OH2', '(OH)2').replace('OH3', '(OH)3'))
    merged_df['Compound_Name'] = merged_df['Compound_Name'].apply(lambda x: x.replace('""', ''))
    merged_df['new_name'] = merged_df.apply(lambda x: f'[BA_core: {x["isomer_label"]}] {x["Compound_Name"]}', axis=1)

    # save
    merged_df.to_csv('out/merged_db_all_metadata.tsv', sep='\t', index=False)

    return


def add_ms2():
    df = pd.read_csv('out/merged_db_all_metadata.tsv', sep='\t')
    df['peaks'] = None

    # For mono spectra
    spectra_from_path = list(load_from_mgf('data/mono_nowaterloss_stage2_result_merged.mgf'))
    for spec in tqdm(spectra_from_path):
        mask = (df['id'].str.startswith('mono_')) & (df['stage2_merged_scan'] == int(spec.metadata['scans']))
        if mask.any():
            first_match_idx = df.index[mask][0]
            peaks = np.column_stack((spec.peaks.mz, spec.peaks.intensities))
            df.at[first_match_idx, 'peaks'] = peaks

    # For di spectra
    spectra_from_path = list(load_from_mgf('data/di_nowaterloss_stage2_result_merged.mgf'))
    for spec in tqdm(spectra_from_path):
        mask = (df['id'].str.startswith('di_')) & (df['stage2_merged_scan'] == int(spec.metadata['scans']))
        if mask.any():
            first_match_idx = df.index[mask][0]
            peaks = np.column_stack((spec.peaks.mz, spec.peaks.intensities))
            df.at[first_match_idx, 'peaks'] = peaks

    # For tri spectra
    spectra_from_path = list(load_from_mgf('data/tri_nowaterloss_stage2_result_merged.mgf'))
    for spec in tqdm(spectra_from_path):
        mask = (df['id'].str.startswith('tri_')) & (df['stage2_merged_scan'] == int(spec.metadata['scans']))
        if mask.any():
            first_match_idx = df.index[mask][0]
            peaks = np.column_stack((spec.peaks.mz, spec.peaks.intensities))
            df.at[first_match_idx, 'peaks'] = peaks

    df.to_csv('out/merged_db_all_metadata_with_ms2.tsv', sep='\t', index=False)
    df.to_pickle('out/merged_db_all_metadata_with_ms2.pkl')


def create_gnps_files():
    df = pd.read_pickle('out/merged_db_all_metadata_with_ms2.pkl')

    out_rows = []
    new_scan = 1
    for i, row in tqdm(df.iterrows(), total=df.shape[0]):

        out_row = {
            'FILENAME': 'ba_isomer.mgf',
            'SEQ': '*..*',
            'COMPOUND_NAME': row['new_name'],
            'MOLECULEMASS': row['precmz'],
            'INSTRUMENT': 'Orbitrap',
            'IONSOURCE': 'LC-ESI',
            'EXTRACTSCAN': new_scan,
            'SMILES': 'N/A',
            'INCHI': 'N/A',
            'INCHIAUX': 'N/A',
            'CHARGE': '1',
            'IONMODE': 'Positive',
            'PUBMED': 'N/A',
            'ACQUISITION': 'Crude',
            'EXACTMASS': row['precmz'] - 1.007276,
            'DATACOLLECTOR': 'Ipsita Mohanty',
            'ADDUCT': '[M+H]+',
            'CASNUMBER': 'N/A',
            'PI': 'Pieter Dorrestein',
            'LIBQUALITY': '4',
            'GENUS': 'N/A',
            'SPECIES': 'N/A',
            'INTEREST': 'N/A',
            'STRAIN': 'N/A'
        }
        out_rows.append(out_row)
        new_scan += 1

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv('out/ba_isomer.tsv', sep='\t', index=False)

    with open('out/ba_isomer.mgf', 'w') as f:
        new_scan = 1
        for i, row in df.iterrows():

            f.write('BEGIN IONS\n')
            f.write(f'TITLE={row["new_name"]}\n')
            f.write(f'PEPMASS={row["precmz"]}\n')
            f.write(f'SCANS={new_scan}\n')
            for mz, intensity in row['peaks']:
                f.write(f'{mz} {intensity}\n')
            f.write('END IONS\n\n')
            new_scan += 1


def load_from_usi(usi):
    """
    Load spectrum from USI
    :param usi: USI
    :return: spectrum
    """
    url = 'https://metabolomics-usi.gnps2.org/json/?usi1=' + usi
    response = get(url, timeout=10)
    json_data = loads(response.text)

    # check if the USI is valid
    if 'error' in json_data:
        raise ValueError

    prec_mz = json_data['precursor_mz']
    peaks = np.asarray(json_data['peaks'])

    return prec_mz, peaks



if __name__ == '__main__':

    # merge_df()

    # add_ms2()

    create_gnps_files()



