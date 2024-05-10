"""
analyze the ms2 library, retrieve the most common frags, NLs, and HNLs
"""

import pandas as pd
import numpy as np
from numba import njit
from tqdm import tqdm


def analyze_all_ms2(df):
    """
    analyze all the ms2 spectra, retrieve most common frags, NLs, and HNLs
    """
    all_frag_ls = np.array([])
    all_nl_ls = np.array([])
    all_hnl_ls = np.array([])
    # n*2 array
    all_frag_pair_ls = np.empty((0, 2), dtype=np.float64)

    for i, row in tqdm(df.iterrows(), total=len(df)):
        mz_ls = np.array(row['mz_ls'])
        int_ls = np.array(row['intensity_ls'])
        prec_mz = float(row['PEPMASS'])
        frag_ls, nl_ls, hnl_ls, frag_pair_ls = calc_ms2(mz_ls, int_ls, prec_mz, top_n=50, rel_int_cutoff=0.05,
                                                        frag_lm=50, frag_um=350,
                                                        nl_lm=30, nl_um=350, hnl_lm=30, hnl_um=350)
        all_frag_ls = np.append(all_frag_ls, frag_ls)
        all_nl_ls = np.append(all_nl_ls, nl_ls)
        all_hnl_ls = np.append(all_hnl_ls, hnl_ls)
        all_frag_pair_ls = np.concatenate((all_frag_pair_ls, frag_pair_ls), axis=0)

    # get the unique frags, nl, hnl and their counts
    frag_ls, frag_counts = np.unique(all_frag_ls, return_counts=True)
    nl_ls, nl_counts = np.unique(all_nl_ls, return_counts=True)
    hnl_ls, hnl_counts = np.unique(all_hnl_ls, return_counts=True)
    # get the unique frag pairs and their counts
    frag_pair_ls, frag_pair_counts = np.unique(all_frag_pair_ls, axis=0, return_counts=True)

    # sort the frags, nl, hnl by counts
    frag_ls = frag_ls[np.argsort(frag_counts)[::-1]]
    nl_ls = nl_ls[np.argsort(nl_counts)[::-1]]
    hnl_ls = hnl_ls[np.argsort(hnl_counts)[::-1]]
    frag_counts = np.sort(frag_counts)[::-1]
    nl_counts = np.sort(nl_counts)[::-1]
    hnl_counts = np.sort(hnl_counts)[::-1]

    # sort the frag pairs by counts
    frag_pair_ls = frag_pair_ls[np.argsort(frag_pair_counts)[::-1]]
    frag_pair_counts = np.sort(frag_pair_counts)[::-1]

    # save
    np.save('data/frag_ls.npy', frag_ls)
    np.save('data/frag_counts.npy', frag_counts)
    np.save('data/nl_ls.npy', nl_ls)
    np.save('data/nl_counts.npy', nl_counts)
    np.save('data/hnl_ls.npy', hnl_ls)
    np.save('data/hnl_counts.npy', hnl_counts)
    np.save('data/frag_pair_ls.npy', frag_pair_ls)
    np.save('data/frag_pair_counts.npy', frag_pair_counts)


def print_stats():
    frag_ls = np.load('data/frag_ls.npy')
    frag_counts = np.load('data/frag_counts.npy')
    nl_ls = np.load('data/nl_ls.npy')
    nl_counts = np.load('data/nl_counts.npy')
    hnl_ls = np.load('data/hnl_ls.npy')
    hnl_counts = np.load('data/hnl_counts.npy')

    # print their lengths
    print(f'Number of unique frags: {len(frag_ls)}')
    print(f'Number of unique NLs: {len(nl_ls)}')
    print(f'Number of unique HNLs: {len(hnl_ls)}')

    print(f'Number of unique frags with counts > 10: {np.sum(frag_counts > 10)}')
    print(f'Number of unique NLs with counts > 10: {np.sum(nl_counts > 10)}')
    print(f'Number of unique HNLs with counts > 10: {np.sum(hnl_counts > 10)}')

    frag_pair_ls = np.load('data/frag_pair_ls.npy')
    frag_pair_counts = np.load('data/frag_pair_counts.npy')
    print(f'Number of unique frag pairs: {len(frag_pair_ls)}')
    print(f'Number of unique frag pairs with counts > 10: {np.sum(frag_pair_counts > 10)}')

    # print top 50 frag pairs
    for i in range(50):
        print(f'{frag_pair_ls[i][0]}-{frag_pair_ls[i][1]}: {frag_pair_counts[i]}')


@njit(cache=True)
def calc_ms2(mz_ls, int_ls, prec_mz, top_n=50, rel_int_cutoff=0.05,
             frag_lm=50, frag_um=350,
             nl_lm=30, nl_um=300, hnl_lm=30, hnl_um=300):
    """
    calculate the frag, nl, hnl array for a single ms2 spectrum
    :param mz_ls: np array, mz list
    :param int_ls: np array, intensity list
    :param prec_mz: float, precursor mz
    :param top_n: int, top n most intense frags
    :param rel_int_cutoff: float, relative intensity cutoff
    :param frag_lm: float, frag lower limit
    :param frag_um: float, frag upper limit
    :param nl_lm: float, nl lower limit
    :param nl_um: float, nl upper limit
    :param hnl_lm: float, hnl lower limit
    :param hnl_um: float, hnl upper limit
    """
    # clean ms2

    # deisotope
    bool_ls = np.zeros(len(mz_ls), dtype=np.bool_)
    bool_ls[0] = True
    for i in range(1, len(mz_ls)):
        if np.min(np.abs(mz_ls[i] - mz_ls[:i] - 1.00335)) < 0.01:
            bool_ls[i] = False
        else:
            bool_ls[i] = True
    mz_ls = mz_ls[bool_ls]
    int_ls = int_ls[bool_ls]

    # relative intensity cutoff
    int_ls = int_ls / np.max(int_ls)
    idx = np.where(int_ls > rel_int_cutoff)[0]
    mz_ls = mz_ls[idx]
    int_ls = int_ls[idx]

    # keep top top_n most intense frags
    if len(mz_ls) > top_n:
        idx = np.argsort(int_ls)[::-1]
        mz_ls = mz_ls[idx][:top_n]

    # de-precursor
    frag_um = min(prec_mz - 1.5, frag_um)

    # frag
    frag_ls = mz_ls[(mz_ls > frag_lm) & (mz_ls < frag_um)]

    # NL
    nl_ls = prec_mz - mz_ls
    nl_ls = nl_ls[(nl_ls > nl_lm) & (nl_ls < nl_um)]

    # hnl: iterate through all the frags, mass diff between each pair
    hnl_size = len(mz_ls) * (len(mz_ls) - 1) // 2
    hnl_ls = np.zeros(hnl_size, dtype=np.float64)
    idx = 0
    for i in range(len(mz_ls)):
        for j in range(i + 1, len(mz_ls)):
            hnl = mz_ls[j] - mz_ls[i]
            if hnl_lm < hnl < hnl_um:
                hnl_ls[idx] = hnl
                idx += 1
    hnl_ls = hnl_ls[:idx]

    frag_ls = np.round(frag_ls, 2)
    nl_ls = np.round(nl_ls, 2)
    hnl_ls = np.round(hnl_ls, 2)

    # get the unique frags, nl, hnl
    frag_ls = np.unique(frag_ls)
    nl_ls = np.unique(nl_ls)
    hnl_ls = np.unique(hnl_ls)

    # frag pairs
    frag_pair_ls = np.zeros((len(frag_ls) * (len(frag_ls) - 1) // 2, 2), dtype=np.float64)
    idx = 0
    for i in range(len(frag_ls) - 1):
        for j in range(i + 1, len(frag_ls)):
            frag_pair_ls[idx] = [frag_ls[i], frag_ls[j]]
            idx += 1

    frag_pair_ls = frag_pair_ls[:idx]

    return frag_ls, nl_ls, hnl_ls, frag_pair_ls


if __name__ == '__main__':
    df = pd.read_pickle('data/all_msms.pkl')
    analyze_all_ms2(df)

    print_stats()
