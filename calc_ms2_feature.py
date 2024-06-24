"""
Calculate the frag, nl, hnl array for all ms2 spectra
"""

import numpy as np
import pandas as pd
from numba import njit
from tqdm import tqdm


def design_ms2_feature(count_cutoff=10):
    """
    design the ms2 feature
    """
    frag_ls = np.load('data/frag_ls.npy')
    nl_ls = np.load('data/nl_ls.npy')
    hnl_ls = np.load('data/hnl_ls.npy')
    frag_pair_ls = np.load('data/frag_pair_ls.npy')

    frag_counts = np.load('data/frag_counts.npy')
    nl_counts = np.load('data/nl_counts.npy')
    hnl_counts = np.load('data/hnl_counts.npy')
    frag_pair_counts = np.load('data/frag_pair_counts.npy')

    # more than count_cutoff counts
    frag_ls = frag_ls[frag_counts > count_cutoff]
    nl_ls = nl_ls[nl_counts > count_cutoff]
    hnl_ls = hnl_ls[hnl_counts > count_cutoff]
    frag_pair_ls = frag_pair_ls[frag_pair_counts > count_cutoff]

    print(f'Count cutoff: {count_cutoff}')

    print(f'Number of unique frags: {len(frag_ls)}\n'
          f'Number of unique NLs: {len(nl_ls)}\n'
          f'Number of unique HNLs: {len(hnl_ls)}\n'
          f'Number of unique frag pairs: {len(frag_pair_ls)}')

    # save
    np.save('data/frag_feature_ls.npy', frag_ls)
    np.save('data/nl_feature_ls.npy', nl_ls)
    np.save('data/hnl_feature_ls.npy', hnl_ls)
    np.save('data/frag_pair_feature_ls.npy', frag_pair_ls)

    feature_no = len(frag_ls) + len(nl_ls) + len(hnl_ls) + len(frag_pair_ls)
    print(f'Number of features: {feature_no}')

    # generate feature names
    frag_feature_names = [f'frag_{frag:.1f}' for frag in frag_ls]
    nl_feature_names = [f'nl_{nl:.1f}' for nl in nl_ls]
    hnl_feature_names = [f'hnl_{hnl:.1f}' for hnl in hnl_ls]
    frag_pair_feature_names = [f'fragIntRatio_{frag1:.1f}_{frag2:.1f}' for frag1, frag2 in frag_pair_ls]

    feature_names = frag_feature_names + nl_feature_names + hnl_feature_names + frag_pair_feature_names
    np.save('data/feature_names.npy', feature_names)

    return


def calc_all_ms2_feature(df):
    """
    calculate the frag, nl, hnl array for all ms2 spectra
    """

    frag_ls = np.load('data/frag_feature_ls.npy')
    nl_ls = np.load('data/nl_feature_ls.npy')
    hnl_ls = np.load('data/hnl_feature_ls.npy')
    frag_pair_ls = np.load('data/frag_pair_feature_ls.npy')

    df['feature'] = [None] * len(df)

    for i, row in tqdm(df.iterrows(), total=len(df)):
        mz_ls = np.array(row['mz_ls'])
        int_ls = np.array(row['intensity_ls'])
        prec_mz = float(row['PEPMASS'])
        this_frag_mz, this_frag_int, this_nl_mz, this_nl_int, this_hnl_ls = calc_ms2(mz_ls, int_ls, prec_mz,
                                                                                     frag_lm=50, frag_um=350,
                                                                                     nl_lm=30, nl_um=350, hnl_lm=30,
                                                                                     hnl_um=350)
        this_feature = calc_ms2_feature(frag_ls, nl_ls, hnl_ls, frag_pair_ls,
                                        this_frag_mz, this_frag_int, this_nl_mz, this_nl_int, this_hnl_ls)

        # fill the feature
        df.at[i, 'feature'] = this_feature

    # save
    df.to_pickle('data/all_msms_feature.pkl')

    return


@njit(cache=True)
def calc_ms2(mz_ls, int_ls, prec_mz, frag_lm=50, frag_um=350,
             nl_lm=30, nl_um=300, hnl_lm=30, hnl_um=300):
    """
    calculate the frag, nl, hnl array for a single ms2 spectrum
    :param mz_ls: np array, mz list
    :param int_ls: np array, intensity list
    :param prec_mz: float, precursor mz
    :param frag_lm: float, frag lower limit
    :param frag_um: float, frag upper limit
    :param nl_lm: float, nl lower limit
    :param nl_um: float, nl upper limit
    :param hnl_lm: float, hnl lower limit
    :param hnl_um: float, hnl upper limit
    """
    int_ls = 100 * int_ls / int_ls.max()  # normalize

    # frag
    frag_idx = (mz_ls > frag_lm) & (mz_ls < frag_um)
    frag_mz = mz_ls[frag_idx]
    frag_int = int_ls[frag_idx]

    # NL
    nl_mz = prec_mz - mz_ls
    nl_idx = (nl_lm < nl_mz) & (nl_mz < nl_um)
    nl_mz = nl_mz[nl_idx]
    nl_int = int_ls[nl_idx]  # use the intensity of the original mz

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

    frag_mz = np.round(frag_mz, 1)
    nl_mz = np.round(nl_mz, 1)
    hnl_ls = np.round(hnl_ls, 1)

    # get the unique frags, nl, hnl
    out_frag_mz = np.unique(frag_mz)
    out_frag_int = np.zeros(len(out_frag_mz), dtype=np.float64)
    for i, frag in enumerate(out_frag_mz):
        out_frag_int[i] = frag_int[frag_mz == frag].max()

    out_nl_mz = np.unique(nl_mz)
    out_nl_int = np.zeros(len(out_nl_mz), dtype=np.float64)
    for i, nl in enumerate(out_nl_mz):
        out_nl_int[i] = nl_int[nl_mz == nl].max()

    out_hnl_ls = np.unique(hnl_ls)

    return out_frag_mz, out_frag_int, out_nl_mz, out_nl_int, out_hnl_ls


@njit(cache=True)
def calc_ms2_feature(frag_ls, nl_ls, hnl_ls, frag_pair_ls,
                     frag_mz, frag_int, nl_mz, nl_int, hnl_mz):
    """
    calculate the ms2 feature
    """
    frag_feature = np.zeros(len(frag_ls), dtype=np.float64)  # intensity
    for i, frag in enumerate(frag_ls):
        matches = frag_int[frag_mz == frag]
        if matches.size > 0:
            frag_feature[i] = matches[0]

    nl_feature = np.zeros(len(nl_ls), dtype=np.float64)  # intensity
    for i, nl in enumerate(nl_ls):
        matches = nl_int[nl_mz == nl]
        if matches.size > 0:
            nl_feature[i] = matches[0]

    hnl_feature = np.zeros(len(hnl_ls), dtype=np.float64)  # 0/1
    for i, hnl in enumerate(hnl_ls):
        hnl_feature[i] = 1.0 if hnl in hnl_mz else 0.0

    frag_pair_feature = np.zeros(len(frag_pair_ls), dtype=np.float64)
    for i, (frag1, frag2) in enumerate(frag_pair_ls):
        matches = frag_int[frag_mz == frag1]
        frag_int_1 = matches[0] if matches.size > 0 else 0.0

        matches = frag_int[frag_mz == frag2]
        frag_int_2 = matches[0] if matches.size > 0 else 0.0

        # ratio, avoid division by zero
        if frag_int_2 == 0:
            frag_pair_feature[i] = 20
        else:
            ratio = frag_int_1 / frag_int_2
            # clip
            ratio = min(ratio, 20)
            frag_pair_feature[i] = ratio

    return np.concatenate((frag_feature, nl_feature, hnl_feature, frag_pair_feature))


if __name__ == '__main__':
    design_ms2_feature(count_cutoff=10)

    df = pd.read_pickle('data/all_msms.pkl')
    calc_all_ms2_feature(df)
