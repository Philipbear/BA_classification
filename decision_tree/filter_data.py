import numpy as np
import pandas as pd


def filter_data(msms_df, label_df, feature_names,
                round_intensity=10, round_intensity_ratio=0.1,
                total_int_pct_50_300=50., amide_only=False):
    """
    prefilter the data used for classification
    """
    # prefilter
    label_df = label_df[label_df['tail_has'] != 'Ignore'].reset_index(drop=True)
    if amide_only:
        label_df = label_df[label_df['tail_has'] == 'Amide'].reset_index(drop=True)

    msms_df = msms_df[msms_df['50_300_int_pct'] >= total_int_pct_50_300].reset_index(drop=True)

    # merge
    df = pd.merge(label_df, msms_df, on='Scans', how='inner')

    # round the intensity
    round_intensity_idx = [i for i, feature_name in enumerate(feature_names) if feature_name.startswith('frag_') or
                           feature_name.startswith('nl_')]
    round_intensity_idx = len(round_intensity_idx)
    round_int_ratio_idx = [i for i, feature_name in enumerate(feature_names) if feature_name.startswith('fragIntRatio_')]
    round_int_ratio_idx = len(round_int_ratio_idx)

    for i, row in df.iterrows():
        feature_arr = row['feature']
        feature_arr[:round_intensity_idx] = np.round(feature_arr[:round_intensity_idx] / round_intensity) * round_intensity
        feature_arr[round_int_ratio_idx:] = np.round(feature_arr[round_int_ratio_idx:] / round_intensity_ratio) * round_intensity_ratio
        df.at[i, 'feature'] = feature_arr

    print('Number of spectra:', len(df))
    print('Number of unique labels:', len(df['label'].unique()))
    print('Number of unique compounds:', len(df['INCHI'].unique()))

    # df OH cnt: how many ',' in the label
    df['OH_cnt'] = 0
    for i, row in df.iterrows():
        _label = row['label']
        if _label == 'No OH':
            _cnt = 0
        elif ',' in _label:
            _cnt = _label.count(',') + 1
        else:
            _cnt = 1
        df.at[i, 'OH_cnt'] = _cnt

    return df


def reshape_data(df):
    """
    reshape the data for classification
    """

    X, y = df['feature'].values.tolist(), df['label'].values.tolist()

    # convert y to numerical
    unique_labels = list(set(y))
    label_dict = {label: i for i, label in enumerate(unique_labels)}
    y = [label_dict[label] for label in y]

    return X, y, label_dict


if __name__ == '__main__':
    label_df = pd.read_pickle('data/label_df.pkl')
    msms_df = pd.read_pickle('data/all_msms_feature.pkl')
    filter_data(msms_df, label_df, np.load('data/feature_names.npy'))
