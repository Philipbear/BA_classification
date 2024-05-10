import pandas as pd


def filter_data(msms_df, label_df, total_int_pct_50_300=50., amide_only=False):
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

    print('Number of spectra:', len(df))
    print('Number of unique labels:', len(df['label'].unique()))
    print('Number of unique compounds:', len(df['INCHI'].unique()))

    return df[['feature', 'label']]


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
    filter_data(msms_df, label_df, total_int_pct_50_300=50., amide_only=False)

