"""
process BA labels, make into groups ready for classification
"""

import pandas as pd


def gen_label(o3, o6, o7, o12):
    """
    generate label for BA
    """
    group_no = 0
    group_no += 1 if o3 > 0 else 0
    group_no += 2 if o6 > 0 else 0
    group_no += 4 if o7 > 0 else 0
    group_no += 8 if o12 > 0 else 0

    if group_no == 0:
        group_name = 'No OH'
    else:
        group_name = ''
        if o3 > 0:
            group_name += '3,'
        if o6 > 0:
            group_name += '6,'
        if o7 > 0:
            group_name += '7,'
        if o12 > 0:
            group_name += '12,'
        group_name = group_name[:-1]
        group_name += '-OH'

    return group_name


def process_ba_labels(label_file):
    """
    process BA labels, make into groups ready for classification
    """
    df = pd.read_csv(label_file)

    # # print value counts
    # group = ['[3*]O', '[6*]O', '[7*]O', '[12*]O']
    # # group = ['[3*]O_alpha', '[3*]O_beta', '[6*]O_alpha', '[6*]O_beta',
    # # '[7*]O_alpha', '[7*]O_beta', '[12*]O_alpha', '[12*]O_beta']
    # for g in group:
    #     print(df[g].value_counts())

    # create labels
    df['label'] = df.apply(lambda x: gen_label(x['[3*]O'], x['[6*]O'], x['[7*]O'], x['[12*]O']), axis=1)
    print(df['label'].value_counts())

    # save as pkl
    df.to_pickle('data/label_df.pkl')

    return df


if __name__ == '__main__':
    process_ba_labels('data/BA_Spectra_for_FDR_names_labelled.csv')
