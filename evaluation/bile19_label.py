import pandas as pd


# [1*]O_beta [2*]O_alpha [2*]O_beta [3*]=O [3*]O_alpha [3*]O_beta  [4*]O_alpha [4*]O_beta
#  [6*]O_alpha	[6*]O_beta [6*]=O	[7*]O_alpha	[7*]O_beta [7*]=O
# [12*]O_alpha [12*]O_beta [12*]=O [14*]O_alpha
# [15*]O_alpha [15*]O_beta
# [16*]O_alpha [16*]=O


def process_unique_smiles():

    df = pd.read_csv('data/BILELIB19_Names_ok_labaled_IM.csv')

    # dereplicate
    df = df.drop_duplicates(subset=['SMILES']).reset_index(drop=True)

    # fill all nan with 0
    df = df.fillna(0)

    # for col [1*]O_beta, fill 1 if 1, else 0
    df['[1*]O_beta'] = df['[1*]O_beta'].apply(lambda x: 1 if x == 1 else 0)

    # create labels
    df['group'] = df.apply(lambda x: gen_label_stereo(x['[1*]O_beta'], x['[2*]O_alpha'], x['[2*]O_beta'], x['[3*]=O'],
                                                      x['[3*]O_alpha'], x['[3*]O_beta'], x['[4*]O_alpha'],
                                                      x['[4*]O_beta'],
                                                      x['[6*]O_alpha'], x['[6*]O_beta'], x['[6*]=O'], x['[7*]O_alpha'],
                                                      x['[7*]O_beta'], x['[7*]=O'], x['[12*]O_alpha'], x['[12*]O_beta'],
                                                      x['[12*]=O'], x['[14*]O_alpha'], x['[15*]O_alpha'],
                                                      x['[15*]O_beta'], x['[16*]O_alpha'], x['[16*]=O']), axis=1)

    df = df[['SMILES', 'group']]

    # save the df
    df.to_csv('data/label/BILELIB19_SMILES_group.tsv', sep='\t', index=False)


def get_label():

    df = pd.read_csv('data/label/BILELIB19_corrected_massql.tsv', sep='\t')

    smiles_df = pd.read_csv('data/label/BILELIB19_SMILES_group.tsv', sep='\t')
    # dictionary of mapping SMILES to group
    smiles_dict = dict(zip(smiles_df['SMILES'], smiles_df['group']))

    df['group'] = df['SMILES'].map(smiles_dict)

    # at least pass mono/di/tri MSQL
    df = df[(df['Monohydroxy'] == 1) | (df['Dihydroxy'] == 1) | (df['Trihydroxy'] == 1)].reset_index(drop=True)

    df.to_csv('data/label/bilelib19_df.tsv', sep='\t', index=False)

    print(df['group'][df['Monohydroxy'] == 1].value_counts())
    '''
    group
3a       8
3a,6a    3
12a      2
12b      2
3keto    1
3a,7a    1
7b       1
    '''
    print(df['group'][df['Dihydroxy'] == 1].value_counts())
    '''
    group
3a,7a        119
3a,12a       103
3a,7b         88
3a,6a         84
3b,12a         3
3b,6a          3
3a             2
3a,6b          1
3keto          1
3a,7a,12a      1
3a,7a,16a      1
    '''

    print(df['group'][df['Trihydroxy'] == 1].value_counts())
    '''
    group
3a,7a,12a     143
3a,6a,7a       97
3a,6b,7b       89
3a,6b,7a       85
3b,6b,7a        7
3a,7b,12a       5
3a,7a,14a       4
3a,7b,12b       4
3a,7a           3
3a,12keto       2
3a,7a,16a       2
3a,12a,16a      2
3keto,7b        1
3a,7keto        1
3keto           1
3a,6a           1
3a,4b,7a        1
3a,6a,7b        1
3a,7a,12b       1
3a,7b           1
3keto,12a       1
3keto,7a        1
    '''


def gen_label_stereo(o1b, o2a, o2b, k3, o3a, o3b, o4a, o4b, o6a, o6b, k6, o7a, o7b, k7, o12a, o12b, k12, o14a,
                     o15a, o15b, o16a, k16):
    """
    generate label for BA
    """

    group_name = ''

    if o1b > 0:
        group_name += '1b,'
    if o2a > 0:
        group_name += '2a,'
    if o2b > 0:
        group_name += '2b,'
    if k3 > 0:
        group_name += '3keto,'
    if o3a > 0:
        group_name += '3a,'
    if o3b > 0:
        group_name += '3b,'
    if o4a > 0:
        group_name += '4a,'
    if o4b > 0:
        group_name += '4b,'
    if o6a > 0:
        group_name += '6a,'
    if o6b > 0:
        group_name += '6b,'
    if k6 > 0:
        group_name += '6keto,'
    if o7a > 0:
        group_name += '7a,'
    if o7b > 0:
        group_name += '7b,'
    if k7 > 0:
        group_name += '7keto,'
    if o12a > 0:
        group_name += '12a,'
    if o12b > 0:
        group_name += '12b,'
    if k12 > 0:
        group_name += '12keto,'
    if o14a > 0:
        group_name += '14a,'
    if o15a > 0:
        group_name += '15a,'
    if o15b > 0:
        group_name += '15b,'
    if o16a > 0:
        group_name += '16a,'
    if k16 > 0:
        group_name += '16keto,'
    if group_name:
        group_name = group_name[:-1]

    return group_name


if __name__ == '__main__':
    # process_unique_smiles()

    get_label()
