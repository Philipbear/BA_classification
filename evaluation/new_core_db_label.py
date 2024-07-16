import pandas as pd


def get_label():
    """
    Get the label of the library
    """
    # read the library
    df = pd.read_csv('data/label/new_core_corrected_massql.tsv', sep='\t')

    # at least pass mono/di/tri MSQL
    df = df[(df['Monohydroxy'] == 1) | (df['Dihydroxy'] == 1) | (df['Trihydroxy'] == 1)].reset_index(drop=True)

    # for the NAME column, romove '_NCE45' if it exists
    df['_NAME'] = df['NAME'].apply(lambda x: x.replace('_NCE45', ''))

    # # if '_' in the name, it is not a free BA
    # df['free_BA'] = df['_NAME'].apply(lambda x: 0 if '_' in x else 1)

    df['group'] = df['_NAME'].apply(lambda x: x.split('_')[0] if '_' in x else x)

    df['group'] = df['group'].apply(lambda x: x.replace('beta', 'b'))
    df['group'] = df['group'].apply(lambda x: x.replace('alpha', 'a'))

    # save the result
    df.to_csv('data/label/new_core_df.tsv', sep='\t', index=False)

    print(df['group'][df['Monohydroxy'] == 1].value_counts())
    '''
    group
3a    47
7b    33
7a    30
3b    24
    '''

    print(df['group'][df['Dihydroxy'] == 1].value_counts())
    '''
    group
3a6a     90
3a12b    80
3b12a    78
3b12b    75
3b7a     67
7a12a    55
3a7b      1
3a        1
    '''

    print(df['group'][df['Trihydroxy'] == 1].value_counts())
    '''
    group
3a6a7a         124
3b7a12a        124
3a6b7b         111
3a7b12a        108
3b7b12a        100
3a7b12b         92
3a7a12b         92
3b7b12b         92
3a7a16a         85
3keta7a         68
3a7keto         28
3a7bDelta22     24
3a7b            22
7a12a            1
    '''


if __name__ == '__main__':
    ##########################################
    # new core library
    get_label()
