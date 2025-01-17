import pandas as pd


def get_label():
    """
    Get the label of the library
    """
    # read the library
    df = pd.read_csv('data/new_core_corrected_massql.tsv', sep='\t')

    # for the NAME column, romove '_NCE45' if it exists
    df['_NAME'] = df['NAME'].apply(lambda x: x.replace('_NCE45', ''))

    df['group'] = df['_NAME'].apply(lambda x: x.split('_')[0] if '_' in x else x)

    # rename 3keta7a to 3keto7a
    df['group'] = df['group'].apply(lambda x: x.replace('keta', 'keto'))

    # label ground truths for mono, di, tri
    def get_oh_gt(name):
        total = 0
        # Convert to lowercase for case-insensitive matching
        text = name.lower()
        # Count occurrences using string methods
        total += text.count('alpha')
        total += text.count('beta')
        total += text.count('delta')
        total += 2 * text.count('keto')
        return total

    df['oh_gt'] = df['group'].apply(lambda x: get_oh_gt(x))

    # label ground truths for mono, di, tri
    df['mono_gt'] = df['oh_gt'].apply(lambda x: 1 if x == 1 else 0)
    df['di_gt'] = df['oh_gt'].apply(lambda x: 1 if x == 2 else 0)
    df['tri_gt'] = df['oh_gt'].apply(lambda x: 1 if x == 3 else 0)

    # should be either mono/di/tri BAs
    df = df[(df['mono_gt'] == 1) | (df['di_gt'] == 1) | (df['tri_gt'] == 1)].reset_index(drop=True)

    df['group'] = df['group'].apply(lambda x: x.replace('beta', 'b'))
    df['group'] = df['group'].apply(lambda x: x.replace('alpha', 'a'))

    # remove 3a7b, trihydroxy
    df = df[~((df['group'] == '3a7b') & (df['Trihydroxy'] == 1))]

    df['di_1_sc_oh'] = 0

    # save the result
    df.to_csv('data/label/new_core_df.tsv', sep='\t', index=False)


if __name__ == '__main__':
    ##########################################
    # new core library
    get_label()
