import pandas as pd

#######
# only use M+H and M+H-H2O spectra ##########

### for mono, di tri FDR, use [*]O (BILELIB19), new core (col name)
### '1-OH-Sidechain; 1-OH-core': Dihydroxy, 1-SC-OH (BILELIB19), No for all (new core)


mono_group_container = {
    # list of massql queries to apply, list of all possible groups
    '3a-OH': [['3-OH', '3a-OH'], ['3a']],
    '7b-OH': [['7b-OH'], ['7b']],
}


di_group_container = {
    # '1-OH-Sidechain; 1-OH-core': ['3a', '3b'],
    '1-ketone': [['1-ketone'], ['3keto', '6keto', '7keto', '12keto']],
    '3,12a-OH; 7,12a-OH': [['3,12a-OH; 7,12a-OH'], ['3a12a', '3b12a', '7a12a', '7b12a']],
    '3,12a-OH': [['3,12a-OH; 7,12a-OH', '3,12a-OH'], ['3a12a', '3b12a']],
    '7,12a-OH': [['3,12a-OH; 7,12a-OH', '7,12a-OH'], ['7a12a', '7b12a']],
    '3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH': [['3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH'],
                                           ['3a7a', '3a7b', '3b7a', '3b7b', '3a6a', '3a6b', '3b6a', '3b6b', '3a12b', '3b12b', '7a12b', '7b12b']],
    '3,7-OH; 3,12b-OH; 7,12b-OH': [['3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH', '3,7-OH; 3,12b-OH; 7,12b-OH'],
                                   ['3a7a', '3a7b', '3b7a', '3b7b', '3a12b', '3b12b', '7a12b', '7b12b']],
    '3,6-OH': [['3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH', '3,6-OH'], ['3a6a', '3a6b', '3b6a', '3b6b']],
    '3,7-OH': [['3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH', '3,7-OH; 3,12b-OH; 7,12b-OH', '3,7-OH'],
               ['3a7a', '3a7b', '3b7a', '3b7b']],
    '3,12b-OH; 7,12b-OH': [['3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH', '3,7-OH; 3,12b-OH; 7,12b-OH', '3,12b-OH; 7,12b-OH'],
                           ['3a12b', '3b12b', '7a12b', '7b12b']],
}


tri_group_container = {
    'All_ketone': [['All_ketone'], ['3keto7a', '3a7keto', '3a7bDelta22', '3a12keto', '3keto7b', '3a7keto']],
    '3k_7or12_OH': [['All_ketone', '3k_7or12_OH'], ['3keto7a', '3keto7b', '3keto12a', '3keto12b']],
    '3a_12k': [['All_ketone', '3a_12k'], ['3a12keto']],

    '3,7a,12a': [['3,7a,12a'], ['3a7a12a', '3b7a12a']],

    '3,6,7-OH': [['3,6,7-OH'], ['3a6a7a', '3a6a7b', '3a6b7a', '3a6b7b', '3b6a7a', '3b6a7b', '3b6b7a', '3b6b7b']],
    '3,6b,7b-OH': [['3,6,7-OH', '3,6b,7b-OH'], ['3a6b7b', '3b6b7b']],
    '3,6a,7a-OH': [['3,6,7-OH', '3,6a,7a-OH'], ['3a6a7a', '3b6a7a']],
    '3,6b,7a-OH; 3,6a,7b-OH': [['3,6,7-OH', '3,6b,7a-OH; 3,6a,7b-OH'], ['3a6b7a', '3b6b7a', '3a6a7b', '3b6a7b']],
    '3,6b,7a-OH': [['3,6,7-OH', '3,6b,7a-OH; 3,6a,7b-OH', '3,6b,7a-OH'], ['3a6b7a', '3b6b7a']],
    '3,6a,7b-OH': [['3,6,7-OH', '3,6b,7a-OH; 3,6a,7b-OH', '3,6a,7b-OH'], ['3a6a7b', '3b6a7b']],
}


def main_evaluation(group='mono', adduct_filter=False):

    if group == 'mono':
        group_container = mono_group_container
    elif group == 'di':
        group_container = di_group_container
    elif group == 'tri':
        group_container = tri_group_container

    bile19_df = pd.read_csv('data/label/bilelib19_df.tsv', sep='\t')
    new_core_df = pd.read_csv('data/label/new_core_df.tsv', sep='\t')

    # ion mode
    bile19_df = bile19_df[bile19_df['IONMODE'] == 'Positive'].reset_index(drop=True)

    # adduct
    if adduct_filter:
        bile19_df = bile19_df[bile19_df['ADDUCT'].isin(['M+H', 'M-H2O+H'])].reset_index(drop=True)
        new_core_df = new_core_df[new_core_df['ADDUCT'].isin(['[M+H]+', '[M-H2O+H]+'])].reset_index(drop=True)

    group_name = f'{group}hydroxy'
    group_name = group_name[0].upper() + group_name[1:]

    out_list = []
    # Mono, Di, Tri queries
    bile19_gt = bile19_df[f'{group}_gt'].values
    bile19_pred = bile19_df[group_name].values
    new_core_gt = new_core_df[f'{group}_gt'].values
    new_core_pred = new_core_df[group_name].values

    bile19_TP = sum((bile19_gt == 1) & (bile19_pred == 1))
    bile19_FP = sum((bile19_gt == 0) & (bile19_pred == 1))
    bile19_TN = sum((bile19_gt == 0) & (bile19_pred == 0))
    bile19_FN = sum((bile19_gt == 1) & (bile19_pred == 0))
    bile19_FN_adduct = bile19_df['ADDUCT'][(bile19_gt == 1) & (bile19_pred == 0)].value_counts()
    bile19_FN_adduct = bile19_FN_adduct.to_dict()

    new_core_TP = sum((new_core_gt == 1) & (new_core_pred == 1))
    new_core_FP = sum((new_core_gt == 0) & (new_core_pred == 1))
    new_core_TN = sum((new_core_gt == 0) & (new_core_pred == 0))
    new_core_FN = sum((new_core_gt == 1) & (new_core_pred == 0))
    new_core_FN_adduct = new_core_df['ADDUCT'][(new_core_gt == 1) & (new_core_pred == 0)].value_counts()
    new_core_FN_adduct = new_core_FN_adduct.to_dict()

    out_list.append([group_name, bile19_TP, bile19_FP, bile19_TN, bile19_FN, bile19_FN_adduct,
                     new_core_TP, new_core_FP, new_core_TN, new_core_FN, new_core_FN_adduct,
                     bile19_TP + new_core_TP, bile19_FP + new_core_FP, bile19_TN + new_core_TN,
                     bile19_FN + new_core_FN])

    # for other MassQL queries
    bile19_df = bile19_df[bile19_df[group_name] == 1].reset_index(drop=True)
    new_core_df = new_core_df[new_core_df[group_name] == 1].reset_index(drop=True)

    # 1-OH-Sidechain; 1-OH-core
    if group == 'di':
        bile19_gt = bile19_df['di_1_sc_oh'].values
        bile19_pred = bile19_df['1-OH-Sidechain; 1-OH-core'].values
        new_core_gt = new_core_df['di_1_sc_oh'].values
        new_core_pred = new_core_df['1-OH-Sidechain; 1-OH-core'].values

        bile19_TP = sum((bile19_gt == 1) & (bile19_pred == 1))
        bile19_FP = sum((bile19_gt == 0) & (bile19_pred == 1))
        bile19_TN = sum((bile19_gt == 0) & (bile19_pred == 0))
        bile19_FN = sum((bile19_gt == 1) & (bile19_pred == 0))
        bile19_FN_adduct = bile19_df['ADDUCT'][(bile19_gt == 1) & (bile19_pred == 0)].value_counts()
        bile19_FN_adduct = bile19_FN_adduct.to_dict()

        new_core_TP = sum((new_core_gt == 1) & (new_core_pred == 1))
        new_core_FP = sum((new_core_gt == 0) & (new_core_pred == 1))
        new_core_TN = sum((new_core_gt == 0) & (new_core_pred == 0))
        new_core_FN = sum((new_core_gt == 1) & (new_core_pred == 0))
        new_core_FN_adduct = new_core_df['ADDUCT'][(new_core_gt == 1) & (new_core_pred == 0)].value_counts()
        new_core_FN_adduct = new_core_FN_adduct.to_dict()

        out_list.append(['1-OH-Sidechain; 1-OH-core', bile19_TP, bile19_FP, bile19_TN, bile19_FN, bile19_FN_adduct,
                         new_core_TP, new_core_FP, new_core_TN, new_core_FN, new_core_FN_adduct,
                         bile19_TP + new_core_TP, bile19_FP + new_core_FP, bile19_TN + new_core_TN,
                         bile19_FN + new_core_FN])

    for _group, _group_items in group_container.items():

        massql_groups, group_ls = _group_items

        if not isinstance(group_ls, list):
            group_ls = [group_ls]

        bile19_ground_truth = bile19_df['group'].apply(
            lambda x: 1 if x and any(g == str(x) for g in group_ls) else 0).values

        bile19_prediction = bile19_df[_group].values
        for msql_group in massql_groups:
            bile19_prediction = bile19_df[msql_group].values & bile19_prediction

        bile19_TP = sum((bile19_ground_truth == 1) & (bile19_prediction == 1))
        bile19_FP = sum((bile19_ground_truth == 0) & (bile19_prediction == 1))
        bile19_TN = sum((bile19_ground_truth == 0) & (bile19_prediction == 0))
        bile19_FN = sum((bile19_ground_truth == 1) & (bile19_prediction == 0))

        # FN adduct forms value_counts
        bile19_FN_adduct = bile19_df['ADDUCT'][(bile19_ground_truth == 1) & (bile19_prediction == 0)].value_counts()
        bile19_FN_adduct = bile19_FN_adduct.to_dict()

        new_core_ground_truth = new_core_df['group'].apply(lambda x: 1 if any(g in str(x) for g in group_ls) else 0)

        new_core_prediction = new_core_df[_group].values
        for msql_group in massql_groups:
            new_core_prediction = new_core_df[msql_group].values & new_core_prediction

        new_core_TP = sum((new_core_ground_truth == 1) & (new_core_prediction == 1))
        new_core_FP = sum((new_core_ground_truth == 0) & (new_core_prediction == 1))
        new_core_TN = sum((new_core_ground_truth == 0) & (new_core_prediction == 0))
        new_core_FN = sum((new_core_ground_truth == 1) & (new_core_prediction == 0))

        # FN adduct forms value_counts
        new_core_FN_adduct = new_core_df['ADDUCT'][(new_core_ground_truth == 1) & (new_core_prediction == 0)].value_counts()
        new_core_FN_adduct = new_core_FN_adduct.to_dict()

        out_list.append([_group, bile19_TP, bile19_FP, bile19_TN, bile19_FN, bile19_FN_adduct,
                         new_core_TP, new_core_FP, new_core_TN, new_core_FN, new_core_FN_adduct,
                         bile19_TP + new_core_TP, bile19_FP + new_core_FP, bile19_TN + new_core_TN,
                         bile19_FN + new_core_FN])

    out_df = pd.DataFrame(out_list, columns=['group', 'bile19_TP', 'bile19_FP', 'bile19_TN', 'bile19_FN', 'bile19_FN_adduct',
                                             'new_core_TP', 'new_core_FP', 'new_core_TN', 'new_core_FN', 'new_core_FN_adduct',
                                             'total_TP', 'total_FP', 'total_TN', 'total_FN'])
    out_df['total_FDR'] = out_df['total_FP'] / (out_df['total_FP'] + out_df['total_TP'])
    out_df['total_FNR'] = out_df['total_FN'] / (out_df['total_FN'] + out_df['total_TP'])

    if adduct_filter:
        out_name = f'data/result/{group}_evaluation_H_adducts.tsv'
    else:
        out_name = f'data/result/{group}_evaluation_all_adducts.tsv'

    out_df.to_csv(out_name, sep='\t', index=False)


if __name__ == '__main__':
    main_evaluation('mono')

    main_evaluation('di')

    main_evaluation('tri')

    main_evaluation('mono', adduct_filter=True)

    main_evaluation('di', adduct_filter=True)

    main_evaluation('tri', adduct_filter=True)
