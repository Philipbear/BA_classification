import graphviz
from sklearn.tree import export_graphviz

import msms_quality
import analyze_ms2_library
import calc_ms2_feature
import process_BA_labels
import filter_data
import train_tree
import joblib
import pandas as pd
import numpy as np


class MLDataset:
    def __init__(self, X, y, feature_names, label_dict):

        self.X = np.array(X)
        self.y = np.array(y)
        self.feature_names = feature_names
        self.label_dict = label_dict  # class label to numerical label

        self.fragment_feature = 0
        self.nl_feature = 0
        self.hnl_feature = 0
        self.frag_pair_feature = 0

        for feature_name in feature_names:
            if feature_name.startswith('frag_'):
                self.fragment_feature += 1
            elif feature_name.startswith('nl_'):
                self.nl_feature += 1
            elif feature_name.startswith('hnl_'):
                self.hnl_feature += 1
            else:
                self.frag_pair_feature += 1


def analyze_data(count_cutoff=10, top_n=50, rel_int_cutoff=0.01, frag_lm=150, frag_um=350):
    """
    analyze data
    """
    # generate msms df, analyze msms quality
    # calculate intensity % of the peaks in mz 50-300
    msms_df = msms_quality.analyze_spec_quality('data/BA_Spectra_for_FDR_2.mgf')

    # analyze the ms2 library, to find common frag, nl, hnl, frag pairs
    analyze_ms2_library.analyze_all_ms2(msms_df, top_n=top_n, rel_int_cutoff=rel_int_cutoff,
                                        frag_lm=frag_lm, frag_um=frag_um)
    analyze_ms2_library.print_stats()

    # design and calculate ms2 features
    calc_ms2_feature.design_ms2_feature(count_cutoff=count_cutoff)
    calc_ms2_feature.calc_all_ms2_feature(msms_df, frag_lm=frag_lm, frag_um=frag_um)  # all_msms_feature.pkl

    # process BA labels
    process_BA_labels.process_ba_labels('data/BA_Spectra_for_FDR_names_labelled.csv')  # label_df.pkl
    process_BA_labels.process_ba_labels_stereo('data/BA_Spectra_for_FDR_names_labelled.csv')  # label_df.pkl

    return


def prepare_data(stereochem=True, round_intensity=10, round_intensity_ratio=0.1,
                 total_int_pct_50_300=40., amide_only=False):
    """
    prepare data for classification
    """
    msms_df = pd.read_pickle('data/all_msms_feature.pkl')
    feature_names = np.load('data/feature_names.npy')

    if stereochem:
        label_df = pd.read_pickle('data/label_df_stereo.pkl')
    else:
        label_df = pd.read_pickle('data/label_df.pkl')

    # filter the data
    df = filter_data.filter_data(msms_df, label_df, feature_names,
                                 round_intensity=round_intensity, round_intensity_ratio=round_intensity_ratio,
                                 total_int_pct_50_300=total_int_pct_50_300,
                                 amide_only=amide_only)

    # save the entire dataset
    X, y, y_label_dict = filter_data.reshape_data(df)
    dataset = MLDataset(X, y, feature_names, y_label_dict)
    dataset_name = 'data/dataset.joblib' if not stereochem else 'data/dataset_stereo.joblib'
    joblib.dump(dataset, dataset_name)

    # split the data into OH specific datasets
    for i in range(1, 5):
        df_specific = df[df['OH_cnt'] == i]
        X, y, y_label_dict = filter_data.reshape_data(df_specific)
        dataset_specific = MLDataset(X, y, feature_names, y_label_dict)
        dataset_specific_name = f'data/dataset_{i}.joblib' if not stereochem else f'data/dataset_{i}_stereo.joblib'
        joblib.dump(dataset_specific, dataset_specific_name)

    return


def train_model(oh_specific=False, stereochem=True,
                use_all_data=True, test_size=0.25,
                use_fragment=True, use_nl=True, use_hnl=True, use_frag_pair_ratio=True,
                max_depth=10, min_samples_split=5, min_samples_leaf=3, random_state=24):
    """
    train the model
    """

    # train the model
    if not oh_specific:
        dataset_name = 'data/dataset.joblib' if not stereochem else 'data/dataset_stereo.joblib'
        dataset = joblib.load(dataset_name)
        model = train_tree.train_tree(dataset,
                                      oh_specific=None,
                                      use_all_data=use_all_data,
                                      test_size=test_size,
                                      use_fragment=use_fragment,
                                      use_nl=use_nl,
                                      use_hnl=use_hnl,
                                      use_frag_pair_ratio=use_frag_pair_ratio,
                                      max_depth=max_depth, min_samples_split=min_samples_split,
                                      min_samples_leaf=min_samples_leaf,
                                      random_state=random_state)

        # save the model
        model_name = 'dtree/stereo/model.joblib' if stereochem else 'dtree/non_stereo/model.joblib'
        joblib.dump(model, model_name)
    else:
        # train OH specific models
        for i in range(1, 5):
            dataset_name = f'data/dataset_{i}.joblib' if not stereochem else f'data/dataset_{i}_stereo.joblib'
            dataset = joblib.load(dataset_name)
            model = train_tree.train_tree(dataset,
                                          oh_specific=i,
                                          use_all_data=use_all_data,
                                          test_size=test_size,
                                          use_fragment=use_fragment,
                                          use_nl=use_nl,
                                          use_hnl=use_hnl,
                                          use_frag_pair_ratio=use_frag_pair_ratio,
                                          max_depth=max_depth, min_samples_split=min_samples_split,
                                          min_samples_leaf=min_samples_leaf,
                                          random_state=random_state)

            # save the model
            model_name = f'dtree/stereo/model_{i}.joblib' if stereochem else f'dtree/non_stereo/model_{i}.joblib'
            joblib.dump(model, model_name)


def plot_tree(oh_specific=False, stereochem=True, max_depth=None):
    if not oh_specific:
        model_name = 'dtree/stereo/model.joblib' if stereochem else 'dtree/non_stereo/model.joblib'
        dataset_name = 'data/dataset.joblib' if not stereochem else 'data/dataset_stereo.joblib'
        dtree = joblib.load(model_name)
        dataset = joblib.load(dataset_name)
        feature_names = np.load('dtree/feature_names.npy')

        # Alternatively, visualize the decision tree using graphviz
        dot_data = export_graphviz(dtree, out_file=None,
                                   max_depth=max_depth,
                                   feature_names=feature_names.tolist(),
                                   class_names=list(dataset.label_dict.keys()),
                                   filled=True, rounded=True,
                                   special_characters=True)

        # Create graph from dot data
        graph = graphviz.Source(dot_data, format="svg")
        out_graph = f"dtree/stereo/decision_tree_graph" if stereochem else f"dtree/non_stereo/decision_tree_graph"
        graph.render(out_graph)  # Saves the tree to a file
        graph.view()  # Displays the tree

    else:
        for i in range(1, 5):
            model_name = f'dtree/stereo/model_{i}.joblib' if stereochem else f'dtree/non_stereo/model_{i}.joblib'
            dataset_name = f'data/dataset_{i}.joblib' if not stereochem else f'data/dataset_{i}_stereo.joblib'
            dtree = joblib.load(model_name)
            dataset = joblib.load(dataset_name)
            feature_names = np.load(f'dtree/feature_names_{i}.npy')

            # Alternatively, visualize the decision tree using graphviz
            dot_data = export_graphviz(dtree, out_file=None,
                                       max_depth=max_depth,
                                       feature_names=feature_names.tolist(),
                                       class_names=list(dataset.label_dict.keys()),
                                       filled=True, rounded=True,
                                       special_characters=True)

            # Create graph from dot data
            graph = graphviz.Source(dot_data, format="svg")
            out_graph = f"dtree/stereo/decision_tree_graph_{i}" if stereochem else f"dtree/non_stereo/decision_tree_graph_{i}"
            graph.render(out_graph)


if __name__ == '__main__':

    OH_specific = True
    stereo_chem = True

    # analyze_data(count_cutoff=10, top_n=50, rel_int_cutoff=0.01, frag_lm=150, frag_um=350)
    #
    # prepare_data(stereochem=stereo_chem,
    #              round_intensity=1, round_intensity_ratio=0.01,
    #              total_int_pct_50_300=40., amide_only=False)

    train_model(oh_specific=OH_specific,  # whether to train OH specific model (mono OH, di OH, tri OH)
                stereochem=stereo_chem,  # whether to consider stereochemistry
                use_all_data=True,  # when False, split the data into training and testing sets
                test_size=0.2,  # test size when use_all_data=False
                use_fragment=True, use_nl=False, use_hnl=False, use_frag_pair_ratio=True,
                max_depth=10, min_samples_split=10, min_samples_leaf=5, random_state=24)

    plot_tree(oh_specific=OH_specific, stereochem=stereo_chem)
