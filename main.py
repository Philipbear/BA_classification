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


def analyze_data(count_cutoff=10):
    """
    analyze data
    """
    # generate msms df, analyze msms quality
    msms_df = msms_quality.analyze_spec_quality('data/BA_Spectra_for_FDR_2.mgf')

    # analyze the ms2 library, to find common frag, nl, hnl, frag pairs
    analyze_ms2_library.analyze_all_ms2(msms_df)
    analyze_ms2_library.print_stats()

    # design and calculate ms2 features
    feature_names = calc_ms2_feature.design_ms2_feature(count_cutoff=count_cutoff)
    msms_df = calc_ms2_feature.calc_all_ms2_feature(msms_df)  # all_msms_feature.pkl

    # process BA labels
    label_df = process_BA_labels.process_ba_labels('data/BA_Spectra_for_FDR_names_labelled.csv')  # label_df.pkl

    return msms_df, label_df, feature_names


def prepare_data(total_int_pct_50_300=40., amide_only=False):
    """
    prepare data for classification
    """
    msms_df = pd.read_pickle('data/all_msms_feature.pkl')
    label_df = pd.read_pickle('data/label_df.pkl')
    feature_names = np.load('data/feature_names.npy')

    # filter the data
    df = filter_data.filter_data(msms_df, label_df, total_int_pct_50_300=total_int_pct_50_300,
                                 amide_only=amide_only)
    X, y, y_label_dict = filter_data.reshape_data(df)

    # save the dataset
    dataset = MLDataset(X, y, feature_names, y_label_dict)
    joblib.dump(dataset, 'data/dataset.joblib')

    return dataset


def train_model(use_fragment=True, use_nl=True, use_hnl=True, use_frag_pair_ratio=True,
                max_depth=10, min_samples_split=5, min_samples_leaf=3, random_state=24):
    """
    train the model
    """
    dataset = joblib.load('data/dataset.joblib')

    # train the model
    model = train_tree.train_tree(dataset,
                                  use_all_data=True,
                                  use_fragment=use_fragment,
                                  use_nl=use_nl,
                                  use_hnl=use_hnl,
                                  use_frag_pair_ratio=use_frag_pair_ratio,
                                  max_depth=max_depth, min_samples_split=min_samples_split,
                                  min_samples_leaf=min_samples_leaf,
                                  random_state=random_state)

    # save the model
    joblib.dump(model, 'dtree/model.joblib')


if __name__ == '__main__':

    analyze_data(count_cutoff=10)

    prepare_data(total_int_pct_50_300=40., amide_only=False)

    train_model(use_fragment=True, use_nl=True, use_hnl=False, use_frag_pair_ratio=True,
                max_depth=15, min_samples_split=6, min_samples_leaf=3, random_state=24)
