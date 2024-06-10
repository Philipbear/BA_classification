import numpy as np
from sklearn.tree import DecisionTreeClassifier, plot_tree, export_text
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.model_selection import cross_val_score, train_test_split
import matplotlib.pyplot as plt
import pandas as pd


def train_tree(dataset, oh_specific=None,
               use_all_data=True, test_size=0.25,
               use_fragment=True, use_nl=True, use_hnl=True, use_frag_pair_ratio=True,
               max_depth=10, min_samples_split=5, min_samples_leaf=3, random_state=24):

    print('OH Specific:', oh_specific)

    # Get the data
    X = calc_X(dataset, use_fragment, use_nl, use_hnl, use_frag_pair_ratio)
    y = dataset.y

    # save dataset.feature_names
    feature_name_file = 'dtree/feature_names.npy' if oh_specific is None else f'dtree/feature_names_{oh_specific}.npy'
    np.save(feature_name_file, dataset.feature_names)

    # Initialize and fit the decision tree classifier
    dtree = DecisionTreeClassifier(random_state=random_state,
                                   max_depth=max_depth,
                                   min_samples_split=min_samples_split,
                                   min_samples_leaf=min_samples_leaf)

    if use_all_data:
        # Use the entire dataset to train the model
        dtree.fit(X, y)
        _y = y  # y_true
        y_pred = dtree.predict(X)
    else:
        # Split the data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                            test_size=test_size, random_state=random_state)
        dtree.fit(X_train, y_train)
        _y = y_test  # y_true
        y_pred = dtree.predict(X_test)

    print('Accuracy:', accuracy_score(_y, y_pred))
    print('Classification Report:')
    print(classification_report(_y, y_pred, target_names=dataset.label_dict.keys()))

    # Confusion matrix
    cm = confusion_matrix(_y, y_pred)
    cm_df = pd.DataFrame(cm, index=dataset.label_dict.keys(), columns=dataset.label_dict.keys())
    print("Confusion Matrix:")
    print(cm_df)
    # save confusion matrix
    file_name = 'dtree/confusion_matrix.npy' if oh_specific is None else f'dtree/confusion_matrix_{oh_specific}.npy'
    cm_df.to_csv(file_name, index=False)

    # # Cross-validate the decision tree
    # cross_validation(dtree, X, y)

    tree_rules = export_text(dtree, feature_names=dataset.feature_names).split('\n')
    print("Decision Tree Rules:")
    for rule in tree_rules:
        print(rule)

    # # Visualize the decision tree using matplotlib
    # plt.figure(figsize=(20, 10))
    # plot_tree(dtree, feature_names=dataset.feature_names.tolist(), max_depth=5,
    #           class_names=list(dataset.label_dict.keys()), filled=True)
    # plt.title('Decision Tree')
    # plt.savefig('dtree/decision_tree_depth5.svg', format='svg', bbox_inches='tight', transparent=True)
    # plt.show()
    return dtree


def calc_X(dataset, use_fragment=True, use_nl=True, use_hnl=True, use_frag_pair_ratio=True):
    X_selected = np.empty((len(dataset.X), 0))
    feature_names = []
    cum_feature_idx = 0
    if use_fragment:
        X_selected = np.hstack((X_selected, dataset.X[:, cum_feature_idx:(cum_feature_idx + dataset.fragment_feature)]))
        feature_names += dataset.feature_names[cum_feature_idx:(cum_feature_idx + dataset.fragment_feature)].tolist()
    cum_feature_idx += dataset.fragment_feature
    if use_nl:
        X_selected = np.hstack((X_selected, dataset.X[:, cum_feature_idx:(cum_feature_idx + dataset.nl_feature)]))
        feature_names += dataset.feature_names[cum_feature_idx:(cum_feature_idx + dataset.nl_feature)].tolist()
    cum_feature_idx += dataset.nl_feature
    if use_hnl:
        X_selected = np.hstack((X_selected, dataset.X[:, cum_feature_idx:(cum_feature_idx + dataset.hnl_feature)]))
        feature_names += dataset.feature_names[cum_feature_idx:(cum_feature_idx + dataset.hnl_feature)].tolist()
    cum_feature_idx += dataset.hnl_feature
    if use_frag_pair_ratio:
        X_selected = np.hstack((X_selected, dataset.X[:, cum_feature_idx:(cum_feature_idx + dataset.frag_pair_feature)]))
        feature_names += dataset.feature_names[cum_feature_idx:(cum_feature_idx + dataset.frag_pair_feature)].tolist()

    dataset.feature_names = np.array(feature_names)
    return X_selected


def cross_validation(dtree, X, y, cv_fold=5):
    scores = cross_val_score(dtree, X, y, cv=cv_fold)
    print("Cross-validated scores:", scores)
    print("Average score:", scores.mean())

