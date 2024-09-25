#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024
"""
import os, glob, shutil
import numpy as np
import pandas as pd
def get_idxs_with_without_val(data, colname, colval_to_split_on):
    train_index = np.array(list(data[data[colname]!=colval_to_split_on].index))
    test_index = np.array(list(data[data[colname]==colval_to_split_on].index))
    return train_index, test_index

def get_nonrandom_splits(data, colname_to_split_on, split_dict={}):
    col_values = list(set(data[colname_to_split_on].tolist()))
    for val in col_values:
        train_index, test_index = get_idxs_with_without_val(data, colname_to_split_on, val)
        split_dict[val] = {'trainval': train_index, 'test': test_index}
    return split_dict

def split_data_to_trainval_test(dataset_fpath_all, n_splits=5):
    split_dict = {}
    data = pd.read_csv(dataset_fpath_all)
    from sklearn.model_selection import KFold
    kf = KFold(n_splits=n_splits, random_state=0, shuffle=True)
    for i, (train_index, test_index) in enumerate(kf.split(data.to_numpy())):
        split_dict.update({i: {'trainval': train_index, 'test': test_index}})
    return split_dict

def get_fpaths_for_precompiled_splits(dataset_dir, dataset_fbase, split_type, get_files_with_substring=None, get_files_without_substring=None):

    # get list of all csv files in directory
    csv_list = [f for f in os.listdir(dataset_dir) if f.endswith('.csv')]
    # filter to get filenames with or without specific substrings
    if get_files_with_substring is not None:
        for substring in get_files_with_substring:
            csv_list = [f for f in csv_list if f.find(substring)>-1]
    if get_files_without_substring is not None:
        for substring in get_files_with_substring:
            csv_list = [f for f in csv_list if f.find(substring)==-1]

    # initialize split_dict
    split_dict = {}

    ## filter files based on split_type
    # RANDOM SPLIT #
    if split_type == 'RandomSplit':
        csv_list = [f for f in csv_list if ((f.find('MutSplit') == -1) & (f.find('EnzSplit') == -1) & (f.find('SubSplit') == -1))]
        # get train sets
        train_fname_list = [f for f in csv_list if f.find('train') > -1]
        train_fname_list.sort()
        if len(train_fname_list) > 0:
            for i, train_fname in enumerate(train_fname_list):
                split_dict[i] = {}
                split_dict[i].update({'trainval': dataset_dir+train_fname})
        # get test sets
        test_fname_list = [f for f in csv_list if f.find('test')>-1]
        test_fname_list.sort()
        if len(test_fname_list) > 0:
            for i, test_fname in enumerate(test_fname_list):
                split_dict[i].update({'test': dataset_dir+test_fname})

    # MUTATION SPLIT #
    elif split_type == 'MutSplit':
        csv_list = [f for f in csv_list if (f.find('MutSplit') > -1)]
        # get unique mutation train/test sets
        mut_list = [s.strip('without_') for s in [f[f.find('MutSplit_')+9:-4] for f in csv_list]]
        mut_list = list(set(mut_list))
        mut_list.sort()
        for mut in mut_list:
            split_dict[mut] = {
                'trainval': f'{dataset_dir}{dataset_fbase}_MutSplit_without_{mut}.csv',
                'test': f'{dataset_dir}{dataset_fbase}_MutSplit_{mut}.csv',
            }
    # Enzyme SPLIT #
    elif split_type == 'EnzSplit':
        csv_list = [f for f in csv_list if (f.find('EnzSplit') > -1)]
        # get unique mutation train/test sets
        enz_list = [s.strip('without_') for s in [f[f.find('EnzSplit_')+9:-4] for f in csv_list]]
        enz_list = list(set(enz_list))
        enz_list.sort()
        for enz in enz_list:
            split_dict[enz] = {
                'trainval': f'{dataset_dir}{dataset_fbase}_EnzSplit_without_{enz}.csv',
                'test': f'{dataset_dir}{dataset_fbase}_EnzSplit_{enz}.csv',
            }
    # SUBSTRATE SPLIT #
    elif split_type == 'SubSplit':
        csv_list = [f for f in csv_list if (f.find('SubSplit') > -1)]
        # get unique mutation train/test sets
        sub_list = [s.strip('without_') for s in [f[f.find('SubSplit_')+9:-4] for f in csv_list]]
        sub_list = list(set(sub_list))
        sub_list.sort()
        for sub in sub_list:
            split_dict[sub] = {
                'trainval': f'{dataset_dir}{dataset_fbase}_SubSplit_without_{sub}.csv',
                'test': f'{dataset_dir}{dataset_fbase}_SubSplit_{sub}.csv',
            }
    # UNSEEN SPLIT #
    elif split_type == 'Unseen':
        csv_list = [f for f in csv_list if (f.find('UnseenSplit') > -1)]
    return split_dict

def load_autogluon_tabular_datasets(split_dict, split_label, y_feature, x_features, dataset_fpath_all=None, save_data_split=None):
    from autogluon.tabular import TabularDataset
    # get train data
    if isinstance(split_dict[split_label]['trainval'], str):
        trainval_fpath = split_dict[split_label]['trainval']
        data = pd.read_csv(trainval_fpath)[[y_feature] + x_features]
        trainval_data = TabularDataset(data)
    else:
        data = pd.read_csv(dataset_fpath_all)[[y_feature] + x_features]
        trainval_idxs = split_dict[split_label]['trainval']
        trainval_data = data.iloc[trainval_idxs, :]
        if save_data_split is not None:
            trainval_data.to_csv(f'{save_data_split}_train.csv')
        trainval_data = TabularDataset(trainval_data)
    trainval_data.head()
    trainval_data[y_feature].describe()
    # get test data, if available
    if 'test' in split_dict[split_label]:
        if isinstance(split_dict[split_label]['test'], str):
            test_fpath = split_dict[split_label]['test']
            data = pd.read_csv(test_fpath)[[y_feature] + x_features]
            test_data = TabularDataset(data)
        else:
            data = pd.read_csv(dataset_fpath_all)[[y_feature] + x_features]
            test_idxs = split_dict[split_label]['test']
            test_data = data.iloc[test_idxs, :]
            if save_data_split is not None:
                test_data.to_csv(f'{save_data_split}_test.csv')
            test_data = TabularDataset(test_data)
    else:
        test_data = None
    return trainval_data, test_data


def get_datasets_for_splittype(split_type, dataset_fbase, data_folder, get_files_with_substring='train', get_files_without_substring=None):

    csv_list = [f for f in os.listdir(f'{data_folder}Input/{dataset_fbase}/') if f.endswith('.csv')]
    # random split
    if split_type == 'RandomSplit':
        csv_list = [f for f in csv_list if ((f.find('MutSplit') == -1) & (f.find('SubSplit') == -1))]
    elif split_type == 'MutSplit':
        csv_list = [f for f in csv_list if (f.find('MutSplit') > -1)]
    elif split_type == 'SubSplit':
        csv_list = [f for f in csv_list if (f.find('SubSplit') > -1)]
    elif split_type == 'Unseen':
        csv_list = [f for f in csv_list if (f.find('UnseenSplit') > -1)]

    # filter to get filenames with or without specific substrings
    if get_files_with_substring is not None:
        for substring in get_files_with_substring:
            csv_list = [f for f in csv_list if f.find(substring)>-1]
    if get_files_without_substring is not None:
        for substring in get_files_with_substring:
            csv_list = [f for f in csv_list if f.find(substring)==-1]

    return csv_list

