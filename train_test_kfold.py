#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def parseArgs():
    from ml_prediction.model_variables import featuresets
    """Parse command line arguments"""
    import argparse
    import traceback
    try:
        parser = argparse.ArgumentParser()
        # dataset_fbase
        parser.add_argument('-f', default='GOh1052mut', help='Base dataset name found in csv filename, also to be used as the folder name for the trained model, e.g. "GOh1052"')
        # data_folder
        parser.add_argument('--data_folder', default='./examples/', help='Relative path of directory in which dataset Input subdirectory is found, e.g. "./examples/"')
        # n_splits
        parser.add_argument('--n_splits', default=4, type=int, help='Number of splits for k-fold testing')
        # num_bag_folds
        parser.add_argument('--num_bag_folds', default=8, type=int, help='Number of folds used for model bagging')
        # num_stack_levels
        parser.add_argument('--num_stack_levels', default=1, type=int, help='Number of levels of model stacking')
        # save_model
        parser.add_argument('-save_model', default='GOh1052mut', help='Folder name for saving the trained model')
    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()
    return parser.parse_args()

def main():

    import pandas as pd
    import matplotlib.pyplot as plt
    from ml_prediction.load_dataset import load_autogluon_tabular_datasets, split_data_to_trainval_test
    from ml_prediction.model_variables import featuresets
    from ml_prediction.model_utils import get_split_metrics_summary, plot_model_metrics_splits, plot_model_metrics_summary
    from ml_prediction.autogluon_models import autogluon_classifier

    # get inputs from commandline
    args = parseArgs()
    data_folder = args.data_folder
    dataset_fbase = args.f
    n_splits = args.n_splits
    num_bag_folds = args.num_bag_folds
    num_stack_levels = args.num_stack_levels
    save_model = f'{data_folder}trained_models/{args.save_model}/'

    # other settings
    load_model = None
    y_feature = 'CategoryV3'
    x_features = featuresets['mut']
    eval_metric_list = ['mcc', 'accuracy']
    save_res = f'{data_folder}Output/{dataset_fbase}'
    split_type = 'RandomSplit'

    # get splits
    dataset_fpath_all = f'{data_folder}Input/{dataset_fbase}.csv'
    split_dict = split_data_to_trainval_test(dataset_fpath_all, n_splits=n_splits)

    # initialize dataframe to store results
    res_cols = ['model', 'train_or_test', 'split_label'] + eval_metric_list
    ensemble_metrics = []
    split_metrics = []

    # iterate through splits and load each dataset
    for i, split_label in enumerate(split_dict):
        print('**************************')
        print(f'Split {i}: {split_label}')
        print('**************************')
        # load data for split
        save_data_split = None # f'{data_folder}Input/{dataset_fbase}_{split_type}{i}' #
        trainval_data, test_data = load_autogluon_tabular_datasets(split_dict, split_label, y_feature, x_features, dataset_fpath_all, save_data_split=save_data_split)
        # train and validate model #
        predictor, y_pred, res = autogluon_classifier(
            trainval_data, test_data, label=y_feature, metrics=eval_metric_list,
            save_model=save_model, load_model=load_model,
            num_bag_folds=num_bag_folds, num_stack_levels=num_stack_levels
        )
        # update overall results
        for train_or_test in ['trainval', 'test']:
            if train_or_test in res:
                leaderboard = res[train_or_test]
                leaderboard['split_label'] = split_label
                leaderboard['train_or_test'] = train_or_test
                split_metrics_selected = leaderboard.iloc[0].to_dict()
                split_metrics.append(split_metrics_selected)
                ensemble_metrics += leaderboard.to_dict('records')

    # save data
    ensemble_metrics = pd.DataFrame(ensemble_metrics)[res_cols]
    split_metrics = pd.DataFrame(split_metrics)[res_cols]
    split_metrics_summary = get_split_metrics_summary(split_metrics, eval_metric_list=eval_metric_list, save_results=None)
    ensemble_metrics.to_csv(f'{save_res}_ensemble_metrics.csv')
    split_metrics.to_csv(f'{save_res}_split_metrics.csv')
    split_metrics_summary.to_csv(f'{save_res}_split_metrics_summary.csv')

    # plot metrics for best selected model
    for eval_metric in eval_metric_list:
        # all splits
        figtitle = f'{eval_metric} for best model from cross-val'
        savefig = f'{save_res}_{eval_metric}_selected_{n_splits}foldCV.png'
        plot_model_metrics_splits(split_metrics, eval_metric=eval_metric, barwidth=0.4, figtitle=figtitle, savefig=savefig, annotate=True)
        plt.show(block=True)
        plt.pause(0.001)
        # average of splits
        figtitle = f'Average {eval_metric} for best model from cross-val'
        savefig = f'{save_res}_{eval_metric}_selected_{n_splits}foldCV_AVG.png'
        plot_model_metrics_summary(split_metrics_summary, eval_metric=eval_metric, figsize=(8,6), barwidth=0.8, figtitle=figtitle, savefig=savefig, plot_errors=True, annotate=True)
        plt.show(block=True)
        plt.pause(0.001)
        plt.show()

    return ensemble_metrics, split_metrics, split_metrics_summary

if __name__ == "__main__":
    ensemble_metrics, split_metrics, split_metrics_summary = main()