#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
def parseArgs():
    """Parse command line arguments"""
    import argparse
    import traceback
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', default='GOh1052mut', help='Base dataset name found in csv filename, also to be used as the folder name for the trained model, e.g. "GOh1052"')
        parser.add_argument('-fsuffix', default='_train', help='Suffix appended to base dataset to get filename for train dataset to load, e.g. "_train"')
        parser.add_argument('--data_folder', default='./examples/', help='Relative path of directory in which dataset Input subdirectory is found, e.g. "./examples/"')
        parser.add_argument('--num_bag_folds', default=8, type=int, help='Number of folds used for model bagging')
        parser.add_argument('--num_stack_levels', default=1, type=int, help='Number of levels of model stacking')
        parser.add_argument('-save_model', default='GOh1052mut', help='Folder name for saving the trained model')
    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()
    return parser.parse_args()
def main():

    import pandas as pd
    from autogluon.tabular import TabularDataset
    from ml_prediction.model_variables import featuresets
    from ml_prediction.autogluon_models import autogluon_classifier

    # get inputs from commandline
    args = parseArgs()
    data_folder = args.data_folder
    dataset_fbase = args.f
    dataset_suffix = args.fsuffix
    num_bag_folds = args.num_bag_folds
    num_stack_levels = args.num_stack_levels
    save_model = f'{data_folder}trained_models/{args.save_model}/'

    # other settings
    load_model = None
    y_feature = 'CategoryV3'
    x_features = featuresets['mut']
    eval_metric_list = ['mcc', 'accuracy']
    save_res = f'{data_folder}Output/{dataset_fbase}'
    train_fpath = f'{data_folder}Input/{dataset_fbase}{dataset_suffix}.csv'

    # load train data
    train_data = TabularDataset(pd.read_csv(train_fpath)[[y_feature] + x_features])

    # perform inference on test data
    predictor, y_pred, res = autogluon_classifier(
        train_data=train_data, test_data=None, label=y_feature, metrics=eval_metric_list,
        save_model=save_model, load_model=load_model,
        num_bag_folds=num_bag_folds, num_stack_levels=num_stack_levels
    )
    # update overall results
    leaderboard = res['trainval']
    trainval_metrics_selected = leaderboard.iloc[0].to_dict()
    print(trainval_metrics_selected)
    leaderboard.to_csv(f'{save_res}_leaderboard_train.csv')

    return trainval_metrics_selected

if __name__ == "__main__":
    trainval_metrics_selected = main()