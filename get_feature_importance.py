import autogluon.tabular as ag
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
from ml_prediction.model_variables import featuresets

def main():
    # Load saved AutoGluon model
    model_name = 'WeightedEnsemble_L3'
    model_dir = './data/ml_prediction/trained_models/GOh1052/'
    y_feature = 'CategoryV3'
    x_features = featuresets['mut']

    # Load new data for prediction
    data = pd.read_csv('./data/ml_prediction/Input/GOh1052.csv')[x_features+[y_feature]]

    # get predictor
    predictor = ag.TabularPredictor.load(model_dir)

    # get feature importance of internally transformed features for the best model
    feature_importance = predictor.feature_importance(data, model=model_name, feature_stage='transformed')
    feature_importance.to_csv(f'./data/ml_prediction/Output/GOh1052_feature_importance_{model}.csv')
    print(feature_importance)

if __name__ == "__main__":
    main()