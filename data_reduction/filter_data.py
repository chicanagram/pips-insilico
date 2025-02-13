import pandas as pd

def filter_by_score(data, filt_by, filt_thres, dataset_fbase):
    n_pre = len(data)
    data_cols = data.columns.tolist()
    print('Dataset size (pre-reduction):', n_pre)
    print(f'Filtering by {filt_by}', filt_thres)
    # get scores to filter by
    scores = pd.read_csv('./data/data_reduction/'+dataset_fbase+'_'+filt_by+'.csv')
    scores = scores.rename(columns={'RealPos':'Position'})
    data = data.merge(scores[['Position',filt_by]], how='left', on='Position')
    # filter data
    data_filt = data[(data[filt_by]>filt_thres[0]) & (data[filt_by]<=filt_thres[1])]
    n_post = len(data_filt)
    print(f'Dataset size (post-reduction):', n_post, f'({round((n_pre-n_post)/n_pre*100,2)}% reduction)')
    # drop columns used for filtration
    data_filt = data_filt[[c for c in data_cols if c not in ['Position', filt_by]]]
    return data_filt
