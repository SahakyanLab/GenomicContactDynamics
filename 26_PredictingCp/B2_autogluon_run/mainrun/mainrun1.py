#!/usr/bin/env python
# coding: utf-8

# In[73]:


import pandas as pd
import re
import numpy as np
import math
from autogluon.tabular import TabularDataset, TabularPredictor


# In[87]:


whorunsit = 'LiezelCluster'  # 'LiezelMac' | 'LiezelCluster'

if whorunsit == 'LiezelMac':
    home_dir = '/Users/ltamon'
elif whorunsit == 'LiezelCluster':
    home_dir = '/project/sahakyanlab/ltamon'
else:
    print("The supplied <whorunsit> option is not created in the script.")

wk_dir = home_dir + '/SahakyanLab/GenomicContactDynamics/26_PredictingCp'
#csv_dir = wk_dir + '/z_ignore_git/out_makeTable'
csv_dir = wk_dir + '/out_makeTable'
out_dir = wk_dir + '/out_autogluon_run'


# In[88]:


generate_train_data_subsample = True
do_subsample = True # Also set to true to load subsampled data for fit_model = True
fit_model = False
num_gpus = 1


# In[89]:


gcb = 'min2Mb'
chrs = ['chr1']
feat_regex = 'grp.compl.|grp.kmer3.|grp.kmer1.|grp.GC.'


# In[90]:


sampling_id = '5percExceptgrthn18' #'0point1percExceptgrthn18' #'5percExceptgrthn18'
percs = np.repeat([5,100], [18,3]) #np.repeat([0.001,0.001], [18,3]) #np.repeat([5,100], [18,3])
groups = list(range(1,22))
seed_val = 234
#[len(percs), len(groups)]
#[percs, groups]


# In[79]:


label_col = 'LABEL'
metric = 'root_mean_squared_error'
problem_type = 'regression' # (options: ‘binary’, ‘multiclass’, ‘regression’, ‘quantile’)

presets = 'best_quality'
time_limit = 5*60*60 # seconds
label_count_threshold = 2 # For multi-class classification problems, this is the minimum number of times a label must appear in dataset in order to be considered an output class.

#holdout_frac
# Baggin/stack ensembling - if enabled, dont provide tuning_data
#num_bag_folds=5 # how many times the k-fold bagging process is repeated to further reduce variance (increasing this may further boost accuracy but will substantially increase training times, inference latency, and memory/disk usage)
#num_bag_sets=1, 
#num_stack_levels=1 
#auto_stack=True # Autogluon will select bagging/stacking numbers
# specifying presets='best_quality' in fit() simply sets auto_stack=True

#num_trials = 5  # try at most 5 different hyperparameter configurations for each type of model
#search_strategy = 'auto'  # to tune hyperparameters using random search routine with a local scheduler
#hyperparameter_tune_kwargs = {  # HPO is not performed unless hyperparameter_tune_kwargs is specified
#    'num_trials': num_trials,
#    'scheduler' : 'local',
#    'searcher': search_strategy,
#}


# In[81]:


chrs_id = ''.join(chrs)
try:
    chrs[chrs.index('chr23')] = 'chrX'
except ValueError:
        print("No chr23 -> chrX.")
        
feat_regex_id = re.compile('[^a-zA-Z0-9]').sub('', feat_regex)
#[chrs_id, feat_regex_id]


# In[82]:


def concat_chrcsv(csv_dir, gcb, chrs):
    
    df = []
    [df.append(pd.read_csv( csv_dir + '/' + gcb + '_' + chr + '_MLtbl.csv' )) for chr in chrs]
    df = pd.concat(df)
    
    return(df)

def choose_features(df, feat_regex, label_col):
    
    final_regex = feat_regex + '|' + label_col
    is_chosen = [ bool(re.search(pattern=final_regex, string=feat)) for feat in df.columns ] 
    is_todrop = [not bool for bool in is_chosen]
    df.drop(columns=list(df.columns[is_todrop]), inplace=True)

    return(df)
    
def switch_contact(df):
    
    features = np.array(df.columns)
    is_icol = ['.i_' in feat for feat in features]
    is_jcol = ['.j_' in feat for feat in features]
    
    features_switch = np.array(features)
    features_switch[is_icol] = features[is_jcol]
    features_switch[is_jcol] = features[is_icol]
    
    df_switch = df.set_axis(features_switch, axis=1, inplace=False)
    df_switch = df_switch[features]
    df = pd.concat([df, df_switch])
    
    return(df)

def sampleGroupInSeries(dfSeries, perc, group, seed_val):
    
    row_ind_group = dfSeries.index
    row_ind_group = row_ind_group[ dfSeries == group ]
    
    np.random.seed(seed_val)
    samp_size = math.ceil( sum(dfSeries == group) * perc / 100 )
    row_ind_group_sampled = list( np.random.choice(a = row_ind_group, size=samp_size, replace=False) )
    
    return(row_ind_group_sampled)

def sampleManyGroupsInSeries(dfSeries, percs, groups, seed_val):
    
    num_groups = len(groups)
    seed_generated = np.random.randint(low=0, high=1000, size=num_groups, dtype=int)
    
    row_ind_manygroups_sampled = []
    
    for i in range(num_groups):
        
        perc = percs[i]
        group = groups[i]
        seed_i = seed_generated[i]
        row_ind_manygroups_sampled = row_ind_manygroups_sampled + sampleGroupInSeries(dfSeries, perc, group, seed_val=seed_i)
        
    return(row_ind_manygroups_sampled)

def processData(df, feat_regex, label_col):
    
    df = choose_features(df, feat_regex, label_col)
    print("Feature columns chosen.")
    df = switch_contact(df)
    print("Contacts switched.")
    
    return(df)


# In[62]:


csv_id = gcb + '_' + chrs_id + '_MLtbl'
if do_subsample:
    csv_id = csv_id + '_' + feat_regex_id + '_seed' + str(seed_val) + '_' + sampling_id 
    
csv_path = csv_dir + '/' + csv_id + '.csv'
[csv_id, csv_path]


# In[83]:


if generate_train_data_subsample:
    
    print("Generating training data subsample...")
    
    train_data = concat_chrcsv(csv_dir, gcb, chrs)
    
    if do_subsample:
        ind = sampleManyGroupsInSeries(train_data['LABEL'], percs, groups, seed_val=seed_val)
        train_data = train_data.iloc[ind,:]
        print('Subsampling percentages: ')
        print(percs)
        print('for each group: ')
        print(groups)
        print("Subsampled " + str(len(train_data)) +  " training datapoints.")
        
        train_data = processData(train_data, feat_regex, label_col)
        
        train_data.to_csv(csv_path)
        print(csv_id + ": Training data saved as CSV.")


# In[86]:


if fit_model:
    
    model_id = gcb + '_' + chrs_id + '_' + feat_regex_id + '_' + presets + '_' + str(time_limit) + 'sec_' + problem_type
    train_data = TabularDataset(csv_path)
    
    if do_subsample:
        train_data.drop(columns=train_data.columns[0], axis=1, inplace=True)
        print(csv_id + ": Loaded subsampled training data from CSV.")
        model_id = model_id + '_subsample_seed' + str(seed_val) + '_' + sampling_id 
    else: # Make function to avoid repetition of code below
        train_data = concat_chrcsv(csv_dir, gcb, chrs)
        print(csv_id + ": Loaded full training data from CSV.")
        train_data = processData(train_data, feat_regex, label_col)
        
    model_id = model_id + '_' + str(len(train_data)) + 'datapoints'
    model_path = out_dir + '/' + 'agModel_' + model_id


# In[ ]:


if fit_model:
    predictor = TabularPredictor(label=label_col, eval_metric=metric, path=model_path, problem_type=problem_type, learner_kwargs={'label_count_threshold': label_count_threshold}).fit(train_data, ag_args_fit={'num_gpus': num_gpus}, presets=presets, time_limit=time_limit)

