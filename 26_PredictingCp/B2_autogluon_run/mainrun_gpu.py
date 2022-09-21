#!/usr/bin/env python
# coding: utf-8

# In[26]:


import pandas as pd
import re
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor


# In[27]:


whorunsit = 'LiezelCluster' # 'LiezelMac' | 'LiezelCluster'

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


# In[28]:


gcb = 'min2Mb'
chrs = ['chr21']

label_col = 'LABEL'
subsample_size = 50  # subsample subset of data for quick run, input '' or any char to skip subsampling

feat_regex = 'grp.compl.|grp.kmer3.|grp.kmer1.|grp.GC.'
model_id = re.compile('[^a-zA-Z]').sub('', feat_regex)
model_id = model_id + '_' + str(subsample_size)

model_path = out_dir + '/' + 'agModel_'  + model_id 
metric = 'accuracy'
problem_type = 'multiclass' # (options: ‘binary’, ‘multiclass’, ‘regression’, ‘quantile’)

# fit()
presets = 'best_quality'
time_limit = 2*60 # seconds
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


# In[29]:


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


# In[31]:


train_data = concat_chrcsv(csv_dir, gcb, chrs)


# In[33]:


try:
    train_data = train_data.sample(n=subsample_size, random_state=0)
except TypeError:
        print("No subsampling of training data.")


# In[35]:


train_data = choose_features(train_data, feat_regex, label_col)


# In[18]:


train_data = switch_contact(train_data)


# In[20]:


train_data = TabularDataset(train_data)


# In[ ]:


predictor = TabularPredictor(label=label_col, eval_metric=metric, path=model_path, problem_type=problem_type).fit(train_data, presets=presets, time_limit=time_limit)

