#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import re
import numpy as np
import math
from autogluon.tabular import TabularDataset, TabularPredictor
from timeit import default_timer as timer


# In[4]:


start = timer()


# In[5]:


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
#out_dir = wk_dir + '/out_autogluon_run_cpu_medium'


# In[6]:


generate_train_data_subsample = False
do_subsample = True # Also set to true to load subsampled data for fit_model = True
fit_model = False
assess_model = True
pred_model = False
num_train_points = 1237618


# In[7]:


num_cpus = 3 # 10G per cpu/gpu
num_gpus = 0

presets = 'best_quality' #'medium_quality'
num_folds_parallel = 2 #3

time_limit = 5*60*60 # seconds

ag_args_fit = {'num_cpus': num_cpus, 'num_gpus': num_gpus}


# In[8]:


gcb = 'min2Mb'
chrs = ['chr1']
chrs_test = ['chr18']
feat_regex = 'grp.compl.|grp.kmer3.|grp.kmer1.|grp.GC.'


# In[9]:


sampling_id = '5percExceptgrthn18' #'0point1percExceptgrthn18' #'5percExceptgrthn18'
percs = np.repeat([5,100], [18,3]) #np.repeat([0.001,0.001], [18,3]) #np.repeat([5,100], [18,3])
groups = list(range(1,22))
seed_val = 234
#[len(percs), len(groups)]
#[percs, groups]


# In[10]:


label_col = 'LABEL'
metric = 'root_mean_squared_error'
problem_type = 'regression' # (options: ‘binary’, ‘multiclass’, ‘regression’, ‘quantile’)
label_count_threshold = 2 # For multi-class classification problems, this is the minimum number of times a label must appear in dataset in order to be considered an output class.


# In[11]:


try:
    chrs[chrs.index('chr23')] = 'chrX'
except ValueError:
        print("No chr23 -> chrX in chrs.")
chrs_id = ''.join(chrs)
        
feat_regex_id = re.compile('[^a-zA-Z0-9]').sub('', feat_regex)
#[chrs_id, feat_regex_id]


# In[12]:


try:
    chrs_test[chrs_test.index('chr23')] = 'chrX'
except ValueError:
        print("No chr23 -> chrX in chrs_test.")
chrs_test_id = ''.join(chrs_test)
#[chrs_test_id]


# In[65]:


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


# In[66]:


csv_id = gcb + '_' + chrs_id + '_MLtbl'
if do_subsample:
    csv_id = csv_id + '_' + feat_regex_id + '_seed' + str(seed_val) + '_' + sampling_id 
    
csv_path = csv_dir + '/' + csv_id + '.csv'
#[csv_id, csv_path]


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


# In[67]:


model_id = 'agModel_' + gcb + '_' + chrs_id + '_' + feat_regex_id + '_' + presets + '_' + str(time_limit) + 'sec_' + problem_type
if do_subsample:
    model_id = model_id + '_subsample_seed' + str(seed_val) + '_' + sampling_id 


# In[68]:


if fit_model:
    
    print("Preparing train data for model fitting...")
    
    train_data = TabularDataset(csv_path)
    
    if do_subsample:
        train_data.drop(columns=train_data.columns[0], axis=1, inplace=True)
        print(csv_id + ": Loaded subsampled training data from CSV.")
    else: # Make function to avoid repetition of code below
        train_data = concat_chrcsv(csv_dir, gcb, chrs)
        print(csv_id + ": Loaded full training data from CSV.")
        train_data = processData(train_data, feat_regex, label_col)
        
    model_id = model_id + '_' + str(len(train_data)) + 'datapoints'

model_path = out_dir + '/' + model_id


# In[ ]:


if fit_model:
    
    print("Model fitting...")
    
    predictor = TabularPredictor(label=label_col, eval_metric=metric, path=model_path, problem_type=problem_type, learner_kwargs={'label_count_threshold': label_count_threshold}).fit(train_data, ag_args_fit=ag_args_fit, ag_args_ensemble={'num_folds_parallel': num_folds_parallel}, presets=presets, time_limit=time_limit)


# In[ ]:


model_path = out_dir + '/' + model_id + '_' + str(num_train_points) + 'datapoints'


# In[ ]:


# Get test data from chr used for training
if assess_model & do_subsample: # by comparing validation and test score; likely overfitted if score_val > score_test
    
    print("Generating test data...")
    
    train_data = concat_chrcsv(csv_dir, gcb, chrs)
    train_ind = sampleManyGroupsInSeries(train_data['LABEL'], percs, groups, seed_val=seed_val)
    test_ind = sampleManyGroupsInSeries(train_data.drop(index=train_ind, inplace=False)['LABEL'], percs / 2, groups, seed_val=seed_val)
    
    if np.intersect1d(np.array(train_ind), np.array(test_ind)).any():
        raise Exception("Overlapping train and test sets.")
    else:
        test_data = train_data.iloc[test_ind,:]


# In[ ]:


if assess_model:
    
    print("Assessing model...")
    
    predictor = TabularPredictor.load(model_path)
    leaderboard = predictor.leaderboard(test_data, silent = True)
    leaderboard.to_csv(model_path + '_leaderboard.csv')


# MODEL ON TEST DATA

# In[74]:


if pred_model:
    test_data = concat_chrcsv(csv_dir, gcb, chrs_test)
    test_data = choose_features(test_data, feat_regex, label_col)
    test_data.drop(columns=label_col, inplace=True)
    test_data = TabularDataset(test_data)
    #test_data = test_data.sample(n=5, random_state=0)
    print('Test data prepared.')


# In[ ]:


if pred_model:
    predictor = TabularPredictor.load(model_path)
    y_pred = predictor.predict(test_data) # No label
    
    csv_pred_path = model_path + '_chrtest_' + chrs_test_id + '_MLtbl_pred.csv'
    print(csv_pred_path)
    y_pred.to_csv(csv_pred_path)


# In[ ]:


end = timer()
print(end - start)  # Time in seconds, e.g. 5.38091952400282


# In[ ]:


# NOTES

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

