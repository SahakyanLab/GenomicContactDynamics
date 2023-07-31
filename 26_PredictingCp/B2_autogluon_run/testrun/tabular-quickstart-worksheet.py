#!/usr/bin/env python
# coding: utf-8

# [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gidler/autogluon-tutorials/blob/main/tutorials/tabular_prediction/tabular-quickstart.ipynb)
# [![Open In SageMaker Studio Lab](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/gidler/autogluon-tutorials/blob/main/tutorials/tabular_prediction/tabular-quickstart.ipynb)
#

# In[ ]:


# Uncomment the code below and run this cell if AutoGluon is not yet installed in the kernel.
# !pip install autogluon==0.5.0  # These tutorials are based on AutoGluon v0.5.0 and might not work with different versions.


# # Predicting Columns in a Table - Quick Start
#
#
#
# Via a simple `fit()` call, AutoGluon can produce highly-accurate models to predict the values in one column of a data table based on the rest of the columns' values. Use AutoGluon with tabular data for both classification and regression problems. This tutorial demonstrates how to use AutoGluon to produce a classification model that predicts whether or not a person's income exceeds $50,000.
#
# To start, import AutoGluon's TabularPredictor and TabularDataset classes:

# In[1]:

from autogluon.tabular import TabularDataset, TabularPredictor
from timeit import default_timer as timer

start = timer()

# ...
#end = timer()
# print(end - start) # Time in seconds, e.g. 5.38091952400282


# Load training data from a [CSV file](https://en.wikipedia.org/wiki/Comma-separated_values) into an AutoGluon Dataset object. This object is essentially equivalent to a [Pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) and the same methods can be applied to both.

# In[2]:


model_id = ''

train_data = TabularDataset(
    'https://autogluon.s3.amazonaws.com/datasets/Inc/train.csv')
# subsample subset of data for faster demo, try setting this to much larger values
subsample_size = 20
train_data = train_data.sample(n=subsample_size, random_state=0)
train_data.head()


# Note that we loaded data from a CSV file stored in the cloud ([AWS s3 bucket](https://aws.amazon.com/s3/)), but you can you specify a local file-path instead if you have already downloaded the CSV file to your own machine (e.g., using [wget](https://www.gnu.org/software/wget/)).
# Each row in the table `train_data` corresponds to a single training example. In this particular dataset, each row corresponds to an individual person, and the columns contain various characteristics reported during a census.
#
# Let's first use these features to predict whether the person's income exceeds $50,000 or not, which is recorded in the `class` column of this table.

# In[3]:


label = 'class'
print("Summary of class variable: \n", train_data[label].describe())


# Now use AutoGluon to train multiple models:

# In[ ]:


# specifies folder to store trained models
save_path = 'agModels-predictClass' + '_' + str(model_id)
# predictor = TabularPredictor(label=label, path=save_path).fit(train_data, ag_args_fit={
#    'num_gpus': 1}, ag_args_ensemble={'num_folds_parallel': 3}, presets='best_quality', time_limit=72*60*60)


# hyperparameters = {
#    'GBM': {},
#    'XGB': {}
#    'CAT': {}
# }

predictor = TabularPredictor(label=label, path=save_path).fit(train_data, ag_args_fit={
    'num_cpus': 1}, presets='medium_quality', time_limit=72*60*60)

# predictor = TabularPredictor(label=label, path=save_path).fit(train_data, hyperparameters=hyperparameters, ag_args_fit={
#    'num_gpus': 1}, ag_args_ensemble={'num_folds_parallel': 3}, presets='best_quality', time_limit=72*60*60)

# Google Colab test run
# label = 'class'
# predictor = TabularPredictor(label=label).fit(train_data, hyperparameters=hyperparameters, ag_args_fit={
#    'num_gpus': 1}, ag_args_ensemble={'num_folds_parallel': 3}, presets='best_quality', time_limit=72*60*60)

end = timer()
print(end - start)  # Time in seconds, e.g. 5.38091952400282
