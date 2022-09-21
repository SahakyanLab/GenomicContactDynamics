#!/usr/bin/env python
# coding: utf-8

from autogluon.tabular import TabularDataset, TabularPredictor
train_data = TabularDataset(
    '/project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/26_PredictingCp/out_makeTable/min2Mb_chr21_MLtbl.csv')
# subsample subset of data for faster demo, try setting this to much larger values
subsample_size = 500
train_data = train_data.sample(n=subsample_size, random_state=0)

label = 'label'
print("Summary of class variable: \n", train_data[label].describe())

# specifies folder to store trained models
save_path = 'agModels-quickrun_' + str(subsample_size)
predictor = TabularPredictor(label=label, path=save_path, eval_metric='accuracy', problem_type='multiclass').fit(
    train_data, ag_args_fit={'num_gpus': 1})
