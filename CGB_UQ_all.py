#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data: rnaseq_FPKM_UQ_all.csv
# Method: Gradient_Boosting based on all gene features
# Date: 12/16/2018
# Name: Donglei Yin

import pandas as pd
import os
import numpy as np
import time


start=time.time()

#os.chdir("/home/donglei/Documents/Thesis/master/temp")

# Part 1: Data preparation

## Import TCGA normalized RNA-seq expression data

df=pd.read_csv("./rnaseq_FPKM_UQ_all.csv", low_memory=False,index_col=1)
df_tumor=df[df['sample_type']=='Tumor']

# gene expression values # 1109 by 56716
df_tumor_gene=df_tumor.iloc[:,1:56717]
df_tumor_gene.index.names=['barcode']

# clinical information   # 1109 by 82
df_tumor_clinical=df_tumor.iloc[:,56717:]
df_tumor_clinical=df_tumor_clinical.set_index('barcode')

## X: all gene features

data_x=df_tumor_gene

## Y: survival outcomes

data_y=df_tumor_clinical[["vital_status","days_to_last_follow_up","days_to_death"]]
time_to_event=[]
status=[]
for index, row in data_y.iterrows():
    if row['vital_status']=='dead':
        time_to_event.append(row['days_to_death'])
        status.append(1)
    elif row['vital_status']=='alive':
        time_to_event.append(row['days_to_last_follow_up'])
        status.append(0)
    else:
        time_to_event.append(None)
        status.append(None)

data_y=data_y.copy()
data_y['status']= pd.Series([x==1 for x in status], index=data_y.index)
data_y['time_to_event']= pd.Series(time_to_event, index=data_y.index)
data_y = data_y.drop(['days_to_last_follow_up','days_to_death','vital_status'], 1)

data_y = data_y.to_records(index=False)


## Check columus with na outcomes

# df=pd.DataFrame(data_y)
# df['index_col'] = df.index
# df_null = df[df.isnull().any(axis=1)]
# df_null

data_x=data_x.drop(data_x.index[[537,944,1070]])
data_y=np.delete(data_y,[537,944,1070])

## Correct the follow uo days less than 0 to 0
for idx, item in enumerate(data_y['time_to_event']):
    if item < 0:
        data_y['time_to_event'][idx] = 0
# data_y
# df.groupby('status').count()

# Part 2: Gradient Boosting for Survival Analysis

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import ShuffleSplit, GridSearchCV


from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.ensemble import GradientBoostingSurvivalAnalysis

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

## create estimator
estimator = GradientBoostingSurvivalAnalysis(loss='coxph', random_state=0,max_depth=3)

## define a function for evaluating the performance of models during grid search using Harrell's concordance index
def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['status'], y['time_to_event'], prediction)
    return result[0]

param_grid = {'learning_rate': [0.001,0.01,0.1,0.5,1], 'n_estimators': [100,200,500,1000]}


cv = ShuffleSplit(n_splits=100, test_size=0.3, random_state=0)

gcv = GridSearchCV(estimator, param_grid, scoring=score_survival_model,
                   n_jobs=4, iid=False, refit=False,
                   cv=cv)


gcv = gcv.fit(data_x, data_y)

print(gcv.best_score_)
print(gcv.best_params_)

end=time.time()

print(end-start)
