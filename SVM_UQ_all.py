
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data: rnaseq_FPKM_UQ_all.csv
# Method: SurvivalSVM based on all gene features(p=56716)
# Date: 12/12/2018
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


# Part 2: Fast Training of Support Vector Machines for Survival Analysis


import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import ShuffleSplit, GridSearchCV


from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastSurvivalSVM

## determine what the amount of censoring for this data is and plot the survival/censoring times
n_censored = data_y.shape[0] - data_y["status"].sum()
print("%.1f%% of records are censored" % (n_censored / data_y.shape[0] * 100))


val, bins, patches = plt.hist((data_y["time_to_event"][data_y["status"]],
                               data_y["time_to_event"][~data_y["status"]]),
                              bins=30, stacked=True)
plt.legend(patches, ["Time of Death", "Time of Censoring"])


plt.savefig('./UQ_all_hist.png')

## create estimator
estimator = FastSurvivalSVM(optimizer="rbtree", max_iter=1000, tol=1e-6, random_state=0)


pd.DataFrame(data_y)['status'].count()

## define a function for evaluating the performance of models during grid search using Harrell's concordance index
def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['status'], y['time_to_event'], prediction)
    return result[0]

param_grid = {'alpha': 2. ** np.arange(-12, 13, 2)}
cv = ShuffleSplit(n_splits=200, train_size=0.5, random_state=0)
gcv = GridSearchCV(estimator, param_grid, scoring=score_survival_model,
                   n_jobs=4, iid=False, refit=False,
                   cv=cv)

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
gcv = gcv.fit(data_x, data_y)

print(gcv.best_score_)
print(gcv.best_params_)

## Finally, we retrieve all 200 test scores for each parameter setting and visualize their distribution by box plots.
def plot_performance(gcv):
    n_splits = gcv.cv.n_splits
    cv_scores = []
    for i, params in enumerate(gcv.cv_results_["params"]):
        validation_scores = np.empty(n_splits, dtype=float)
        for j in range(n_splits):
            validation_scores[j] = gcv.cv_results_["split%d_test_score" % j][i]
        name = "%.5f" % params["alpha"]
        cv_scores.append((name, validation_scores))

    sns.boxplot(pd.DataFrame.from_items(cv_scores))
    _, xtext = plt.xticks()
    for t in xtext:
        t.set_rotation("vertical")

plot_performance(gcv)

plt.savefig('./FTSVM_UQ_all_box.png')


end=time.time()

print(end-start)

