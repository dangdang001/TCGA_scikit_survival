
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data: rnaseq_FPKM_UQ_all.csv
# Method: FastKernelSurvivalSVM based on PCA transformed features(n=300)
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


## Apply PCA to all gene features:

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# standardize gene features before applying PCA, make mean 0 and unit variance
sc=StandardScaler()
X_norm=sc.fit_transform(df_tumor_gene)

# apply PCA to gene features, select the first 300 PCs and apply the transformation to data:

pca=PCA(n_components=300)
data_x=pca.fit_transform(X_norm)

data_x = pd.DataFrame(data_x)


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


# Part 2: FastKernelSurvivalSVM


from sklearn.model_selection import ShuffleSplit, GridSearchCV
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastKernelSurvivalSVM
from sksurv.kernels import clinical_kernel

kernel_matrix = clinical_kernel(data_x)
kssvm = FastKernelSurvivalSVM(optimizer="rbtree", kernel="precomputed", random_state=0)

## define a function for evaluating the performance of models during grid search using Harrell's concordance index
def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['status'], y['time_to_event'], prediction)
    return result[0]

param_grid = {'alpha': [0.001,0.01,0.1,0.5,1,10,100,1000]}
cv = ShuffleSplit(n_splits=200, test_size=0.3, random_state=0)

kgcv = GridSearchCV(kssvm, param_grid, score_survival_model,
                    n_jobs=4, iid=False, refit=False,
                    cv=cv)
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
kgcv = kgcv.fit(kernel_matrix, data_y)


print(kgcv.best_score_)
print(kgcv.best_params_)

## Finally, we retrieve all 200 test scores for each parameter setting and visualize their distribution by box plots.
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")


def plot_performance(gcv):
    n_splits = gcv.cv.n_splits
    cv_scores = {"alpha": [], "test_score": [], "split": []}
    order = []
    for i, params in enumerate(gcv.cv_results_["params"]):
        name = "%.5f" % params["alpha"]
        order.append(name)
        for j in range(n_splits):
            vs = gcv.cv_results_["split%d_test_score" % j][i]
            cv_scores["alpha"].append(name)
            cv_scores["test_score"].append(vs)
            cv_scores["split"].append(j)
    df = pd.DataFrame.from_dict(cv_scores)
    _, ax = plt.subplots(figsize=(11, 6))
    sns.boxplot(x="alpha", y="test_score", data=df, order=order, ax=ax)
    _, xtext = plt.xticks()
    for t in xtext:
        t.set_rotation("vertical")

plot_performance(kgcv)

plt.savefig('KSVM_UQ_PCA.png',dpi=600,bbox_inches = "tight")


end=time.time()

print(end-start)
