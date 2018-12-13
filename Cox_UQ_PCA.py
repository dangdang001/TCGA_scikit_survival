
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data: rnaseq_FPKM_UQ_all.csv
# Method: Apply Cox-PH model (with Ridge penalty) based on PCA transformed features(n=300)
# Date: 12/13/2018
# Name: Donglei Yin


### Part 1: Preprocessing

import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


df=pd.read_csv("./rnaseq_FPKM_UQ_all.csv", low_memory=False,index_col=1)
df_tumor=df[df['sample_type']=='Tumor']

# gene expression values # 1109 by 56716
df_tumor_gene=df_tumor.iloc[:,1:56717]
df_tumor_gene.index.names=['barcode']

# clinical information   # 1109 by 82
df_tumor_clinical=df_tumor.iloc[:,56717:]
df_tumor_clinical=df_tumor_clinical.set_index('barcode')


# apply PCA to gene features

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# standardize gene features before applying PCA, make mean 0 and unit variance
sc=StandardScaler()
X_norm=sc.fit_transform(df_tumor_gene)

# apply PCA to gene features
# pca=PCA()
# X_pca=pca.fit_transform(X_norm)

# temp=pca.explained_variance_ratio_.cumsum()

# plt.plot(temp)
# plt.xlabel('Number of Principle Components')
# plt.ylabel('Proportion of Variance Explained')
# plt.show();
 
# select the first 300 PCs and apply the transformation to data:

pca=PCA(n_components=300)
X_pca=pca.fit_transform(X_norm)

temp=pca.explained_variance_ratio_.cumsum()

data_x=X_pca

# survival outcomes
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


df=pd.DataFrame(data_y)
df['index_col'] = df.index
df_null = df[df.isnull().any(axis=1)]
df_null

## Check columus with na outcomes

df=pd.DataFrame(data_y)
df['index_col'] = df.index
df_null = df[df.isnull().any(axis=1)]
df_null

data_x=np.delete(data_x,[537,944,1070],axis=0)
data_y=np.delete(data_y,[537,944,1070])

## Correct the follow uo days less than 0 to 0
for idx, item in enumerate(data_y['time_to_event']):
    if item < 0:
        data_y['time_to_event'][idx] = 0
data_y

df.groupby('status').count()

### Part 2: Cox-PH model with Ridge Penalty

# tuning parameter alpha over grid search according to c-index

from sksurv.linear_model import CoxPHSurvivalAnalysis
from sklearn.model_selection import GridSearchCV

coxph = CoxPHSurvivalAnalysis()
grid_values={'alpha': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]}

grid_c = GridSearchCV(coxph, param_grid = grid_values, scoring = None)
grid_c.fit(data_x, data_y) 

print('Grid best parameter (max c-index): ', grid_c.best_params_)
print('Grid best score (c-index): ', grid_c.best_score_)


# Apply Cox-PH model based on 5-fold 10-repeated CV using optimal alpha selected from grid search:

from sklearn.model_selection import RepeatedKFold

rkf = RepeatedKFold(n_splits=5, n_repeats=10, random_state=0) # 5-fold 10-repeated CV

c_index_train,c_index_test=[],[]

for train_index, test_index in rkf.split(data_x):
    x_train, x_test = data_x[train_index], data_x[test_index]
    y_train, y_test = data_y[train_index], data_y[test_index]
    coxph = CoxPHSurvivalAnalysis(alpha=float(grid_c.best_params_['alpha'])).fit(x_train, y_train)
    c_index_train.append(coxph.score(x_train,y_train))
    c_index_test.append(coxph.score(x_test,y_test))
    
print("Averaged c-index from 5-fold 10 repeated CV(training): {:.3f}".format(np.mean(c_index_train)))
print("Averaged c-index from 5-fold 10 repeated CV(test): {:.3f}".format(np.mean(c_index_test)))




