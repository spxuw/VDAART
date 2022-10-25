import copy
import datetime
import logging
import numpy as np
import os
import sys
import pdb
import torch

from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier

from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.linear_model import SGDClassifier

from sklearn.multiclass import OneVsOneClassifier 
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OutputCodeClassifier

from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB 

from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import RadiusNeighborsClassifier

from sklearn.svm import SVC

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import ExtraTreeClassifier
from sklearn.feature_selection import SelectKBest, f_classif

from sklearn.neural_network import MLPClassifier

from sklearn.model_selection import GridSearchCV,ParameterGrid


from pytorch_tabnet.tab_model import TabNetClassifier

# used for normalization
from sklearn.preprocessing import StandardScaler,MinMaxScaler

# used for cross-validation
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline

# this is an incredibly useful function
from pandas import read_csv
import pandas as pd 
import itertools

from train_test import train_test


def create_tabnet(n_d=32, n_a=32,n_steps=5, lr=0.02, gamma=1.5, 
                  n_independent=2, n_shared=2, lambda_sparse=1e-4, 
                  momentum=0.3, clip_value=2.):
    model = TabNetClassifier(
        n_d=n_d, n_a=n_a, n_steps=n_steps,
        gamma=gamma, n_independent=n_independent, n_shared=n_shared,
        lambda_sparse=lambda_sparse, momentum=momentum, clip_value=clip_value,
        optimizer_fn=torch.optim.Adam,
        scheduler_params = {"gamma": 0.95,"step_size": 20},
        scheduler_fn=torch.optim.lr_scheduler.StepLR, epsilon=1e-15, verbose = 0,
        mask_type = 'entmax'
    )
    return model


# classifiers
models = []
models.append(('LR', LogisticRegression(),
	{'penalty' : ['l2'],'C' : [100,10,1.0,0.1,0.01],'solver' : ['lbfgs','newton-cg','liblinear']}))
models.append(('LRCV', LogisticRegressionCV(),
	{'penalty' : ['l2'],'solver' : ['lbfgs','newton-cg','liblinear']}))
models.append(('KNN', KNeighborsClassifier(),
	{'leaf_size': [1,2,5,10,50], 'n_neighbors': [1,2,5,10,50], 'p': np.linspace(1,5,5)}))
models.append(('SVC', SVC(probability=True),{'C': [0.1,1,10,100], 'gamma': [1,0.1,0.01,0.001],'kernel': ['rbf', 'poly', 'sigmoid']}))
models.append(('RF', RandomForestClassifier(),
	{'bootstrap': [True, False],'max_depth': [10, 20, 50, 100, 500, None],'min_samples_leaf': [1, 2, 4],'min_samples_split': [2, 5, 10],
	'n_estimators': [200, 500, 1000, 2000]}))
models.append(('Boosting', GradientBoostingClassifier(),
	{'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[200, 500, 1000, 2000], 'max_depth': [10, 20, 50, 100, 500, None]}))
models.append(('AdaBoost', AdaBoostClassifier(),
	{'learning_rate':[0.15,0.1,0.05,0.01,0.005,0.001], 'n_estimators':[200, 500, 1000, 2000]}))
models.append(('Bagging', BaggingClassifier(), {'n_estimators':[200, 500, 1000, 2000]}))
models.append(('BernoulliNB', BernoulliNB(),{'alpha': [0.01, 0.1, 0.5, 1.0, 10.0]}))
models.append(('GaussianNB', GaussianNB(),{'var_smoothing': [1e-09,1e-10,1e-11]}))
models.append(('DecisionTree', DecisionTreeClassifier(),
	{'max_depth': [10, 20, 50, 100, 500, None], 'min_samples_split': [2, 5, 10, 15, 100], 'min_samples_leaf': [1, 2, 5, 10]}))
models.append(('ExtraTree', ExtraTreeClassifier(),
	{'max_depth': [10, 20, 50, 100, 500, None], 'min_samples_split': [2, 5, 10, 15, 100], 'min_samples_leaf': [1, 2, 5, 10]}))
models.append(('MLP', MLPClassifier(),
	{'hidden_layer_sizes': [(50,50,50), (50,100,50), (100,)],'activation': ['tanh', 'relu'],'solver': ['sgd', 'adam'],'alpha': [0.0001, 0.05],'learning_rate': ['constant','adaptive']}))

without_probability = ['Lasso','ElasticNet','ElasticNetCV','Ridge','RidgeCV','SGD','NearestCentroid','']

omics_combination_results = pd.DataFrame()

classifierTopFeatures = dict()
for fold in range(1):
	Train_dat =  read_csv('../filtered/external_all/Median_Imputed_Train_external.csv', header=0,sep=',')
	Test_dat =  read_csv('../filtered/external_all/Median_Imputed_Test_external.csv', sep=',')
	Train_label =  read_csv('../filtered/external_all/Imputed_Train_Label_external.csv', sep=',')
	Test_label =  read_csv('../filtered/external_all/Imputed_Test_Label_external.csv', sep=',')
	#feature_num = read_csv('../filtered/internal/Imputed_feature_num_Fold_'+str(fold)+'.csv', sep=',')

	Train_value = Train_dat.values
	Test_value = Test_dat.values
	y_train = Train_label['asthmawhz'].values
	y_test = Test_label['asthmawhz'].values
	y_train = np.multiply(y_train, 1)
	y_test = np.multiply(y_test, 1)
	#feature_num = np.cumsum(feature_num.values)
	biomarkers = Train_dat.columns

	gwas_train = Train_value[:,0:6]
	mirna_train = Train_value[:,6:306]
	mrna_train = Train_value[:,306:606]
	microbiome_train = Train_value[:,606:906]
	metabolomics_train = Train_value[:,906:1206]
	methylation_train = Train_value[:,1206:1506]


	gwas_test = Test_value[:,0:6]
	mirna_test = Test_value[:,6:306]
	mrna_test = Test_value[:,306:606]
	microbiome_test = Test_value[:,606:906]
	metabolomics_test = Test_value[:,906:1206]
	methylation_test = Test_value[:,1206:1506]

	Dict_train = {1:gwas_train, 2: mirna_train, 3: mrna_train, 4: microbiome_train, 5: metabolomics_train, 6: methylation_train}
	Dict_test = {1:gwas_test, 2: mirna_test, 3: mrna_test, 4: microbiome_test, 5: metabolomics_test, 6: methylation_test}

	N_test = len(y_test)
	N_classifier =  15

	a = [1,2,3,4,5,6]
	for nchoose in range(1,len(a)+1):
		comb = list(itertools.combinations(a,nchoose))
		for j in range(len(comb)):
			X_train = np.array([])
			X_test = np.array([])
			for k in range(nchoose):
				X_train = np.concatenate([np.transpose(Dict_train[comb[j][x]]) for x in range(nchoose)])
				X_train = np.transpose(X_train)
				X_test = np.concatenate([np.transpose(Dict_test[comb[j][x]]) for x in range(nchoose)])
				X_test = np.transpose(X_test)
			# end the feature combine
			pred_binary = np.array([])
			pred_probability = np.array([])
			for name, model,parameter_space in models:
				print(name)
				# transform data and fit the model
				scaler = StandardScaler()
				X_train = scaler.fit_transform(X_train)
				X_test = scaler.transform(X_test)
				classifier = copy.deepcopy(model)

				classifier = GridSearchCV(classifier, parameter_space, n_jobs=5, cv=3, scoring="roc_auc")
				classifier.fit(X_train,y_train)
				y_pred = classifier.predict(X_test)
				y_pred_prob = classifier.predict_proba(X_test)[:,1]

				# merge the classifiers
				pred_binary = np.concatenate((pred_binary,y_pred),axis=0)
				pred_probability = np.concatenate((pred_probability,y_pred_prob),axis=0)
				precision, recall, thresholds = precision_recall_curve(y_test, y_pred_prob)
				performance = [{'Method': name, 'Omics used': comb[j],'Fold': fold,'Accuracy': accuracy_score(y_test, y_pred), 'F1': f1_score(y_test, y_pred), 'AUROC': roc_auc_score(y_test, y_pred_prob), 'AUPRC': auc(recall, precision)}]
				results_sub = pd.DataFrame(performance)
				omics_combination_results = omics_combination_results.append(results_sub)

			
			# Tabnet
			M = X_train.shape[0]
			random_indices = np.random.choice(M, size=int(0.2*M), replace=False)
			X_val_tab = X_train[random_indices,:]
			X_train_tab =  X_train[np.setdiff1d(range(0,M),random_indices),:]
			y_val_tab = y_train[random_indices]
			y_train_tab =  y_train[np.setdiff1d(range(0,M),random_indices)]


			param_grid = dict(n_d = [8,10,16,32],
                  n_a = [24],
                  gamma = [0.5,1,1,5,2],
                  n_independent = [2]
                  )

			grid = ParameterGrid(param_grid)
			search_results = pd.DataFrame()
			val_AUC = 0  
			for params in grid:
				print(params)
				params['n_a'] = params['n_d']
				tabnet = create_tabnet()
				tabnet.set_params(**params)
				tabnet.fit(X_train=X_train_tab, y_train=y_train_tab,
    				eval_set=[(X_train_tab, y_train_tab),(X_val_tab,y_val_tab)],
    				eval_name = ['train','valid'],
    				eval_metric = ['auc'],
    				max_epochs=200, patience=20,
    				batch_size=20, virtual_batch_size=128,
    				drop_last = False)

				y_prob_val = tabnet.predict_proba(X_val_tab)  
				AUC_val = roc_auc_score(y_val_tab, y_prob_val[:,1])
				if val_AUC<AUC_val:
					y_pred_prob = tabnet.predict_proba(X_test)[:,1]
					y_pred = tabnet.predict(X_test)
					val_AUC = AUC_val

			pred_binary = np.concatenate((pred_binary,y_pred),axis=0)
			pred_probability = np.concatenate((pred_probability,y_pred_prob),axis=0)
			precision, recall, thresholds = precision_recall_curve(y_test, y_pred_prob)
			performance = [{'Method': 'Tabnet', 'Omics used': comb[j],'Fold': fold,'Accuracy': accuracy_score(y_test, y_pred), 'F1': f1_score(y_test, y_pred), 'AUROC': roc_auc_score(y_test, y_pred_prob), 'AUPRC': auc(recall, precision)}]
			results_sub = pd.DataFrame(performance)
			omics_combination_results = omics_combination_results.append(results_sub)

			# MONGNET
			num_epoch_pretrain = 500
			num_epoch = 2500
			lr_e_pretrain = 1e-3
			lr_e = 5e-4
			lr_c = 1e-3
			num_class = 2

			# prepare the input for MONGNET
			labels_tr = y_train
			labels_te = y_test
			data_tr_list = []
			data_te_list = []
			for k in range(nchoose):
				dat_training = Dict_train[comb[j][k]]
				scaler = MinMaxScaler()
				scaler.fit(X_train)
				dat_training = scaler.transform(X_train)
				dat_test = scaler.transform(X_test)
				data_tr_list.append(dat_training)
				data_te_list.append(dat_test)

			num_tr = data_tr_list[0].shape[0]
			num_te = data_te_list[0].shape[0]
			data_mat_list = []
			for i in range(nchoose):
				data_mat_list.append(np.concatenate((data_tr_list[i], data_te_list[i]), axis=0))
			data_tensor_list = []
			for i in range(len(data_mat_list)):
				data_tensor_list.append(torch.FloatTensor(data_mat_list[i]))
			idx_dict = {}
			idx_dict["tr"] = list(range(num_tr))
			idx_dict["te"] = list(range(num_tr, (num_tr+num_te)))
			data_train_list = []
			data_all_list = []
			for i in range(len(data_tensor_list)):
				data_train_list.append(data_tensor_list[i][idx_dict["tr"]].clone())
				data_all_list.append(torch.cat((data_tensor_list[i][idx_dict["tr"]].clone(),data_tensor_list[i][idx_dict["te"]].clone()),0))
			labels = np.concatenate((labels_tr, labels_te))

			data_tr_list = data_train_list
			data_trte_list = data_all_list
			trte_idx = idx_dict
			labels_trte = labels

			y_pred, y_pred_prob= train_test(data_tr_list, data_trte_list, trte_idx, labels_trte, nchoose, num_class,lr_e_pretrain, lr_e, lr_c,num_epoch_pretrain, num_epoch)
			precision, recall, thresholds = precision_recall_curve(y_test, y_pred_prob)
			performance = [{'Method': 'MONGNET', 'Omics used': comb[j],'Fold': fold,'Accuracy': accuracy_score(y_test, y_pred), 'F1': f1_score(y_test, y_pred), 'AUROC': roc_auc_score(y_test, y_pred_prob), 'AUPRC': auc(recall, precision)}]   
			results_sub = pd.DataFrame(performance)
			omics_combination_results = omics_combination_results.append(results_sub)
			
			# ensemble
			pred_binary = np.concatenate((pred_binary,y_pred),axis=0)
			pred_probability = np.concatenate((pred_probability,y_pred_prob),axis=0)
			pred_binary = np.reshape(pred_binary,(N_test,N_classifier))
			pred_probability = np.reshape(pred_probability,(N_test,N_classifier))
			pred_probability = np.mean(pred_probability,axis=1)
			pred_binary = np.sum(pred_binary,axis=1)
			pred_binary[pred_binary<7] = 0
			pred_binary[pred_binary>=7] = 1
			
			precision, recall, thresholds = precision_recall_curve(y_test, pred_probability)
			performance = [{'Method': 'Ensemble', 'Omics used': comb[j],'Fold': fold,'Accuracy': accuracy_score(y_test, pred_binary), 'F1': f1_score(y_test, pred_binary), 'AUROC': roc_auc_score(y_test, pred_probability), 'AUPRC': auc(recall, precision)}]
			results_sub = pd.DataFrame(performance)
			omics_combination_results = omics_combination_results.append(results_sub)

			omics_combination_results.to_csv('../results/Ensemble/Median_Python_Metrics_external.csv', index=False)	


Metric = pd.DataFrame({'AUROC': AUROC, 'AUPRC': AUPRC, 'F1': F_1, 'Accuracy': Accuracy},'Method': Method)
