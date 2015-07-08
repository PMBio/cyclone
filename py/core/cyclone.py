# Cell cycle classification and feature selection - perform feature selection and prediction with internal and external cross-validation; results are saved as figures and in a hdf5 file.
# 

########
######## PACKAGE LOADING
########
from __future__ import print_function
import h5py
import sys
import os
from os.path import exists
sys.path.append('./')
import scipy as SP
from numpy import inf
import pdb
import matplotlib as mpl
if 0:
	mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab as PL
from sklearn.decomposition import PCA

from sklearn.naive_bayes import GaussianNB
from sklearn import metrics, svm, ensemble, linear_model, preprocessing
from sklearn.grid_search import GridSearchCV
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.svm import SVC, l1_min_c
from sklearn.feature_selection import RFECV, SelectKBest

from sklearn import metrics, ensemble
from sklearn.cross_validation import StratifiedKFold, LeaveOneOut
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.datasets import make_classification
import warnings
warnings.filterwarnings('ignore')
# hacky hack we need for Windows users
#os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
	
class cyclone:
	"""
	cyclone
	This class takes care of assigning cell cycle stage based on cell cycle annotated genes in the  transcriptome.  
	"""
	def __init__(self,Y,row_namesY = None,cc_geneNames = None,labels = None, Y_tst = None, row_namesY_tst = None, labels_tst = None, learn_reverse = False, norm = 'rank'):	
		print("Initialise model...")
		assert row_namesY!=None, 'cyclone: provide gene names for Y'
		assert cc_geneNames!=None, 'cyclone: provide cell cycle genes'
		assert labels!=None, 'cyclone: provide labels for training'
		assert norm in ['rank', 'tc', 'median'] , 'normalisation has to be either "rank" (rank normalisation) or "tc" (total count) or "median" (log2 fold change to median)' 

		cc_ens = cc_geneNames
		genes = row_namesY 
		genesT = row_namesY_tst
		ind_cc_tr = SP.zeros((len(cc_ens),))
		if Y_tst!=None and row_namesY_tst!=None:
			ind_cc_tr = SP.zeros((len(cc_ens),))
			ind_cc_ts = SP.zeros((len(cc_ens),))
			for i in range(len(cc_ens)):
				ind_cc_ts[i] = SP.where(map(lambda x:x==cc_ens[i], genesT))[0][0]
				ind_cc_tr[i] = SP.where(map(lambda x:x==cc_ens[i], genes))[0][0]
			inds_tr = ind_cc_tr.astype('int32')
			inds_ts = ind_cc_ts.astype('int32')
			Y_tst = Y_tst[:,inds_ts]
			Y = Y[:,inds_tr]
			if norm == 'rank':
				for ir in range(Y_tst.shape[0]):
					Y_tst[ir,:] = SP.stats.rankdata(Y_tst[ir,:], method='average')		
				for ir in range(Y.shape[0]):
					Y[ir,:] = SP.stats.rankdata(Y[ir,:], method='average')
			elif norm=='tc':
		   		Y_tst = SP.transpose(SP.transpose(Y_tst)/SP.sum(Y_tst,1))	 
				Y = SP.transpose(SP.transpose(Y)/SP.sum(Y,1))	
			else:	
				for ir in range(Y_tst.shape[0]):
					med_ = SP.median(Y_tst[ir,:])
					Y_tst[ir,:] = SP.log2(Y_tst[ir,:]/med_)		
				for ir in range(Y.shape[0]):
					med_ = SP.median(Y[ir,:])
					Y[ir,:] = SP.log2(Y[ir,:]/med_)		
				Y[Y==-inf] = min(Y[Y>-inf])-1 
				Y_tst[Y_tst==-inf] = min(Y_tst[Y_tst>-inf])-1 
			if learn_reverse == True:	
				labels_ = labels_tst
				Y_ = Y_tst
				Y_tst = Y
				labels_tst = labels
				Y = Y_
				labels = labels_
			self.Y = Y
			self.Y_tst = Y_tst
			self.scores = None
			self.labels = labels
			self.labels_tst = labels_tst
			self.numClasses_tst = len(SP.unique((labels_tst)))
			self.inds_tr = inds_tr
			self.inds_tst = inds_ts
		else:
			ind_cc_tr = SP.zeros((len(cc_ens),))
			for i in range(len(cc_ens)):
				ind_cc_tr[i] = SP.where(map(lambda x:x==cc_ens[i], cc_genes))[0][0]
			inds_tr = ind_cc_tr.astype('int32')
			for ir in range(Y.shape[0]):
				Y[ir,:] = SP.stats.rankdata(Y[ir,:])
			self.Y = Y[:,inds_tr]
			slef.Y_tst=None
			
		self.geneNames = cc_geneNames
		self.numClasses = max(labels)#+1

	
	def trainModel(self, do_pca = False,out_dir='./cache', rftop = 40, class_labels = SP.array(['G1','S','G2M']), cv=10, npc=3, is_SVM=1, is_RFE=0 , scale=False):
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)	
		CFG = {}
		CFG['is_RFE'] = is_RFE # use recursive feature selection (can be slow for large datasets)
		CFG['is_SVM'] = is_SVM # use SVM with univariate feature selection (faster than RFE)
		CFG['CV_inner'] = cv #inner CV for RFE_CV: either an int or 'LOOCV'
		CFG['out_dir'] = out_dir
		CFG['do_pca'] = do_pca
		CFG['lassotop'] = 20
		self.cv = cv
		Y = self.Y
		labels = self.labels
		var_names = self.geneNames
		numClasses = self.numClasses
		predRF = SP.zeros((len(labels),numClasses))
		predSVM = SP.zeros((len(labels),numClasses)) 
		predSVMrbf = SP.zeros((len(labels),numClasses))
		predGNB = SP.zeros((len(labels),numClasses))
		predLR = SP.zeros((len(labels),numClasses))
		predLRall = SP.zeros((len(labels),numClasses))
		names_dict={}
		if self.cv == 'LOOCV':
			loo = LeaveOneOut(len(labels))
			CV_list = (list(iter(loo)))
			CV_list.append((SP.array(range(Y.shape[0])), SP.array(range(Y.shape[0]))))#all data...
		else:
			skf = StratifiedKFold(labels, n_folds=self.cv)
			CV_list = (list(iter(skf)))
			CV_list.append((SP.array(range(Y.shape[0])), SP.array(range(Y.shape[0]))))#all data...
		lambda_best = SP.zeros((1,len(CV_list))).ravel()	
		print("Performing cross validation ...")
		for i in range(len(CV_list)):
			if i<len(CV_list)-1:
				print("Fold " + str(i+1) + " of " + str(len(CV_list)-1))
			else:
				print("Final model")
			# string label for this fold
			#get data of a CV run
			cv_tr = CV_list[i][0]
			cv_tst = CV_list[i][1]
			lab_tr = labels[cv_tr]
			Ytr = Y[cv_tr,:]
			Ytst = Y[cv_tst,:]
			lab_tst = labels[cv_tst]
			if (i==len(CV_list)-1):
				foldlabel = 'full'
				if (self.Y_tst==None):
					Ytst = Y[cv_tst,:]
					lab_tst = labels[cv_tst]
				else:
					foldlabel = 'Test'
					Ytst = self.Y_tst
					lab_tst = self.labels_tst
			else:
				foldlabel = str(i)	
			if do_pca>=1:
				npc = npc#3	
				#do PCA to get features
				pcaCC = PCA(n_components=npc, whiten=False)
				pcaCC.fit(Ytr)
				pcaTst=pcaCC.transform(Ytst)
				pcaTr=pcaCC.transform(Ytr)
				#selection = SelectKBest(k=1)
				#combined_features = FeatureUnion([("pca", pcaCC), ("univ_select", selection)])
				combined_features = FeatureUnion([("pca", pcaCC)])
				gnb = GaussianNB()
				y_pred = gnb.fit(pcaTr, lab_tr).predict_proba(pcaTst)
				if i<len(CV_list)-1:
					predGNB[cv_tst,:] =y_pred#[:,1]
				else:
					predGNB_ts = y_pred#[:,1]
			if do_pca==2:
				Ytr = SP.concatenate((Ytr, pcaTr),1)
				Ytst = SP.concatenate((Ytst, pcaTst),1)
				pcnames = []
				for pci in range(npc):
					pcnames.append('PC'+str(pci+1))
				var_names = SP.concatenate((var_names, SP.array(pcnames)),1)				
			print("  Computing random forest ...")
			
			if CFG['is_RFE']==1:#Recursive feature selection with SVM
				print("  Computing RFE with SVM ...")
				svc = SVC(kernel="linear", probability=False, class_weight='auto')#use linear SVM for selection
				rfecv = RFECV(estimator=svc, step=1,scoring='f1')
				param_grid = dict(estimator__C=[0.1, 1, 10, 100, 1000])
				clf_rfe = GridSearchCV(rfecv, param_grid=param_grid, cv=3, scoring='f1')#GridSearch to find optimal parameters
				clf_rfe.fit(Ytr, lab_tr)
				svc = SVC(kernel="linear", probability=False,C=clf_rfe.best_estimator_.estimator.C, class_weight='auto')#use linear SVM for selection
				if CFG['CV_inner']=='':
					rfecv = RFECV(estimator=svc, step=1,scoring='f1')
				elif CFG['CV_inner']=='LOOCV':
					rfecv = RFECV(estimator=svc, step=1,scoring='f1', cv=LeaveOneOut(len(lab_tr)))
				else:
					rfecv = RFECV(estimator=svc, step=1,scoring='f1', cv=StratifiedKFold(lab_tr, n_folds=CFG['CV_inner']))
				clf_rfe.best_estimator_.fit(Ytr, lab_tr)
				predicted = clf_rfe.best_estimator_.predict(Ytst)
				if i<len(CV_list)-1:
					predSVM[cv_tst,:] = predicted
				else:
					predSVM_ts[cv_tst] = predicted
				classifier = svm.SVC(kernel='rbf', gamma=0.05, class_weight='auto', probability=True)#rbf kernel for prediction
				param_grid = dict(C=[0.1, 1], gamma=[1e-1,1e-2,1e-3])
				clf_rbf = GridSearchCV(classifier, param_grid=param_grid, cv=3, scoring='f1')
				clf_rbf.fit(Ytr[:,clf_rfe.best_estimator_.ranking_==1], lab_tr)
				clf_rbf.best_estimator_.fit(Ytr[:,clf_rfe.best_estimator_.ranking_==1], lab_tr)
				predicted = clf_rbf.best_estimator_.predict_proba(Ytst[:,clf_rfe.best_estimator_.ranking_==1])
				if i<len(CV_list)-1:
					predSVMrbf[cv_tst,:] = predicted
				fpr, tpr, thresholds = metrics.roc_curve(lab_tst, predicted[:,1])
				if (i==len(CV_list)-1) | CFG["CV_plots"]>0:
					PL.figure()
					PL.plot(fpr, tpr)
					PL.savefig(CFG['out_dir']+'/RF_SVM_'+foldlabel+'.pdf')
					names_dict[foldlabel+'_SVM']=self.geneNames[clf_rfe.best_estimator_.ranking_==1]			
			elif CFG['is_SVM']==1:#univariate FS with rbf SVM; choose this if you hava a large data set (many features, eg RNAseq)
				print("  SVM feature selection ...")
				classifier = svm.SVC(kernel='rbf', gamma=0.05, class_weight='auto', probability=True)
				selection = SelectKBest(k=1)
				combined_features = FeatureUnion([("univ_select", selection)])				
				
				X_features = combined_features.fit(Ytr, lab_tr).transform(Ytr)
				scaler = preprocessing.StandardScaler().fit(Ytr)
				YtrS = scaler.transform(Ytr)
				YtstS = scaler.transform(Ytst)
				
				classifier.fit(X_features, lab_tr)
				pipeline = Pipeline([("features", combined_features), ("svm", classifier)])
				if CFG['do_pca']==3:
					param_grid = dict(features__pca__n_components=SP.unique(SP.round_(SP.logspace(1.0,max(SP.log2(Ytr.shape[1]), SP.log2(10)),num=min(5,Ytr.shape[1]),base=2.0))),
									  features__univ_select__k=SP.unique(SP.round_(SP.logspace(3.0,SP.log2(Ytr.shape[1]),num=min(10,Ytr.shape[1]),base=2.0))),
									  svm__C=[0.1, 1, 10], svm__gamma=[1e-1,1e-2,1e-3])
				else:
					C_range = 10. ** SP.arange(0, 2)
					gamma_range = 10. ** SP.arange(-5, 1)
					param_grid = dict(features__univ_select__k=SP.unique(SP.round_(SP.logspace(3.0,SP.log2(Ytr.shape[1]),num=min(10,Ytr.shape[1]),base=2.0))),
									  svm__C=C_range, svm__gamma=gamma_range)
				clf = GridSearchCV(pipeline, param_grid=param_grid, cv=5, scoring='f1')
				clf.fit(YtrS, lab_tr)
				print("The best classifier is: ", clf.best_estimator_)
				select_best=clf.best_estimator_.get_params()['features__univ_select']
				#names_dict[foldlabel+'_SVM']=self.geneNames[SP.argsort(-1.0*select_best.scores_)[0:(select_best.k-1)]]
				expected = lab_tst
				predicted = clf.best_estimator_.predict_proba(YtstS)
				if i<len(CV_list)-1:
					predSVM[cv_tst,:] = predicted
				else:
					predSVM_ts = predicted
				#print(clf.best_estimator_)

				classifier = svm.SVC(kernel='rbf', gamma=0.05, class_weight='auto', probability=True)#rbf kernel for prediction
				param_grid = dict(C=[1,10], gamma=[ 1e-1,1e-2,1e-3])
				clf_rbf = GridSearchCV(classifier, param_grid=param_grid, cv=5, scoring='f1')
				clf_rbf.fit(Ytr, lab_tr)
				clf_rbf.best_estimator_.fit(Ytr, lab_tr)
				predicted = clf_rbf.best_estimator_.predict_proba(Ytst)
				if i<len(CV_list)-1:
					predSVMrbf[cv_tst,:] = predicted
				else:
					predSVMrbf_ts = predicted
		
			#do lasso with regularisation path
			cs = l1_min_c(Ytr, lab_tr, loss='log') * SP.logspace(0, 3)
			print("  Computing regularization path ...")
		
			lasso = linear_model.LogisticRegression(C=cs[0]*10.0, penalty='l1', tol=1e-6)
			param_grid = dict(C=cs)
			clf_lr = GridSearchCV(lasso, param_grid=param_grid, cv=5, scoring='f1')
			clf_lr.fit(Ytr, lab_tr)
			clf_lr.best_estimator_.fit(Ytr, lab_tr)
			lambda_best[i] = clf_lr.best_params_.get('C')
			predicted = clf_lr.best_estimator_.predict_proba(Ytst)

			clf = linear_model.LogisticRegression(C=cs[0]*10.0, penalty='l1', tol=1e-6)
			coefs_ = []
			for c in cs:
				clf.set_params(C=c)
				clf.fit(Ytr, lab_tr)
				coefs_.append(clf.coef_.ravel().copy())
		
		
			if i<len(CV_list)-1:
				predLR[cv_tst,:] = predicted
			else:
				predLR_ts = predicted
			coefs_ = SP.array(coefs_)
			# get ordering by importance (how many times they appear)
			order=(coefs_!=0).sum(axis=0).argsort()
			order=order[::-1] # descending
			# store this order
			featrank_lasso = order
			showtop= min(Ytr.shape[1], CFG['lassotop'])

			clfAll = linear_model.LogisticRegression(C=1e5, penalty='l2', tol=1e-6)
			clfAll.fit(Ytr, lab_tr)
			predicted = clfAll.predict_proba(Ytst)
			if i<len(CV_list)-1:
				predLRall[cv_tst,:] = predicted
			else:
				predLRall_ts = predicted
			forest = ExtraTreesClassifier(n_estimators=500,
										  random_state=0, criterion="entropy", bootstrap=False)
			#forest = RandomForestClassifier(n_estimators=500,
			#							 random_state=0, criterion="entropy")
			forest.fit(Ytr, lab_tr)
			pred = forest.predict_proba(Ytst)
			#pdb.set_trace()
			if i<len(CV_list)-1:
				predRF[cv_tst,:] = pred#[:,1]
			else:
				predRF_ts = pred#[:,1]
			importances = forest.feature_importances_
			std = SP.std([tree.feature_importances_ for tree in forest.estimators_],
						 axis=0)
		
			topfeat=min(Ytr.shape[1], rftop)
			indices = SP.argsort(importances)[::-1][0:topfeat]
			# store full feature ranking
			featrank_rf = SP.argsort(importances)[::-1]
			# Plot the feature importances of the forest
			if (i==len(CV_list)-1):
				PL.figure()
				#PL.title("Feature importances, Fold "+foldddPPlabel+', AUC='+str(SP.round_(metrics.auc(fpr, tpr),3)))
				PL.title("Feature importances")
				#PL.bar(range(topfeat), importances[indices],color="r", yerr=std[indices], align="center")
				PL.bar(range(topfeat), importances[indices],color="r", align="center")
				PL.xticks(range(topfeat), indices, rotation=70)
				PL.gca().set_xticklabels(var_names[indices])
				PL.setp(PL.gca().get_xticklabels(), fontsize=8)
				PL.xlim([-1, topfeat])
				PL.savefig(out_dir+'/RF_featureimportance_'+foldlabel+'.pdf')
	
	
		f2 = open(os.path.join(out_dir,'classification_reportCV.txt')  ,'w')
		
		predRFv = SP.argmax(predRF_ts,axis=1)+1
		predRF_trv = SP.argmax(predRF,axis=1)+1
		self.scores = predRF
		self.scores_tst = predRF_ts
		self.ranking = var_names[indices]

		predLRv = SP.argmax(predLR_ts,axis=1)+1
		predLR_trv = SP.argmax(predLR,axis=1)+1
		self.scoresLR = predLR
		self.scoresLR_tst = predLR_ts
		
		predLRallv = SP.argmax(predLRall_ts,axis=1)+1
		predLRall_trv = SP.argmax(predLRall,axis=1)+1
		self.scoresLRall = predLRall
		self.scoresLRall_tst = predLRall_ts
		if CFG['is_SVM']==1:
			predSVMv = SP.argmax(predSVM_ts,axis=1)+1
			predSVM_trv = SP.argmax(predSVM,axis=1)+1
			self.scoresSVM = predSVM
			self.scoresSVM_tst = predSVM_ts
			
			predSVMrbfv = SP.argmax(predSVMrbf_ts,axis=1)+1
			predSVMrbf_trv = SP.argmax(predSVMrbf,axis=1)+1
			self.scoresSVMrbf = predSVMrbf
			self.scoresSVMrbf_tst = predSVMrbf_ts

		predGNBv = SP.argmax(predGNB_ts,axis=1)+1
		predGNB_trv = SP.argmax(predGNB,axis=1)+1
		self.scoresGNB = predGNB
		self.scoresGNB_tst = predGNB_ts

		print("Classification report for classifier %s:\n%s\n" % ('Gaussian Naive Bayes', metrics.classification_report(self.labels, predGNB_trv)))
		print("Classification report for classifier %s:\n%s\n" % ('Random Forest', metrics.classification_report(self.labels, predRF_trv)))
		print("Classification report for classifier %s:\n%s\n" % ('LR', metrics.classification_report(self.labels, predLR_trv)))
		print("Classification report for classifier %s:\n%s\n" % ('LRall', metrics.classification_report(self.labels, predLRall_trv)))
		if CFG['is_RFE']==1:
			print("Classification report for classifier %s:\n%s\n" % ('SVM ', metrics.classification_report(labels, predSVM>0.5)),file=f2)
		elif CFG['is_SVM']==1:
			print("Classification report for classifier %s:\n%s\n" % ('SVM', metrics.classification_report(self.labels, predSVM_trv)))
			print("Classification report for classifier %s:\n%s\n" % ('SVMrbf', metrics.classification_report(self.labels, predSVMrbf_trv)))
		f2.close()
		

	def plotHistograms(self, class_labels = SP.array(['G1','S','G2M']),out_dir='./cache', plot_tst = True, method = 'RF',do_h=False): 
		plparams = {'backend': 'pdf',
		  'axes.labelsize': 14,
		  'text.fontsize': 14,
		  'legend.fontsize': 13,
		  'xtick.labelsize': 14,
		  'ytick.labelsize': 14,
		  'text.usetex': False}
		PL.rcParams.update(plparams)	
		assert self.scores != None, 'cyclone: first train the model before attempting to plot'
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		width =0.75
		ind_b = SP.arange(3)
		if method=='RF':
			predRF = self.scores
		elif  method=='SVM':
			predRF = self.scoresSVM
		elif  method=='SVMrbf':
			predRF = self.scoresSVMrbf
		elif  method=='GNB':
			predRF = self.scoresGNB
		elif  method=='LR':
			predRF = self.scoresLR
		elif  method=='LRall':
			predRF = self.scoresLRall
		predRF_trv = SP.argmax(predRF,axis=1)+1
		labels = self.labels
		if plot_tst ==True:
			if method=='RF':
				predRF_ts = self.scores_tst 
			elif method=='SVM':
				predRF_ts = self.scoresSVM_tst
			elif method=='SVMrbf':
				predRF_ts = self.scoresSVMrbf_tst
			elif method=='LR':
				predRF_ts = self.scoresLR_tst
			elif method=='LRall':
				predRF_ts = self.scoresLRall_tst
			elif method=='GNB':
				predRF_ts = self.scoresGNB_tst
			predRFv = SP.argmax(predRF_ts,axis=1)+1
			lab_tst = self.labels_tst
			u_classes = SP.unique(lab_tst)
			ind_b = SP.arange(len(u_classes))
			G1pred = SP.array([SP.sum(predRFv[lab_tst==iclass]==1) for iclass in u_classes])
			Spred = SP.array([SP.sum(predRFv[lab_tst==iclass]==2) for iclass in u_classes])
			G2Mpred = SP.array([SP.sum(predRFv[lab_tst==iclass]==3) for iclass in u_classes])

			nClass = SP.array([SP.sum(lab_tst==iclass).astype('double') for iclass in u_classes])

			G1pred_rel = SP.array([SP.sum(predRFv[lab_tst==iclass]==1)/nClass[iclass-1] for iclass in u_classes])
			Spred_rel = SP.array([SP.sum(predRFv[lab_tst==iclass]==2)/nClass[iclass-1] for iclass in u_classes])
			G2Mpred_rel = SP.array([SP.sum(predRFv[lab_tst==iclass]==3)/nClass[iclass-1] for iclass in u_classes])

			cols = ['#1b9e77', '#d95f02', '#7570b3']
			if do_h == True:
				fig_s = 2.0/3.0*len(class_labels)
				PL.figure(figsize=(5.0,fig_s))
				p1 = plt.barh(ind_b, G1pred, height = 0.9, color=cols[0])
				p2 = plt.barh(ind_b, Spred, color=cols[1], left=G1pred, height = 0.9)
				p3 = plt.barh(ind_b, G2Mpred, color=cols[2], left=G1pred+Spred, height = 0.9)
				PL.xlabel('# cells')
				PL.yticks(ind_b+width/2., class_labels )
			else:
				PL.figure(figsize=(5,5))
				p1 = plt.bar(ind_b, G1pred,   width, color=cols[0])
				p2 = plt.bar(ind_b, Spred,width, color=cols[1], bottom=G1pred)
				p3 = plt.bar(ind_b, G2Mpred,width, color=cols[2], bottom=G1pred+Spred)
				PL.ylabel('# cells')
				PL.title('True class')
				PL.xticks(ind_b+width/2., class_labels )
			#PL.yticks(SP.arange(0,20,5))
				# remove top and right axis
			ax = PL.gca()
			ax.spines["right"].set_visible(False)
			ax.spines["top"].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().tick_left()
			lgd = PL.legend( (p1[0], p2[0],p3[0]), ('Pred. G1', 'Pred. S','Pred. G2M') ,loc='center left', bbox_to_anchor=(1, 0.5))
			PL.savefig(out_dir+'/barchart_h'+str(do_h)+method+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
			if do_h == True:
				fig_s = 2.0/3.0*len(class_labels)
				PL.figure(figsize=(5.0,fig_s))
				p11 = plt.barh(ind_b, G1pred_rel, color=cols[0], height = 0.9)
				p22 = plt.barh(ind_b, Spred_rel, color=cols[1], left=G1pred_rel, height = 0.9)
				p33 = plt.barh(ind_b, G2Mpred_rel, color=cols[2], left=G1pred_rel+Spred_rel, height = 0.9)
				PL.xlabel('Fraction')
				PL.yticks(ind_b+width/2., class_labels)
			else:
				PL.figure(figsize=(5,5))
				p11 = plt.bar(ind_b, G1pred_rel,   width, color=cols[0])
				p22 = plt.bar(ind_b, Spred_rel,width, color=cols[1], bottom=G1pred_rel)
				p33 = plt.bar(ind_b, G2Mpred_rel,width, color=cols[2], bottom=G1pred_rel+Spred_rel)
				PL.ylabel('Fraction')
				PL.title('True cell cycle phase')
				PL.xticks(ind_b+width/2., class_labels)
			#PL.yticks(SP.arange(0,20,5))
			ax = PL.gca()
			ax.spines["right"].set_visible(False)
			ax.spines["top"].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().tick_left()
			lgd = PL.legend( (p11[0], p22[0],p33[0]), ('Pred. G1', 'Pred. S','Pred. G2M') ,loc='center left', bbox_to_anchor=(1, 0.5))
			PL.savefig(out_dir+'/barchart_rel_h'+str(do_h)+method+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')

		G1pred_tr = SP.array((SP.sum(predRF_trv[labels==1]==1),SP.sum(predRF_trv[labels==2]==1),SP.sum(predRF_trv[labels==3]==1)))
		Spred_tr = SP.array((SP.sum(predRF_trv[labels==1]==2),SP.sum(predRF_trv[labels==2]==2),SP.sum(predRF_trv[labels==3]==2)))
		G2Mpred_tr = SP.array((SP.sum(predRF_trv[labels==1]==3),SP.sum(predRF_trv[labels==2]==3),SP.sum(predRF_trv[labels==3]==3)))
		class_labels = ['G1', 'S', 'G2M']
		ind_b = SP.arange(3)
		if do_h == True:
			PL.figure(figsize=(5,2))
			p1cv = plt.barh(ind_b, G1pred_tr,   height=0.85, color=cols[0])
			p2cv = plt.barh(ind_b, Spred_tr,height=0.85, color=cols[1], left=G1pred_tr)
			p3cv = plt.barh(ind_b, G2Mpred_tr,height=0.85, color=cols[2], left=G1pred_tr+Spred_tr)
			PL.xlabel('# cells')
			PL.yticks(ind_b+width/2., class_labels )
		else:
			PL.figure(figsize=(5,5))
			p1cv = plt.barh(ind_b, G1pred_tr,   width, color=cols[0])
			p2cv = plt.barh(ind_b, Spred_tr,width, color=cols[1], left=G1pred_tr)
			p3cv = plt.barh(ind_b, G2Mpred_tr,width, color=cols[2], left=G1pred_tr+Spred_tr)
			PL.ylabel('# cells')
			PL.title('True cell cycle phase')
			PL.xticks(ind_b+width/2., class_labels )
		#PL.yticks(SP.arange(0,20,5))
		ax = PL.gca()
		ax.spines["right"].set_visible(False)
		ax.spines["top"].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		lgd = PL.legend( (p1cv[0], p2cv[0],p3cv[0]), ('Pred. G1', 'Pred. S','Pred. G2M') ,loc='center left', bbox_to_anchor=(1, 0.5))
		PL.savefig(out_dir+'/barchart_cv'+method+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')	
	
	def plotScatter(self, file_name = 'scatter_ts',out_dir = './cache', plot_test = True, xaxis = 0, yaxis = 2, xlab = 'G1 score', ylab = 'G2M score', class_labels = ['G1 phase','S phase','G2M phase'], method = 'RF', decision_lines=True):	
		plparams = {'backend': 'pdf',
		  'axes.labelsize': 14,
		  'text.fontsize': 14,
		  'legend.fontsize': 13,
		  'xtick.labelsize': 14,
		  'ytick.labelsize': 14,
		  'text.usetex': False}
		PL.rcParams.update(plparams)
		assert self.scores != None, 'cyclone: first train the model before attempting to plot'
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		file_name = file_name+method+'.pdf'
		if plot_test==False:
			file_name = file_name+method+'_cv.pdf'
		if plot_test ==True:
			if method=='RF':
				labs = self.labels_tst
				scores = self.scores_tst
			elif method=='LR':
				labs = self.labels_tst
				scores = self.scoresLR_tst
			elif method=='LRall':
				labs = self.labels_tst
				scores = self.scoresLRall_tst
			elif method=='GNB':
				labs = self.labels_tst
				scores = self.scoresGNB_tst
			elif method=='SVM':
				labs = self.labels_tst
				scores = self.scoresSVM_tst
			elif method=='SVMrbf':
				labs = self.labels_tst
				scores = self.scoresSVMrbf_tst
		else:	
			if method=='RF':
				labs = self.labels
				scores = self.scores
			elif method=='LR':
				labs = self.labels
				scores = self.scoresLR
			elif method=='LRall':
				labs = self.labels
				scores = self.scoresLRall
			elif method=='GNB':
				labs = self.labels
				scores = self.scoresGNB
			elif method=='SVM':
				labs = self.labels
				scores = self.scoresSVM
			elif method=='SVMrbf':
				labs = self.labels
				scores = self.scoresSVMrbf
		cols = ['r', 'b', 'g', 'y', 'Crimson', 'DeepPink','LightSalmon','Lime', 'Olive']
		cols_d = {}
		cols_d['G1'] = '#1b9e77'
		cols_d['S'] = '#d95f02'
		cols_d['G2M'] = '#7570b3'
		cols_d['mid'] = '#e7298a'
		cols_d['early'] = '#e6ab02'
		cols_d['late'] = '#66a61e'

		labs = labs.astype('int')
		lab_col = list()
		fig = PL.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		#ax=PL.axes([0.25,0.25,0.65,0.65])
		ax.set_position([0.1,0.1,0.7,0.7])
		hList =list()
		for iplot in range(len(labs)):
			#hList.append(PL.plot(scores[iplot,xaxis],scores[iplot,yaxis],'.',markersize=15,c=cols[labs[iplot]-1], alpha=0.75))
			if class_labels[labs[iplot]-1] in cols_d.keys():
				hList.append(PL.plot(scores[iplot,xaxis],scores[iplot,yaxis],'.',markersize=15,c=cols_d[class_labels[labs[iplot]-1]], alpha=0.75))
			else:
				hList.append(PL.plot(scores[iplot,xaxis],scores[iplot,yaxis],'.',markersize=15,c='#8c510a', alpha=0.75))
		PL.xlabel(xlab)
		PL.ylabel(ylab)
		if xlab == 'G1 score':
			PL.xlim(xmin = 0.0, xmax = .95)
			x_max = 0.75
		else:
			x_max = scores[:,xaxis].max() #+ 0.05
			PL.xlim(xmax = x_max+0.05)
		if ylab == 'G2M score':
			PL.ylim(ymin = 0.0, ymax = .95)
			y_max = 0.75
		else:
			y_max = scores[:,yaxis].max() #+ 0.05
			PL.ylim(ymax = y_max+0.05)

		if decision_lines ==True:
			x_min = 0.0
			y_min = 0.0
			h = 0.001
			xx, yy = SP.meshgrid(SP.arange(x_min, x_max, h),
							 SP.arange(y_min, y_max, h))
			zz = 1 - (xx + yy)
			Z = SP.argmax(SP.dstack((xx,yy,zz)), 2)
			PL.contour(xx, yy, Z, levels = [0,1])

		legH=list()
		u_classes = SP.unique(labs)
		for ileg in u_classes:
			legH.append(hList[SP.where(labs==ileg)[0][0]][0])
		lh=PL.legend(legH,class_labels,loc='upper center',bbox_to_anchor=(0.5, 1.15),ncol=3, numpoints=1,scatterpoints=1)
		lh.set_frame_on(False)
		ax.spines["right"].set_visible(False)
		ax.spines["top"].set_visible(False)
		ax.get_xaxis().tick_bottom()
		PL.savefig(out_dir+'/'+file_name,bbox_extra_artists=[lh])#,bbox_inches='tight')

		#print("Generating HDF5 result file ...")
		
		
		
		#fi = h5py.File(CFG['out_file'],'w')
		#fi['X'] = Y
		#fi['labels'] = labels
		#fi['predLR'] = predLR
		#fi['predLR_ts'] = predLR_ts
		#fi['predRF'] = predRF
		#fi['scores'] = scores
		#fi.close()
	def plotPerformance(self,perfType = 'ROC',plot_test=True, out_dir = './cache',class_labels = ['G1 phase','S phase','G2M phase'], method = 'RF'):
		plparams = {'backend': 'pdf',
		  'axes.labelsize': 14,
		  'text.fontsize': 14,
		  'legend.fontsize': 13,
		  'xtick.labelsize': 14,
		  'ytick.labelsize': 14,
		  'text.usetex': False}
		PL.rcParams.update(plparams)
		assert(perfType in ['ROC', 'PR'])
		if plot_test ==True:
			if method=='RF':
				labs = self.labels_tst
				scores = self.scores_tst
			elif method=='LR':
				labs = self.labels_tst
				scores = self.scoresLR_tst
			elif method=='LRall':
				labs = self.labels_tst
				scores = self.scoresLRall_tst
			elif method=='GNB':
				labs = self.labels_tst
				scores = self.scoresGNB_tst
			elif method=='SVM':
				labs = self.labels_tst
				scores = self.scoresSVM_tst
			elif method=='SVMrbf':
				labs = self.labels_tst
				scores = self.scoresSVMrbf_tst
		else:
			if method=='RF':
				labs = self.labels
				scores = self.scores
			elif method=='LR':
				labs = self.labels
				scores = self.scoresLR
			elif method=='LRall':
				labs = self.labels
				scores = self.scoresLRall
			elif method=='GNB':
				labs = self.labels
				scores = self.scoresGNB
			elif method=='SVM':
				labs = self.labels
				scores = self.scoresSVM
			elif method=='SVMrbf':
				labs = self.labels
				scores = self.scoresSVMrbf	
		PL.figure()
		#col = ['r', 'b', 'g']
		col = ['#1b9e77', '#d95f02', '#7570b3']
		aucList = SP.zeros((scores.shape[1],))
		for ind in range(scores.shape[1]):
			labels_i = labs.copy()
			scores_i = scores[:,ind]
			labels_i[labs==ind+1] = 1
			labels_i[labs!=ind+1] = 0
			if perfType == 'ROC':
				fpr, tpr, thresholds = metrics.roc_curve(labels_i, scores_i)
				aucList[ind] = metrics.auc(fpr, tpr)
			elif perfType == 'PR':
				fpr, tpr, thresholds = metrics.precision_recall_curve(labels_i, scores_i)
			PL.plot(fpr, tpr, '-', c=col[ind])
		if perfType=='ROC':
			PL.title("ROC curve ccClassify")
			PL.xlabel('FPR')
			PL.ylabel('TPR')
			leg_str = list()
			for i in range(len(class_labels)):
				leg_str.append(class_labels[i]+', AUC = '+str(SP.round_(aucList[i],3)))
			PL.legend(leg_str, loc='lower-right')
			ax = plt.gca()
			ax.spines["right"].set_visible(False)
			ax.spines["top"].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().tick_left()
			PL.savefig(out_dir+'/ROC_test'+str(plot_test)+'_'+method+'.pdf',bbox_inches='tight')
		else:
			PL.title("PR curve ccClassify")
			PL.xlabel('Precision')
			PL.ylabel('Recall')
			leg_str = list()
			for i in range(len(class_labels)):
				leg_str.append(class_labels[i])#+', AUC = '+str(SP.round_(aucList[i],3)))
			PL.legend(leg_str, loc='lower-right')
			ax = plt.gca()
			ax.spines["right"].set_visible(False)
			ax.spines["top"].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().tick_left()
			PL.savefig(out_dir+'/ROC_test'+str(plot_test)+'_'+method+'.pdf',bbox_inches='tight')

	def plotF1(self,plot_test=True, out_dir = './cache',class_labels = ['G1','S','G2M']):
		plparams = {'backend': 'pdf',
		  'axes.labelsize': 14,
		  'text.fontsize': 14,
		  'legend.fontsize': 13,
		  'xtick.labelsize': 14,
		  'ytick.labelsize': 14,
		  'text.usetex': False}
		PL.rcParams.update(plparams)
		fig = PL.figure(figsize=(6,6))
		ax = fig.add_subplot(111)
		#ax=PL.axes([0.25,0.25,0.65,0.65])
		ax.set_position([0.1,0.1,0.7,0.7])
		hList =list()
		cols_d = {}
		cols_d['G1'] = '#1b9e77'
		cols_d['S'] = '#d95f02'
		cols_d['G2M'] = '#7570b3'
		cols_d['mid'] = '#e7298a'
		cols_d['early'] = '#e6ab02'
		cols_d['late'] = '#66a61e'
		m_i = 0
		all_meth = ['GNB','RF', 'LR', 'LRall','SVM']
		for method in all_meth: 
			if plot_test ==True:
				if method=='RF':
					labs = self.labels_tst
					scores = self.scores_tst
				elif method=='LR':
					labs = self.labels_tst
					scores = self.scoresLR_tst
				elif method=='LRall':
					labs = self.labels_tst
					scores = self.scoresLRall_tst
				elif method=='GNB':
					labs = self.labels_tst
					scores = self.scoresGNB_tst
				elif method=='SVM':
					labs = self.labels_tst
					scores = self.scoresSVM_tst
				elif method=='SVMrbf':
					labs = self.labels_tst
					scores = self.scoresSVMrbf_tst
			else:
				if method=='RF':
					labs = self.labels
					scores = self.scores
				elif method=='LR':
					labs = self.labels
					scores = self.scoresLR
				elif method=='LRall':
					labs = self.labels
					scores = self.scoresLRall
				elif method=='GNB':
					labs = self.labels
					scores = self.scoresGNB
				elif method=='SVM':
					labs = self.labels
					scores = self.scoresSVM
				elif method=='SVMrbf':
					labs = self.labels
					scores = self.scoresSVMrbf	
			pred = SP.argmax(scores,axis=1)+1
			
			f1_scores = metrics.f1_score(labs, pred, average=None)
			f1_score_av = metrics.f1_score(labs, pred, average="macro")
	
			for iplot in range(len(f1_scores)):
				#hList.append(PL.plot(scores[iplot,xaxis],scores[iplot,yaxis],'.',markersize=15,c=cols[labs[iplot]-1], alpha=0.75))
				if class_labels[iplot] in cols_d.keys():
					hList.append(PL.plot(m_i,f1_scores[iplot],'^',markersize=15,c=cols_d[class_labels[iplot]], alpha=0.65))		
			PL.plot([m_i-0.25, m_i+0.25], [f1_score_av,f1_score_av], linewidth=2, color='r', hold=True)
			m_i+=1.0
		PL.xticks(range(len(all_meth)), all_meth, rotation=46)
		PL.ylabel('F1 score')
		PL.ylim(ymin = -0.05, ymax = 1.05)
		ax = plt.gca()
		ax.spines["right"].set_visible(False)
		ax.spines["top"].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()
		PL.savefig(out_dir+'/F1_test'+str(plot_test)+'_'+method+'.pdf',bbox_inches='tight')
