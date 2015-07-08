from __future__ import print_function
import h5py
import sys
import os
sys.path.append('./')
sys.path.append('./CFG')
import scipy as SP

def load_data(CFG, is_Ens=True, gene_set='GOCB', het_only = True, het_onlyCB=True, pairs=False, filter_median = True, combine=False, filter_expressed = 0):
	f = h5py.File(CFG['train_file'],'r')
	Y = f['LogNcountsMmus'][:]
	labels = f['labels'][:].ravel()
	
	futil = h5py.File(CFG['util_file'],'r')
	Y_util = futil['LogNcountsQuartz'][:]
	
	ftst = h5py.File(CFG['test_file'],'r')
	if is_Ens ==True:
		genes = f['EnsIds'][:]
		genes_util = futil['gene_names_all'][:]
	else:
		genes = SP.char.lower(f['sym_names'][:])
		genes_util = SP.char.lower(futil['sym_namesQ'][:])

	#test file
	labels_util = futil['phase_vecS'][:]*2+futil['phase_vecG2M'][:]*3+futil['phase_vecG1'][:]
	if CFG['util_file']==CFG['test_file']:
		genes_tst = genes_util 
		YT = ftst['LogNcountsQuartz'][:]
		labels_tst = ftst['phase_vecS'][:]*2+ftst['phase_vecG2M'][:]*3+ftst['phase_vecG1'][:]
	elif is_Ens == False:
		ftst = h5py.File(CFG['test_file'],'r')
		YT = ftst['counts'][:]
		genes_tst = SP.char.lower(ftst['sym_names'][:])
		#genes_tst = ftst['ensIds'][:]
		#labels_tst = SP.array([1,1,1,1,1])#ftst['labels'][:].ravel() 
		labels_tst = ftst['labels'][:].ravel()
	elif is_Ens == True:
		ftst = h5py.File(CFG['test_file'],'r')
		YT = ftst['counts'][:]
		#genes_tst = ftst['sym_names'][:]
		genes_tst = ftst['ensIds'][:]
		#labels_tst = SP.array([1,1,1,1,1])#ftst['labels'][:].ravel() 
		labels_tst = ftst['labels'][:].ravel() 
	
	if 'class_labels' in ftst.keys():
		class_labels = ftst['class_labels'][:]
	else:
		class_labels = [i.astype('str') for i in labels_tst]
		class_labels = SP.sort(SP.unique(class_labels))
	heterogen_util = genes_util[SP.intersect1d(SP.where(Y_util.mean(0)>0)[0],SP.where(futil['genes_heterogen'][:]==1)[0])]
	heterogen_train = genes[SP.intersect1d(SP.where(Y.mean(0)>0)[0],SP.where(f['genes_heterogen'][:]==1)[0])]
	

	cellcyclegenes_GO = genes[SP.unique(f['cellcyclegenes_filter'][:].ravel() -1)] # idx of cell cycle genes
	cellcyclegenes_CB = genes[f['ccCBall_gene_indices'][:].ravel() -1]		# idxof cell cycle genes ...
	


	if SP.any(gene_set=='GOCB'):	
		cc_ens = SP.union1d(cellcyclegenes_GO,cellcyclegenes_CB)
	elif SP.any(gene_set=='GO'):
		cc_ens = cellcyclegenes_GO 
	elif SP.any(gene_set=='CB'):
		cc_ens = cellcyclegenes_CB 
	elif SP.any(gene_set=='all'):
		cc_ens = genes 
	else:
		#assert(gene_set in CFG.keys()), str(gene_set+' does not exist. Chose different gene set.')
		cc_ens = gene_set 

	
	if het_only==True:
		cc_ens = SP.intersect1d(cc_ens, heterogen_train)
		if pairs==True:
			Y = Y[:,SP.where(f['genes_heterogen'][:]==1)[0]]
			genes = genes[SP.where(f['genes_heterogen'][:]==1)[0]]
	if het_onlyCB==True:
		cc_ens = SP.intersect1d(cc_ens, heterogen_util)
	
	#filter_expressed = .2
	lod = 0
	if filter_expressed>0: 
		medY = SP.sum(Y>lod,0)*1.0
		idx_filter = (medY/SP.float_(Y.shape[0]))>filter_expressed
		Y = Y[:,idx_filter]
		genes = genes[idx_filter]
		
		#medY_tst = SP.sum(Y_tst>lod,0)
		#Y_tst = Y_tst[:,medY_tst>filter_expressed]
		#genes_tst = genes_tst[medY_tst>filter_expressed]		
		
		medY_util = SP.sum(Y_util>lod,0)
		idx_filter = (medY_util/SP.float_(Y_util.shape[0]))>filter_expressed
		Y_util = Y_util[:,idx_filter]
		genes_util = genes_util[idx_filter]		
	
	cc_ens = SP.intersect1d(cc_ens, genes)
	cc_ens = SP.intersect1d(cc_ens, genes_tst)
	cc_ens = SP.intersect1d(cc_ens, genes_util)
		
	if combine==True:
		genes = list(genes)
		genes_util = list(genes_util)
		genes_intersect = SP.intersect1d(genes,genes_util)
		cidx_tr = [ genes.index(x) for x in genes_intersect ]
		cidx_util = [genes_util.index(x) for x in genes_intersect]	
		genes = SP.array(genes)[cidx_tr]
		genes_util = SP.array(genes_util)[cidx_util]
		Y = SP.vstack([Y[:,cidx_tr],Y_util[:,cidx_util]])
		genes = genes_intersect
		labels = SP.hstack([labels, labels_util])				


	Y_tst = YT
	cc_data = {}
	cc_data['cc_ens'] = cc_ens
	cc_data['labels_tst'] = labels_tst	
	cc_data['labels'] = labels
	cc_data['genes_tst'] = genes_tst 
	cc_data['genes'] = genes 
	cc_data['Y'] = Y 
	cc_data['Y_test'] = Y_tst 
	cc_data['class_labels'] = class_labels 
	return cc_data
	


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()	
