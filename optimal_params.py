import sys
import filereader as fr
import lsdr
from scipy.io import loadmat
from scipy.sparse import csr_matrix
import sklearn.model_selection as ms
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import average_precision_score
import numpy as np
import pickle
import general_routines as gr
import os
import sklearn.metrics as skm

#prints the method parameters that had the best performance in the inner loop,
#for each of the 3 folds 


directory = sys.argv[1]

outer_folds = range(3)

nr_out_folds = 3
nr_in_folds = 3

chooseWithProtein = 1


directory_suffix = directory.split('/')[1]

[similarity, ontology, classifier, lsdr_method, nrPrototypes, label_fraction] = gr.decomposeName(directory_suffix)

organism = 'thalia'

files = sorted(os.listdir(directory))

nrFiles = len(files)
totalFolds = nr_in_folds * nr_out_folds

assert (nrFiles % totalFolds) == 0

nrCombos = len(files) / totalFolds

combos = []

for fil in files:
	if fil[:3] == '0_1':
		break
	
	combos.append(fil[3:])
	
assert len(combos) == nrCombos

#print combos

aucs = np.zeros((2, nrCombos, nr_in_folds), float)


for outer_fold in outer_folds:
	for i in range(len(combos)):
		for j in range(nr_in_folds):

			try:
				with open(directory + str(outer_fold) + '_' + str(j) + combos[i]) as f:

					for k, line in enumerate(f):
				               	#k=0 protein, k=1 term
						aucs[k, i, j] = float(line)
		                                
			except IOError:
			 	print 'exception'
				aucs[0, i, j] = np.nan
				aucs[1, i, j] = np.nan


	if chooseWithProtein:
		optimalHP_ind = np.argmax(np.nanmean(aucs[0,:,:],1))
	else:
		optimalHP_ind = np.argmax(np.nanmean(aucs[1,:,:],1))	

	bestFile = files[optimalHP_ind]

	print bestFile
