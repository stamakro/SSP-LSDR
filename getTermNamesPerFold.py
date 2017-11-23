import filereader as fr
from scipy.io import loadmat
import numpy as np
import general_routines as gr
import os
import sklearn.model_selection as ms
import pickle

print 'Loading data...'
[X, Y, termNames, geneNames, parentsCoord] = gr.loadData('sequence', 'P', 'thalia')


np.random.seed(104)
outer_seed = 4294967295 * np.random.rand(1)

outer_cv = ms.KFold(n_splits=3, random_state=int(outer_seed[0]), shuffle=True)

i = -1


label_fraction = 0.1

for outer_train, outer_test in outer_cv.split(X):
	i += 1

	#if i != outer_fold:
	#	continue

	Xtrain = X[outer_train,:][:,outer_train]
	Ytrain = Y[outer_train,:]

	[_, Ytrain, tokeep, _] = gr.removeRareTerms(Xtrain, Ytrain, label_fraction)

	Ytest = Y[outer_test,:][:, tokeep]

	newTermNames = []

	for trmInd in tokeep:
		newTermNames.append(termNames[trmInd])



	[_, Ytest, tokeep, _] = gr.removeRareTerms(Xtrain, Ytest, label_fraction)

	newTermNames2 = []

	for trmInd in tokeep:
		newTermNames2.append(newTermNames[trmInd])

	with open('out_scores/gt/terms' + str(i) + '_0.1.pkl', 'w') as f:
		pickle.dump(newTermNames2, f)
