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


directory = sys.argv[1]
outer_fold = int(sys.argv[2])

nr_out_folds = int(sys.argv[3])
nr_in_folds = int(sys.argv[4])

chooseWithProtein = int(sys.argv[5])


directory_suffix = directory.split('/')[1]

print directory_suffix

[similarity, ontology, classifier, lsdr_method, nrPrototypes, label_fraction] = gr.decomposeName(directory_suffix)

print 'My parameters are: ',
print similarity,
print ontology,
print nr_out_folds,
print nr_in_folds,
print classifier,
print lsdr_method,
print nrPrototypes,
print label_fraction,
print outer_fold

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

fields = bestFile.split('_')

if lsdr_method == "None":
	k = int(fields[2].split('.')[0])


else:
	if lsdr_method == "CPLST":
		k = int(fields[2])
		var_scale = int(fields[3])
		nr_dims = int(fields[4])

		lamda = float(fields[5].split('.')[0])

	else:
		k = int(fields[2])
		var_scale = int(fields[3])
		nr_dims = int(fields[4].split('.')[0])



print 'Loading data...'
[X, Y, termNames, geneNames, parentsCoord] = gr.loadData(similarity, ontology, organism)


np.random.seed(104)
outer_seed = 4294967295 * np.random.rand(1)
inner_seed = 4294967295 * np.random.rand(nr_out_folds)


outer_cv = ms.KFold(n_splits=nr_out_folds, random_state=int(outer_seed[0]), shuffle=True)

i = -1


for outer_train, outer_test in outer_cv.split(X):
	i += 1

	if i != outer_fold:
		continue

	Xtrain = X[outer_train,:][:,outer_train]
	Ytrain = Y[outer_train,:]


	[Xtrain, Ytrain, tokeep, _] = gr.removeRareTerms(Xtrain, Ytrain, label_fraction)



	Xtest = X[outer_test,:][:,outer_train]
	Ytest = Y[outer_test,:][:, tokeep]

	if nrPrototypes > 0:
		np.random.seed(892374980)
		indices = np.random.permutation(Xtrain.shape[0])[:nrPrototypes]
		Xtrain = Xtrain[:, indices]
		Xtest = Xtest[:, indices]


	newTermNames = []

	for trmInd in tokeep:
		newTermNames.append(termNames[trmInd])

	parentsCoord = gr.getParentsCoord(newTermNames, ontology)

	print 'ic calc...'
	ic = gr.calculateIC(Ytrain, parentsCoord)


	print 'dim red...'

	if lsdr_method == 'None':
		Ytrain_prime = Ytrain

	elif lsdr_method == 'PLST':
		mapping = lsdr.PLST()
		Ytrain_prime = mapping.fit(Ytrain, ndims = nr_dims, var_correction = var_scale)

	elif lsdr_method == 'CPLST':
		mapping = lsdr.CPLST()
		Ytrain_prime = mapping.fit(Xtrain, Ytrain, ndims = nr_dims, var_correction = var_scale, lamda=lamda)

	elif lsdr_method == 'GOAT':
		mapping = lsdr.GOAT(ontology)
		Ytrain_prime = mapping.fit(Ytrain, ndims = nr_dims, var_correction = var_scale, cols=tokeep)

	elif lsdr_method == 'SEM':
		mapping = lsdr.SEM(ontology)
		Ytrain_prime = mapping.fit(Ytrain, newTermNames, ndims = nr_dims, var_correction = var_scale)

	elif lsdr_method == 'DAG':
		mapping = lsdr.DAG(ontology)
		Ytrain_prime = mapping.fit(Ytrain, newTermNames, ndims = nr_dims, var_correction = var_scale)


	else:
		print 'Invalid LSDR method'
		exit()


	print 'classification'
	#classification
	Ytest_pred_prime = gr.predict(Xtrain, Ytrain_prime, Xtest, classifier, k)


	if lsdr_method == 'None':
		Ytest_pred = Ytest_pred_prime

	else:
		Ytest_pred = mapping.inverseMap(Ytest_pred_prime)


	gr.respectGO(Ytest_pred, parentsCoord)


print 'eval...'


[_, Ytest, tokeep, _] = gr.removeRareTerms(Xtest, Ytest, label_fraction)


ypred = gr.transform(Ytest_pred)

ypred = ypred[:, tokeep]


ic = ic[tokeep]



(nrProteins, nrTerms) = Ytest.shape

thresholds = np.linspace(0.0, 1.0, 51)

prot_fm = np.zeros((thresholds.shape[0], nrProteins), float)
prot_pr = np.zeros((thresholds.shape[0], nrProteins), float)
prot_rc = np.zeros((thresholds.shape[0], nrProteins), float)

wprot_fm = np.zeros((thresholds.shape[0], nrProteins), float)
wprot_pr = np.zeros((thresholds.shape[0], nrProteins), float)
wprot_rc = np.zeros((thresholds.shape[0], nrProteins), float)



prot_sd = np.zeros((thresholds.shape[0], nrProteins), float)

term_fm = np.zeros((thresholds.shape[0], nrTerms), float)
term_pr = np.zeros((thresholds.shape[0], nrTerms), float)
term_rc = np.zeros((thresholds.shape[0], nrTerms), float)

wterm_fm = np.zeros((thresholds.shape[0], nrTerms), float)



print 'term_auc'
term_auc = skm.average_precision_score(Ytest, ypred, average=None)


print 'term_roc'
term_roc = skm.roc_auc_score(Ytest, ypred, average=None)

#term_auc = np.nanmean(term_auc)

print 'prot_auc'
#prot_auc = skm.average_precision_score(Ytest, ypred, average='samples')
prot_auc = skm.average_precision_score(Ytest.T, ypred.T, average=None)

for j, th in enumerate(thresholds):

	print th

	ythreshold = (ypred >= th).astype(int)

	#[prot_pr[j], prot_rc[j], prot_fm[j], _] = skm.precision_recall_fscore_support(Ytest, ythreshold, average='samples', pos_label=1)
	[prot_pr[j], prot_rc[j], prot_fm[j], _] = skm.precision_recall_fscore_support(Ytest.T, ythreshold.T, average=None, pos_label=1)
	[wprot_pr[j], wprot_rc[j], wprot_fm[j], _] = skm.precision_recall_fscore_support(Ytest.T, ythreshold.T, average=None, pos_label=1, sample_weight=ic)


	[term_pr[j, :], term_rc[j, :], term_fm[j, :], _] = skm.precision_recall_fscore_support(Ytest, ythreshold, average=None, pos_label=1)


	[_, _, prot_sd[j,:]] = gr.semanticDistance(Ytest, ythreshold, ic)

	#prot_sd[j] = np.mean(sd)


prefix = '/fold' + str(outer_fold)

print 'saving...'

if lsdr_method == 'None':
	with open('out_scores/gt/gt' + str(outer_fold) + '_' + str(label_fraction) + '.pkl', 'wb') as f:
		pickle.dump(Ytest, f)

	with open('out_scores/gt/ic' + str(outer_fold) + '_' + str(label_fraction) + '.pkl', 'wb') as f:
		pickle.dump(ic, f)


with open('out_scores/' + directory_suffix + prefix + '_predictions.pkl', 'wb') as f:
	pickle.dump(ypred, f)

#print term_fm.shape

with open('out_scores/' + directory_suffix + prefix + '_pfm.pkl', 'wb') as f:
	pickle.dump(prot_fm, f)

with open('out_scores/' + directory_suffix + prefix + '_ppr.pkl', 'wb') as f:
	pickle.dump(prot_pr, f)

with open('out_scores/' + directory_suffix + prefix + '_prc.pkl', 'wb') as f:
	pickle.dump(prot_rc, f)

with open('out_scores/' + directory_suffix + prefix + '_pwfm.pkl', 'wb') as f:
	pickle.dump(wprot_fm, f)

with open('out_scores/' + directory_suffix + prefix + '_pwpr.pkl', 'wb') as f:
	pickle.dump(wprot_pr, f)

with open('out_scores/' + directory_suffix + prefix + '_pwrc.pkl', 'wb') as f:
	pickle.dump(wprot_rc, f)


with open('out_scores/' + directory_suffix + prefix + '_tfm.pkl', 'wb') as f:
	pickle.dump(term_fm, f)

with open('out_scores/' + directory_suffix + prefix + '_tpr.pkl', 'wb') as f:
	pickle.dump(term_pr, f)

with open('out_scores/' + directory_suffix + prefix + '_trc.pkl', 'wb') as f:
	pickle.dump(term_rc, f)

with open('out_scores/' + directory_suffix + prefix + '_psd.pkl', 'wb') as f:
	pickle.dump(prot_sd, f)

with open('out_scores/' + directory_suffix + prefix + '_pauc.pkl', 'wb') as f:
	pickle.dump(prot_auc, f)

with open('out_scores/' + directory_suffix + prefix + '_tauc.pkl', 'wb') as f:
	pickle.dump(term_auc, f)

with open('out_scores/' + directory_suffix + prefix + '_troc.pkl', 'wb') as f:
	pickle.dump(term_roc, f)
