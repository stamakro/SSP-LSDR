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
from sklearn.utils import resample

def getIndices(termsTrain, termsTest):

	ttr = set(termsTrain)
	tts = set(termsTest)

	trInd = []
	tsInd = []

	train2test = dict()

	for t in ttr.intersection(tts):
		trInd.append(termsTrain.index(t))
		tsInd.append(termsTest.index(t))

	for termindTrain, termindTest in zip(trInd, tsInd):
		train2test[termindTrain] = termindTest


	return [np.array(trInd), np.array(tsInd), train2test]



directory = sys.argv[1]

nr_in_folds = int(sys.argv[2])


directory_suffix = directory.split('/')[1]

[similarity, ontology, classifier, lsdr_method, _, _] = gr.decomposeName(directory_suffix)


organism = 'thalia'

files = sorted(os.listdir(directory))

nrFiles = len(files)


assert (nrFiles % nr_in_folds) == 0

nrCombos = len(files) / nr_in_folds

combos = []

for fil in files:
	if fil[:2] == '1_':
		break

	combos.append(fil[2:])

assert len(combos) == nrCombos

#print combos

aucs = np.zeros((nrCombos, nr_in_folds), float)


for i in range(len(combos)):
        for j in range(nr_in_folds):

		try:

			with open(directory + str(j) + '_' + combos[i]) as f:

				for k, line in enumerate(f):
		                       	aucs[i, j] = float(line)

		except IOError:
	         	print 'exception'
			aucs[i, j] = np.nan
			aucs[i, j] = np.nan


optimalHP_ind = np.argmax(np.nanmean(aucs,1))


bestFile = files[optimalHP_ind]

print bestFile


fields = bestFile.split('_')

print fields

if lsdr_method == "None":
	k = int(fields[1].split('.')[0])


else:
	if lsdr_method == "CPLST":
		k = int(fields[1])
		var_scale = int(fields[2])
		nr_dims = int(fields[3])

		lamda = float(fields[4].split('.')[0])

	else:
		k = int(fields[1])
		var_scale = int(fields[2])
		nr_dims = int(fields[3].split('.')[0])



print 'Loading data...'

X = loadmat('../cafa/analysis/X.mat')['identityMatrix']

Xtrain = X[:,:6077][:6077, :]

Xtest = X[:,:6077][6077:, :]


with open('../cafa/benchmark20170605/labelMatrixTest.pkl', 'rb') as f:
	Ytest = pickle.load(f)

with open('../cafa/benchmark20170605/proteinsTest.pkl', 'rb') as f:
	proteinsTest = pickle.load(f)

with open('../cafa/benchmark20170605/termsTest.pkl', 'rb') as f:
	termsTest = pickle.load(f)



with open('../cafa/training_data/labelMatrixTrain.pkl', 'rb') as f:
	Ytrain = pickle.load(f)

with open('../cafa/training_data/proteinsTrain.pkl', 'rb') as f:
	proteinsTrain = pickle.load(f)

with open('../cafa/training_data/termsTrain.pkl', 'rb') as f:
	termsTrain = pickle.load(f)


[trainIndices, testIndices, train2test] = getIndices(termsTrain, termsTest)


parentsCoord = gr.getParentsCoord(termsTrain, ontology)

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
	Ytrain_prime = mapping.fit(Ytrain, ndims = nr_dims, var_correction = var_scale)

elif lsdr_method == 'SEM':
	mapping = lsdr.SEM(ontology)
	Ytrain_prime = mapping.fit(Ytrain, termsTrain, ndims = nr_dims, var_correction = var_scale)


elif lsdr_method == 'DAG':
	mapping = lsdr.DAG(ontology)
	Ytrain_prime = mapping.fit(Ytrain, termsTrain, ndims = nr_dims, var_correction = var_scale)


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

print Ytest_pred.shape

print 'eval...'


#[_, Ytest, tokeep, _] = gr.removeRareTerms(Xtest, Ytest, label_fraction)


ypred_temp = gr.transform(Ytest_pred)


ypred = np.zeros((ypred_temp.shape[0], Ytest.shape[1]), float)

ic2 = np.zeros((Ytest.shape[1],), float)

for i in range(ypred_temp.shape[1]):
	if i in train2test:
		ypred[:, train2test[i]] = ypred_temp[:, i]
		ic2[train2test[i]] = ic[i]


print ypred.shape
print Ytest.shape

ypred = ypred[:, testIndices]
Ytest = Ytest[:, testIndices]
ic = ic2[testIndices]


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

wterm_fm = np.zeros((thresholds.shape[0], ), float)



print 'term_auc'
term_auc = skm.average_precision_score(Ytest, ypred, average=None)

wterm_auc = np.average(term_auc, weights=ic)


print 'term_roc'
term_roc = skm.roc_auc_score(Ytest, ypred, average=None)


wterm_roc = np.average(term_roc, weights=ic)

#term_auc = np.nanmean(term_auc)

print 'prot_auc'
#prot_auc = skm.average_precision_score(Ytest, ypred, average='samples')
prot_auc = skm.average_precision_score(Ytest.T, ypred.T, average=None)



for j, th in enumerate(thresholds):

	#print th

	ythreshold = (ypred >= th).astype(int)


	[prot_pr[j], prot_rc[j], prot_fm[j], _] = skm.precision_recall_fscore_support(Ytest.T, ythreshold.T, average=None, pos_label=1)
	[wprot_pr[j], wprot_rc[j], wprot_fm[j], _] = skm.precision_recall_fscore_support(Ytest.T, ythreshold.T, average=None, pos_label=1, sample_weight=ic)


	[term_pr[j, :], term_rc[j, :], term_fm[j, :], _] = skm.precision_recall_fscore_support(Ytest, ythreshold, average=None, pos_label=1)

	wterm_fm[j] = np.average(term_fm[j, :], weights=ic)

	[_, _, prot_sd[j,:]] = gr.semanticDistance(Ytest, ythreshold, ic)



nrBootstraps = 1000
np.random.seed(1002003445)
seedonia = np.random.randint(low=0, high=4294967295, size=nrBootstraps)


bootstraps_pfm = np.zeros((nrBootstraps,), float)
bootstraps_ppr = np.zeros((nrBootstraps,), float)
bootstraps_prc = np.zeros((nrBootstraps,), float)

bootstraps_pwfm = np.zeros((nrBootstraps,), float)
bootstraps_pwpr = np.zeros((nrBootstraps,), float)
bootstraps_pwrc = np.zeros((nrBootstraps,), float)

bootstraps_psd = np.zeros((nrBootstraps,), float)
bootstraps_pauc = np.zeros((nrBootstraps,), float)

bootstraps_tfm = np.zeros((nrBootstraps,), float)
bootstraps_tpr = np.zeros((nrBootstraps,), float)
bootstraps_trc = np.zeros((nrBootstraps,), float)

bootstraps_tauc = np.zeros((nrBootstraps,), float)
bootstraps_troc = np.zeros((nrBootstraps,), float)

bootstraps_twfm = np.zeros((nrBootstraps,), float)
bootstraps_twauc = np.zeros((nrBootstraps,), float)
bootstraps_twroc = np.zeros((nrBootstraps,), float)



bestThresPfm = thresholds[np.argmax(np.mean(prot_fm, axis=1))]
bestThresPwfm = thresholds[np.argmax(np.mean(wprot_fm, axis=1))]
bestThresTfm = thresholds[np.argmax(np.mean(term_fm, axis=1))]
bestThresTwfm = thresholds[np.argmax(wterm_fm)]
bestThresPsd = thresholds[np.argmin(np.mean(prot_sd, axis=1))]


for m in range(nrBootstraps):

	print 'Bootstraping ', m
	sys.stdout.flush()
	[newYtest, newYpost] = resample(Ytest, ypred, random_state=seedonia[m])

	tokeep = np.where(np.sum(newYtest, 0) > 0)[0]

	newYtest = newYtest[:, tokeep]
	newYpost = newYpost[:, tokeep]

	newYpred = (newYpost >= bestThresPfm).astype(int)


	[bootstraps_ppr[m], bootstraps_prc[m], bootstraps_pfm[m], _] = skm.precision_recall_fscore_support(newYtest, newYpred, average='samples')

	newYpred = (newYpost >= bestThresPwfm).astype(int)

	[bootstraps_pwpr[m], bootstraps_pwrc[m], bootstraps_pwfm[m], _] = skm.precision_recall_fscore_support(newYtest.T, newYpred.T, average='macro', sample_weight=ic[tokeep])

	newYpred = (newYpost >= bestThresPsd).astype(int)

	[_, _, sd] = gr.semanticDistance(newYtest, newYpred, ic[tokeep])

	bootstraps_psd[m] = np.mean(sd)

	newYpred = (newYpost >= bestThresTfm).astype(int)

	[bootstraps_tpr[m], bootstraps_trc[m], bootstraps_tfm[m], _] = skm.precision_recall_fscore_support(newYtest, newYpred, average='macro')

	newYpred = (newYpost >= bestThresTwfm).astype(int)

	[_, _, temp_twfm, _] = skm.precision_recall_fscore_support(newYtest, newYpred, average=None)



	ind = np.where(~np.isnan(temp_twfm))[0]


	bootstraps_twfm[m] = np.average(temp_twfm[ind], weights=ic[ind])

	tempauc = skm.average_precision_score(newYtest, newYpost, average=None)
	bootstraps_tauc[m] = np.nanmean(tempauc)


	ind = np.where(~np.isnan(tempauc))[0]
	bootstraps_twauc[m] = np.average(tempauc[ind], weights=ic[ind])


	temproc = skm.roc_auc_score(newYtest, newYpost, average=None)
	bootstraps_troc[m] = np.nanmean(temproc)

	ind = np.where(~np.isnan(temproc))[0]
	bootstraps_twroc[m] = np.average(temproc[ind], weights=ic[ind])



	bootstraps_pauc[m] = skm.average_precision_score(newYtest, newYpost, average='samples')




print 'saving...'

if lsdr_method == 'None':
	with open('cafaout_scores/gt/ic.pkl', 'wb') as f:
		pickle.dump(ic, f)


with open('cafaout_scores/' + directory_suffix + '/predictions.pkl', 'wb') as f:
	pickle.dump(ypred, f)

#print term_fm.shape

with open('cafaout_scores/' + directory_suffix + '/pfm.pkl', 'wb') as f:
	pickle.dump(prot_fm, f)

with open('cafaout_scores/' + directory_suffix + '/ppr.pkl', 'wb') as f:
	pickle.dump(prot_pr, f)

with open('cafaout_scores/' + directory_suffix + '/prc.pkl', 'wb') as f:
	pickle.dump(prot_rc, f)

with open('cafaout_scores/' + directory_suffix + '/pwfm.pkl', 'wb') as f:
	pickle.dump(wprot_fm, f)

with open('cafaout_scores/' + directory_suffix + '/pwpr.pkl', 'wb') as f:
	pickle.dump(wprot_pr, f)

with open('cafaout_scores/' + directory_suffix + '/pwrc.pkl', 'wb') as f:
	pickle.dump(wprot_rc, f)


with open('cafaout_scores/' + directory_suffix + '/tfm.pkl', 'wb') as f:
	pickle.dump(term_fm, f)

with open('cafaout_scores/' + directory_suffix + '/tpr.pkl', 'wb') as f:
	pickle.dump(term_pr, f)

with open('cafaout_scores/' + directory_suffix + '/trc.pkl', 'wb') as f:
	pickle.dump(term_rc, f)

with open('cafaout_scores/' + directory_suffix + '/psd.pkl', 'wb') as f:
	pickle.dump(prot_sd, f)

with open('cafaout_scores/' + directory_suffix + '/pauc.pkl', 'wb') as f:
	pickle.dump(prot_auc, f)

with open('cafaout_scores/' + directory_suffix + '/tauc.pkl', 'wb') as f:
	pickle.dump(term_auc, f)

with open('cafaout_scores/' + directory_suffix + '/troc.pkl', 'wb') as f:
	pickle.dump(term_roc, f)

with open('cafaout_scores/' + directory_suffix + '/twauc.pkl', 'wb') as f:
	pickle.dump(wterm_auc, f)

with open('cafaout_scores/' + directory_suffix + '/twroc.pkl', 'wb') as f:
	pickle.dump(wterm_roc, f)

with open('cafaout_scores/' + directory_suffix + '/twfm.pkl', 'wb') as f:
	pickle.dump(wterm_fm, f)



#bootstraps

with open('cafaout_scores/' + directory_suffix + '/btstrp_pfm.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_pfm, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_ppr.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_ppr, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_prc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_prc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_pwfm.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_pwfm, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_pwpr.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_pwpr, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_pwrc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_pwrc, [2.5, 97.5]), f)


with open('cafaout_scores/' + directory_suffix + '/btstrp_tfm.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_tfm, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_tpr.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_tpr, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_trc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_trc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_psd.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_psd, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_pauc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_pauc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_tauc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_tauc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_troc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_troc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_twauc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_twauc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_twroc.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_twroc, [2.5, 97.5]), f)

with open('cafaout_scores/' + directory_suffix + '/btstrp_twfm.pkl', 'wb') as f:
	pickle.dump(np.percentile(bootstraps_twfm, [2.5, 97.5]), f)
