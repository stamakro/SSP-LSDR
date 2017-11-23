import sys
import general_routines as gr
import numpy as np
import sklearn.model_selection as ms
from sklearn.metrics import average_precision_score
import lsdr
import pickle
from scipy.io import loadmat


def getIndices(termsTrain, termsTest):

	ttr = set(termsTrain)
	tts = set(termsTest)

	trInd = []
	tsInd = []

	for t in ttr.intersection(tts):
		trInd.append(termsTrain.index(t))
		tsInd.append(termsTest.index(t))

	return [np.array(trInd), np.array(tsInd)]






similarity = 'sequence'
ontology = 'P'

nr_in_folds = 3

classifier = sys.argv[1]
lsdr_method = sys.argv[2]


inner_fold = int(sys.argv[3])

k = int(sys.argv[4])


if lsdr_method != "None":
	var_scale = int(sys.argv[5])
	nr_dims = int(sys.argv[6])

	if lsdr_method == "CPLST":
		lamda = float(sys.argv[7])


organism = 'thalia'

myname = 'cafascores/' + gr.composeName(similarity, ontology, classifier, lsdr_method, 0, 0.0) + '/'


print 'Loading data...'


X = loadmat('files/X.mat')['identityMatrix']

#training set only
Xtrain = X[:,:6077][:6077, :]



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


[trainIndices, testIndices] = getIndices(termsTrain, termsTest)


np.random.seed(17021991)
seed = 4294967295 * np.random.rand()


cv = ms.KFold(n_splits=nr_in_folds, random_state=int(seed), shuffle=True)

i = -1


for train, test in cv.split(Xtrain):
	i += 1

	if i != inner_fold:
		continue

	print 'In loop..'
	Xtrain_cv = Xtrain[train,:][:,train]

	Ytrain_cv = Ytrain[train,:]

	#[Xtrain_cv, Ytrain_cv, tokeep, removedProteins] = gr.removeRareTerms(Xtrain_cv, Ytrain_cv, 0.0)

	Xtest_cv = Xtrain[test,:][:, train]

	#Xtest_cv = np.delete(Xtest_cv, removedProteins, axis=1)

	Ytest_cv = Ytrain[test,:][:, trainIndices]


	[Xtest_cv, Ytest_cv] = gr.removeEmptyTestProteins(Xtest_cv, Ytest_cv)


	parentsCoord = gr.getParentsCoord(termsTrain, ontology)


	if lsdr_method != 'None' and nr_dims > Ytrain_cv.shape[1]:
		print 'too many dims'
		filename = myname + str(inner_fold) + '_' + str(k)


		if lsdr_method != "None":
			filename = filename + '_' + str(var_scale) + '_' + str(nr_dims)

			if lsdr_method == "CPLST":
				filename = filename + '_' + str(lamda)


		with open(filename + '.txt', 'w') as f:
			f.write(str(-10000) + '\n')
			f.write(str(-10000) + '\n')


		sys.exit(0)

	print 'dim red...'

	if lsdr_method == 'None':
		Ytrain_cv_prime = Ytrain_cv

	elif lsdr_method == 'PLST':
		mapping = lsdr.PLST()
		Ytrain_cv_prime = mapping.fit(Ytrain_cv, ndims = nr_dims, var_correction = var_scale)

	elif lsdr_method == 'CPLST':
		mapping = lsdr.CPLST()
		Ytrain_cv_prime = mapping.fit(Xtrain_cv, Ytrain_cv, ndims = nr_dims, var_correction = var_scale, lamda=lamda)

	elif lsdr_method == 'GOAT':
		mapping = lsdr.GOAT(ontology)
		Ytrain_cv_prime = mapping.fit(Ytrain_cv, ndims = nr_dims, var_correction = var_scale)

	elif lsdr_method == 'SEM':
		mapping = lsdr.SEM(ontology)
		Ytrain_cv_prime = mapping.fit(Ytrain_cv, termsTrain, ndims = nr_dims, var_correction = var_scale)

	elif lsdr_method == 'DAG':
		mapping = lsdr.DAG(ontology)
		Ytrain_cv_prime = mapping.fit(Ytrain_cv, termsTrain, ndims = nr_dims, var_correction = var_scale)

	else:
		print 'Invalid LSDR method'
		exit()


	print 'classification'
	#classification
	Ytest_cv_pred_prime = gr.predict(Xtrain_cv, Ytrain_cv_prime, Xtest_cv, classifier, k)


	if lsdr_method == 'None':
		Ytest_cv_pred = Ytest_cv_pred_prime

	else:
		Ytest_cv_pred = mapping.inverseMap(Ytest_cv_pred_prime)


	gr.respectGO(Ytest_cv_pred, parentsCoord)
	print 'eval...'

	Ytest_cv_pred = Ytest_cv_pred[:, trainIndices]


	prot_auc = average_precision_score(Ytest_cv, Ytest_cv_pred, average='samples')
	#term_auc = average_precision_score(Ytest_cv, Ytest_cv_pred, average=None)

	#term_auc = np.nanmean(term_auc)


	filename = myname + str(inner_fold) + '_' + str(k)


	if lsdr_method != "None":
		filename = filename + '_' + str(var_scale) + '_' + str(nr_dims)

		if lsdr_method == "CPLST":
			filename = filename + '_' + str(lamda)


	with open(filename + '.txt', 'w') as f:
		f.write(str(prot_auc) + '\n')
