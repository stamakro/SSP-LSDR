import sys
import filereader as fr
import lsdr
from scipy.io import loadmat
from scipy.sparse import csr_matrix
import sklearn.model_selection as ms
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neighbors import NearestCentroid
from sklearn.metrics import average_precision_score
from sklearn.metrics.pairwise import cosine_distances
from sklearn.multioutput import MultiOutputClassifier
import numpy as np
import pickle
import copy

#routines for loading datasets, training methods and analyzing data



def composeName(similarity, ontology, classifier, lsdr, prototypes, label_fraction):

	return similarity + '_' + ontology + '_' + classifier + '_' + lsdr + '_' + str(prototypes) + '_' + str(label_fraction)



def decomposeName(name):

	fields = name.split('_')


	if fields[4] != 'allanc' and fields[4] != 'CAS' and fields[4] != 'LIN' and fields[4] != 'Resnik' and fields[4] != 'rel':
		fields[4] = int(fields[4])
		fields[5] = float(fields[5])

	else:
		fields[3] += '_' + fields[4]
		fields[4] = int(fields[5])
		fields[5] = float(fields[6])

		del fields[6]

	return fields



def respectGO_1term(Y, termColumn, parentsCoord):
	#helper function of respectGO

        yterm = Y[:,termColumn]

        for pcolumn in parentsCoord[termColumn]:
                Y[:, pcolumn] = np.maximum(Y[:, pcolumn], yterm)

                respectGO_1term(Y, pcolumn, parentsCoord)


def respectGO(Y, parentsCoord):
	#ensure that prediction respect DAG structure

        for i in range(Y.shape[1]):
                respectGO_1term(Y, i, parentsCoord)


def countInconsistentPredictions(Y, parentsCoord):
	countInconsistent = 0

	parents = set()
	for p in parentsCoord:
		parents = parents.union(set(parentsCoord[p]))

	countTotal = Y.shape[0] * len(parents)

	Y2 = copy.deepcopy(Y)

	respectGO(Y2, parentsCoord)

	countInconsistent = len(np.where(Y2 > Y)[0])



	return [countInconsistent, countTotal]


def countInconsistentPairs(Y, parentsCoord):
	countTotal = 0
	countInconsistent = 0

	nrProteins = Y.shape[0]

	for termColumn in xrange(Y.shape[1]):
	        yterm = Y[:,termColumn]


		countTotal += len(parentsCoord[termColumn]) * nrProteins

        	for pcolumn in parentsCoord[termColumn]:
			countInconsistent += len(np.where(yterm > Y[:, pcolumn])[0])



	return [countInconsistent, countTotal]


def predict(Xtrain, Ytrain, Xtest, classifier, k):

	assert Xtrain.shape[1] == Xtest.shape[1]

	if classifier == 'ssp':
		return simknnPredict(Xtrain, Ytrain, Xtest, k, 'euclidean')

	elif classifier == 'ms':
		return msknnPredict(Ytrain, Xtest, k)

	elif classifier == 'cafablast':
		return blastPredict(Ytrain, Xtest)

	elif classifier == 'transferblast':
		assert k == 1
		return transferPredict(Ytrain, Xtest, 1)
	else:
		print 'wrong classifier'
		sys.exit(1)



def simknnPredict(Xtrain, Ytrain, Xtest, k, dist, dd=dict()):

	clf = KNeighborsRegressor(n_neighbors=k, metric=dist, metric_params=dd)

	clf.fit(Xtrain, Ytrain)

	return clf.predict(Xtest)


def msknnPredict(Ytrain, Xtest, k):
	(nrTestProteins, nrTrainProteins) = Xtest.shape

	nrTerms = Ytrain.shape[1]

	Ypred = np.zeros((nrTestProteins, nrTerms), float)

	for i in range(nrTestProteins):
		x = Xtest[i,:]

		sortInd = np.argsort(x)[::-1][:k]

		Ypred[i, :] = x[sortInd].dot(Ytrain[sortInd[:k], :]) / k

 	return Ypred


def transferPredict(Ytrain, Xtest, k):
	(nrTestProteins, nrTrainProteins) = Xtest.shape

	nrTerms = Ytrain.shape[1]

	Ypred = np.zeros((nrTestProteins, nrTerms), float)

	for i in range(nrTestProteins):
		x = Xtest[i,:]

		sortInd = np.argmax(x)

		Ypred[i, :] = Ytrain[sortInd, :]

 	return Ypred




def blastPredict(Ytrain, Xtest):
	(nrTestProteins, nrTrainProteins) = Xtest.shape

	nrTerms = Ytrain.shape[1]



	posteriors = np.zeros((nrTestProteins, nrTerms), float)

	for i in range(nrTestProteins):
		posteriors[i,:] = np.amax(Ytrain.T * Xtest[i,:], axis=1)


 	return posteriors




def loadData(similarity, ontology, organism):

	assert ontology == 'P'
	root  = 'GO:0008150'


	[dag, mapping] = fr.read_ontology_from_file('files/terms.obo')



	[Y, termNames, geneNames] = fr.getDataLabels(organism, True, ontology, 8000)

	data = loadmat('files/identities.mat')

	X = data['identityMatrix']

	X = np.array(X)

	rootInd = termNames.index(root)


	Y = np.array(Y)
	assert X.shape[0] == Y.shape[0]

	Y = np.delete(Y, rootInd, axis=1)
	del termNames[rootInd]


	[X, Y, deletedProteins] = removeEmptyProteins(X, Y)

	parentsCoord = getParentsCoord(termNames, ontology, dag, mapping)

	for i in sorted(deletedProteins, reverse=True):

		del geneNames[i]


	return [X, Y, termNames, deletedProteins, parentsCoord]







def getParentsCoord(termNames, ontology, dag=None, mapping=None):

	if dag is None:
		[dag, mapping] = fr.read_ontology_from_file('files/terms.obo')
		#[dag, mapping] = fr.read_ontology_from_file('files/terms_09_17.obo')

	assert ontology == 'P'
	root  = 'GO:0008150'

	parentsCoord = dict()

	for i, tname in enumerate(termNames):
		dagInd = mapping[tname]

		parentsCoord[i] = []

		for pInd in dag[dagInd].parents:

		        if dag[pInd].ontology != ontology or dag[pInd].ID == root:
		                continue

		        parLoc = termNames.index(dag[pInd].ID)

		        parentsCoord[i].append(parLoc)

	return parentsCoord


def removeEmptyProteins(X, Y):

	nrTerms = np.sum(Y, 1)

	zeros = np.where(nrTerms == 0)

	X = np.delete(X, zeros, axis=0)
	X = np.delete(X, zeros, axis=1)
	Y = np.delete(Y, zeros, axis=0)

	return [X, Y, zeros[0]]


def removeEmptyTestProteins(X, Y):

	nrTerms = np.sum(Y, 1)

	zeros = np.where(nrTerms == 0)

	X = np.delete(X, zeros, axis=0)
	Y = np.delete(Y, zeros, axis=0)

	return [X, Y]


def removeRareTermsOnly(Y, G, label_fraction):
	label_fraction /= 100.0

	[nrProteins, nrTerms] = Y.shape

	fractions = np.sum(Y, 0) / float(nrProteins)

	tokeep = np.where(fractions > label_fraction)[0]

	Y = Y[:, tokeep]

	G = G[:, tokeep][tokeep, :]

	return [Y, G]



def removeRareTerms(X, Y, label_fraction):
	label_fraction /= 100.0

	[nrProteins, nrTerms] = Y.shape

	fractions = np.sum(Y, 0) / float(nrProteins)

	tokeep = np.where(fractions > label_fraction)[0]

	Y = Y[:, tokeep]

	[X, Y, removedProteins] = removeEmptyProteins(X, Y)

	return [X, Y, tokeep, removedProteins]



def calculateIC(Ytrain, parentsCoord):

        (nrProteins, nrTerms) = Ytrain.shape

        ic = np.zeros((nrTerms,), float)

        for i in range(nrTerms):
                numerator = Ytrain[:,i]

                denominator = np.ones(nrProteins,)

                for p in parentsCoord[i]:
                        denominator = np.minimum(denominator, Ytrain[:,p])

                if np.sum(denominator) > 0 and np.sum(numerator) > 0:

                        p = float(np.sum(numerator)) / float(np.sum(denominator))

                        if p < 0.0 or p > 1.0:
                                print 'wtf'

                        ic[i] = -1.0 * np.log2(p)

                else:
                        ic[i] = 0.0

        return ic



def remainingUncertainty(Ytrue, Ypred, termIC, avg=False):
	(nrProteins, nrTerms) = Ytrue.shape

	ru = np.zeros((nrProteins,), float)

	for i in range(nrProteins):
		pr = Ypred[i,:]
		gt = Ytrue[i,:]

		fn =  np.intersect1d(np.where(gt == 1), np.where(pr == 0))

		for fni in fn:
			ru[i] += termIC[fni]

	if avg:
		ru = np.mean(ru)

	return ru


def misInformation(Ytrue, Ypred, termIC, avg=False):
	(nrProteins, nrTerms) = Ytrue.shape

	mi = np.zeros((nrProteins,), float)

	for i in range(nrProteins):
		pr = Ypred[i,:]
		gt = Ytrue[i,:]

		fp =  np.intersect1d(np.where(gt == 0), np.where(pr == 1))

		for fpi in fp:
			mi[i] += termIC[fpi]


	if avg:
		mi = np.mean(mi)

	return mi


def semanticDistance(Ytrue, Ypred, termIC):
	ru = remainingUncertainty(Ytrue, Ypred, termIC, False)
	mi = misInformation(Ytrue, Ypred, termIC, False)

	sd = np.sqrt(ru ** 2 + mi ** 2)

	return [ru, mi, sd]


def transform(X):

	m = float(np.amin(X))
	M = float(np.amax(X))

	return (X - m) / (M - m)

def calculateGOseq(ontology):
	organism = 'thalia'

	[dag, mapping] = fr.read_ontology_from_file('files/terms.obo')


	[Y, termNames, _] = fr.getDataLabels(organism, True, ontology, 8000)

	go = np.eye(len(termNames))


	for i, tname in enumerate(termNames):
		dagInd = mapping[tname]

		for pInd in dag[dagInd].parents:

			if dag[pInd].ontology != ontology:
				continue

			parLoc = termNames.index(dag[pInd].ID)

			go[i, parLoc] = 1.0
			go[parLoc, i] = 1.0


	#[Y, go] = removeRareTermsOnly(Y, go, label_fraction)

	go = csr_matrix(go)

	with open('go' + ontology + '.pkl', 'wb') as f:
		pickle.dump(go, f)

 	return go


def calculateGOallseq(ontology):
	organism = 'thalia'


	[dag, mapping] = fr.read_ontology_from_file('files/terms.obo')


	[Y, termNames, _] = fr.getDataLabels(organism, True, ontology, 8000)

	go = DFS(dag, termNames, mapping, ontology)

	go = csr_matrix(go)

	with open('go_allancestors' + ontology + '.pkl', 'wb') as f:
		pickle.dump(go, f)

 	return go


def DFS(dag, termNames, mapping, ontology):

	if ontology == 'C':
		root = 5895
	elif ontology == 'F':
		root = 1114
	else:
		root = 6

	gomat = np.eye(len(termNames))

	[gomat2, uselesspath] = visitNode(gomat, dag, termNames, mapping, [], termNames.index(dag[root].ID))

	return gomat2


def visitNode(gomat, dag, termNames, mapping, path, node):

	#print 'In node: ', node
	path.append(node)

	#print 'Path: ', path

	for cc in path:
		#print 'Updating entry ' + str(node) + ', ' + str(cc)
		gomat[node, cc] = 1
		gomat[cc, node] = 1

	for ch in dag[mapping[termNames[node]]].children:

		if dag[ch].ID not in termNames:
			continue



		[gomat, path] = visitNode(gomat, dag, termNames, mapping, path, termNames.index(dag[ch].ID))

	del path[-1]

	return [gomat, path]

def getLCAs(termNames, ontology, dag=None, mapping=None):

	if dag is None:
		[dag, mapping] = fr.read_ontology_from_file('files/terms.obo')

	assert ontology == 'P'
	#find all ancestors of each term
	ancestors = dict()

	for i, t in enumerate(termNames):
		ancestors[t] = set()

		for pInd in dag[mapping[t]].parents:
			getAncestorsRecursively(dag, pInd, ancestors[t], ontology)

		ancestors[t].remove(root)


	lcas = dict()
	#get common ancestors
	for i, t1 in enumerate(termNames):
		print i
		for j, t2 in enumerate(termNames):

			#print t1, t2

			if i == j:
				lcas[t1, t2] = t1
			else:
				cas = ancestors[t1].intersection(ancestors[t2])
				#print cas

				for t in cas:

					others = cas.difference(set([t]))


					if others == ancestors[t]:
						lcas[t1, t2] = t

					break

	with open('lowestCommonAncestors.pkl', 'wb') as f:
		pickle.dump(lcas, f)

	with open('ancestors.pkl', 'wb') as f:
		pickle.dump(ancestors, f)


def getAncestorsRecursively(dag, current, ancestors, ontology):
	if dag[current].ontology == ontology:
		ancestors.add(dag[current].ID)

	for pInd in dag[current].parents:
		getAncestorsRecursively(dag, pInd, ancestors, ontology)


def getLCAsCAFA(termNames, ontology, dag=None, mapping=None):

	if dag is None:
		[dag, mapping] = fr.read_ontology_from_file('files/terms_09_17.obo')

	if ontology == 'P':
		root  = 'GO:0008150'
	elif ontology == 'F':
		root = 'GO:0003674'
	else:
		print 'wrong ontology'
		sys.exit(1)


	#find all ancestors of each term

	ancestors = dict()

	for i, t in enumerate(termNames):
		ancestors[t] = set()

		for pInd in dag[mapping[t]].parents:
			getAncestorsRecursively(dag, pInd, ancestors[t], ontology)

		ancestors[t].remove(root)


	lcas = dict()
	#get common ancestors
	for i, t1 in enumerate(termNames):
		print i
		for j, t2 in enumerate(termNames):

			#print t1, t2

			if i == j:
				lcas[t1, t2] = t1
			else:
				cas = ancestors[t1].intersection(ancestors[t2])
				#print cas

				for t in cas:

					others = cas.difference(set([t]))


					if others == ancestors[t]:
						lcas[t1, t2] = t

					break

	with open('lowestCommonAncestors.pkl', 'wb') as f:
		pickle.dump(lcas, f)

	with open('ancestors.pkl', 'wb') as f:
		pickle.dump(ancestors, f)


def calculateGOseqCAFA(ontology, termNames):
	organism = 'thalia'

	[dag, mapping] = fr.read_ontology_from_file('files/terms_09_17.obo')

	if ontology == 'P':
		root  = 'GO:0008150'
	elif ontology == 'F':
		root = 'GO:0003674'
	else:
		print 'wrong ontology'
		sys.exit(1)



	go = np.eye(len(termNames))


	for i, tname in enumerate(termNames):
		dagInd = mapping[tname]

		for pInd in dag[dagInd].parents:

			if dag[pInd].ontology != ontology or dag[pInd] == root:
				continue


			parLoc = termNames.index(dag[pInd].ID)

			go[i, parLoc] = 1.0
			go[parLoc, i] = 1.0


	#[Y, go] = removeRareTermsOnly(Y, go, label_fraction)

	go = csr_matrix(go)

	with open('go' + ontology + '.pkl', 'wb') as f:
		pickle.dump(go, f)

 	return go
