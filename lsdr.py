import numpy as np
import scipy.linalg as la
from sklearn.decomposition import FastICA
from sklearn.decomposition import KernelPCA
import pickle



def regularizedPinv(X, lamda):

	return la.inv((X.T).dot(X) + lamda * np.eye(X.shape[1])).dot(X.T)


'''-----------------------------------------------------------------'''
class PLST:

	#Principle label space transformation, Tai & Lin 2010
	def __init__(self):
		#transformation matrix
		self.w = np.array([])
		#inverse map
		self.winv = np.array([])

		self.bias = np.array([])
		self.sigma = np.array([])
		self.varianceExplained = 0.0


	def fit(self, Y, **kwargs):
		#Y: #objects x #labels

		try:
			var_correction = kwargs.pop('var_correction')

		except KeyError:
			var_correction = False


		try:
			ndims =  kwargs.pop('ndims')

			dim = True

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)

		except KeyError:

			var =  kwargs.pop('var')

			dim = False

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)


		if dim and ndims > Y.shape[1]:
			raise ValueError('Invalid number of dimensions specified')
		if not dim and (var < 0.0 or var > 1.0):
			raise ValueError('Invalid variance fraction specified')


		self.bias = np.mean(Y, axis=0)

		Y = Y - np.tile(self.bias, [Y.shape[0], 1])

		if var_correction:
			self.sigma = np.std(Y, axis=0)

			for i in range(self.sigma.shape[0]):
				if self.sigma[i] == 0:
					self.sigma[i] = 1


			Y = Y / np.tile(self.sigma, (Y.shape[0],1))

		else:
			self.sigma = np.ones((Y.shape[1],))


		#print np.sum(np.isnan(Y).astype(int))

		[U, s, VT] = la.svd(Y)


		totalVariance = np.sum(s)

		if dim:
			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		else:

			variances = np.cumsum(s) / totalVariance

			found = False
			i = 0

			ndims = -1

			while not found and i != variances.shape[0]:
				if variances[i] > var:
					found = True
					ndims = i + 1
				i += 1

			if not found:
				ndims = Y.shape[1]

			print ndims

			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		return Y.dot(self.w)


	def inverseMap(self, Yprime):

		Yrec = Yprime.dot(self.winv) * np.tile(self.sigma, [Yprime.shape[0], 1])

		Yrec = Yrec + np.tile(self.bias, [Yprime.shape[0], 1])

		return Yrec

'''-----------------------------------------------------------------'''


class CPLST:

	#Conditional principle label space transformation, Chen & Lin, 2012
	def __init__(self):
		#transformation matrix
		self.w = np.array([])
		#inverse map
		self.winv = np.array([])

		self.bias = np.array([])
		self.sigma = np.array([])
		self.varianceExplained = 0.0


	def fit(self, X, Y, **kwargs):
		#Y: #objects x #labels

		try:
			var_correction = kwargs.pop('var_correction')

		except KeyError:
			var_correction = False


		try:
			lamda = kwargs.pop('lamda')
		except KeyError:
			lamda = 0.0

		try:
			ndims =  kwargs.pop('ndims')

			dim = True

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)

		except KeyError:

			var =  kwargs.pop('var')

			dim = False

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)


		if dim and ndims > Y.shape[1]:
			raise ValueError('Invalid number of dimensions specified')
		if not dim and (var < 0.0 or var > 1.0):
			raise ValueError('Invalid variance fraction specified')


		self.bias = np.mean(Y, axis=0)

		Y = Y - np.tile(self.bias, [Y.shape[0], 1])

		if var_correction:
			self.sigma = np.std(Y, axis=0)

			for i in range(self.sigma.shape[0]):
				if self.sigma[i] == 0:
					self.sigma[i] = 1

			Y = Y / np.tile(self.sigma, (Y.shape[0],1))

		else:
			self.sigma = np.ones((Y.shape[1],))


		ZtHZ = (((Y.T).dot(X)).dot(regularizedPinv(X, lamda))).dot(Y)


		[U, s, VT] = la.svd(ZtHZ)


		totalVariance = np.sum(s)

		if dim:
			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		else:

			variances = np.cumsum(s) / totalVariance

			found = False
			i = 0

			ndims = -1

			while not found and i != variances.shape[0]:
				if variances[i] > var:
					found = True
					ndims = i + 1
				i += 1

			if not found:
				ndims = Y.shape[1]

			print ndims

			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		return Y.dot(self.w)


	def inverseMap(self, Yprime):

		Yrec = Yprime.dot(self.winv) * np.tile(self.sigma, [Yprime.shape[0], 1])

		Yrec = Yrec + np.tile(self.bias, [Yprime.shape[0], 1])

		return Yrec

'''-----------------------------------------------------------------'''


class GOAT:

	#Label Space Dimensionality Reduction with DAG-aware transformation

	def __init__(self, ontology):
		#transformation matrix
		self.w = np.array([])
		#inverse map
		self.winv = np.array([])

		self.bias = np.array([])
		self.sigma = np.array([])
		self.varianceExplained = 0.0

		with open('files/go' + ontology + '.pkl', 'rb') as fgo:
			#file generated by general_routines.calculateGOseq() or general_routines.calculateGOseqCAFA()
			self.gomat = pickle.load(fgo)

		self.gomat = self.gomat.toarray()
		print self.gomat.shape

	def fit(self, Y, **kwargs):
		#Y: #objects x #labels

		try:
			columns = kwargs.pop('cols')

		except KeyError:
			columns = None

		try:
			var_correction = kwargs.pop('var_correction')

		except KeyError:
			var_correction = False


		try:
			ndims =  kwargs.pop('ndims')

			dim = True

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)

		except KeyError:

			var =  kwargs.pop('var')

			dim = False

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)


		if dim and ndims > Y.shape[1]:
			raise ValueError('Invalid number of dimensions specified')
		if not dim and (var < 0.0 or var > 1.0):
			raise ValueError('Invalid variance fraction specified')


		self.bias = np.mean(Y, axis=0)

		Y = Y - np.tile(self.bias, [Y.shape[0], 1])

		if var_correction:
			self.sigma = np.std(Y, axis=0)

			for i in range(self.sigma.shape[0]):
				if self.sigma[i] == 0:
					self.sigma[i] = 1


			Y = Y / np.tile(self.sigma, (Y.shape[0],1))

		else:
			self.sigma = np.ones((Y.shape[1],))



		if columns is not None:
			self.gomat = self.gomat[:, columns][columns,:]

		[U, s, VT] = la.svd(Y.dot(self.gomat))


		totalVariance = np.sum(s)

		if dim:
			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		else:

			variances = np.cumsum(s) / totalVariance

			found = False
			i = 0

			ndims = -1

			while not found and i != variances.shape[0]:
				if variances[i] > var:
					found = True
					ndims = i + 1
				i += 1

			if not found:
				ndims = Y.shape[1]

			print ndims

			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		return Y.dot(self.w)


	def inverseMap(self, Yprime):

		Yrec = Yprime.dot(self.winv) * np.tile(self.sigma, [Yprime.shape[0], 1])

		Yrec = Yrec + np.tile(self.bias, [Yprime.shape[0], 1])

		return Yrec




'''-----------------------------------------------------------------'''


class SEM:

	#Label Space Dimensionality Reduction using Semantic Similarity between GO terms

	def __init__(self, ontology):
		#transformation matrix
		self.w = np.array([])
		#inverse map
		self.winv = np.array([])

		self.bias = np.array([])
		self.sigma = np.array([])
		self.varianceExplained = 0.0



	def fit(self, Y, termNames, **kwargs):
		#Y: #objects x #labels

		try:
			var_correction = kwargs.pop('var_correction')

		except KeyError:
			var_correction = False


		try:
			ndims =  kwargs.pop('ndims')

			dim = True

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)

		except KeyError:

			var =  kwargs.pop('var')

			dim = False

			if kwargs:
				raise TypeError('Unexpected **kwargs: %r' % kwargs)


		if dim and ndims > Y.shape[1]:
			raise ValueError('Invalid number of dimensions specified')
		if not dim and (var < 0.0 or var > 1.0):
			raise ValueError('Invalid variance fraction specified')

		#calculate resnik similarity

		with open('files/lowestCommonAncestors.pkl') as f:
			#file generated by general_routines.getLCAs
			lcas = pickle.load(f)


		[nrProteins, nrTerms] = Y.shape

		nrProteins = float(nrProteins)

		#this is minus the IC!!!!!!!!!!
		termICs = np.log2(np.sum(Y, axis=0) / nrProteins)

		termICs[np.where(np.isinf(termICs))[0]] = 0.0


		resnik = np.zeros((nrTerms, nrTerms), float)


		for ii in range(nrTerms):
			#just for speed
			yyy = Y[:,ii]
			t1 = termNames[ii]


			for jj in range(ii, nrTerms):
				t2 = termNames[jj]

				try:
					c = lcas[t1, t2]
				except KeyError:
					try:
						c = lcas[t2, t1]
					except KeyError:
						c = -1

				if c == -1:
					resnik[ii, jj] = 0.0
					resnik[jj, ii] = 0.0

				else:
					resnik[ii, jj] = -1.0 * termICs[termNames.index(c)]
					resnik[jj, ii] = resnik[ii, jj]
		#resnik calculated


		self.bias = np.mean(Y, axis=0)

		Y = Y - np.tile(self.bias, [Y.shape[0], 1])

		if var_correction:
			self.sigma = np.std(Y, axis=0)

			for i in range(self.sigma.shape[0]):
				if self.sigma[i] == 0:
					self.sigma[i] = 1


			Y = Y / np.tile(self.sigma, (Y.shape[0],1))

		else:
			self.sigma = np.ones((Y.shape[1],))


		[U, s, VT] = la.svd(Y.dot(resnik))

		totalVariance = np.sum(s)

		if dim:
			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		else:

			variances = np.cumsum(s) / totalVariance

			found = False
			i = 0

			ndims = -1

			while not found and i != variances.shape[0]:
				if variances[i] > var:
					found = True
					ndims = i + 1
				i += 1

			if not found:
				ndims = Y.shape[1]

			print ndims

			self.w = VT[:ndims,:].T

			self.winv = self.w.T

			self.varianceExplained = np.sum(s[:ndims]) / totalVariance

		return Y.dot(self.w)


	def inverseMap(self, Yprime):

		Yrec = Yprime.dot(self.winv) * np.tile(self.sigma, [Yprime.shape[0], 1])

		Yrec = Yrec + np.tile(self.bias, [Yprime.shape[0], 1])

		return Yrec


'''-----------------------------------------------------------------'''

class DAG:
	#Bi & Kwok


	def __init__(self, ontology):
		#transformation matrix
		self.w = np.array([])
		#inverse map
		self.winv = np.array([])

		self.bias = np.array([])
		self.sigma = np.array([])
		self.varianceExplained = 0.0



	def fit(self, Y, termNames, **kwargs):
		#Y: #objects x #labels

		ndims =  kwargs.pop('ndims')

		try:
			var_correction = kwargs.pop('var_correction')

		except KeyError:
			var_correction = False


		if kwargs:
			raise TypeError('Unexpected **kwargs: %r' % kwargs)


		if ndims > Y.shape[1]:
			raise ValueError('Invalid number of dimensions specified')

		#load ancestor vectors


		self.bias = np.mean(Y, axis=0)

		Y = Y - np.tile(self.bias, [Y.shape[0], 1])

		if var_correction:
			self.sigma = np.std(Y, axis=0)

			for i in range(self.sigma.shape[0]):
				if self.sigma[i] == 0:
					self.sigma[i] = 1


			Y = Y / np.tile(self.sigma, (Y.shape[0],1))

		else:
			self.sigma = np.ones((Y.shape[1],))

		with open('files/ancestors.pkl') as f:
			#file generated by general_routines.getLCAs()
			ancs = pickle.load(f)


		#calculate DAG ancestor vectors
		K = np.zeros((Y.shape[1], Y.shape[1]), float)

		for i in range(Y.shape[1]):

			currentAncs = ancs[termNames[i]]

			for tt in currentAncs:
				K[i, termNames.index(tt)] = 1


		[U, s, VT] = la.svd(K)

		self.w = VT[:ndims,:].T

		self.winv = self.w.T

		totalVariance = np.sum(s)

		self.varianceExplained = np.sum(s[:ndims]) / totalVariance


		return Y.dot(self.w)




	def inverseMap(self, Yprime):

		Yrec = Yprime.dot(self.winv) * np.tile(self.sigma, [Yprime.shape[0], 1])

		Yrec = Yrec + np.tile(self.bias, [Yprime.shape[0], 1])

		return Yrec
