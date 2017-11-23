import sys
import general_routines as gr
import numpy as np
import sklearn.model_selection as ms
from sklearn.metrics import average_precision_score
import lsdr

similarity = sys.argv[1]
ontology = sys.argv[2]

nr_out_folds = int(sys.argv[3])
nr_in_folds = int(sys.argv[4])

classifier = sys.argv[5]
lsdr_method = sys.argv[6]

nrPrototypes = int(sys.argv[7])
label_fraction = float(sys.argv[8])

outer_fold = int(sys.argv[9])
inner_fold = int(sys.argv[10])

k = int(sys.argv[11])


if lsdr_method != "None":
	var_scale = int(sys.argv[12])
	nr_dims = int(sys.argv[13])

	if lsdr_method == "CPLST":
		lamda = float(sys.argv[14])


organism = 'thalia'


myname = 'scores/' + gr.composeName(similarity, ontology, classifier, lsdr_method, nrPrototypes, label_fraction) + '/'


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


	inner_cv = ms.KFold(n_splits=nr_in_folds, random_state=int(inner_seed[i]), shuffle=True)

	j = -1

	for inner_train, inner_test in inner_cv.split(Xtrain):
		j += 1

		if j != inner_fold:
			continue

		print 'In loop..'
		Xtrain_in = Xtrain[inner_train,:][:,inner_train]
		Ytrain_in = Ytrain[inner_train,:]

		[Xtrain_in, Ytrain_in, tokeep, removedProteins] = gr.removeRareTerms(Xtrain_in, Ytrain_in, label_fraction)

		Xtest_in = Xtrain[inner_test,:][:,inner_train]
		
		Xtest_in = np.delete(Xtest_in, removedProteins, axis=1)
		
		Ytest_in = Ytrain[inner_test,:][:, tokeep]

		[Xtest_in, Ytest_in] = gr.removeEmptyTestProteins(Xtest_in, Ytest_in)


		if nrPrototypes > 0:
			np.random.seed(892374980)
			indices = np.random.permutation(Xtrain_in.shape[0])[:nrPrototypes]
			Xtrain_in = Xtrain_in[:, indices]
			Xtest_in = Xtest_in[:, indices]


		newTermNames = []
		
		for trmInd in tokeep:
			newTermNames.append(termNames[trmInd])


		parentsCoord = gr.getParentsCoord(newTermNames, ontology)



		if lsdr_method != 'None' and nr_dims > Ytrain_in.shape[1]:
			print 'too many dims'
			filename = myname + str(outer_fold) + '_' + str(inner_fold) + '_' + str(k)


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
			Ytrain_in_prime = Ytrain_in
			
		elif lsdr_method == 'PLST':
			mapping = lsdr.PLST()
			Ytrain_in_prime = mapping.fit(Ytrain_in, ndims = nr_dims, var_correction = var_scale)

		elif lsdr_method == 'CPLST':
			mapping = lsdr.CPLST()
			Ytrain_in_prime = mapping.fit(Xtrain_in, Ytrain_in, ndims = nr_dims, var_correction = var_scale, lamda=lamda)

		elif lsdr_method == 'GOAT':
			mapping = lsdr.GOAT(ontology)
			Ytrain_in_prime = mapping.fit(Ytrain_in, ndims = nr_dims, var_correction = var_scale, cols=tokeep) 

		elif lsdr_method == 'SEM':
			mapping = lsdr.SEM(ontology)
			Ytrain_in_prime = mapping.fit(Ytrain_in, newTermNames, ndims = nr_dims, var_correction = var_scale) 

		elif lsdr_method == 'DAG':
			mapping = lsdr.DAG(ontology)
			Ytrain_in_prime = mapping.fit(Ytrain_in, newTermNames, ndims = nr_dims, var_correction = var_scale) 			

		else:
			print 'Invalid LSDR method'
			exit()


		print 'classification'
		#classification
		Ytest_in_pred_prime = gr.predict(Xtrain_in, Ytrain_in_prime, Xtest_in, classifier, k)

	
		if lsdr_method == 'None':
			Ytest_in_pred = Ytest_in_pred_prime
			
		else: 
			Ytest_in_pred = mapping.inverseMap(Ytest_in_pred_prime)

	
		gr.respectGO(Ytest_in_pred, parentsCoord)
		print 'eval...'


		prot_auc = average_precision_score(Ytest_in, Ytest_in_pred, average='samples') 
		term_auc = average_precision_score(Ytest_in, Ytest_in_pred, average=None) 
		
		term_auc = np.nanmean(term_auc)
		

		filename = myname + str(outer_fold) + '_' + str(inner_fold) + '_' + str(k)


		if lsdr_method != "None":
			filename = filename + '_' + str(var_scale) + '_' + str(nr_dims)

			if lsdr_method == "CPLST":
				filename = filename + '_' + str(lamda)
		
		
		with open(filename + '.txt', 'w') as f:
			f.write(str(prot_auc) + '\n')
			f.write(str(term_auc) + '\n')			


