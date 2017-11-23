import filereader as fr
import numpy as np
import pickle


def addAnnotation(protein, term, annotations):
	if protein in annotations:
		annotations[protein].add(term)
		
	else:
		annotations[protein] = set([term])

def getAllTerms(leaves, dag, mapping, ontology):

	allterms = set()
		
	while len(leaves) > 0:
		
		termInd = mapping[leaves.pop()]
		
		
		if dag[termInd].ontology == ontology:		
			allterms.add(dag[termInd].ID)
		
		
		for pInd in dag[termInd].parents:
			if dag[pInd].ontology == ontology:
				leaves.add(dag[pInd].ID)	
	
	return allterms
	
	
def createLabelMatrix(annotations, dag, mapping, ontology):	

	#big binary label matrix with all ontologies
	Y = np.zeros((len(annotations), len(dag)), int)

	proteins = [p for p in annotations]

	print len(proteins)

	for i, protein in enumerate(proteins):

		#if (i % 10) == 0:
		#	print i
						
		anno = annotations[protein]
		
		while len(anno) > 0:
		
			#print i, len(anno)
			#print anno
		
			termInd = mapping[anno.pop()]

			#print 'after'
			#print anno

			if dag[termInd].ontology == ontology:
				Y[i, termInd] = 1
				
			for pInd in dag[termInd].parents:
				anno.add(dag[pInd].ID)


	print Y.shape

	ontologyIndices = np.array([i for i in range(len(dag)) if dag[i].ontology == ontology])

	termNames = [dag[i].ID for i in ontologyIndices]

	Y = Y[:, ontologyIndices]
	
	print Y.shape	
	
	nonEmptyTerms = np.where(np.sum(Y, 0) > 0)[0]
	
	Y = Y[:, nonEmptyTerms]
	
	print Y.shape
	
	
	termNames = [termNames[i] for i in nonEmptyTerms]

	if 'GO:0008150' in termNames:
	
		root = termNames.index('GO:0008150')
	
		del termNames[root]
	
		Y = np.delete(Y, root, 1)

	
	return [Y, proteins, termNames]
	



trainProteins = set()

leaves = set()



with open('../training_data/thalia_proteins.list') as f:
	for line in f:
		trainProteins.add(line.rsplit('\n')[0])		


print 'Proteins in training set:', len(trainProteins)


[dag, mapping] = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms_09_17.obo')

anno = dict()

absentTerms = set()

with open('../training_data/uniprot_sprot_exp.txt') as f:
	for line in f:
	
		fields = line.split('\t')

		protein = fields[0]
	
		if protein not in trainProteins:
			continue
	
		term = fields[1].rsplit('\n')[0]
	
	
		if term not in mapping:
			absentTerms.add(term)
			
			continue
			
		if term not in leaves:
			leaves.add(term)
			
			
		if dag[mapping[term]].ontology == 'P':	
			#fw.write(line)
			addAnnotation(protein, term, anno)



print 'Different Go terms in file:', len(leaves)


print 

print 'Absent GO terms:', len(absentTerms)

terms = getAllTerms(leaves, dag, mapping, 'P')

proteins = [t for t in trainProteins]

[Y, pp, tt] = createLabelMatrix(anno, dag, mapping, 'P')	


permission2Write = False

with open('../training_data/uniprot_sprot_exp.fasta') as f:
	for i, line in enumerate(f):

		if line[0] == '>':
						
			if i > 0 and permission2Write:
				fw.close()


			permission2Write = False
				
			name = line[1:].rsplit('\n')[0]
			
			
			if name in pp:
				#print name
				fw = open('../training_data/train_fasta/' + name + '.fasta', 'w')
				permission2Write = True


		if permission2Write:
			fw.write(line)


with open('../training_data/labelMatrixTrain.pkl', 'wb') as f:
	pickle.dump(Y, f)


with open('../training_data/proteinsTrain.pkl', 'wb') as f:
	pickle.dump(pp, f)


with open('../training_data/termsTrain.pkl', 'wb') as f:
	pickle.dump(tt, f)


