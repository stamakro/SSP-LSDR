import numpy as np
import matplotlib
matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import sklearn.metrics as skm
import sys
import pickle

#a lot of helper/parser routines, mostly related to the GO dag


class Gene:

        def __init__(self, name_):

                self.name = name_

                self.functions = set()

                self.predictedFunctions = set()


        def add_function(self, newF):

                if newF not in self.functions:
                        self.functions.add(newF)



        def getFunctions(self):
                return self.functions


        def getName(self):
                return self.name





class Term:
	def __init__(self):
		self.ID = ''
		self.manyIDs = False
		self.altIDs = []

		self.ontology = None

		self.isObsolete = False

		self.isReplaced = False
		self.replacedBy = None
		self.isLeaf = True
		self.parents = []
        	self.children = []
		self.associatedGenes = set()
		self.predictedGenes = set()
		self.IC = 0.0
		self.freq = 0.0


	def getCount(self):
		return len(self.associatedGenes)



def merge_annotation_datasets(dag1, dag2):
	new_dag = dag1

	if len(dag1) != len(dag2):
		print 'Different sizes'
		exit()

	for i in range(len(dag1)):
		if dag1[i].ID != dag2[i].ID:
			print 'Different order'
			exit()

		new_dag[i].associatedGenes = (dag1[i].associatedGenes).union(dag2[i].associatedGenes)
		
		
	return new_dag



def computePath2Root(dag, ontology):
	path2root = dict()

	if ontology == 'C':
		root = 5895
	elif ontology == 'F':
		root = 1114
	else:
		root = 6

	computePath2RootAux(dag, root, 0, path2root)

	return path2root



def computePath2RootAux(dag, ind, pos, path2root):
			
		if ind not in path2root:
			path2root[ind] = pos

		elif path2root[ind] < pos:
			path2root[ind] = pos

		for c in dag[ind].children:
			computePath2RootAux(dag, c, pos+1, path2root)




def calculate_IC(dag, ontology):
	new_dag = dag

	for i,t in enumerate(dag):

		if t.ontology != ontology:
			continue

		if len(t.parents) == 0 or t.getCount() == 0:
			t.IC = 0.0

		else:
			parentSet = dag[t.parents[0]].associatedGenes
			for pInd in t.parents:
				
				parentSet = parentSet.intersection(dag[pInd].associatedGenes)
				#print '%d %d' % (pInd, len(parentSet)) 

			if len(parentSet) > 0:
				prob = float(t.getCount()) / float(len(parentSet))		

				if prob > 1.0:

					#print 'current term'
					print t.ID

					prob = 1.0
					'''
					print prob

					print parentSet
					print t.associatedGenes

					print len(parentSet)
					print len(t.associatedGenes)


					sys.exit()
					'''

				t.IC = (-1.0) * np.log2(prob)
					
			else:
				#print 'Divide by 0'	
				#print t.ID
	
				t.IC = 0.0

		new_dag[i] = t


	return new_dag



def calculate_IC_all_ontologies(dag):
	new_dag = dag

	for i,t in enumerate(dag):

		#if t.ontology != ontology:
		#	continue

		if len(t.parents) == 0 or t.getCount() == 0:
			t.IC = 0.0

		else:
			parentSet = dag[t.parents[0]].associatedGenes
			for pInd in t.parents:
				
				parentSet = parentSet.intersection(dag[pInd].associatedGenes)
				#print '%d %d' % (pInd, len(parentSet)) 

			if len(parentSet) > 0:
				prob = float(t.getCount()) / float(len(parentSet))		

				if prob > 1.0:

					print 'current term'

					exit()


				t.IC = (-1.0) * np.log2(prob)
					
			else:
				#print 'Divide by 0'	
				#print t.ID
	
				t.IC = 0.0

		new_dag[i] = t


	return new_dag


def find_all_ancestors(dag, tInd):
    ancestors = set(dag[tInd].parents)

    for pInd in dag[tInd].parents:
        ancestors = ancestors.union(find_all_ancestors(dag, pInd))

    return ancestors


def find_ontology_ancestors(dag, tInd, ontology):
    ancestors = set([ind for ind in dag[tInd].parents if dag[ind].ontology == ontology])
    #ancestors = set(dag[tInd].parents)

    for pInd in dag[tInd].parents:
        ancestors = ancestors.union(find_ontology_ancestors(dag, pInd, ontology))

    return ancestors



def find_common_ancestors_denovo(dag, t1Ind, t2Ind, ontology, all):

    if not all:
        t1ancestors = find_ontology_ancestors(dag, t1Ind, ontology)
        t2ancestors = find_ontology_ancestors(dag, t2Ind, ontology)
    else:
        t1ancestors = find_all_ancestors(dag, t1Ind)
        t2ancestors = find_all_ancestors(dag, t2Ind)


    return t1ancestors.intersection(t2ancestors)


def get_all_frequencies(dag, ontology):

    individual_probability = dict()
    for i, t in enumerate(dag):
        if t.ontology == ontology:
          individual_probability[i] = float(t.getCount()) / float(nrPresentGenes)

    return individual_probability



def dyn_ancestors(dag, ontology, anc):
	
	leaves = [i for i in range(len(dag)) if dag[i].ontology == ontology and len(dag[i].children) == 0]
	
	ontoTerms = [i for i in range(len(dag)) if dag[i].ontology == ontology]

	isDone = np.zeros((len(dag),1), int)
	
	for leaf in leaves:
		recAncestors(dag, anc, leaf, isDone)
	
	return anc
	
	
def recAncestors(dag, anc, i, isDone):
	
	if isDone[i]:
		return anc[i]
	
	
	anc[i] = set(dag[i].parents)
		
		
		
	for pInd in dag[i].parents:		
			anc[i] = anc[i].union(recAncestors(dag, anc, pInd, isDone))
	
	isDone[i] = 1
	#print '!'

	return anc[i]


def getLCA(dag, ontology, ancestors, ti, tj):
	
	anci = ancestors[ti]
	
	ancj = ancestors[tj]
	
	if ti in ancj:
		ms = ti
		
	elif tj in anci:
		ms = tj
	
	else:
		common = anci.intersection(ancj)

		if len(common) == 1:
			ms = common.pop()
		
		else:
	
			tmpMs = -1.0
			ind = -1
		
			for t in common:
				if dag[t].IC > tmpMs:
					tmpMs = dag[t].IC
					ind = t
				
			ms = ind
		
			
	return ms	




def lin(dag, ontology):
	#dag should be deep-copied at call
	
	
	ontoTerms = [i for i in range(len(dag)) if dag[i].ontology == ontology]

	setTerms = set(ontoTerms)

	ancestors = [set() for i in range(len(dag))]
	
	'''
	if ontology == 'C':
		root = 5895
	elif ontology == 'F':
		root = 1114
	else:
		root = 6
	'''
		
	dyn_ancestors(dag, ontology, ancestors)	


	#return ancestors

	eps = 1e-5
	#test for 0

	lowestCommonAnc = np.zeros((len(ontoTerms), len(ontoTerms)), int)

	linMatrix = np.zeros((len(ontoTerms), len(ontoTerms)), float)

	lll = len(ontoTerms)


	for i in range(len(ontoTerms)):
	
		#if i % 100 == 0:
		#	print '%d out of %d' % (i, lll)

		lowestCommonAnc[i,i] = i


		pi = dag[ontoTerms[i]].IC
		
		if pi < eps:
			continue

		linMatrix[i,i] = 1.0		
		
		for j in range(i + 1, len(ontoTerms)):

			pj = dag[ontoTerms[j]].IC

			if pj < eps:
				continue


			lca = getLCA(dag, ontology, ancestors, ontoTerms[i], ontoTerms[j])
					
			lowestCommonAnc[i, j] = lca

			lsij = 2.0 * dag[lca].IC / (pi + pj)
			
			linMatrix[i,j] = lsij
			
			linMatrix[j,i] = lsij
						
			
	return [ontoTerms, linMatrix]



def annotate_gene_with_term(gene, dag, termInd):

	if dag[termInd].ID not in gene.getFunctions():
	
		gene.add_function(dag[termInd].ID)
	
		dag[termInd].associatedGenes.add(gene)
	
	for pInd in dag[termInd].parents:

		l =  annotate_gene_with_term(gene, dag, pInd)	
		dag = l[0]
		gene = l[1]


	return [dag, gene]



def read_go_per_gene(fil):

	ontology = read_ontology_from_file('files/terms.obo')

	tree = ontology[0]
	mapping = ontology[1]

	
	dupl = [i for i in range(len(tree)) if tree[i].manyIDs]


        allGenes = []

	currentGene = None

	with open(fil) as f:
        	for line in f:
                	if line[0]== '!':
                        	continue

                        fields = line.split('\t')


                        name = fields[1]
                        term = fields [4]

			stt = term[:3]
			if stt != 'GO:':
            			print 'fuck'

			if term not in mapping:
				found = False
				i = 0

				while i < len(dupl) and not found:
				
					if term in tree[dupl[i]].altIDs:
						term = tree[dupl[i]].ID
						found = True
					
					i += 1


				if not found:
					print 'still not found: ',
					print term


			if currentGene is None:
				
				currentGene = Gene(name)
				
				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]



                        elif name == currentGene.getName():
                                #print 'Continuation: ' + name + '\t' + term
				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]


                        else:
				allGenes.append(currentGene)

                                currentGene = Gene(name)

				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]




        return [allGenes, tree, mapping]


def read_go(organism, includeExclude, onto):

	ontology = read_ontology_from_file('files/terms.obo')

	fil = 'files/goa-allOrganisms-noIEA-reviewed-unique.gaf'

	if organism == 'human':
		organismCode = "9606"

	elif organism == 'thalia':
		organismCode = "3702"

	targetTaxon = "taxon:" + organismCode

	tree = ontology[0]
	mapping = ontology[1]

	
	dupl = [i for i in range(len(tree)) if tree[i].manyIDs]


        allGenes = []

	currentGene = None

	with open(fil) as f:
        	for line in f:
                	if line[0]== '!':
                        	continue

                        fields = line.split('\t')

                        if fields[6] == 'IEA' or fields[8] != onto:
                        	continue

                        currentTaxon = fields[12]

                        if includeExclude:
                        	#i want those that are from targetTaxon
                        	if currentTaxon != targetTaxon:
                        		continue

                        else:
                        	#i want everything but targetTaxon
                        	if currentTaxon == targetTaxon:
                        		continue

                        name = fields[1]
                        term = fields [4]

                        stt = term[:3]

                        if stt != 'GO:':
                        	print 'fuck'

                        if term not in mapping:
                        	found = False

                        	i = 0

                        	while i < len(dupl) and not found:
                        		if term in tree[dupl[i]].altIDs:
                        			term = tree[dupl[i]].ID
                        			found = True

                        		i += 1

                        	if not found:
                        		print 'still not found: ',
                        		print term


			if currentGene is None:
				
				currentGene = Gene(name)
				
				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]



                        elif name == currentGene.getName():
                                #print 'Continuation: ' + name + '\t' + term
				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]


                        else:
				allGenes.append(currentGene)

                                currentGene = Gene(name)

				
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]


        return [allGenes, tree, mapping]


def read_ontology_from_file(fil):

	originalIDs = set()
	obsolete = set()

	mapping = dict()

	useless = set(['\n', 'disjoint_from', 'comment', 'consider', 'created_by', 'creation_date', 'def', 'name', 'property_value', 'subset', 'synonym', 'xref', 'intersection_of'])

	usefullRelations = set(['part_of', 'regulates', 'positively_regulates', 'negatively_regulates', 'occurs_in'])

	dag = []

        with open(fil) as f:
                for line in f:
                        if line == '\n':
                                if dag[currentInd].isObsolete == True:
                                        obsolete.add(dag[currentInd].ID)
                                        del dag[currentInd]



                        if line.rstrip('\n') == '[Term]':
                                t = Term()

                                dag.append(t)

                                currentInd = len(dag) - 1

                        else:
                                fields = line.split(':')

                                pr = fields[0]


                                if pr in useless:
                                        continue


                                elif pr == 'id':

                                        ID = 'GO:' + fields[2].rstrip('\n')

                                        originalIDs.add(ID)

                                        if ID not in mapping:
                                                '''new term'''
                                                dag[currentInd].ID = ID
                                                mapping[t.ID] = currentInd

                                        else:
                                                '''seen before as someone's parent or replacement'''
                                                '''remove newly added object and get right position in the list'''
                                                del dag[currentInd]

                                                currentInd = mapping[ID]



                                elif pr == 'namespace':
                                        if fields[1].strip() == 'cellular_component':
                                                dag[currentInd].ontology = 'C'

                                        elif fields[1].strip() == 'molecular_function':
                                                dag[currentInd].ontology = 'F'


                                        elif fields[1].strip() == 'biological_process':
                                                dag[currentInd].ontology = 'P'

                                        else:
                                                print 'DANGER'




                                elif pr == 'alt_id':
                                        dag[currentInd].manyIDs = True
                                        newID = 'GO:' + fields[2].rstrip('\n')

                                        dag[currentInd].altIDs.append(newID)
					#print newID
					mapping[newID] = currentInd


                                elif pr == 'is_obsolete':
                                        if fields[1].strip() == 'true':
                                                dag[currentInd].isObsolete = True
                                        else:
                                                print fields[1]




                                elif pr == 'replaced_by':
                                        continue
                                        ''' if obsolete, just delete it. obsolete terms do not have parents   '''


                                elif pr == 'is_a':
                                        parentID = 'GO:' + fields[2].split('!')[0].strip()

                                        if parentID not in mapping:
                                                tt = Term()
                                                tt.ID = parentID
                                                dag.append(tt)

                                                mapping[parentID] = len(dag) - 1

                                        if mapping[parentID] not in dag[currentInd].parents:
                                                '''this is because some are and regulate'''
                                                dag[currentInd].parents.append( mapping[parentID])
                                                dag[mapping[parentID]].children.append(currentInd)

                                        dag[mapping[parentID]].isLeaf = False


                                elif pr == 'relationship':

                                        relType = fields[1].split()[0]

                                        if relType not in usefullRelations:
                                                continue

                                        parentID = 'GO:' + fields[2].split('!')[0].strip()

                                        if parentID not in mapping:
                                                tt = Term()
                                                tt.ID = parentID
                                                dag.append(tt)

                                                mapping[parentID] = len(dag) - 1

                                        if mapping[parentID] not in dag[currentInd].parents:
											'''this is because some are and regulate'''
											dag[currentInd].parents.append( mapping[parentID])
											dag[mapping[parentID]].children.append(currentInd)
												
                                        dag[mapping[parentID]].isLeaf = False



                                else:
                                        print pr

        return [dag, mapping]
                                                                                        








def getBlastData(genes, organism, number):

	hits = set()

	with open('files/' + organism + '-proteins-with-' + organism + '-hits.txt') as f:
		for line in f:
			hits.add(line.rsplit('\n')[0])


	tobedel = [i for i in range(len(genes)) if genes[i].name not in hits]



	for ind in sorted(tobedel, reverse=True):
        	del genes[ind]

	'''
	Learn how to program or quit!
	If you ask for all the genes, nothing is permuted
	'''
	np.random.seed(31578)
	tobedel = np.random.permutation(len(genes))[number:]


	for ind in sorted(tobedel, reverse=True):
		del genes[ind]

	#print len(genes)
	assert len(genes) == number
	

	return genes


#does not work for stuff before lsdr 
def getDataLabels(organism, includeExclude, ontology, nrObjects):

	[thaliaGenes, dag_thalia, map_thalia] = read_go(organism, True, ontology)	


	#print 'Keeping genes with blast hists...'
	thaliaGenes = getBlastData(thaliaGenes, organism, nrObjects)


	#print 'Removing unused terms and creating data labels...'
	[dataLabels, usedTerms, usedGenes] =  getLabelMatrix(thaliaGenes, dag_thalia, ontology)


	#used to have only 2 return arguments, should not work for older code
	return [dataLabels, usedTerms, usedGenes]


def getRNAData(genes, organism, number):
	hits = set()

	with open('../rna-seq/genes.pkl', 'rb') as f:
		gg = pickle.load(f)

	for g in gg:
		hits.add(g)


	tobedel = [i for i in range(len(genes)) if genes[i].name not in hits]



	for ind in sorted(tobedel, reverse=True):
        	del genes[ind]


	np.random.seed(31578)
	tobedel = np.random.permutation(len(genes))[number:]


	for ind in sorted(tobedel, reverse=True):
		del genes[ind]

	#print len(genes)
	assert len(genes) == number
	

	return genes



def getDataLabelsRNA(organism, includeExclude, ontology, nrObjects):

	[thaliaGenes, dag_thalia, map_thalia] = read_go(organism, True, ontology)	


	#print 'Keeping genes with blast hists...'
	thaliaGenes = getRNAData(thaliaGenes, organism, nrObjects)


	#print 'Removing unused terms and creating data labels...'
	[dataLabels, usedTerms, usedGenes] =  getLabelMatrix(thaliaGenes, dag_thalia, ontology)

	return [dataLabels, usedTerms, usedGenes]



def getLabelMatrix(genes, dag, ontology):

	#this function removes terms with no annotated proteins
	#CAREFUL! After calling it, the dag[x].parents, dag[x].children are no longer valid!

	[usedTermsIndices, usedTermsNames] = getUsedTerms(dag, ontology, set([g.name for g in genes]))
		
	usedGenesNames = [g.name for g in genes]	
		
	geneDict = dict()
	
	for i, t in enumerate(usedGenesNames):
		geneDict[t] = i
		
	data = np.zeros((len(genes), len(usedTermsIndices)), int)
	

	for j, tInd in enumerate(usedTermsIndices):
		assocGenes = [g.name for g in dag[tInd].associatedGenes]
	
		for gg in assocGenes:
			try:
	
				data[geneDict[gg], j] = 1

			except KeyError:
				continue
	
	return [data, usedTermsNames, usedGenesNames]



def getUsedTerms(dag, ontology, targetGenes):

	used = []
	names = []
	for i, t in enumerate(dag):
	
		if t.ontology == ontology:
		
			assocG = set([g.name for g in t.associatedGenes])
			
			if len(assocG.intersection(targetGenes)) > 0:		
				used.append(i)
				names.append(t.ID)
			
	return [used, names]
