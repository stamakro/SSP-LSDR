import general_routines as gr
import subprocess
import os



similarity = 'sequence'

ontology = 'P'

classifiers = ['ssp', 'ms', 'transferblast', 'cafablast']

print 'Choose classifier:'

for i,s in enumerate(classifiers):
	print 'For ' + s + ': ' + str(i)
	
clfInd = int(raw_input())	

classifier = classifiers[clfInd]


lsdr_methods = ['None', 'PLST', 'CPLST', 'GOAT', 'SEM', 'DAG']


print
print
print


print 'Choose LSDR method:'

for i,s in enumerate(lsdr_methods):
	print 'For ' + s + ': ' + str(i)
	
lsdrInd = int(raw_input())	
lsdr = lsdr_methods[lsdrInd]


print
print
print


if classifier == 'sim':

	prototypes = int(raw_input('How many prototypes? (0 for no prototype selection):'))

	print
	print
	print

else:
	prototypes = 0


	
label_fraction = float(raw_input('Minimum percentage of positive samples:'))


dirName = gr.composeName(similarity, ontology, classifier, lsdr, prototypes, label_fraction)


command = './dir_test.sh scores/' + dirName
subprocess.call(command, shell=True)


command = './runFullDatasetInnerLoop.sh ' + similarity + ' ' + ontology + ' ' + classifier + ' ' + lsdr + ' ' + str(prototypes) + ' ' + str(label_fraction)
subprocess.call(command, shell=True)


