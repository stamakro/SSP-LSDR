import general_routines as gr
import subprocess
import os



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


dirName = gr.composeName('sequence', 'P', classifier, lsdr, 0, 0.0)


command = './dir_test.sh cafascores/' + dirName
subprocess.call(command, shell=True)


command = './runCafaTrain.sh ' + classifier + ' ' + lsdr  
subprocess.call(command, shell=True)


