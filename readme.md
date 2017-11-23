# SSP-LSDR
Code for experiments

## Full Dataset
* Calculate pairwise similarities with *alignments/fulldataset/allVsAllAlign.m*
* For the inner loop, use *runFullDatasetInnerLoop.py*
* For the outer loop, use *runFullDatasetOuterLoop.sh scores/directory*


## CAFA Dataset
* Download the CAFA3 data from *http://biofunctionprediction.org/cafa/*
* Parse the data with *createCafaTrainingSet.py* and *createCafaTestSet.py*
* Calculate pairwise similarities with *alignments/fulldataset/cafaAlign.m*
* For training a method, use *runCafaTrain.py*
* For evaluating, use *runCafaTest.sh*


## Dependencies
python2, numpy, scipy, scikit-learn, matplotlib, pickle, matlab Bioinformatics toolbox
