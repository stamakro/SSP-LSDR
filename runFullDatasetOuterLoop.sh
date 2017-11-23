
outFolds=3
inFolds=3

str=$1
name=${str#*/}
mkdir "out_scores/"$name


for ((fold = 0; fold < $outFolds; fold++ ))
do
	for pt in {0..1}
	do
		python fulldataset_out.py $1 $fold $outFolds $inFolds $pt

	done
done
