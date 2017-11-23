#!/bin/bash
declare -a ssp_ks=(1 3 5 7 11 15 17 19 21 25 30 35 40 45 50 55 60 80 100 150 200 300 400 500 1000)
declare -a msknn_ks=(1 3 5 9 13 17 25 30)

declare -a var_scales=(0 1)
declare -a cplst_lamdas=(0.0001 0.1 0.2 0.5 1.0 2.0 5.0 10.0)
declare -a lsdr_dims=(50 100 200 500 1000 2000)


outFolds=3
inFolds=3

similarity=$1
ontology=$2
classifier=$3
lsdr=$4
prototypes=$5
labelFraction=$6


if [ "$classifier" == "ssp" ]
then
	ks=${ssp_ks[@]}
	
elif [ "$classifier" == "ms" ]
then
	ks=${msknn_ks[@]}
	
elif [ "$classifier" == "transferblast" ]
then
	ks=(1)		
	
elif [ "$classifier" == "cafablast" ]
then
	ks=(1)		

elif [ "$classifier" == "nmc" ]
then
	ks=(1)		

else
	echo "Invalid classifier"
	exit
fi

if [ "$similarity" != "sequence" ] 
then
	echo "Invalid similarity measure"
	exit
fi



for ((i = 0; i < $outFolds; i++ ))
do
	for ((j = 0; j < $inFolds; j++ ))
	do
		if [ "$lsdr" == "None" ]
		then
			for k in $ks
			do
			
				if  ! [ -f "scores/"$similarity"_"$ontology"_"$classifier"_"$lsdr"_"$prototypes"_"$labelFraction"/"$i"_"$j"_"$k".txt" ] 
				then
					python fulldataset.py $similarity $ontology $outFolds $inFolds $classifier $lsdr $prototypes $labelFraction $i $j $k
			
				fi
				
	
			done

		elif [ "$lsdr" == "CPLST" ]
		then
			for k in $ks
			do
				for var in ${var_scales[@]}
				do
					for dim in ${lsdr_dims[@]}
					do
						for l in ${cplst_lamdas[@]}
						do


							if  ! [ -f "scores/"$similarity"_"$ontology"_"$classifier"_"$lsdr"_"$prototypes"_"$labelFraction"/"$i"_"$j"_"$k"_"$var"_"$dim"_"$l".txt" ] 
							then
								python fulldataset.py $similarity $ontology $outFolds $inFolds $classifier $lsdr $prototypes $labelFraction $i $j $k $var $dim $l
							fi
						done							
					done
				done
			done

		else
			for k in $ks
			do
				for var in ${var_scales[@]}
				do
					for dim in ${lsdr_dims[@]}
					do					

						if  ! [ -f "scores/"$similarity"_"$ontology"_"$classifier"_"$lsdr"_"$prototypes"_"$labelFraction"/"$i"_"$j"_"$k"_"$var"_"$dim".txt" ] 
						then

							python fulldataset.py $similarity $ontology $outFolds $inFolds $classifier $lsdr $prototypes $labelFraction $i $j $k $var $dim
						fi
					done
				done
			done
		fi
	done
done


