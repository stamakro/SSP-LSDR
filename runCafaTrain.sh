#!/bin/bash
declare -a simknn_ks=(1 3 5 7 11 15 17 19 21 25 30)
declare -a msknn_ks=(1 3 5 9 13 17 25 30)

declare -a var_scales=(0 1)
declare -a cplst_lamdas=(0.0001 0.1 0.2 0.5 1.0 2.0 5.0 10.0)
declare -a lsdr_dims=(50 100 200 500 1000 2000)



inFolds=3


classifier=$1
lsdr=$2


if [ "$classifier" == "sim" ]
then
	ks=${simknn_ks[@]}
	
elif [ "$classifier" == "ms" ]
then
	ks=${msknn_ks[@]}
	
elif [ "$classifier" == "blast" ]
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




for ((j = 0; j < $inFolds; j++ ))
do
	if [ "$lsdr" == "None" ]
	then
		for k in $ks
		do
			
			if  ! [ -f "cafascores/sequence_P_"$classifier"_"$lsdr"_0_0.0/"$i"_"$j"_"$k".txt" ] 
			then
				python cafaTrain.py $classifier $lsdr $j $k
			
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


						if  ! [ -f "cafascores/sequence_P_"$classifier"_"$lsdr"_0_0.0/"$i"_"$j"_"$k"_"$var"_"$dim"_"$l".txt" ] 
						then
							python cafaTrain.py $classifier $lsdr $j $k $var $dim $l
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

					if  ! [ -f "cafascores/sequence_P_"$classifier"_"$lsdr"_0_0.0/"$i"_"$j"_"$k"_"$var"_"$dim".txt" ] 
					then

						python cafaTrain.py $classifier $lsdr $j $k $var $dim
					fi
				done
			done
		done
	fi
done



