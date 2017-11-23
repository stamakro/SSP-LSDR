
str=$1
name=${str#*/}
mkdir "cafaout_scores/"$name


python cafaTest.py $1 3
