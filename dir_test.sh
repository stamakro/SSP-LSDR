#check if a directory is present. If not, create it

if ! [ -d $1 ]
then
	echo "Creating directory"
	mkdir $1
else
	echo "Directory exists"
fi
	
