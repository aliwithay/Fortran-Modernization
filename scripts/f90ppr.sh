#! /bin/sh
DIRECTORY="convertoutput"
if [ -d $DIRECTORY ]; then
	rm -r convertoutput
fi
mkdir $DIRECTORY
LOG=$DIRECTORY/convert_log.txt
module load ftools
n=0
echo "Converting files..."
for f in *.f; do
	n=$((n+1))
	f=${f%.f}
	echo "    "
	echo $f.f
	f90ppr <$f.f >$f.f90
done>$LOG
echo "Converted" $n "files."
mv *.f90 $DIRECTORY
CHECK=$DIRECTORY/syntax_log.txt
scount=0
pcount=0
echo "Analyzing converted files..."
for f in $DIRECTORY/*.f90; do
	#read something
	#clear
	scount=$((scount+1))
	#echo "   "
	#echo $f
	var=$(gfortran -std=f2008 -fsyntax-only $f 2>&1)
	if [ -z "$var" ]; then
		pcount=$((pcount+1))
	else
		echo "$var"
		echo "   "
	fi
done>$CHECK
ecount=$((scount-pcount))
echo "Analyzed" $scount "files."
echo $ecount "files failed."
echo "Check syntax_log.txt for details on the files."


