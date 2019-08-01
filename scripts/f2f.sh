#! /bin/sh

# Developed by Aly Ammar at NCAR
# as a SIParCS intern summer 2019

DIRECTORY="f2foutput"
if [ ! -d $DIRECTORY ]; then
	mkdir $DIRECTORY
fi
LOG=$DIRECTORY/convert_log.txt
module load ftools
n=0
echo "Converting files..."
for f in *.f; do
	n=$((n+1))
	f=${f%.f}
	echo "    "
	echo $f.f
	perl /glade/work/alyammar/ff/f2f.pl $f.f $f.f90 
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


