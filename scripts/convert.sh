#! /bin/sh

# Developed by Aly Ammar at NCAR
# as a SIParCS intern, summer 2019

DIRECTORY="convertoutput"
if [ -d $DIRECTORY ]; then
	rm -r convertoutput
fi
mkdir $DIRECTORY
LOG=$DIRECTORY/convert_log.txt
module load ftools
for f in *.f; do
	f=${f%.f}
	echo "    "
	echo $f.f
	/glade/u/apps/ch/opt/ftools/bin/convert <<<$f,3,8,T,F/ 
done > $LOG
mv *.f90 $DIRECTORY
for f in $DIRECTORY/*.f90; do
	read something
	clear
	echo $f
	gfortran -std=f2008 -fsyntax-only $f
	echo done
done
