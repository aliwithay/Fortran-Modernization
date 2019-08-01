#! /bin/bash
rm *.f90
rm -r tof90output
module load ftools
for f in *.f; do
	echo $f
done | to_f90
mkdir tof90output
mv *.f90 tof90output
for f in tof90output/*.f90; do
	read something
	clear
	echo $f
	gfortran -std=f95 -fsyntax-only $f
	echo done
done
