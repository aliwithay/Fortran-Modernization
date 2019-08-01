# Developed by Aly Ammmar at NCAR
# as a SIParCS intern, summer 2019.

import os
currdir = os.getcwd()
outputdir = currdir + "/convertoutput"
convertlog = outputdir + "/convert_log.txt"
syntaxlog = outputdir + "/syntax_log.txt"
if not os.path.exists(outputdir):
    os.mkdir(outputdir)
if os.path.exists(convertlog):
    os.remove(convertlog)
if os.path.exists(syntaxlog):
    os.remove(syntaxlog)
os.system("module load ftools")
files = os.listdir(currdir)
ccount = 0
print ("Converting files...")
for f in files:
    if ".f" in f:
        outputfile = open(convertlog,'a+')
        outputfile.write('\n' + f + '\n')
        outputfile.close()
        f = f.replace(".f", "")
        ccount = ccount + 1
        os.system("/glade/u/apps/ch/opt/ftools/bin/convert <<<" + f + "/ 1>>" + convertlog)
print("Converted " + str(ccount) + " files.")
os.system("mv *.f90 " + outputdir)
files = os.listdir(outputdir)
acount = 0
print("Analyzing files...")
for f in files:
    if ".f90" in f:
        acount = acount + 1
        os.system("gfortran -std=f2008 -fsyntax-only " + outputdir + "/" + f + " 2>>" + syntaxlog)
print("Finished analyzing " + str(acount) + " files.")
