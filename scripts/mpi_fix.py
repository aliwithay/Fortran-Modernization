import os
import re
import sys
import fileinput

def remove_mpif(filename):
    mpiexists = False
    for line in fileinput.input(filename, inplace=1):
        replaceexp = '!' + line
        if "include 'mpif.h'" in line:
            mpiexists = True
            if "!" not in line:
                line = line.replace(line, replaceexp)
        if "integer :: mype" in line and "!" not in line:
            sys.stdout.write(replaceexp)
            line = re.sub(r"mype[,]?", "", line)
        if "common" in line and "!" not in line and "mype" in line:
            line = line.replace(line, replaceexp)
        sys.stdout.write(line)
    return mpiexists

def add_mpi_mod(filename):
    skip = False
    for line in fileinput.input(filename, inplace=1):
        if "use mpi" in line:
            skip = True
        replacewith = "use mpi\nuse mpi_mod\n" + line
        if "implicit none" in line and "!" not in line and not skip:
            line = replacewith
            skip = False
        sys.stdout.write(line)

path = sys.argv[1]
if os.path.isdir(path):
    flist = os.listdir(path)
else:
    flist.append(path)
for f in flist:
    if "mod_mpi.F90" not in f:
        mpiexists = remove_mpif(dir + f)
        if mpiexists:
            add_mpi_mod(dir + f)
