import sys
import os
import subprocess
import fileinput

# Find common block
# Comment out common
# Create new modue file (template)
# copy common block to module
# Replace common block with use statement

def replacecommon (filename, name):
    #open filename for inplace editing
    #Find implicit none statement
    #Replace with use 'name'
    #Add implicit none
    for line in fileinput.input(filename, inplace=1):
        if "implicit none" in line:
            line = line.replace(line, "use " + name + "_mod\n" + line)
        sys.stdout.write(line)
        sys.stdout.flush()
    fileinput.close()
    return

def match (filename, common):
    #Check if both common declarations are the same size
    #if in-file shorter, replace with common
    for line in fileinput.input(filename, inplace=1):
        if "common" in line and common not in line:
            line = common
        sys.stdout.write(line)
        sys.stdout.flush()
    fileinput.close()
    return

#def  declarations ():
    #ensure proper variable declarations
    
#main
    #Copy name of file
    #Rename to filename_bak
    #Open filename_bak for reading
    #Create new file filename for writing
    #copy content from filename_bak to filename
    #if common block, comment out statement.
    #parse statement to find name
    #close filename
    #call replace common, pass name
    #Check module file with name mod_name exists
    #if module file does not exists,
        #Create a new module file called mod_name
        #Template:
        #module mod_name
            #common/name/[variables]
        #end mod_name
    #if module file exits:
        #match common mode statements

filename = sys.argv[1]
if not os.path.exists(filename):
    print("File " + filename + " does not exist")
    exit(1)
bfilename, ext = os.path.splitext(filename)
bfilename = bfilename + "_bak" + ext
os.rename(filename, bfilename)
if '/' in filename:
    path, filen = os.path.split(filename)
else:
    path = '.'
filein = open(bfilename, 'r')
fileout = open(filename, 'w')
for line in filein:
    if "common/" in line and "!" not in line:
        name = line[line.find("/")+1:line.rfind("/")]
        line = "!" + line
        fileout.write(line)
        fileout.flush()
        #fileout.close()
        replacecommon(filename,name)
        modfilename = path + "/mod_" + name + ".f90"
        #fileout = open(filename, 'w')
        if os.path.exists(modfilename):
            match(filename, line)
        else:
            modfile = open(modfilename, 'w')
            modfile.write("subroutine " + name + "_mod\n")
            modfile.write(line)
            modfile.write("end module " + name + "_mod")
            modfile.close()
    else:
        fileout.write(line)
        fileout.flush()
fileout.close()
filein.close()
exit(0)
fileout = open(filename,'r')
for line in fileout:
    if "common/" + name in line:
        variables = line[line.rfind('/'):len(line) - 1]
        varlist = variables.split(variables, ',')