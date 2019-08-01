import sys

flist = os.listdir(".")
searchexp = "include 'mpif.h'"
for f in flist:
    fopen = open(f, 'rw')
    for line in fopen:
        print(line)
        if line.contain(searchexp) and not line.contain('!'):
            print("found include mpi!")
            replaceexp = '!' + line
            line = line.replace(searchexp, replaceexp)
            print("Replacing with: " + line)
        sys.stdout.write(line)