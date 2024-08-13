import numpy as np
import argparse

parser = argparse.ArgumentParser()
#Cheat sheet of arguments to be used with the script
parser.add_argument("-inp", "--input", dest = "input", default = "POSCAR", help="Name of input file")
parser.add_argument("-out", "--output", dest = "output", default = "POSCAR.xyz", help="Name of output file")
args = parser.parse_args()

#Read in POSCAR or CONTCAR file
f = open(args.input,"r")
content = f.readlines()

#Get all atoms
Atoms = content[5].split()
AtomTot = content[6].split()

#Get total number of atoms
sumAtoms = 0
for i in AtomTot:
    sumAtoms += int(i)

#Get the coordinates of each atom
AtomCoords = []
for i in range(9,sumAtoms+9):
    AtomCoords.append([float(content[i].split()[0]),float(content[i].split()[1]),float(content[i].split()[2])])

#Write output file
f = open(args.output,"w")
#Write number of atoms and a title. Needed for XYZ file
print(sumAtoms,file=f)
print("POSCAR to XYZ",file=f)

#Print the coordinates
for i in range(len(AtomTot)):
    if i==0:
        for j in range(int(AtomTot[i])):
            MatrixMult = AtomCoords[j]
            print(Atoms[i],MatrixMult[0],MatrixMult[1],MatrixMult[2],file=f)
    else:
        sumVals = 0
        for j in range(i):
            sumVals+=int(AtomTot[j])
        for j in range(sumVals,int(AtomTot[i])+sumVals):
            MatrixMult = AtomCoords[j]
            print(Atoms[i],MatrixMult[0],MatrixMult[1],MatrixMult[2],file=f)
f.close()
