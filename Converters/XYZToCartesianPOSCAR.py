import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser()
#Cheat sheet of arguments to be used with the script
parser.add_argument("-inp1", "--input1", dest = "input1", default = "VECTORS", help="Name of input file that contains only the scaling factor and the cell vectors")
parser.add_argument("-inp2", "--input2", dest = "input2", default = "0.xyz", help="Name of .xyz file that contains the cartesian coordinates")
parser.add_argument("-out", "--output", dest = "output", default = "POSCAR", help="Name of your output file")
args = parser.parse_args()

#Read in the input file containing the vector
f = open(args.input1,"r")
content = f.readlines()
#Get the scaling factor/lattice constants
ScaleFactor = float(content[0])
#Get the Vector
Vectors = []
for i in range(1,4):
    Vectors.append(content[i].split())
VectorsAdj = []
#Multiply the Scaling Factor with the Vector
for i in Vectors:
    VectorsAdj.append([float(i[0])*ScaleFactor,float(i[1])*ScaleFactor,float(i[2])*ScaleFactor])
#Find all atoms that must be fixed
FixedVals = content[4]
commapos = []
#Check if the list of fixed values has a comma in-between. If so, append it to list commapos.
for i in range(len(FixedVals)):
    if FixedVals[i]==',':
        commapos.append(i)

FRange = []

#If there are no commas in the list of fixed values, then check if there is a dash (range of values)
if len(commapos)==0:
    x = FixedVals[:len(FixedVals)]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            FRange.append(i)
#If there is no comma or dash in the list, then it is either only one value or none that must be frozen. We check first if it is none.
    else:
    	if x!='None':
             FRange.append(int(x)-1)
#If there are commas in the list of fixed values, here is where we get all the values.
else:
    for i in range(len(commapos)):
        if i==0:
            x = FixedVals[:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    FRange.append(i)
            else:
                FRange.append(int(x)-1)
        else:
            x = FixedVals[commapos[i-1]+1:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    FRange.append(i)
            else:
                FRange.append(int(x)-1)
#This is to get values after the last comma
if len(commapos)>0:
    x =  FixedVals[commapos[len(commapos)-1]+1:]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            FRange.append(i)
    else:
        FRange.append(int(x)-1)
#We read in the XYZ file
g = open(args.input2,"r")
content = g.readlines()

AtomCoords = []
Atoms = []

for i in range(2,len(content)):
    Atoms.append(content[i].split()[0])
    AtomCoords.append([float(content[i].split()[1]),float(content[i].split()[2]),float(content[i].split()[3])])

#We write the output file
h = open(args.output,"w")
#We add a title to the top of the POSCAR file
print("XYZ to POSCAR",file=h)
#We add the scaling factor next
print("%0.16f" % ScaleFactor,file=h)
#We add the vectors after this
for i in Vectors:
    print("%0.16f" % float(i[0]),"%0.16f" % float(i[1]),"%0.16f" % float(i[2]),file=h)
#Now we find the order of atoms and how many of them there are
AtomList = []
AtomNum = []
for i in range(len(Atoms)):
	#We take the first atom and we have it take value of 1
    if i==0:
        AtomList.append(Atoms[i])
        AtomNum.append(1)
    #If the next atom is different, assign a new atom and give it value of 1
    else:
        if Atoms[i] != AtomList[len(AtomList)-1]:
            AtomList.append(Atoms[i])
            AtomNum.append(1)
        #But if the next atom is the same, increment the number of atoms there are by 1
        else:
            AtomNum[len(AtomNum)-1] += 1
#Print each of the atom numbers
for i in AtomList:
    print(i, end=" ",file=h)
#Print each of the atoms (in order)
print("",file=h)
for i in AtomNum:
    print(i, end=" ",file=h)
print("",file=h)
print("Selective dynamics",file=h)
print("Cartesian",file=h)

#Print the cartesian coordinates
AtomTot = 0
for i in AtomNum:
    AtomTot+=i

for i in range(AtomTot):
    MatrixMult = AtomCoords[i]
    if i in FRange:
        print('%0.16f' % MatrixMult[0], '%0.16f' % MatrixMult[1],'%0.16f' % MatrixMult[2],'F','F','F',file=h)  
    else:
        print('%0.16f' % MatrixMult[0], '%0.16f' % MatrixMult[1],'%0.16f' % MatrixMult[2],'T','T','T',file=h)  
h.close()
