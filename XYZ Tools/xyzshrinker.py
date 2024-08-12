import os.path
import numpy as np
from ase.io import read
from ase.io import write
import sys

#Load in structure file, whether single-framed or multi-framed
class openStruct:
    def __init__(self,file,start=None,last=None):
        if file[-3:] != 'xyz':
            try:
                # Attempt to read the file with ASE and write to a temporary file if successful
                temp = read(file, index=':')
                write("temp.xyz", temp)
                file = "temp.xyz"
            except Exception as e:
                print(f"Error: Unable to read the file with ASE. {e}")
                sys.exit(1)  # Exit the program with an error status

        g = open(file,"r")
        self.content = g.readlines()
        self.numberofatoms = int(int(self.content[0].split()[0]))
        self.numberofrows = int(len(self.content)/(self.numberofatoms+2))
        allcoordinates = []
        if start==None:
            start = 0
        if last==None:
            last = self.numberofrows
        for i in range(start,last):
            coordinates = []
            for j in range(i*(self.numberofatoms+2)+2,i*(self.numberofatoms+2)+self.numberofatoms+2):
                coordinates.append([self.content[j].split()[0],float(self.content[j].split()[1]),float(self.content[j].split()[2]),float(self.content[j].split()[3])])
            allcoordinates.append(coordinates)
        self.coordinates = allcoordinates
        g.close()

        if file=="temp.xyz":
            os.remove(file)

#This class deals with multi-frame files
class shrinker:
    def __init__(self, structure_file = "temp.xyz", Reduction_X = 1, Reduction_Y = 1, Reduction_Z = 1, Aa = 10, Bb = 10, Cc = 10, Alpha = 90, Beta = 90, Gamma = 90):

        #Checks if the structure file exists
        def file_exists(filename):
            return os.path.exists(filename)

        if not file_exists(structure_file):
            print("Structure file does not exist where specified!")
            sys.exit(1)

        readStructure = openStruct(structure_file)

        #Float check
        def isFloat(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        #Checks if the reduction parameters used are correct
        def isInteger(number):
            return isinstance(number, int) and number >= 1 

        if not isInteger(Reduction_X) or not isInteger(Reduction_Y) or not isInteger(Reduction_Z):
            print("Your reduction factors are incorrect!")
            sys.exit(1)

        #Place the altered side lengths and angles into lists
        Sides = [[Aa/Reduction_X,Bb/Reduction_Y,Cc/Reduction_Z]]
        Angles = [[Alpha, Beta, Gamma]]

        #Create a cell vector from the new lattice parameters
        A = [Sides[0][0], 0, 0]
        B = [Sides[0][1]*np.cos((np.pi/180)*Angles[0][2]), Sides[0][1]*np.sin((np.pi/180)*Angles[0][2]), 0]
        C1 = Sides[0][2]*np.cos((np.pi/180)*Angles[0][1])
        C2 = Sides[0][2]*(np.cos((np.pi/180)*Angles[0][0]) - np.cos((np.pi/180)*Angles[0][1])*np.cos((np.pi/180)*Angles[0][2]))/(np.sin((np.pi/180)*Angles[0][2]))
        C3 = np.sqrt(Sides[0][2]**2 - C1**2 - C2**2)
        C = [C1, C2, C3]
        Cell = [A,B,C]

        #Invert the old cell vectors
        Inverted_Cell = np.linalg.inv([Reduction_X*np.array(A), Reduction_Y*np.array(B), Reduction_Z*np.array(C)])

        TotalAtoms = []
        TotalCoordinates = []

        for p in range(len(readStructure.coordinates)):
            #Adds the atoms and coordinates of the reduced cell into lists
            Atoms = []
            Coordinates = []
            for i in range(len(readStructure.coordinates[p])):

                #converts cartesian to fractional coordinates
                content = np.matmul([readStructure.coordinates[p][i][1], readStructure.coordinates[p][i][2], readStructure.coordinates[p][i][3]],Inverted_Cell)

                #Make sure all atoms are within the cell
                X = float(content[0]) % 1
                if X < 0:
                    X += 1
                Y = float(content[1]) % 1
                if Y < 0:
                    Y += 1
                Z = float(content[2]) % 1
                if Z < 0:
                    Z += 1    
             
                #Check if the atoms are in the reduced cell, and if they are, append them to Atoms and Coordinates
                if (X >= 0 and X <= 1/Reduction_X) and (Y >= 0 and Y <= 1/Reduction_Y) and (Z >=0 and Z <= 1/Reduction_Z):
                    Atoms.append(readStructure.coordinates[p][i][0])
                    Coordinates.append([X*Reduction_X,Y*Reduction_Y,Z*Reduction_Z])

            #In the event that the Atoms list is empty, it means that the reduced cell is empty
            if len(Atoms) == 0:
                print("The reduced cell contains no atoms, try to decrease the reduction parameters")
                sys.exit(1)

            #Calculate the distance between the coordinates to ensure no periodicity/duplicates
            def distance(X,Y,Z,Cell):
                Y = [Y[0] - Z[0], Y[1] - Z[1], Y[2] - Z[2]]
                X = np.matmul(X, Cell)
                Y = np.matmul(Y, Cell)
                return np.sqrt((X[0] - Y[0])**2 + (X[1] - Y[1])**2 + (X[2] - Y[2])**2)

            #Remove any duplicates
            removalList = []
            for i in range(len(Coordinates)):
                for j in range(i + 1, len(Coordinates)):
                    for k in [0, 1, -1]:
                        for l in [0, 1, -1]:
                            for m in [0, 1, -1]:
                                dist = distance(Coordinates[i], Coordinates[j], [k,l,m], Cell)
                                if dist < 0.1:
                                    removalList.append(i)

            if len(removalList) > 0:
                # Find unique elements by converting to a set
                unique_elements_set = set(removalList)

                # Convert the set back to a list
                unique_elements_list = list(unique_elements_set)

                # Sort the list in reverse order (i.e., descending)
                unique_elements_list.sort(reverse=True)

                for i in unique_elements_list:
                    Atoms.pop(i)
                    Coordinates.pop(i)

            TotalAtoms.append(Atoms)
            TotalCoordinates.append(Coordinates)

        #Saved as a multi-framed .xyz file
        XYZ_Output = "%s/Multi_Shrunk.xyz" % (os.getcwd())
        with open(XYZ_Output, "w") as file_save:
            for p in range(len(TotalCoordinates)):
                print(len(TotalCoordinates[p]), file=file_save)
                print("", file=file_save)
                for i in range(len(TotalCoordinates[p])):
                    cartesianCoordinate = np.matmul(TotalCoordinates[p][i],[A, B, C])
                    print(TotalAtoms[p][i], cartesianCoordinate[0], cartesianCoordinate[1], cartesianCoordinate[2], file=file_save)
