import os.path
import numpy as np
from ase.io import read
from ase.io import write
import sys
import copy

#Reads in single-framed or multi-framed structure files
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

#Interpolate N frames between each frame in a .xyz file. N = Images in the class on how many extra frames to add between frame M and frame M-1, while the lattice parameters need to be provided as well
class interpXYZ:
    def __init__(self,structure_file_1, Images = 6, Aa = 10, Bb = 10, Cc = 10, Alpha = 90, Beta = 90, Gamma = 90):

        #Checks if the structure file exists
        def file_exists(filename):
            return os.path.exists(filename)

        if not file_exists(structure_file_1):
            print("Structure file #1 does not exist where specified!")
            sys.exit(1)

        readStructure1 = openStruct(structure_file_1)

        #Float check
        def isFloat(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        #Checks if the number of frames to be interpolated is correct
        def isInteger(number):
            return isinstance(number, int) and number >= 1 

        if not isInteger(Images):
            print("Your number of images value is incorrect!")
            sys.exit(1)

        #Place the altered side lengths and angles into lists
        Sides = [[Aa,Bb,Cc]]
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
        Inverted_Cell = np.linalg.inv(Cell)

        #Adds the atoms and coordinates into lists
        TotalAtoms_1 = []
        TotalCoordinates_1 = []

        for j in range(len(readStructure1.coordinates)):
            Atoms_1 = []
            Coordinates_1 = []
            for i in range(len(readStructure1.coordinates[j])):

                #converts cartesian to fractional coordinates
                content = np.matmul([readStructure1.coordinates[j][i][1], readStructure1.coordinates[j][i][2], readStructure1.coordinates[j][i][3]],Inverted_Cell)

                Atoms_1.append(readStructure1.coordinates[j][i][0])
                Coordinates_1.append(content)

            #In the event that the Atoms list is empty, it means that the cell is empty
            if len(Atoms_1) == 0:
                print("One of your structure's contain no atoms")
                sys.exit(1)
            TotalAtoms_1.append(Atoms_1)
            TotalCoordinates_1.append(Coordinates_1)

        def difference(Coordinates1, Coordinates2):
            return [Coordinates1[0] - Coordinates2[0], (Coordinates1[1] - Coordinates2[1]), (Coordinates1[2] - Coordinates2[2])]

        def dist(Coordinate):
            return np.sqrt(np.sum(np.array(Coordinate)**2))

        #Interpolate coordinates between two sets of coordinates
        def interpolate(Atoms1, Atoms2, Coordinates1, Coordinates2, Images, Cell):

            Atoms = []
            Coordinates = []

            for j in range(0, Images + 2):
                Atoms.append(Atoms1)
                Coordinates.append(copy.deepcopy(Coordinates1))

            for i in range(len(Atoms1)):
                Temp_Coordinate = Coordinates2[i]
                mindist = 10000000
                for x in [-1, 0, 1]:
                    for y in [-1, 0, 1]:
                        for z in [-1, 0, 1]:
                            coordinate2 = Coordinates2[i] + np.array([x,y,z])
                            diff = difference(Coordinates1[i], coordinate2)
                            distance = dist(np.matmul(diff, Cell))
                            if distance < mindist:
                                mindist = distance
                                Temp_Coordinate = coordinate2

                diff = np.array(np.matmul(Coordinates1[i], Cell)) - np.array(np.matmul(Temp_Coordinate, Cell))

                for j in range(0, Images + 2):
                    frame = np.matmul(np.array(Coordinates1[i]), Cell) + diff * j/(Images+1)
                    Cell_inv = np.linalg.inv(Cell)
                    frame_fractional = np.matmul(frame, Cell_inv)
                    Coordinates[j][i] = frame_fractional


            return [Atoms, Coordinates]

        i = 1

        N = len(TotalAtoms_1)-1

        #Save the interpolated frames and the original frames as an .XYZ file called Interpolated.xyz
        XYZ_Output = "%s/Interpolated.xyz" % (os.getcwd())
        with open(XYZ_Output, "w") as file_save:
            while i <= N:
                coord1 = TotalCoordinates_1[i]
                coord2 = TotalCoordinates_1[i-1]
                Atoms1 = TotalAtoms_1[i]
                Atoms2 = TotalAtoms_1[i-1]

                intermediate_step = interpolate(Atoms1, Atoms2, coord1, coord2, Images, Cell)

                Atoms = intermediate_step[0]
                Coordinates = intermediate_step[1]

                if i < N:
                    Atoms = Atoms[:len(Atoms)-1]
                    Coordinates = Coordinates[:len(Coordinates)-1]
                else:
                    Atoms = Atoms[:]
                    Coordinates = Coordinates[:]

                for p in range(len(Atoms)):
                    print(len(Atoms[p]), file=file_save)
                    print("", file=file_save)
                    for j in range(len(Atoms[p])):
                        cartesianCoordinate = np.matmul(Coordinates[p][j],[A, B, C])
                        print(Atoms[p][j], cartesianCoordinate[0], cartesianCoordinate[1], cartesianCoordinate[2], file=file_save)

                i += 1
