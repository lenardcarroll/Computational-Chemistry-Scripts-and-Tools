import os.path
import numpy as np
from ase.io import read
from ase.io import write
import sys

#Reads in the ASE supported structure files for creating a supercell.
class openStruct:
    def __init__(self,file,start=None,last=None):
        file_made = 0

        if file[-3:] != 'xyz': # Checks if the structure file is not of the XYZ format
            try:
                # Attempt to read the file with ASE and write to a temporary file if successful
                temp = read(file, index=':') #If it isn't, the script tries to ASE to read the file
                write("temp.xyz", temp) #Then the script converts the structure file to  the temporary file temp.xyz
                file = "temp.xyz" #variable file is set to "temp.xyz"
                file_made += 1
            except Exception as e: #If the file cannot be reead, the script will be quit.
                print(f"Error: Unable to read the file with ASE. {e}")
                sys.exit(1)  # Exit the program with an error status

        self.coordinates = []

        with open(file,"r") as g: #Now that we have the structure file definitely in the XYZ format we can open it
            self.content = g.readlines() #Read the content
            self.numberofatoms = int(int(self.content[0].split()[0])) #get the number of atoms from the first element in the file
            self.numberofrows = int(len(self.content)/(self.numberofatoms+2)) #this essentially calculates the number of frames in the file, by dividing the total number of lines by the number of atoms + 2
            if start==None: #if the user does not choose a start value, it will be set to start from 0
                start = 0
            if last==None: #if the user does not choose a ending/last value, it will be set to the total number of rows in the file
                last = self.numberofrows
            for i in range(start,last): #We then go through the range of start to last, so we load in the start=x frame to last=y frame
                coordinates = []
                '''
                The XYZ format looks like:
                3
                Comment Line
                A 0 0 0
                B 1 1 1
                C 4 5 7

                So we add the +2 to handle the info on the number of atoms per frame and the comment line. The range of range(i*(self.numberofatoms+2)+2,i*(self.numberofatoms+2)+self.numberofatoms+2) 
                means that for every frame, we take the frame number and read the atoms of that frame. For example, if we have 10 atoms and want to read the 1st frame, the range would be:
                0*(10+2)+2 -> 0*(10+2)+10+2
                basically
                0*(10+2)+2 -> 1*(10+2)
                OR
                2 -> 12
                So in this case we are reading the info from lines 2 to 12. For frame 2 this would be from 14 to 24.
                The coordinates are then appeneded to the coordinates list
                '''
                for j in range(i*(self.numberofatoms+2)+2,i*(self.numberofatoms+2)+self.numberofatoms+2):
                    coordinates.append([self.content[j].split()[0],float(self.content[j].split()[1]),float(self.content[j].split()[2]),float(self.content[j].split()[3])])
                self.coordinates.append(coordinates)

        #Here we check that the filename is temp.xyz, but also confirm that it was created with this script
        if file=="temp.xyz" and file_made == 1:
            os.remove(file)

#This creates a supercell for a multi-frame XYZ file, based on structure file, atoms selected to be duplicated, scaling factors and lattice parameters.
#Many information is duplicated from the other class, so only changes are included with comments
class superCell:
    def __init__(self,structure_file = "coordinates.xyz", selection = None, Scale_X = 1, Scale_Y = 1, Scale_Z = 1, Aa = 10, Bb = 10, Cc = 10, Alpha = 90, Beta = 90, Gamma = 90):

        def file_exists(filename):
            return os.path.exists(filename)

        if not file_exists(structure_file):
            print("Structure file does not exist where specified!")
            sys.exit(1)

        readStructure = openStruct(structure_file)

        def isInteger(number):
            return isinstance(number, int) and number >= 1 

        if not isInteger(Scale_X) or not isInteger(Scale_Y) or not isInteger(Scale_Z):
            print("Your scaling factors are incorrect!")
            sys.exit(1)

        A = [Aa, 0, 0]
        B = [Bb*np.cos((np.pi/180)*Gamma), Bb*np.sin((np.pi/180)*Gamma), 0]
        C1 = Cc*np.cos((np.pi/180)*Beta)
        C2 = Cc*(np.cos((np.pi/180)*Alpha) - np.cos((np.pi/180)*Beta)*np.cos((np.pi/180)*Gamma))/(np.sin((np.pi/180)*Gamma))
        C3 = np.sqrt(Cc**2 - C1**2 - C2**2)
        C = [C1, C2, C3]
        Cell = [A,B,C]

        Inverted_Cell = np.linalg.inv(Cell)

        TotalAtoms = []
        TotalCoordinates = []

        if selection == None:
            selection = range(len(readStructure.coordinates[-1]))

        #Instead of the last frame, here we go through each frame now.
        for p in range(len(readStructure.coordinates)):
            #Adds the atoms and coordinates of the cell into lists
            Atoms1 = []
            Coordinates1 = []

            for i in selection:

                #Reads the coordinates i of a particular frame p and converts it to fractional coordinates
                content = np.matmul([readStructure.coordinates[p][i][1], readStructure.coordinates[p][i][2], readStructure.coordinates[p][i][3]],Inverted_Cell)

                Atoms1.append(readStructure.coordinates[p][i][0])
                Coordinates1.append(content)

            if len(Atoms1) == 0:
                print("Your cell contains no atoms")
                sys.exit(1)

            superCellAtoms = []
            superCellCoordinates = []

            def scale(Param,Direction,Atoms,Coordinates):

                scaledAtoms = []
                scaledCoordinates = []

                for j in range(len(Coordinates)):
                    for i in np.arange(0, Param, 1):
                        if Direction == 0:
                            coord = np.array(Coordinates[j]) + np.array([i, 0, 0])
                        elif Direction == 1:
                            coord = np.array(Coordinates[j]) + np.array([0, i, 0])
                        elif Direction == 2:
                            coord = np.array(Coordinates[j]) + np.array([0, 0, i])
                        scaledAtoms.append(Atoms[j])
                        scaledCoordinates.append(coord)

                return [scaledAtoms, scaledCoordinates]

            if len(Coordinates1) > 0:   
                runX = scale(Scale_X, 0, Atoms1, Coordinates1)
                Atoms1 = runX[0]
                Coordinates1 = runX[1]
                runY = scale(Scale_Y, 1, Atoms1, Coordinates1)
                Atoms1 = runY[0]
                Coordinates1 = runY[1]
                runZ = scale(Scale_Z, 2, Atoms1, Coordinates1)
                Atoms1 = runZ[0]
                Coordinates1 = runZ[1]

                Coordinates1 = (np.array(Coordinates1)) / [Scale_X, Scale_Y, Scale_Z]

            Atoms = []
            Coordinates = []

            #Since we are dealing with coordinates we add it to two different lists
            for i in range(len(Atoms1)):
                Atoms.append(Atoms1[i])
            for i in range(len(Coordinates1)):
                Coordinates.append(Coordinates1[i])

            #Since we are dealing with many frames, we add the coordinates and atoms of frame N to different lists
            TotalAtoms.append(Atoms)
            TotalCoordinates.append(Coordinates)

        A = np.array(A) * Scale_X
        B = np.array(B) * Scale_Y
        C = np.array(C) * Scale_Z

        #Write the coordinates out as a multi-framed .xyz file
        XYZ_Output = "%s/Multi_Super.xyz" % (os.getcwd())
        with open(XYZ_Output, "w") as file_save:
            #We go through each frame
            for p in range(len(TotalCoordinates)):
                #Each frame has to have print the total number of atoms followed by a comment like, which is kept blank here
                print(len(TotalCoordinates[p]), file=file_save)
                print("", file=file_save)
                #Then we go through the coordinates of each frame
                for i in range(len(TotalCoordinates[p])):
                    #We multiply the fractional coordinates with the supercell vectors to get the cartesian coordinates, which are then saved as well
                    cartesianCoordinate = np.matmul(TotalCoordinates[p][i],[A, B, C])
                    print(TotalAtoms[p][i], cartesianCoordinate[0], cartesianCoordinate[1], cartesianCoordinate[2], file=file_save)
