import numpy as np
import os
import sys
from mathutils import Vector
import math

'''
start = Initial frame to use from the .xyz file
end = Final frame to use from the .xyz file
skip = Choose to skip N frames. If skip = 1, then each frame is used.
'''
class read:
    def __init__(self, file=None, start=None, end=None, skip=1):

        if file==None:
            file = "temp.xyz"

        if not os.path.exists(file):
            print("E: Structure file not found at path!")
            sys.exit()

        if start == None:
            start = 0

        Atoms1 = []
        X1 = []
        Y1 = []
        Z1 = []

        try:
            with open(file, "r") as f:
                fread = f.readlines()

                num_atoms = int(fread[0])
                num_lines = len(fread)
                num_frames = int(num_lines/(num_atoms + 2)) # The extra 2 comes from the number of atoms line and the comment line

                if end == None:
                    end = num_frames

                for i in np.arange(start + 1, end + 1, skip):
                    Atoms = []
                    X = []
                    Y = []
                    Z = []

                    for j in range(num_atoms*(i-1) + 2*i, num_atoms*(i) + 2*i): # The extra 2 comes from the number of atoms line and the comment line
                        content = fread[j].split()

                        Atoms.append(content[0])
                        X.append(float(content[1]))
                        Y.append(float(content[2]))
                        Z.append(float(content[3]))

                    Atoms1.append(Atoms)
                    X1.append(X)
                    Y1.append(Y)
                    Z1.append(Z)
                    
        except FileNotFoundError:
            print("E: File not found at path!")
            sys.exit()

        self.coordinates = [Atoms1, X1, Y1, Z1]
 
class distance:
    def __init__(self, coordinates, frame=None, item1=None, item2=None):

        if frame == None:
            frame = 0

        if item1 == None:
            item1 = 0

        if item2 == None:
            item2 = 1

        coord1 = [coordinates[1][frame][item1], coordinates[2][frame][item1], coordinates[3][frame][item1]]
        coord2 = [coordinates[1][frame][item2], coordinates[2][frame][item2], coordinates[3][frame][item2]] 

        '''
        vector 1 is subtracted from vector 2
        '''
        dist = np.array(coord2) - np.array(coord1)

        '''
        Vector distance is squared
        '''
        dist_squared = dist**2

        '''
        Vector is added together and the sum is square-rooted
        '''
        self.distance = np.sqrt(np.sum(dist_squared))

class bondangle:
    def __init__(self, coordinates, frame=None, item1=None, item2=None, item3=None):

        if frame == None:
            frame = 0

        if item1 == None:
            item1 = 0

        if item2 == None:
            item2 = 1

        if item3 == None:
            item3 = 2

        coord1 = [coordinates[1][frame][item1], coordinates[2][frame][item1], coordinates[3][frame][item1]]
        coord2 = [coordinates[1][frame][item2], coordinates[2][frame][item2], coordinates[3][frame][item2]] 
        coord3 = [coordinates[1][frame][item3], coordinates[2][frame][item3], coordinates[3][frame][item3]] 

        '''
        Calculates the vector difference between coordinates 1 and 2, and 3 and 2.
        '''
        BA = (np.array(coord1) - np.array(coord2))
        BC = (np.array(coord3) - np.array(coord2))

        '''
        Confirm that the atoms are not overlapping, to avoid divide by 0
        '''
        if np.linalg.norm(BA) == 0 or np.linalg.norm(BC) == 0:
            print("E: Atoms are overlapping") 

        else:
            cosine_angle = BA.dot(BC) / (np.linalg.norm(BA) * np.linalg.norm(BC))
            cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
            self.angle = math.degrees(math.acos(cosine_angle))

class dihedralangle:
    def __init__(self, coordinates, frame=None, item1=None, item2=None, item3=None, item4=None):

        if frame == None:
            frame = 0

        if item1 == None:
            item1 = 0

        if item2 == None:
            item2 = 1

        if item3 == None:
            item3 = 2

        if item4 == None:
            item4 = 3

        coord1 = [coordinates[1][frame][item1], coordinates[2][frame][item1], coordinates[3][frame][item1]]
        coord2 = [coordinates[1][frame][item2], coordinates[2][frame][item2], coordinates[3][frame][item2]] 
        coord3 = [coordinates[1][frame][item3], coordinates[2][frame][item3], coordinates[3][frame][item3]] 
        coord4 = [coordinates[1][frame][item4], coordinates[2][frame][item4], coordinates[3][frame][item4]] 


        '''
        Calculates the three vector differences first
        '''
        BA = (np.array(coord1) - np.array(coord2))
        CB = (np.array(coord3) - np.array(coord2))
        DC = (np.array(coord4) - np.array(coord3))

        '''
        Confirm that the vector differences are non-zero, that is, there is no atom overlap
        '''
        if np.linalg.norm(BA) == 0 or np.linalg.norm(CB) == 0 or np.linalg.norm(DC) == 0:
            print("E: Atoms are overlapping") 

        else:
            n1 = np.cross(BA, CB) / np.linalg.norm(np.cross(BA, CB))  # Normalized cross product of BA and CB
            n2 = np.cross(CB, DC) / np.linalg.norm(np.cross(CB, DC))  # Normalized cross product of CB and DC
            m1 = np.cross(n1, CB / np.linalg.norm(CB))  # Normalized cross product of n1 and normalized CB

            x = np.dot(n1, n2) # dot product of n1 and n2
            y = np.dot(m1, n2) # dot product of m1 and n2
            angle = math.degrees(math.atan2(y, x)) # arctam of y and x

            self.angle = 180 -1 * angle # subtract the angle from 180
