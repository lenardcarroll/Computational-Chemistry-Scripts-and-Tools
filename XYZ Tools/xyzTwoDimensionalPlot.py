#Load in All The Modules Needed
import argparse
import pandas as pd
import numpy as np
import csv
import os
from os import path
import shutil
import matplotlib
matplotlib.use('Agg') # You might want to remove this or use a different backend, depending on the errors you get
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText
import multiprocessing
from multiprocessing import Pool
import cv2
import sys

'''
This class reads in the the structure file, converts it to the .xyz file format (if not already) and puts the coordinates into lists
'''
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

#Function to get the color of the atom
def atom_color(X, colorlist, ptable):
    return colorlist[ptable.index(X)]

#Function to get the color and label for the atom plot legend
def atom_labels(X, colorlist, ptable):
    label = X
    color = atom_color(X, colorlist, ptable)
    return "mpatches.Patch(color='%s', label='%s')" % (color,label)

#Function to determine the size of the atom in the plot
def atom_size(X, radiuslist, ptable):
    return radiuslist[ptable.index(X)]

'''
This class plots the structure file structure_file into a folder output_folder/Images, from the direction "X", "Y" or "Z" (corresponding to the axis direction), skipping every skip=N frames
'''
class plotFrames:
    def __init__(self,structure_file, skip_value, direction, output_folder):
    
        #Checks if the structure file exists
        def file_exists(structure_file):
            return os.path.exists(structure_file)

        if not file_exists(structure_file):
            print("Structure file does not exist where specified!")
            sys.exit(1)

        readStructure = openStruct(structure_file)
        
        #Checks if the skip_value is used are correct
        def isInteger(number):
            return isinstance(number, int) and number >= 1 

        if not isInteger(skip_value):
            print("Your skip value is incorrect!")
            sys.exit(1)

        '''
        Creates the directory to where the plots need to be saved
        '''
            
        directory = os.path.join(output_folder, "Images")
        parent_dir = "."
        if os.path.exists(directory) == False:
            os.mkdir(directory)
        else:
            shutil.rmtree(directory)
            os.mkdir(directory)

        '''
        Here the script goes through each frame, skipping every Nth frame (as specified by you)
        '''

        min_X = "A"
        min_Y = "B"
        min_Z = "C"
        max_X = "D"
        max_Y = "E"
        max_Z = "F"
        for p in np.arange(0,len(readStructure.coordinates),skip_value):
            #Adds the atoms and coordinates of the cell into lists
            Atoms1 = []
            x = []
            y = []
            z = []

            #Here the coordinates of the frames are stored into separate lists

            for i in range(len(readStructure.coordinates[p])):
                Atoms1.append(readStructure.coordinates[p][i][0])
                temp = [readStructure.coordinates[p][i][1], readStructure.coordinates[p][i][2], readStructure.coordinates[p][i][3]]
                x.append(temp[0])
                y.append(temp[1])
                z.append(temp[2])

            #In the event that the Atoms list is empty, it means that the cell is empty
            if len(Atoms1) == 0:
                print("Your cell contains no atoms")
                sys.exit(1)
        
            #We setup the plot here
            fig,ax = plt.subplots()
            plt.rc('grid', linestyle=':', color='gray', linewidth=1)

            #The lists below are storing all the elements of the periodic table, their size for the plot and the color assigned to them. The color was randomly chosen.
            #The radius was carefully tested based on various situations
            ptable = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Un']
            colorlist = ['#ADADAD', '#0578c6', '#838931', '#e82e6a', '#e58918', '#404040', '#004A7F', '#FF0000', '#c70102', '#bfdf85', '#013e99', '#98f7c6', '#ee383d', '#d5a4ca', '#b3c8be', '#960892', '#fb7a9d', '#26d389', '#5ebc33', '#89263e', '#10f9bd', '#bee760', '#205df5', '#b10bab', '#4163da', '#9ca948', '#31a9f3', '#b33d74', '#FF6A00', '#30ff79', '#b7d1b8', '#a092bc', '#8139f6', '#57ffd2', '#21685d', '#9eb570', '#4ccf48', '#41be51', '#2f1861', '#c90799', '#508a1d', '#7c1464', '#48a64b', '#31fa8e', '#9e5fbd', '#c77094', '#7f60b7', '#f66f74', '#715856', '#d221e4', '#de99e2', '#29c2a9', '#6069c2', '#0cc49d', '#642786', '#e405b6', '#560dd3', '#7006d7', '#f04d9d', '#5c1b63', '#cb920f', '#61521a', '#e6da8a', '#da25bc', '#bd4b95', '#94a8ec', '#f054a6', '#bc7424', '#b6582d', '#494d80', '#239f1f', '#3d738d', '#cf9c04', '#77b7c7', '#e6a211', '#c6f2cb', '#246d90', '#054298', '#FFD800', '#2e0e29', '#b03d4b', '#612d01', '#347408', '#65d366', '#057faf', '#69cf90', '#587960', '#24d8c0', '#ac6a12', '#58226c', '#ac851a', '#97adf3', '#da8f0c', '#726883', '#bb7b8b', '#40150a', '#a6c2b6', '#790946', '#38232e', '#41e8c0', '#41d468', '#4edd68', '#3f2fb9', '#3082a8', '#3cbd53', '#4b9395', '#29bdf3', '#d40c42', '#2c219c', '#6de5ef', '#d7cd83', '#6ef13a', '#e03729', '#677b60', '#632e7e', '#9a265a', '#741b0c', '#261741','#000000']
            radiuslist = [1000, 600, 5800, 4200, 3400, 2800, 2600, 2400, 2000, 1600, 7200, 6000, 5000, 4400, 4000, 4000, 4000, 2840, 8800, 7200, 6400, 5600, 5400, 5600, 5600, 5600, 5400, 5400, 5400, 5400, 5200, 5000, 4600, 4600, 4600, 4400, 9400, 8000, 7200, 6200, 5800, 5800, 5400, 5200, 5400, 5600, 6400, 6200, 6200, 5800, 5800, 5600, 5600, 5400, 10400, 8600, 7800, 7400, 7400, 7400, 7400, 7400, 7800, 7200, 7000, 7000, 7000, 7000, 7000, 7000, 7000, 6200, 5800, 5400, 5400, 5200, 5400, 5400, 5400, 6000, 7600, 7200, 6400, 7600, 8000, 8400, 11200, 8600, 7800, 7200, 7200, 7000, 7000, 7000, 7000, 7000, 7000, 6800, 6800, 6800, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]

            ax.clear() #clears the plots as to make sure multiple plots aren't made on the same file

            col = []

            #If the direction of "Z" is chosen, the plots are generated from the z-axis    
            if direction == "Z":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in combined_lists
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[2])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)

                #Here we go through each atom, and append its colour from the function atom_color to the list col
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []

                #Here the minimum and maximum z-value of the atoms from the specific frame is calculated.
                #Each of the z-values are subtracted from the maximum value, as to ensure that the top z-avalue is at z=0
                minZ = min(z_sorted)
                maxZ = max(z_sorted)
                z_alt = -maxZ + np.array(z_sorted)

                #The minimum and maximum X and Y values for the first frame is calculated, with a padding of 4 placed around it
                if min_X == "A":
                    min_X = min(x_sorted) -4
                if min_Y == "B":
                    min_Y = min(y_sorted) -4
                if max_X == "D":
                    max_X = max(x_sorted) + 4
                if max_Y == "E":
                    max_Y = max(y_sorted) + 4

                #Here we go through each atom, and find its size from the atom_size function.
                #This size is then divided by (1+z_alt[i]/10). Essentially we take the altered z-axis values from earlier
                #divide it by 10 and add it to 1. This is adjust the radius size based on its z-axis value, making
                #atoms further away look smaller, but we only want a small adjustment of the size.

                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+z_alt[i]/10))

                #The atoms are plotted based on the x and y values, their colour from the col list and their size, which is further manipulated by:
                #multiplying the size value from earlier with 10*((17/(max_Y-min_Y))*(17/(max_X-min_X))) adjusts the size of the atom more
                #This is my way of adding a perspective to the size of the atoms
                #When putting the sizes together, I tested it out at box lengths of 17, that's where the 17 above comes from. So if the user chooses
                #box lenghts of 17, the size  is not further adjusted. If the user chooses box lengths of 34, it means that the size is
                #divided by (0.5)*(0.5) = 0.25. So doubling the box size, makes the size of the atom decreases by a quarter.
                for i in range(len(col)):
                    ax.scatter(x_sorted[i], y_sorted[i], c=col[i], s=size[i]*((17/(max_Y-min_Y))*(17/(max_X-min_X))), edgecolors='black', marker='o', lw=2)

                ax.set_ylim(min_Y, max_Y)
                ax.set_xlim(min_X, max_X)    
                
                atom_types = []
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    

                #Here we make the plots
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("x-coordinates (in Å)",fontsize=20)
                plt.ylabel("y-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min$_Z$ = %0.5f\nmax$_Z$ = %0.5f" % (minZ,max(z_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')

            elif direction == "Y":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in list3
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[1])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)
                
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []
                minY = min(y_sorted)
                maxY = max(y_sorted)
                y_alt = -maxY + np.array(y_sorted)

                if min_X == "A":
                    min_X = min(x_sorted) -4
                if min_Z == "C":
                    min_Z = min(y_sorted) -4
                if max_X == "D":
                    max_X = max(x_sorted) + 4
                if max_Z == "F":
                    max_Z = max(y_sorted) + 4              
                
                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+y_alt[i]/10))
                for i in range(len(col)):
                    ax.scatter(x_sorted[i], z_sorted[i], c=col[i], s=size[i]*((17/(max_Z-min_Z))*(17/(max_X-min_X))), edgecolors='black', marker='o', lw=2)

                ax.set_ylim(min_Z, max_Z)
                ax.set_xlim(min_X, max_X)    
                
                atom_types = []
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    
                
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("x-coordinates (in Å)",fontsize=20)
                plt.ylabel("z-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min$_Y$ = %0.5f\nmax$_Y$ = %0.5f" % (minY,max(y_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')

            elif direction == "X":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in list3
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[0])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)
                
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []
                minX = min(x_sorted)
                maxX = max(x_sorted)
                x_alt = maxX - np.array(x_sorted)

                if min_Y == "B":
                    min_Y = min(x_sorted) -4
                if min_Z == "C":
                    min_Z = min(y_sorted) -4
                if max_Y == "E":
                    max_Y = max(x_sorted) + 4
                if max_Z == "F":
                    max_Z = max(y_sorted) + 4               
                
                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+x_alt[i]/10))
                for i in range(len(col)):
                    ax.scatter(y_sorted[i], z_sorted[i], c=col[i], s=size[i]*((17/(max_Z-min_Z))*(17/(max_Y-min_Y))), edgecolors='black', marker='o', lw=2)

                ax.set_ylim(min_Z, max_Z)
                ax.set_xlim(min_Y, max_Y) 
                
                atom_types = []   
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    
                
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("y-coordinates (in Å)",fontsize=20)
                plt.ylabel("z-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min$_X$ = %0.5f\nmax$_X$ = %0.5f" % (minX,max(x_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')


'''
This class lets you specify the bounds of the plot. If direction = "Z", Min1 = Min_X, Min2 = Min_Y. If direction = "X", then Min1 = Min_Y, etc.
'''
class plotFramesCustom:
    def __init__(self,structure_file, skip_value, direction, output_folder, Min1, Min2, Max1, Max2):
    
        #Checks if the structure file exists
        def file_exists(structure_file):
            return os.path.exists(structure_file)

        if not file_exists(structure_file):
            print("Structure file does not exist where specified!")
            sys.exit(1)

        readStructure = openStruct(structure_file)
        
        #Checks if the skip_value is used are correct
        def isInteger(number):
            return isinstance(number, int) and number >= 1 

        if not isInteger(skip_value):
            print("Your skip value is incorrect!")
            sys.exit(1)
            
        directory = os.path.join(output_folder, "Images")
        parent_dir = "."
        if os.path.exists(directory) == False:
            os.mkdir(directory)
        else:
            shutil.rmtree(directory)
            os.mkdir(directory)
            
        for p in np.arange(0,len(readStructure.coordinates),skip_value):
            #Adds the atoms and coordinates of the cell into lists
            Atoms1 = []
            x = []
            y = []
            z = []

            for i in range(len(readStructure.coordinates[p])):

                #converts cartesian to fractional coordinates
                content = [readStructure.coordinates[p][i][1], readStructure.coordinates[p][i][2], readStructure.coordinates[p][i][3]]
                
                if direction == "Z":

                    #Make sure all atoms are within the cell
                    X = float(content[0])
                    Y = float(content[1])
                    Z = float(content[2])   
                    
                    if X >= Min1-4 and X <= Max1+4:
                        if Y >= Min2-4 and Y<=Max2+4:
                            Atoms1.append(readStructure.coordinates[p][i][0])
                            x.append(X)
                            y.append(Y)
                            z.append(Z)
                elif direction == "Y":

                    #Make sure all atoms are within the cell
                    X = float(content[0])
                    Y = float(content[1])
                    Z = float(content[2])   
                    
                    if X >= Min1-4 and X <= Max1+4:
                        if Z >= Min2-4 and Z<=Max2+4:
                            Atoms1.append(readStructure.coordinates[p][i][0])
                            x.append(X)
                            y.append(Y)
                            z.append(Z)
                            
                elif direction == "X":

                    #Make sure all atoms are within the cell
                    X = float(content[0])
                    Y = float(content[1])
                    Z = float(content[2])   
                    
                    if Y >= Min1-4 and Y <= Max1+4:
                        if Z >= Min2-4 and Z<=Max2+4:
                            Atoms1.append(readStructure.coordinates[p][i][0])
                            x.append(X)
                            y.append(Y)
                            z.append(Z)

            #In the event that the Atoms list is empty, it means that the cell is empty
            if len(Atoms1) == 0:
                print("Your cell contains no atoms")
                sys.exit(1)
                        
            #We setup the plot here
            fig,ax = plt.subplots()
            plt.rc('grid', linestyle=':', color='gray', linewidth=1)

            #The lists below are storing all the elements of the periodic table, their size for the plot and the color assigned to them. The color was randomly chosen.
            ptable = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Un']
            colorlist = ['#ADADAD', '#0578c6', '#838931', '#e82e6a', '#e58918', '#404040', '#004A7F', '#FF0000', '#c70102', '#bfdf85', '#013e99', '#98f7c6', '#ee383d', '#d5a4ca', '#b3c8be', '#960892', '#fb7a9d', '#26d389', '#5ebc33', '#89263e', '#10f9bd', '#bee760', '#205df5', '#b10bab', '#4163da', '#9ca948', '#31a9f3', '#b33d74', '#FF6A00', '#30ff79', '#b7d1b8', '#a092bc', '#8139f6', '#57ffd2', '#21685d', '#9eb570', '#4ccf48', '#41be51', '#2f1861', '#c90799', '#508a1d', '#7c1464', '#48a64b', '#31fa8e', '#9e5fbd', '#c77094', '#7f60b7', '#f66f74', '#715856', '#d221e4', '#de99e2', '#29c2a9', '#6069c2', '#0cc49d', '#642786', '#e405b6', '#560dd3', '#7006d7', '#f04d9d', '#5c1b63', '#cb920f', '#61521a', '#e6da8a', '#da25bc', '#bd4b95', '#94a8ec', '#f054a6', '#bc7424', '#b6582d', '#494d80', '#239f1f', '#3d738d', '#cf9c04', '#77b7c7', '#e6a211', '#c6f2cb', '#246d90', '#054298', '#FFD800', '#2e0e29', '#b03d4b', '#612d01', '#347408', '#65d366', '#057faf', '#69cf90', '#587960', '#24d8c0', '#ac6a12', '#58226c', '#ac851a', '#97adf3', '#da8f0c', '#726883', '#bb7b8b', '#40150a', '#a6c2b6', '#790946', '#38232e', '#41e8c0', '#41d468', '#4edd68', '#3f2fb9', '#3082a8', '#3cbd53', '#4b9395', '#29bdf3', '#d40c42', '#2c219c', '#6de5ef', '#d7cd83', '#6ef13a', '#e03729', '#677b60', '#632e7e', '#9a265a', '#741b0c', '#261741','#000000']
            radiuslist = [1000, 600, 5800, 4200, 3400, 2800, 2600, 2400, 2000, 1600, 7200, 6000, 5000, 4400, 4000, 4000, 4000, 2840, 8800, 7200, 6400, 5600, 5400, 5600, 5600, 5600, 5400, 5400, 5400, 5400, 5200, 5000, 4600, 4600, 4600, 4400, 9400, 8000, 7200, 6200, 5800, 5800, 5400, 5200, 5400, 5600, 6400, 6200, 6200, 5800, 5800, 5600, 5600, 5400, 10400, 8600, 7800, 7400, 7400, 7400, 7400, 7400, 7800, 7200, 7000, 7000, 7000, 7000, 7000, 7000, 7000, 6200, 5800, 5400, 5400, 5200, 5400, 5400, 5400, 6000, 7600, 7200, 6400, 7600, 8000, 8400, 11200, 8600, 7800, 7200, 7200, 7000, 7000, 7000, 7000, 7000, 7000, 6800, 6800, 6800, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600, 6600]

            ax.clear() #clears the plots as to make sure multiple plots aren't made on the same file

            col = []
                
            if direction == "Z":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in list3
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[2])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)
                
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []
                maxZ = max(z_sorted)
                minZ = min(z_sorted)
                z_alt = -maxZ + np.array(z_sorted)
                
                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+z_alt[i]/10))
                for i in range(len(col)):
                    ax.scatter(x_sorted[i], y_sorted[i], c=col[i], s=size[i]*((17/(Max1-Min1))*(17/(Max2-Min2))), edgecolors='black', marker='o', lw=2)
                ax.set_ylim(Min2, Max2)
                ax.set_xlim(Min1, Max1)    
                
                atom_types = []
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    
                
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("x-coordinates (in Å)",fontsize=20)
                plt.ylabel("y-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min_Z = %0.5f\nmax_Z = %0.5f" % (minZ,max(z_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')

            elif direction == "Y":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in list3
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[1])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)
                
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []
                minY = min(y_sorted)
                maxY = max(y_sorted)
                y_alt = -maxY + np.array(y_sorted)
                
                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+y_alt[i]/10))
                for i in range(len(col)):
                    ax.scatter(x_sorted[i], z_sorted[i], c=col[i], s=size[i]*((17/(Max1-Min1))*(17/(Max2-Min2))), edgecolors='black', marker='o', lw=2)
                ax.set_ylim(Min2, Max2)
                ax.set_xlim(Min1, Max1)    
                
                atom_types = []
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    
                
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("x-coordinates (in Å)",fontsize=20)
                plt.ylabel("z-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min_Y = %0.5f\nmax_Y = %0.5f" % (minY,max(y_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')

            elif direction == "X":
                
                # Zip the lists together
                combined_lists = zip(x, y, z, Atoms1)

                # Sort the zipped list based on the values in list3
                sorted_combined_lists = sorted(combined_lists, key=lambda item: item[0])

                # Unzip the sorted list back into separate lists
                x_sorted, y_sorted, z_sorted, Atoms1 = zip(*sorted_combined_lists)
                
                for i in range(len(Atoms1)):
                    col.append(atom_color(Atoms1[i], colorlist, ptable))
                    
                size = []
                minX = min(x_sorted)
                maxX = max(x_sorted)
                x_alt = maxX - np.array(x_sorted)
                
                for i in range(len(Atoms1)):
                    size.append(atom_size(Atoms1[i], radiuslist, ptable)/(1+x_alt[i]/10))
                for i in range(len(col)):
                    ax.scatter(y_sorted[i], z_sorted[i], c=col[i], s=size[i]*((17/(Max1-Min1))*(17/(Max2-Min2))), edgecolors='black', marker='o', lw=2)
                ax.set_ylim(Min2, Max2)
                ax.set_xlim(Min1, Max1) 
                
                atom_types = []   
                
                for i in range(len(Atoms1)):
                    if Atoms1[i] not in atom_types:
                        atom_types.append(Atoms1[i])    
                
                plt.legend(handles=[eval(atom_labels(i, colorlist, ptable)) for i in atom_types],bbox_to_anchor=(1,1), loc="upper left",fontsize=20)
                plt.grid(True)
                plt.xlabel("y-coordinates (in Å)",fontsize=20)
                plt.ylabel("z-coordinates (in Å)",fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(fontsize=20)
                plt.gcf().set_size_inches(19.6, 14.4)
                #Minimum and maxmum z values are plotted as to keep track how the z component is changing
                anchored_text = AnchoredText("min_X = %0.5f\nmax_X = %0.5f" % (minX,max(x_sorted)), loc="lower left")
                ax.add_artist(anchored_text)
                plt.draw()
                fig.savefig(path.join(directory,f"frame_{p:05d}.jpg"),bbox_inches ="tight")
                plt.close('all')

class makeVideo:
    def __init__(self, output_folder, video_name, fps):
        image_folder = os.path.join(output_folder, "Images")

        # Get list of image files and sort them
        images = sorted([img for img in os.listdir(image_folder) if img.endswith(".jpg")])

        # Check if there are images in the folder
        if not images:
            raise ValueError("No images found in the directory!")

        # Read the first frame to get dimensions
        frame = cv2.imread(os.path.join(image_folder, images[0]))
        height, width, layers = frame.shape
        
        # Initialize the video writer
        video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'MJPG'), fps, (width, height))
        
        # Write each image to the video
        for image in images:
            video.write(cv2.imread(os.path.join(image_folder, image)))
        
        # Release the video writer and close all windows
        video.release()
        cv2.destroyAllWindows()
