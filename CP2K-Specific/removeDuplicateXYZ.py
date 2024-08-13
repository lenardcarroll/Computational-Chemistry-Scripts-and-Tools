'''
Author: Dr. Lenard Leslie Carroll
Group: AbinitSim Unit | IFF-CSIC, Madrid, Spain

This script is used for removing duplicate MD frames from your .xyz file.
Upon restarting a file in CP2K, CP2K will start off with the last saved MD frame (energy, velocity, position, etc.). This might not align with the settings you have.

For example, if your restart file is saved every 50 MD frames and your .xyz file appended to every 10 MD frames, if the simulation has stopped or failed at frame 480, it will have only saved
the restart file at frame 450, but the position at 480. When you restart the file, CP2K will start over from 450 and append the restarted run's coordinates to the main coordinate file,
but you will still have the first run's MD frames of 460, 470 and 480. This will cause duplicate frames in your trajectory file.

The below script removes the older duplicate frames (replaces them with the new ones) by using:

from removeDuplicateXYZ import removeFrames
removeFrames("YOUR_FILE_NAME", "YOUR_OUPUT_FILE_NAME")
'''

class removeFrames:
    def __init__(self, filename, outputname):

        full_coordinates = []
        time = []
        index = []
        energy = []

        with open(filename,"r") as f:
            fread = f.readlines()
            atom_count = int(fread[0])
            total_lines = int(len(fread))
            frames = int(total_lines/(atom_count+2))
            for i in range(frames):
                frame_coordinates = []
                print()
                time_frame = float(fread[i*(atom_count +2) + 1].split(',')[1].split('=')[-1])
                index_frame = int(float(fread[i*(atom_count +2) + 1].split(',')[0].split('=')[-1]))
                energy_frame = float(fread[i*(atom_count +2) + 1].split(',')[2].split('=')[-1])
                for j in range(atom_count):
                    coord = fread[i*(atom_count +2) + 2 + j].split()
                    frame_coordinates.append([coord[0], float(coord[1]), float(coord[2]), float(coord[3])])
                if time_frame not in time:
                    full_coordinates.append(frame_coordinates)
                    time.append(time_frame)
                    index.append(index_frame)
                    energy.append(energy_frame)
                else:
                    for j in range(len(time)):
                        if time_frame == time[j]:
                            full_coordinates[j] = frame_coordinates
                            index[j] = index_frame
                            energy[j] = energy_frame

        with open(outputname,"w") as f:
            for i in range(len(full_coordinates)):
                print(len(full_coordinates[i]), file=f)
                print("i = %d, time = %f, E = %f" %(index[i], time[i], energy[i]), file=f)
                for j in range(len(full_coordinates[i])):
                    content = full_coordinates[i][j] 
                    print(content[0], content[1], content[2], content[3], file=f)
