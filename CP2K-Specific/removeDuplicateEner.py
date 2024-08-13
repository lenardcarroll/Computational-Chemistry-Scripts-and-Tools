'''
Author: Dr. Lenard Leslie Carroll
Group: AbinitSim Unit | IFF-CSIC, Madrid, Spain

This script is used for removing duplicate MD information from your .ener file.
Upon restarting a file in CP2K, CP2K will start off with the last saved MD frame (energy, velocity, position, etc.). This might not align with the settings you have.

For example, if your restart file is saved every 50 MD frames and your .ener file appended to every 2 MD frames, if the simulation has stopped or failed at frame 480, it will have only saved
the restart file at frame 450, but the .ener information at 480. When you restart the simulation, CP2K will start over from 450 and append the restarted run's .ener info to the main .enerfile,
but you will still have the first run's MD info from 450 to 480. This will cause duplicate information in your .ener file.

The below script removes the older duplicate info (replaces them with the new ones) by using:

from removeDuplicateEner import removeFramesEner
removeFramesEner("YOUR_FILE_NAME", "YOUR_OUPUT_FILE_NAME")
'''

class removeFramesEner:
    def __init__(self, filename, outputname):

        step = []
        time = []
        kinetic = []
        temperature = []
        potential = []
        cons = []
        usedtime = []

        with open(filename,"r") as f:
            fread = f.readlines()
            for i in range(1,len(fread)):
                content = fread[i].split()
                step_frame = int(content[0])
                time_frame = float(content[1])
                kinetic_frame = float(content[2])
                temperature_frame = float(content[3])
                potential_frame = float(content[4])
                cons_frame = float(content[5])
                usedtime_frame = float(content[6])

                if time_frame in time:
                    for j in range(len(time)):
                        if time[j] == time_frame:
                            step[j] = step_frame
                            kinetic[j] = kinetic_frame
                            temperature[j] = temperature_frame
                            potential[j] = potential_frame
                            cons[j] = cons_frame
                            usedtime[j ] = usedtime_frame
                else:
                    step.append(step_frame)
                    time.append(time_frame)
                    kinetic.append(kinetic_frame)
                    temperature.append(temperature_frame)
                    potential.append(potential_frame)
                    cons.append(cons_frame)
                    usedtime.append(usedtime_frame)

        with open(outputname,"w") as f:
            print("#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            "
    "Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]", file=f)
            for i in range(len(time)):
                f.write("{:10d} {:18.6f} {:18.9f} {:18.9f} {:18.9f} {:18.9f} {:18.9f}\n".format(step[i], time[i], kinetic[i], temperature[i], potential[i], cons[i], usedtime[i]))
