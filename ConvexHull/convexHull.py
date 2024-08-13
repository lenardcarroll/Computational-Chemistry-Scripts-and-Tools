#All modules/packages needed
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#Arguments Cheat/Help Sheet
parser = argparse.ArgumentParser()
#parser.add_argument("-in", "--input", dest = "input", default = "input.xyz", help="Name of the input file")
parser.add_argument("-print", "--print", dest = "print", default = "N", help="Choose if you want to print the convex points to the console. Y for print, N for don't print.")
parser.add_argument("-save", "--save", dest = "save", default = "N", help="Choose if you want to save the convex points. Y for save, N for don't save.")
parser.add_argument("-plot", "--plot", dest = "plot", default = "Y", help="Specify if you want to plot your coordinates, with the convex points highlighted.")
parser.add_argument("-conf", "--confim", dest = "confirm", default = "N", help="If Y, the script will confirm that no points lie outside the convex, if it does, it will let you know. N is used to not do the check.")
args = parser.parse_args()

#Insert your points in here in the form [[x1,y1].[x2,y2],...,[xn,yn]]
coordinates = []

#If no coordinates are chosen, random ones will be chosen
if len(coordinates) == 0:
    import random
    ra = np.arange(-10,10,0.01)
    for i in range(100):
        coordinates.append([random.choice(ra),random.choice(ra)])

#Coordinates split up into dictionary
coordinates_dict = {"X":[],"Y":[]};

#Coordinates from dictionary to added to dataframe
#Dataframe sorted according to X-axis minimum value to maximum
for i in coordinates:
    coordinates_dict["X"].append(i[0])
    coordinates_dict["Y"].append(i[1])
df = pd.DataFrame.from_dict(coordinates_dict)
df = df.sort_values(by ='X', ascending = True)

#Coordinates list and dictionary cleared
coordinates_dict.clear()
coordinates.clear()

#Added sorted X and Y coordinates to list points
points = []
for i in range(len(df)):
    points.append([df['X'].iloc[i],df['Y'].iloc[i]])

#Add first point into list convex_points
convex_points = [[points[0][0],points[0][1]]]

#Give the first point an angle 90
angles = [90]

def sqrtdist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

#Algorithm starts, check https://github.com/lenardcarroll/myConvexHull.py/blob/main/README.md for an explanation
j = 0
while j>-1:
    altered_angles = []
    k_points = []
    unaltered_angles = []
    i = convex_points[len(convex_points)-1]
    for k in range(len(points)):
        if i!=points[k]:
            altx1y1 = [i[0] - i[0],i[1] - i[1]]
            altx2y2 = [points[k][0] - i[0],points[k][1] - i[1]]
            altx2y2 = altx2y2/np.linalg.norm(altx2y2)
            ang = np.arctan2(altx2y2[1],altx2y2[0]) - np.arctan2(altx1y1[1],altx1y1[0])
            altang = ang - angles[j] + 2*np.pi
            if altang > 2*np.pi:
                altang = altang - 2*np.pi
            altered_angles.append(altang)
            unaltered_angles.append(ang)
            k_points.append(k)
    max_altered_angle = max(altered_angles)
    convex_points_max = points[k_points[altered_angles.index(max_altered_angle)]]
    convex_points.append(convex_points_max)
    convex_points_unaltered_angle = unaltered_angles[altered_angles.index(max_altered_angle)]
    angles.append(convex_points_unaltered_angle)

    if sqrtdist(convex_points_max,convex_points[0])<0.001:
        j=-1
    if j > 0:
        points.remove(convex_points_max)
        j+=1
    elif j == 0:
        j+=1

#This is a function used to check if point is on the lines connecting the convex/polygon
def dist2(a,b,c):
    if np.abs(np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)-np.sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)-np.sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2))<0.001:
        return True

#Script confirms that point is either in or on the polygon
if args.confirm == 'Y':
    polygonsurface = []
    for i in convex_points:
        polygonsurface.append((i[0],i[1]))
    polygon = Polygon(polygonsurface)
    for i in df.index:
        if polygon.contains(Point(df['X'].iloc[i],df['Y'].iloc[i])) == False:
            convex_close = []
            for j in convex_points:
                dist = np.abs((df['X'].iloc[i]-j[0])**2+(df['Y'].iloc[i]-j[1])**2)
                if dist<0.001:
                    convex_close.append(dist)
            if len(convex_close)==0:
                dotp = []
                for k in range(len(convex_points)):
                    for l in range(len(convex_points)):
                        if k!=l:
                            if dist2(convex_points[k],convex_points[l],[df['X'].iloc[i],df['Y'].iloc[i]]) == True:
                                dotp.append(1)
                if len(dotp)==0:
                    print("Convex Hull Script Failed! Sorry :(")
                    print(df['X'].iloc[i],df['Y'].iloc[i])

#Convex points are printed to console
if args.print == 'Y':
    for i in convex_points:
        print(i[0], i[1])

#Convex points are saved to output.xyz
if args.save == 'Y':
    f = open('output.xyz', 'w')
    for i in convex_points:
        print(i[0], i[1], file=f)
    f.close()

#A basic plot is made chosing the convex points.
if args.plot == 'Y':
    X1 = []
    Y1 = []
    for i in convex_points:
        X1.append(i[0])
        Y1.append(i[1])
    plt.plot(X1, Y1, color='blue',lw=3)
    plt.scatter(df['X'], df['Y'], color='red', edgecolors='black',s=250,lw=2)
    plt.scatter(X1,Y1, color='yellow', edgecolors='black',s=250,lw=2)
    plt.gcf().set_size_inches(19.2, 14.4)
    plt.grid(False)
    plt.xlabel("x-coordinates",fontsize=20)
    plt.ylabel("y-coordinates",fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.show()
