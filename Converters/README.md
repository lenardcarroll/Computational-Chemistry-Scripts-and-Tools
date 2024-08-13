To convert a .xyz file to POSCAR format, but using the cartesian format use:
```
python XYZToCartesianPOSCAR.py -inp1 VECTORS -inp2 Structure.xyz -out POSCAR
```
The VECTORS file has the format of:
```
L
X1 X2 X3
Y1 Y2 Y3
Z1 Z2 Z3
FIXED_LIST
```
L is the scaling parameter for the cell vectors. X1, X2, X3, Y1, Y2, Y3, Z1, Z2 and Z3 are the cell vectors. FIXED_LIST is a list of atoms that needs to be constrained. Use: 1, 2, 3, 4, 5 OR 1-5. You can combine and use just ```None``` if no atoms need to be fixed.

To convert a .xyz file to POSCAR format, but using the fractional format use:
```
python XYZToFractionalPOSCAR.py -inp1 VECTORS -inp2 Structure.xyz -out POSCAR
```

To convert a fractional POSCAR file to the .xyz format use:
```
python FractionalPOSCARToXYZ.py -inp POSCAR -out POSCAR.xyz
```

To convert a cartesian POSCAR file to the .xyz format use:
```
python CartesianPOSCARToXYZ.py -inp POSCAR -out POSCAR.xyz
```
