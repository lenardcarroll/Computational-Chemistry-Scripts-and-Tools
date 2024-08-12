Create a supercell of a single-frame or multi-frame .xyz file using:
```
from xyzsupercell import superCell
superCell(filename, list_of_atoms, scale_X, scale_Y, scale_Z, A_length, B_length, C_length, A_angle, B_angle, C_angle)
```

To shrink the unit cell and atoms of a single-frame or multi-frame .xyz file use:
```
from xyzshrinker import shrinker
shrinker(filename, reduction_factor_X, reduction_factor_Y, reduction_factor_Z, A_length, B_length, C_length, A_angle, B_angle, C_angle)
```

To interpolate additional N frames between the frames in a multi-frame .xyz file, use:
```
from xyzinterpolation import interpXYZ
interpXYZ(filename, Interpolated_Frames, A_length, B_length, C_length, A_angle, B_angle, C_angle)
```
