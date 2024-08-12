# XYZ-File-Analysis

To calculate the distance between atom X and Y from frame Z, use:
```
from xyzanalysis import read, distance
c = read(filename, start_frame_num, end_frame_num, skip_frame_num).coordinates
dist = distance(c, frame_num, atom_1_num, atom_2_num).distance
```

To calculate the bond angle between atoms W, X and Y from frame Z, use:
```
from xyzanalysis import read, bondangle
c = read(filename, start_frame_num, end_frame_num, skip_frame_num).coordinates
ba = bondangle(c, frame_num, atom_1_num, atom_2_num, atom_3_num).angle
```

To calculate the dihedral angle between atoms V, W, X and Y from frame Z, use:
```
from xyzanalysis import read, dihedralangle
c = read(filename, start_frame_num, end_frame_num, skip_frame_num).coordinates
da = dihedralangle(c, frame_num, atom_1_num, atom_2_num, atom_3_num, atom_4_num).angle
```

