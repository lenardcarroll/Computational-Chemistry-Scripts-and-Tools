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

To make the two dimensional python plots, use:
```
from xyzTwoDimensionalPlot import plotFrames, makeVideo
plotFrames(filename, skip_frame_by, axis_direction, output_folder)
makeVideo(output_folder, output_filename, frames_per_second)

#OR

from xyzTwoDimensionalPlot import plotFramesCustom, makeVideo
plotFrames(filename, skip_frame_by, direction, output_folder, min_Axis_1, min_Axis_2, max_Axis_1, max_Axis_2)
makeVideo(output_folder, output_filename, frames_per_second)
```

direction = "X", "Y" or "Z". skip_frame_by means that after the current frame is done we move on by skip_frame_by. So if the value is 1, we move to the next frame. If the value is 2, we skip the next frame and move to the one after that.
