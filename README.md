# Interactive-Frustratometer
Interfaces between the Python Frustratometer and tools such as VMD

## Talking to VMD
VMD TCL scripts can pass coordinates to Python through a pipe and receive the results of the frustration calculation as a formatted block of text. This allows the frustratogram to be visualized in the same way as any precomputed result, but without the latency of loading the data from disk or the need to pick the conformers of interest in advance and run them through the python frustratometer outside of the VMD session.

Usage: `vmd traj.pdb -e interactive_frustration.tcl`

## To Do
- Update the Python script to actually calculate frustration instead of passing back a predetermined result
- Add a function to the TCL script to save upon user request the frustratogram currently shown on the screen
