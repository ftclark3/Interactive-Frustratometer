# Interactive-Frustratometer
Interfaces between the Python Frustratometer and molecular visualization tools. At the moment, this includes only VMD, but I hope to add others in the future.

## VMD
VMD TCL scripts can pass coordinates to Python through a pipe and receive the results of the frustration calculation as a formatted block of text. This allows the frustratogram to be visualized in the same way as any precomputed result, but without the latency of loading the data from disk or the need to pick the conformers of interest in advance and run them through the python frustratometer outside of the VMD session.

The `interactive_frustratometry.tcl` script and the `update_frustration.py` script should be in the working directory along with the protein whose structure and frustratogram is to be visualized.

The TCL script `interactive_frustratometry.tcl` can be launched immediately after VMD startup using the `-e` command line switch, or sourced in the `tk console` after loading a molecule. A frustratogram is displayed only for the "top" molecule. If the "top" molecule is changed, the old frustratogram will be deleted from the screen and a new frustratogram will be drawn for the new "top" molecule. 

Starting the TCL script or changing the "top" molecule can take ~10 s because the TCL script must start up a python script, which must import the `Frustratometer` and initialize a few objects. However, switching through multiple conformers within the "top" entry (e.g., different structures of an NMR ensemble represented in a single PDB file or snapshots of a molecular dynamics trajectory) is very fast. For a small-to-medium-sized protein (I have tested up to ~600 residues), changes in the frustratogram as a function of a conformational transition are rendered with the same smoothness and clarity as the protein structure itself, with no noticeable latency.

While it is not necessary, the user may modify `update_frustration.py` to control how frustration is calculated or which draw commands will be sent back to VMD.

### Usage Examples
- `vmd traj.pdb -e interactive_frustratometry.tcl`, where `interactive_frustratometry.tcl` and `update_frustration.py` are in the working directory alongside `traj.pdb`
- `vmd top.pdb traj.dcd -e interactive_frustratometry.tcl`, where `interactive_frustratometry.tcl` and `update_frustration.py` are in the working directory alongside `top.pdb` and `traj.dcd`

## AI Use Disclosure
These scripts are under active development, though I hope to finalize them soon. The development of these scripts has been aided by a large language model (LLM). I am thoroughly reviewing, documenting, and testing all code, including the portions written by me and the portions written by the LLM.

## Acknowledgement
I would like to thank Steven Luo and Seonghui Kim for using/testing draft versions of these scripts.

## See Also
- I haven't tried it, but you may also be interested in [frustrapy](https://github.com/engelberger/frustrapy/tree/main)
