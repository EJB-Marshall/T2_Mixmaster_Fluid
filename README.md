# T2_Mixmaster_Fluid
Fortran code to evolve the Einstein-Euler equations in T2 Symmetry. (Requires HDF5 library to save simulation results!)

After compiling the code, a simulation can be run using the python file write_namelist.py.
This script creates a namelist file with the simulation parameters and takes the following the command line arguments:

N - The number of grid cells<br/>
K - The sound speed parameter<br/>
f - The name of the file to be saved (NB: Need to include the .hdf5 extension in this argument! <br/>
e.g. -f results.hdf5 will save an hdf5 file called results)

Other simulation parameters, e.g. simulation length, can be modified inside the python file.

Initial data is specified in the evolve_system.f90 file. 

Output files can be plotted using plot_sim.py (requires numpy, matplotlib etc.)





