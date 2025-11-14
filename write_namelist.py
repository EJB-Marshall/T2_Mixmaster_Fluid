import subprocess
import numpy as np
import argparse
import sys



#######################################
#Parser Settings
#######################################

# Initialise parser
parser = argparse.ArgumentParser(description=\
"""This program numerically solves an IVP for the T2-symmetric 
    Einstein-Euler equations.""")

# Parse files
parser.add_argument('-K', type=float,help=\
"""The value of the sound speed parameter K.""")
parser.add_argument('-N', type=int,help=\
"""The number of grid points.""")
parser.add_argument('-f','-file', help=\
"""The name of the hdf file to be produced.""")
args = parser.parse_args()


#Check Inputs
if args.K is None:
    print("Error: argument -K is required.")
    sys.exit(1)
if args.N is None:
    print("Error: argument -N is required.")
    sys.exit(1)
if args.f is None:
    print("Error: argument -f is required.")
    sys.exit(1)


### Define Grid Parameters

nml_1 = "&Grid_Params"
x_min = 0.0
x_max = 2.0*np.pi
Nx = args.N
Ngz = 2
Nvar = 12
K = args.K


nml_2 = "&Sim_Params"
t0 = 0.0
tend = -300.0
filename = args.f


### Create namelist file
f = open("params.nml","w")

f.write(nml_1+"\n")
f.write("x_min = " + str(x_min) + "\n")
f.write("x_max = " + str(x_max) + "\n")
f.write("Nx = " + str(Nx) + "\n")
f.write("Ngz = " + str(Ngz) + "\n")
f.write("Nvar = " + str(Nvar) + "\n")
f.write("K= " + str(K) + "\n")
f.write("/\n")
f.write("\n")
f.write(nml_2+"\n")
f.write("t0 = " + str(t0) + "\n")
f.write("tend = " + str(tend) + "\n")
f.write("filename = '" + str(filename) + "'" + "\n")
f.write("/\n")
f.close() # Need to close namelist file before running Fortran!


### Run our Fortran script
subprocess.run("./main")
