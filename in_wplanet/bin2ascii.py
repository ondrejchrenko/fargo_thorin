#!/usr/bin/env python

### An auxiliary python script to convert binary
### output files from the first example calculation into ascii
### input files for the second example calculation.
###
### Author: Ondrej Chrenko
### chrenko@sirrah.troja.mff.cuni.cz
### 2017

import numpy as np

### Set the grid parameters
Nrad=256
Nsec=512
### Select the files' id
nout=40

filetype="gasdens"
outfile = open(filetype + ".cfg", 'w')
fargogrid1D = np.fromfile("../out_relax/" + filetype + str(nout) + ".dat", dtype=np.double)
for i in range (0,Nrad):
  for j in range (0,Nsec): 
    l =  j + i*Nsec
    outfile.write ("%#.18g\n" % fargogrid1D[l])
outfile.close()

filetype="gasvrad"
outfile = open(filetype + ".cfg", 'w')
fargogrid1D = np.fromfile("../out_relax/" + filetype + str(nout) + ".dat", dtype=np.double)
for i in range (0,Nrad):
  for j in range (0,Nsec): 
    l =  j + i*Nsec
    outfile.write ("%#.18g\n" % fargogrid1D[l])
outfile.close()

filetype="gasvtheta"
outfile = open(filetype + ".cfg", 'w')
fargogrid1D = np.fromfile("../out_relax/" + filetype + str(nout) + ".dat", dtype=np.double)
for i in range (0,Nrad):
  for j in range (0,Nsec): 
    l =  j + i*Nsec
    outfile.write ("%#.18g\n" % fargogrid1D[l])
outfile.close()

filetype="gastemper"
outfile = open(filetype + ".cfg", 'w')
fargogrid1D = np.fromfile("../out_relax/" + filetype + str(nout) + ".dat", dtype=np.double)
for i in range (0,Nrad):
  for j in range (0,Nsec): 
    l =  j + i*Nsec
    outfile.write ("%#.18g\n" % fargogrid1D[l])
outfile.close()
