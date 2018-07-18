#!/usr/bin/env python

### An auxiliary python script to convert binary
### output files from the first example calculation into ascii
### input files for the second example calculation.
### Expands the azimuthal span of the grid.
###
### Author: Ondrej Chrenko
### chrenko@sirrah.troja.mff.cuni.cz
### 2017

import numpy as np

### Set the grid parameters
Nrad=256
NsecRelax=4
NsecFull=512
### Select the files' id
nout=40

for filetype in ["gasdens", "gastemper", "gasvrad", "gasvtheta"]:
  outfile = open(filetype + ".cfg", 'w')
  fargogrid1D = np.fromfile("../out_relax/" + filetype + str(nout) + ".dat", dtype=np.double)
  for i in range (0,Nrad):
    avrg = 0.0
    for j in range (0,NsecRelax): 
      l =  j + i*NsecRelax
      avrg = avrg + fargogrid1D[l]
    avrg = avrg / NsecRelax
    for j in range (0,NsecFull):
      outfile.write ("%#.18g\n" % avrg)
  outfile.close()
