#!/bin/csh

# FARGO_THORIN compilation C-shell script
# for a parallel MPI build.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# Note: uses the icc compiler

setenv CC icc
setenv OPENMP 
make -f makefile.reb
setenv FARGO_ARCH THORIN
make BUILD=parallel
