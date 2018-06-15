#!/bin/bash

# FARGO_THORIN compilation shell script
# for a parallel MPI build.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# Similar to make.sh, but setting BUILD=parallel
# enables the MPI support.

export CC=gcc
export OPENMP=
make -f makefile.reb
export FARGO_ARCH=THORIN
make BUILD=parallel
