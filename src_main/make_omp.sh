#!/bin/bash

# FARGO_THORIN compilation shell script
# for a single-machine build with multithread support.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# Similar to make.sh, but setting OPENMP=1 enables
# the OpenMP support.

export CC=gcc
export OPENMP=1
make -f makefile.reb
export FARGO_ARCH=THORIN
make BUILD=sequential
