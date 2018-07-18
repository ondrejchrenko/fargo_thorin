#!/bin/csh

# FARGO_THORIN compilation C-shell script
# for a single-CPU sequential build.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# Similar to make.sh, but uses the icc compiler
# instead and runs under C-shell.

setenv CC icc
setenv OPENMP
make -f makefile.reb
setenv FARGO_ARCH THORIN
make BUILD=sequential
