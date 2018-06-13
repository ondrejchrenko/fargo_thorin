#!/bin/bash

# FARGO_THORIN shell script to clear the previous build.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# The script works with 'makefile' and 'makefile.reb'.
# It cleans the object files and the REBOUND shared library.

make clean -f makefile.reb
make clean
