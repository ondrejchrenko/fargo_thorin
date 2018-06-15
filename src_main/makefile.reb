# #THORIN: this makefile is new here; it is a modification
# of the makefile which is used in one of the example
# directories in the REBOUND package distribution.
# Purpose: compile the REBOUND code and put it into
# a shared library.

include ../src_reb/Makefile.defs

all: librebound

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../src_reb/
	@-rm -f librebound.so
	@ln -s ../src_reb/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../src_reb/ clean
