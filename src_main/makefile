# FARGO_THORIN makefile
#
# THORIN modifications by
# Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
# original code by Frédéric Masset
#

# Standard Linux sequential or MPI build
ifeq ($(CC), gcc)
	CC_THORIN = gcc
	OPT_THORIN = -O3 -march=native -Wall -Wno-unknown-pragmas -Wno-unused-result\
              -I../src_reb/ -std=c99 -Wl,-rpath,src_main/ -Wpointer-arith
	PARACC_THORIN = mpicc
endif
ifeq ($(CC), icc)
	CC_THORIN = icc
	OPT_THORIN = -O3 -I../src_reb/ -std=c99 -Wl,-rpath,src_main/ -axCORE-AVX2 -xSSE4.2 -ip -Wno-unknown-pragmas
	PARACC_THORIN = icc
endif
PARAOPT_THORIN =
OPTSEQ_THORIN =

# Generic options valid for any platfom
# which has 'cc' in the path
# These are used if you leave FARGO_ARCH undefined
CC_  = cc
OPT_ = -O
OPTSEQ_ = 
PARAOPT_ =
PARACC_ = cc

#
#
#
#--------------------No Changes Needed after this line------------------------
#
#
#
SHELL		=  /bin/sh

# #THORIN: adds new object files
MAINOBJ         = LowTasks.o SideEuler.o Output.o Init.o main.o Theo.o\
		  Interpret.o SourceEuler.o TransportEuler.o\
		  Planet.o Viscosity.o Psys.o Force.o\
		  var.o Pframeforce.o split.o merge.o commbound.o fpe.o rebin.o\
		  mpiTasks.o\
                  EnergySources.o\
                  ReboundInterface.o\
		  Pebbles.o\

DUMMY		= mpi_dummy.o

COMP        = $(CC_$(FARGO_ARCH))
OPT         = $(OPT_$(FARGO_ARCH))
# #THORIN: OpenMP is enabled by setting the environment variable
# OPENMP=1. This adds the respective compilation switch as follows.
ifeq ($(OPENMP), 1)
ifeq ($(COMP), icc)
	OPT += -qopenmp
endif
ifeq ($(COMP), gcc)
	OPT += -fopenmp
endif
endif
OPTSEQ      = $(OPTSEQ_$(FARGO_ARCH))
PARAOPT     = $(PARAOPT_$(FARGO_ARCH)) -D_PARALLEL
PARACOMP    = $(PARACC_$(FARGO_ARCH))
# #THORIN: shared library librebound.so added through -L. -lrebound
LIBS        = -L. -lrebound -lm
PARALIBS    = -lmpi
AUTOINCL    = param.h param_noex.h global_ex.h

include	.config

EXENAME        = ../thorin

ARCHIVE		= $(EXENAME:../%=%.tar)
ARCHIVECOMP	= $(EXENAME:../%=%.tar.gz)

SRC = *.c
INCLUDE = *.h

ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ)
else
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(DUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(OBJ) $(EXENAME) archive
	@echo "" 
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)."
else
	@echo "This built is SEQUENTIAL."
endif
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo ""
ifeq ($(BUILD),parallel)
	@echo "gmake BUILD=sequential"
else
	@echo "gmake BUILD=parallel"
endif
	@echo ""
# #THORIN: dump note in case of a multi-threaded build
ifeq ($(OPENMP),1)
	@echo "This build also supports MULTI-THREADING (OpenMP)"
endif

$(EXENAME): $(OBJ)
	$(COMPILER) $(OBJ) $(OPTIONS) -o $(EXENAME) $(LIBRARIES)

.PHONY: conditionalrebuild
ifneq ($(BUILD),$(OLDBUILD))
conditionalrebuild: clean
	@echo "OLDBUILD = $(BUILD)" > .config
	@echo "BUILD = $(BUILD)" >> .config
else
conditionalrebuild:
endif

.oldconfig:
.config:

archive : $(SRC) $(INCL) makefile varparser.pl	
	@echo "Creating ../.source.tar.bz2"
	@tar cf ../.source.tar *.c
	@tar rf ../.source.tar *.h
	@tar rf ../.source.tar makefile
	@tar rf ../.source.tar varparser.pl
	@bzip2 -9 -f ../.source.tar

para:
	@gmake BUILD=parallel

seq:
	@gmake BUILD=sequential

$(AUTOINCL) : var.c global.h makefile varparser.pl
	@./varparser.pl

$(OBJ): fargo.h fondam.h param.h param_noex.h types.h makefile

.PHONY: clean mrproper package

#mrproper:
#	rm -f *.o *~ *.s *.il $(AUTOINCL) $(EXENAME) ../core.*\
#	*.tex *.dvi *.pdf *.ps *.log *.aux *.lint $(ARCHIVE)\
#	$(ARCHIVECOMP)

clean:
	rm -f *.o *~ *.s *.il

# package: $(ARCHIVECOMP)

# The Doxyfile provided here is relatively general and should allow
# you to build your own documentation if you modify the source. You
# should have the required executable files in your path
# (doxygen, dot, latex, epstopdf, dvips, gs). On Mac Os X I made soft
# links in /usr/bin to these executables (either in /sw/bin or
# /Applications/Doxygen.app/Content/Resources).
doc:
	doxygen Doxyfile

#$(ARCHIVECOMP) : $(SRC) $(INCL) makefile varparser.pl	Doxyfile
#	@echo "Creating $(ARCHIVECOMP)"
#	@cd ../..;tar cf $(ARCHIVE) fargo/src/*.c
#	@cd ../..;tar rf $(ARCHIVE) fargo/src/*.h
#	@cd ../..;tar rf $(ARCHIVE) fargo/src/varparser.pl
#	@cd ../..;tar rf $(ARCHIVE) fargo/src/makefile
#	@cd ../..;tar rf $(ARCHIVE) fargo/src/Doxyfile
#	@cd ../..;tar rf $(ARCHIVE) fargo/in/template.par
#	@cd ../..;tar rf $(ARCHIVE) fargo/in/Jup.cfg
#	@cd ../..;tar rf $(ARCHIVE) fargo/in/Nep.cfg
#	@cd ../..;tar rf $(ARCHIVE) fargo/in/stockholm*par
#	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out1/
#	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out2/
#	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out3/
#	@cd ../..;tar rf $(ARCHIVE) --exclude="*.*" fargo/out4/
#	@cd ../..;tar rf $(ARCHIVE) fargo/mp.pro
#	@cd ../..;tar rf $(ARCHIVE) fargo/INSTALL
#	@cd ../..;gzip -9 $(ARCHIVE)
#	@cd ../..;mv $(ARCHIVECOMP) fargo/

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS)
