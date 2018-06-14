/** \file Init.c

Contains the functions needed to initialize the hydrodynamics arrays.
These can be initialized by reading a given output (in the case of a
restart) or by calling a function, InitEuler (), which contains
analytic prescription for the different hydrodynamics fields. Note
that this function InitEuler() is located in SourceEuler.c, which
itself calls InitGas(), in the file Pframeforce.c.
Also, note that the present file contains InitLabel(), which sets
the initial value of a passive scalar.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
PolarGrid *array;
char *fileprefix;
int filenumber;
{
  int nr,ns,c, foo=0;
  real *field;
  char name[256];
  FILE *input;
/* Simultaneous read access to the same file have been observed to give wrong results. */
/* A sequential reading is imposed below. */
				/* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &fargostat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything meanwhile */
}

void InitLabel (array)
PolarGrid *array;
{
  int nr,ns,i,j,l;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      field[l] = (Rmed[i]-RMIN)/(RMAX-RMIN);
    }
  }
}

void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label)	/* #THORIN */
PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_energy, *gas_label;
{
  int i, j, l, nr, ns;
  real *rho, *energy;
  /* ----- */
  ReadPrevDim ();
  InitEuler (gas_density, gas_v_rad, gas_v_theta, gas_energy);
  InitLabel (gas_label);
  if (Restart || InitFromFile) {
    nr = gas_density->Nrad;
    ns = gas_density->Nsec;
    rho = gas_density->Field;
    energy = gas_energy->Field;
    if (InitFromFile) {
      if (!Restart) {
	mastererr ("Initializing from files...\n");
      } else {
	mastererr ("Reading initial files to refill boundary zones for the damping condition (if used)...\n");
      }
      fflush (stderr);
      ReadfromAsciiFile (gas_density, DENSINFILE);
      ReadfromAsciiFile (gas_v_rad, VRADINFILE);
      ReadfromAsciiFile (gas_v_theta, VTHETAINFILE);
      if (EnergyEq) {
	ReadfromAsciiFile (gas_energy, TEMPERINFILE);	// 2DO option for gas label should be included
        for (i=0; i<nr; i++) {				/* !At this point, we have temperature values in the energy HD field */
  	  for (j=0; j<ns; j++) {
            l = j + i*ns;
  	    energy[l] = energy[l]*rho[l]/(ADIABIND-1.0); /* E=(T*\rho)/(\gamma-1) */
	  }
        }	
      }
    }
    RefillSigma (gas_density);			/* this is needed for BC only - even when restarting, these field reflect the INITIAL conditions */
    if (EnergyEq) RefillEnergy (gas_energy);	/* see Theo.c */
    if (DampInit) FillVtheta (gas_v_theta);     /* needed for damping BC */
    if (Restart) {	/* if Restart==YES, it overrides InitFromFile option */
      CheckRebin (NbRestart);
      MPI_Barrier (MPI_COMM_WORLD); /* Don't start reading before master has finished rebining... */
  				  /* It shouldn't be a problem though since a sequential read is */
                                  /* imposed in the ReadfromFile function below */
      mastererr ("Reading restart files...\n");
      fflush (stderr);
      ReadfromFile (gas_density, "gasdens", NbRestart);
      ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
      ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
      if (EnergyEq) {		/* !At this point, we have temperature values in the energy HD field */
	ReadfromFile (gas_energy, "gastemper", NbRestart); 
        for (i=0; i<nr; i++) {
	  for (j=0; j<ns; j++) {
            l = j + i*ns;
	    energy[l] = energy[l]*rho[l]/(ADIABIND-1.0); /* E=(T*\rho)/(\gamma-1) */
	  }
        }
      }
      ReadfromFile (gas_label, "gaslabel", NbRestart);
    } 
    /* #THORIN: following lines (in similar order to InitEuler()) are new - we need to update cs,P,T using rho and energy read at restart/or from files */
    ComputeSoundSpeed (gas_density, gas_energy);
    ComputePressureField (gas_density, gas_energy);
    ComputeTemperatureField (gas_density, gas_energy);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  WriteDim (); 
  if (Pebbles) {
    if (Restart) {
      RestartPebbleDisk (gas_density, NbRestart);	// done in case of restart
    } else {
      EquilPebbleDisk (gas_density, gas_v_rad, gas_v_theta);	// done in standard runs or if the HD fields are initialized from files
    }
  }
}

/** Enables to read a polar grid array from an ascii file */
void ReadfromAsciiFile (array, path)	/* #THORIN */
PolarGrid *array;
char *path;
{
  char ignore[512], *err;
  int foo=0, nr, ns, i, j, l, chck;
  real *field, tmp;
  FILE *input;
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &fargostat);
  input = fopen (path, "r");
  if (input == NULL) {
    fprintf (stderr, "\nError! Can't find the initialization file %s \n", path);
    prs_exit (1);
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (i = 0; i < IMIN; i++) {
    for (j = 0; j < ns; j++) {
      err = fgets (ignore, sizeof(ignore), input);
      if (err == NULL) {
        fprintf (stderr, "\nError! Number of values in the initialization file %s\n", path);
        fprintf (stderr, "is smaller than the size of hydrodynamic arrays.\n");
	prs_exit(1);
      }      
    }
  }
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      chck = fscanf (input, "%lf", &field[l]);
      if (chck == EOF) {
        fprintf (stderr, "\nError! Number of values in the initialization file %s\n", path);
        fprintf (stderr, "is smaller than the size of hydrodynamic arrays.\n");
	prs_exit(1);
      }
    }
  }
  if (CPU_Rank == CPU_Number-1) {
    chck = fscanf (input, "%lf", &tmp);
    if (chck != EOF) {
      fprintf (stderr, "\nError! Number of values in the initialization file %s\n", path);
      fprintf (stderr, "is larger than the size of hydrodynamic arrays.\n");
      prs_exit(1);
    }
  }
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD); /* previous CPUs do not touch anything meanwhile */
}
