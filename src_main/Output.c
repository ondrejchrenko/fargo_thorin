/** \file Output.c

Contains most of the functions that write the output files.
In addition to the writing of hydrodynamics files (handled by
SendOutput ()), this file also contains the functions that update
the planet.dat and bigplanet.dat files, and the functions that
seek information about the planets at a restart.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

static real     Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual;
extern real     LostMass, OmegaFrame;
extern boolean  Write_Density, Write_Velocity, IsDisk;

void EmptyPlanetSystemFile (sys)
PlanetarySystem *sys;
{
  FILE *output;
  char name[256];
  int i, n;
  n = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < n; i++) {
    sprintf (name, "%splanet%d.dat", OUTPUTDIR, i);
    output = fopenp (name, "w"); /* This empties the file */
    fclose (output);
  }
}

void WritePlanetFile (TimeStep, n)
int TimeStep;
int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Updating 'planet%d.dat'...", n);
  fflush (stdout);
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  output = fopenp (name, "a");
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",\
	   TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

void WritePlanetSystemFile (sys, t)
PlanetarySystem *sys;
int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WritePlanetFile (t, i);
  }
}
   

void WriteBigPlanetFile (TimeStep, n)
int TimeStep;
int n;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  sprintf (name, "%sbigplanet%d.dat", OUTPUTDIR, n);
  output = fopenp (name, "a");
  fprintf (output, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",\
	   TimeStep, Xplanet, Yplanet, VXplanet, VYplanet, MplanetVirtual, LostMass, PhysicalTime, OmegaFrame);
  fclose (output);
}

void WriteBigPlanetSystemFile (sys, t)
PlanetarySystem *sys;
int t;
{
  int i, n;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    Xplanet = sys->x[i];
    Yplanet = sys->y[i];
    VXplanet = sys->vx[i];
    VYplanet = sys->vy[i];
    MplanetVirtual = sys->mass[i];
    WriteBigPlanetFile (t, i);
  }
}

real GetfromPlanetFile (TimeStep, column, n)
int TimeStep, column, n;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%splanet%d.dat", OUTPUTDIR, n);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'planet%d.dat' file. Aborting restart.\n",n);
    prs_exit (1);
  }
  if (column < 2) {
    mastererr ("Invalid column number in 'planet%d.dat'. Aborting restart.\n",n);
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'planet%d.dat' file. Aborting restart.\n", TimeStep,n);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  while (column > 1) {
    pt += strspn(pt, "eE0123456789-.");
    pt += strspn(pt, "\t :=>_");
    column--;
  }
  sscanf (pt, "%lf", &value);
  return (real)value;
}

void RestartPlanetarySystem (timestep, sys)
PlanetarySystem *sys;
int timestep;
{
  int k;
  for (k = 0; k < sys->nb; k++) {
    sys->x[k] = GetfromPlanetFile (timestep, 2, k);
    sys->y[k] = GetfromPlanetFile (timestep, 3, k);
    sys->vx[k] = GetfromPlanetFile (timestep, 4, k);
    sys->vy[k] = GetfromPlanetFile (timestep, 5, k);
    sys->mass[k] = GetfromPlanetFile (timestep, 6, k);
  }
}

void WriteDiskPolar(array, number)
PolarGrid 	*array;
int 	         number;
{
  int             Nr, Ns;
  FILE           *dump;
  char 		name[80];
  real 		*ptr;
  ptr = array->Field;
  if (CPU_Master)
    sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  else
    sprintf (name, "%s%s%d.dat.%05d", OUTPUTDIR, array->Name, number, CPU_Rank);
  Nr = array->Nrad;
  Ns = array->Nsec;
  dump = fopenp (name, "w");
  masterprint ("Writing '%s%d.dat'...", array->Name, number);
  fflush (stdout);
  MPI_Barrier (MPI_COMM_WORLD);
/* We  strip  the first CPUOVERLAP rings  if  the  current  CPU  is not  the
   innermost one */
  if (CPU_Rank > 0) {
    ptr += CPUOVERLAP*Ns;
    Nr -=CPUOVERLAP ;
  }
/* We strip the last CPUOVERLAP rings if the current CPU is not the outermost
   one */
  if (CPU_Rank < CPU_Number-1) {
    Nr -=CPUOVERLAP;
  }
  fwrite (ptr, sizeof(real), Nr*Ns,dump);
  fclose(dump);
  fprintf(stdout, "%d/", CPU_Rank);  
  fflush(stdout);
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint("done\n");
}

void WriteDim () {	  
  char filename[200];
  FILE 	*dim;
  if (!CPU_Master) return;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  dim = fopenp (filename, "w");
  fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",\
	   0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  fclose (dim);
}

void SendOutput (index, dens, gasvr, gasvt, gasenerg, label)	/* #THORIN */
int          index;
PolarGrid   *dens, *gasvr, *gasvt, *gasenerg, *label;
{
  int local_index;
  if (CPU_Master)
    printf ("\n*** OUTPUT %d ***\n", index);
  local_index = index;
  if (IsDisk == YES) {
    if (Write_Density == NO) local_index = 0;
    WriteDiskPolar (dens, local_index);
    local_index = index;
    if (Write_Velocity == NO) local_index = 0;
    WriteDiskPolar (gasvr, local_index);
    WriteDiskPolar (gasvt, local_index);
    local_index = index;
    if (AdvecteLabel == YES) WriteDiskPolar (label, index);
    /* #THORIN: Following options can create files: gasenergy.dat, gastemper.dat,
     * gasdivv.dat, gasqplus.dat, gasqbalance.dat, gaspebbledens.dat,
     * gaspebblevrad.dat, gaspebblevtheta.dat */
    if (Write_Energy == YES) WriteDiskPolar (gasenerg, index);
    if (Write_Temperature == YES) WriteDiskPolar (Temperature, index);
    if (Write_Divergence == YES) WriteDiskPolar (DivergenceVelocity, index);
    if (Write_Qplus == YES) WriteDiskPolar (Qplus, index);
    if (Write_Qbalance == YES) {
      ActualizeQbalance ();
      WriteDiskPolar (Qbalance, index);
    }
    if (Pebbles == YES) {
      WritePebbles (index);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    if (Merge && (CPU_Number > 1)) merge (local_index);
  }
}

/* #THORIN: 2DO Deprecated and should be discarded */
void ActualizeQbalance ()	/* #THORIN */
{
  int nr, ns, i, j, l;
  real *qp, *qm, *qb;
  /* ----- */
  qp = Qplus->Field;
  qm = Qminus->Field;
  qb = Qbalance->Field;
  nr = Qplus->Nrad;
  ns = Qminus->Nsec;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      qb[l] = qp[l] - qm[l];
    }
  }
}

/** Writes the angular velocity of the
 * coordinate system */
void DumpOmegaFrame (TimeStep)	/* #THORIN */
int TimeStep;
{
  FILE *output;
  char name[256];
  if (!CPU_Master) return;
  printf ("Dumping 'OmegaFrame' ... ");
  fflush (stdout);
  sprintf (name, "%somegaframe.dat", OUTPUTDIR);
  output = fopenp (name, "a");
  fprintf (output, "%d\t%#.18g\n",\
	   TimeStep, OmegaFrame);
  fclose (output);
  printf ("done\n");
  fflush (stdout);
}

/** Finds the angular velocity of the
 * coordinate system in 'omegaframe.dat'
 * at a given TimeStep */
real GetOmegaFrame (TimeStep)	/* #THORIN */
int TimeStep;
{
  FILE *input;
  char name[256];
  char testline[256];
  int time;
  char *pt;
  double value;
  sprintf (name, "%somegaframe.dat", OUTPUTDIR);
  input = fopen (name, "r");
  if (input == NULL) {
    mastererr ("Can't read 'omegaframe.dat' file. Aborting restart.\n");
    prs_exit (1);
  }
  do {
    pt = fgets (testline, 255, input);
    sscanf (testline, "%d", &time);
  } while ((time != TimeStep) && (pt != NULL));
  if (pt == NULL) {
    mastererr ("Can't read entry %d in 'omegaframe.dat' file. Aborting restart.\n", TimeStep);
    prs_exit (1);
  }
  fclose (input);
  pt = testline;
  pt += strspn(pt, "eE0123456789-.");
  pt += strspn(pt, "\t :=>_");
  sscanf (pt, "%lf", &value);
  return (real)value;
}

