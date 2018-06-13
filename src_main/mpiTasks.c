/** \file mpiTasks.c
 *
 * Contains the function to create a 1D azimuthally-averaged
 * array from a polar array. Works also for a MPI-split
 * grid.
 *
 * @author Taken from FARGO-ADSG by Clément Baruteau.
 * */

#include "fargo.h"

void mpi_make1Dprofile (gridfield, axifield)
     real* gridfield;
     real* axifield;
{
  MPI_Request req1;
  int i, j, l;
  real *localaxifield;
  localaxifield = (real*) malloc(sizeof(real) * NRAD);
  if ( localaxifield == NULL ) {
    mastererr ("Error: Not enough memory in mpi_make1Dprofile().");
    prs_exit (1);
  }
  for ( i = 0; i < NRAD; i++ )
    localaxifield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      localaxifield[i] += gridfield[l];
    }
    localaxifield[i] /= (real)NSEC;
  }
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++ )
      axifield[i] = localaxifield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for ( i = 0; i < GLOBALNRAD; i++ ) {
	if ( i < Max_or_active )
	  axifield[i] = localaxifield[i];
	else
	  axifield[i] = 0.;
      }
      MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++)
	axifield[i+IMIN] = localaxifield[i];
      if ( CPU_Rank != (CPU_Number - 1) ) {
	MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD, &req1);
	MPI_Wait (&req1, &fargostat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  }
  free (localaxifield);
}
