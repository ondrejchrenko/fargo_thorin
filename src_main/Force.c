/**
 * @file Force.c
 *
 * @brief Contains the function to evaluate and write the disk torques acting
 * on planets and also the function to get the thickness smoothing parameter.
 *
 * @author Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>
 *
 * @details The original ComputeForce() function was discarded, the torque computation
 * follows from the formalism of the vertical averaging which directly
 * provides the planet accelerations. The normalised torque is now part
 * of the output. UpdateLog() is called from AdvanceSystemFromDisk() when needed.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Ondřej Chrenko. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

extern boolean OpenInner, NonReflecting;

/** Calculates and writes the disk torques (both specific and normalized) 
 * acting on the planets. Uses the accelerations provided by the vertical
 * averaging approach. */
void UpdateLog (Rho, psys)	/* #THORIN */
PolarGrid *Rho;
PlanetarySystem *psys;
{
  int i, nb, ii;
  real x, y, r, m, ax, ay;
  real *rho, *cs, rho1D[MAX1D], cs1D[MAX1D];
  real r2, ifrac, frac, rhopl, cspl, omega, h, norm;
  real tq, tq_norm;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  rho = Rho->Field;
  cs = SoundSpeed->Field;
  mpi_make1Dprofile (rho, rho1D);	/* #THORIN: needed for the torque normalisation */
  mpi_make1Dprofile (cs, cs1D);
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    ax = psys->ax[i];
    ay = psys->ay[i];
    m = psys->mass[i];
    if (CPU_Rank == CPU_Number-1) {
      r2 = x*x+y*y;
      r = sqrt(r2);
      ifrac = GetGlobalIFrac (r);              /* get radial index w.respect to the global grid */
      frac = ifrac-floor(ifrac);
      ii = (int)ifrac;
      rhopl = rho1D[ii]*(1.0-frac) + rho1D[ii+1]*frac;
      cspl = cs1D[ii]*(1.0-frac) + cs1D[ii+1]*frac;
      omega = pow(r, -1.5);
      h = cspl/(omega*r)*SQRT_ADIABIND_INV;
      norm = m*omega*r2/h;
      norm *= norm;
      norm *= rhopl;
      tq = x*ay - y*ax;
      tq_norm = tq*m*ADIABIND/norm;
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopenp (filename, "a");
      fprintf (out, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", PhysicalTime,tq,tq_norm,x,y,ax,ay);
      fclose (out);
    }
  }
}

/** Computes the local thickness from the sound speed and applies
 * the thickness smoothing parameter. */
/* 2DO add an option for the standard isothermal version!!! */
real ThicknessSmoothing (x,y)		/* #THORIN */
real x,y;
{
  real smlength=0.0, globsmlength, r, ang, H, csavrg;
  int i, j, l, ns;
  real *cs;
  /* ----- */
  cs = SoundSpeed->Field;
  ns = SoundSpeed->Nsec;
  smlength = 0.0;		/* set the smoothing length 0 on all processes */
  r = sqrt(x*x + y*y);		/* midplane projection of the planet */
  if (r >= Rinf[Zero_or_active] && r <= Rsup[Max_or_active-1]) {	/* condition satisfied only for a CPU which contains the planet within its active zone */
    ang = atan2(y,x);		/* get the smoothing on this cpu */
    if (ang < 0.0) ang += 2.0*PI;
    i = Zero_or_active;
    while (Rsup[i] < r) i++;
    csavrg = 0.0;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      csavrg += cs[l];
    }
    csavrg /= (real)ns;
    H = csavrg*OmegaInv[i]*SQRT_ADIABIND_INV;
    smlength = THICKNESSSMOOTHING*H;
  }
  MPI_Allreduce (&smlength, &globsmlength, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); /* other CPUs have smlength = 0.0, thus the correct number (>0) calculated above is distributed */
  return globsmlength;
}
/* <--- */
