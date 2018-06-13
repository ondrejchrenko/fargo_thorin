/** \file SideEuler.c

Total mass and angular momentum monitoring, and boundary conditions.
In addition, this file contains a few low-level functions that
manipulate PolarGrid 's or initialize the forces evaluation.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

extern boolean OpenInner, NonReflecting, OuterSourceMass;


real GasTotalMass (array)
PolarGrid *array;
{
   int i,j,ns;
   real *density, total = 0.0, fulltotal=0.0;
   ns = array->Nsec;
   density = array->Field;
   if (FakeSequential && (CPU_Rank > 0)) 
     MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
   for (i = Zero_or_active; i < Max_or_active; i++) {
     for (j = 0; j < ns; j++) {
       total += Surf[i]*density[j+i*ns];
     }
   }
   if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
   }
   else
     MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (FakeSequential) {
     MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
     fulltotal = total;
   }
   return fulltotal;
}

real GasMomentum (Density, Vtheta)
PolarGrid *Density, *Vtheta;
{
   int i,j,ns;
   real *density, *vtheta, total = 0.0, fulltotal=0.0;
   ns = Density->Nsec;
   density = Density->Field;
   vtheta = Vtheta->Field;
   if (FakeSequential && (CPU_Rank > 0)) 
     MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
   for (i = Zero_or_active; i < Max_or_active; i++) {
     for (j = 1; j < ns; j++) {
       total += Surf[i]*(density[j+i*ns]+density[j-1+i*ns])*Rmed[i]*(vtheta[j+i*ns]+OmegaFrame*Rmed[i]);
     }
     total += Surf[i]*(density[i*ns]+density[i*ns+ns-1])*Rmed[i]*(vtheta[i*ns]+OmegaFrame*Rmed[i]);
   }
   if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
   }
   else
     MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (FakeSequential) {
     MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
     fulltotal = total;
   }
   return 0.5*fulltotal;
}

void DivisePolarGrid (Num, Denom, Res)
PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
       l = j+i*ns;
       abs[l] = Rmed[i] * cos(2.0*PI*(real)j/(real)ns);
       ord[l] = Rmed[i] * sin(2.0*PI*(real)j/(real)ns);
    }
  }
}
  
void OpenBoundary (Vrad, Rho, Energy)	/* #THORIN */
PolarGrid *Vrad, *Rho, *Energy;
{
  int i,j,l,ns;
  real *rho, *vr, *energy;
  if (CPU_Rank != 0) return;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vr  = Vrad->Field;
  energy = Energy->Field;
  i = 1;
#pragma omp parallel for private(l)
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l-ns] = rho[l];		/* copy first ring into ghost ring */
    energy[l-ns] = energy[l];	/* #THORIN */
    if ((vr[l+ns] > 0.0) || (rho[l] < SigmaMed[0]))
      vr[l] = 0.0; /* we just allow outflow [inwards] */
    else
      vr[l] = vr[l+ns];
  }
}

void NonReflectingBoundary (Vrad, Rho, Energy)	/* #THORIN */
PolarGrid *Vrad, *Rho, *Energy;
{
  int i,j,l,ns,nr,jp,lp,i_angle;
  real *rho, *vr, *energy, *cs;
  real dangle, mean, vr_med, csin, csout;
  /* ----- */
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  energy = Energy->Field;
  cs = SoundSpeed->Field;
  if (CPU_Rank == 0) {
    csin = 0.0;
    csout = 0.0;
    for (j=0; j<ns; j++) {
      csin += cs[j];
      csout += cs[ns+j];
    }
    csin /= (real)ns;
    csout /= (real)ns;
    i = 1;			/* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[1],-1.5)-1.0)/(.5*(csin + csout));
    dangle *= (Rmed[1]-Rmed[0]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j+i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp;
      rho[lp] = rho[l];		/* copy first ring into ghost ring */
      energy[lp] = energy[l];
      vr_med = -cs[l]*(rho[l]-SigmaMed[1])/SigmaMed[1];
      vr[l] = 2.*vr_med-vr[l+ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j] += SigmaMed[0]-mean;
    } 
    /* #THORIN ---> */
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += energy[j];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      energy[j] += EnergyMed[0]-mean;
    }
    /* <--- */
  }
  if (CPU_Rank == CPU_Number-1) {
    csin = 0.0;
    csout = 0.0;
    for (j=0; j<ns; j++) {
       csin += cs[(nr-2)*ns+j];
       csout += cs[(nr-1)*ns+j];
    }
    csin /= (real)ns;
    csout /= (real)ns;
    i = nr-1;			/* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[nr-2],-1.5)-1.0)/(.5*(csin + csout));
    dangle *= (Rmed[nr-1]-Rmed[nr-2]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j-i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp+(i-1)*ns;
      rho[l] = rho[lp];		/* copy first ring into ghost ring */
      energy[l] = energy[lp];   /* #THORIN */
      vr_med = cs[l]*(rho[l-ns]-SigmaMed[nr-2])/SigmaMed[nr-2];
      vr[l] = 2.*vr_med-vr[l-ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j+(nr-1)*ns] += SigmaMed[nr-1]-mean;
    }
    /* #THORIN ---> */
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += energy[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      energy[j+(nr-1)*ns] += EnergyMed[nr-1]-mean;
    }
    /* <--- */
  }
}

/** Sets the wave-killing factors within the damping zones;
 * inspired by de Val-Borro et al. (2006). */
void SetWaveKillingZones ()	/* #THORIN */
{
  int i;
  real DRin, DRout, tauin, tauout, damp;
  /* ----- */
  DRin = RMIN*DAMPINGRMINFRAC;     /* Boundaries of killling zones */
  DRout = RMAX*DAMPINGRMAXFRAC;
  tauin = DAMPINGPERIODFRAC*2.0*PI*pow(RMIN,1.5);
  tauout = DAMPINGPERIODFRAC*2.0*PI*pow(DRout, 1.5);	/* DRout (shortest frequency) used in outer zone */
  /* compute damping coefficients */
  for (i=Zero_or_active; i<Max_or_active; i++) {
    WaveKiller[i] = 0.0;
    if (Rmed[i] < DRin) {
      damp = (DRin - Rmed[i])/(DRin - RMIN);
      WaveKiller[i] = damp*damp/tauin;
    }
    if (Rmed[i] > DRout) {
      damp = (Rmed[i] - DRout)/(RMAX - DRout);
      WaveKiller[i] = damp*damp/tauout;
      WaveKiller[i] = damp/tauout;
    }
  }
}

/** Imposes the wave-killing boundary condition. Currently,
 * the condition is set to always damp the radial velocity
 * to zero. Additionaly, the density, azimuthal velocity
 * and energy can be damped to their initial values. */
void DampingBoundary (Vrad, Vtheta, Rho, Energy, step)	/* #THORIN */
PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
real step;
{
  int i, j, l, ns;
  real *vrad, *vtheta, *rho, *energy;
  real vrad0, vtheta0=0.0, rho0, energy0;
  real DRin, DRout, lambda;
  /* ----- */
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  energy = Energy->Field;
  DRin = RMIN*DAMPINGRMINFRAC;
  DRout = RMAX*DAMPINGRMAXFRAC;
#pragma omp parallel default(none) \
  shared(Zero_or_active,Max_or_active,Rmed,DRin,DRout,DampInit,\
         ns,vrad,rho,energy,vtheta,WaveKiller,OmegaFrame,SigmaMed,EnergyMed,VthetaMed,step) \
  private(i,vrad0,rho0,energy0,vtheta0,lambda,j,l)
  {
#pragma omp for
  for (i=Zero_or_active; i<Max_or_active; i++) {
    if (Rmed[i] < DRin || Rmed[i] > DRout) {
      vrad0 = 0.0;
      if (DampInit) {
        rho0 = SigmaMed[i];
	energy0 = EnergyMed[i];
	vtheta0 = VthetaMed[i] - Rmed[i]*OmegaFrame;	/* from inertial to corot */
      }
      lambda = WaveKiller[i]*step;
      for (j=0; j<ns; j++) {
        l = j + i*ns;
        vrad[l] = (vrad[l] + lambda*vrad0)/(1.0 + lambda);
	if (DampInit) {
  	  rho[l] = (rho[l] + lambda*rho0)/(1.0 + lambda);
          energy[l] = (energy[l] + lambda*energy0)/(1.0 + lambda);
	  vtheta[l] = (vtheta[l] + lambda*vtheta0)/(1.0 + lambda);
	}
      }
    }
  }
  }  // end of the omp parallel region
}

/** Damps the pebble disk inside the wave-killing zones
 * towards its equilibrium state. */
void DampPebbles (PebbleDens, PebbleVrad, PebbleVtheta, dt)	/* #THORIN */
PolarGrid *PebbleDens, *PebbleVrad, *PebbleVtheta;
real dt;
{
  int i, j, l, ns;
  real DRin, DRout, lambda;
  real pdens0, pvr0, pvt0;
  real *pdens, *pvr, *pvt;
  /* ----- */
  ns = PebbleDens->Nsec;
  pdens = PebbleDens->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  DRin = RMIN*DAMPINGRMINFRAC;
  DRout = RMAX*DAMPINGRMAXFRAC;
#pragma omp parallel for default(none) \
  shared (Zero_or_active,Max_or_active,Rmed,DRin,DRout, \
	  WaveKiller,PebDensInit,PebVradInit,PebVthetaInit,OmegaFrame,dt,pdens,pvr,pvt,ns) \
  private (i,j,l,lambda,pdens0,pvr0,pvt0)
  for (i=Zero_or_active; i<Max_or_active; i++) {
    if (Rmed[i] < DRin || Rmed[i] > DRout) {
      lambda = WaveKiller[i]*dt;
      pdens0 = PebDensInit[i];
      pvr0 = PebVradInit[i];
      pvt0 = PebVthetaInit[i] - Rmed[i]*OmegaFrame;
      for (j=0; j<ns; j++) {
        l = j + i*ns;
	pdens[l] = (pdens[l] + lambda*pdens0)/(1.0 + lambda);
        pvr[l] = (pvr[l] + lambda*pvr0)/(1.0 + lambda);
        pvt[l] = (pvt[l] + lambda*pvt0)/(1.0 + lambda);
      }
    }
  }
}

/* <--- */

void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i,j,l,nr,ns;
  real *rho, average_rho = 0.0, *vr, penul_vr;
  if (CPU_Rank != CPU_Number-1) return;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    average_rho += rho[l];
  }
  average_rho /= (real)ns;
  average_rho = SigmaMed[nr-1]-average_rho;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l] += average_rho;
  }
  i = nr-1;
  penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),-SIGMASLOPE);
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    vr[l] = penul_vr;
  }
}


void ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt)	/* #THORIN */
PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
real dt;
{
  /* Note: Stockholm boundary was discarded, it is now refined in DampingBoundary() */
  if (OpenInner == YES) OpenBoundary (Vrad, Rho, Energy);	/* #THORIN */
  if (NonReflecting == YES) {
    if (EnergyEq) ComputeSoundSpeed (Rho, Energy);	/* #THORIN */
    NonReflectingBoundary (Vrad, Rho, Energy);		/* #THORIN */	
  }
  if (Damping) DampingBoundary (Vrad, Vtheta, Rho, Energy, dt);	/* #THORIN */
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
PolarGrid *vtheta;
real domega;
{
  int i,j,l,nr,ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}
 
/* see e.g. GasTotalMass() to compare the MPI implementation */
real GasTotalEnergy (Rho, Vrad, Vtheta, Energy)
PolarGrid *Rho, *Vrad, *Vtheta, *Energy;
{
  int i, j, l, ns;
  real vr_cent, vt_cent;
  real total = 0.0, fulltotal = 0.0;
  real *density, *vrad, *vtheta, *energy;
  /* ----- */
  ns = Rho->Nsec;
  density = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  energy = Energy->Field;
  if (FakeSequential && (CPU_Rank > 0))
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j + i*ns;
      /* centered-in-cell radial velocity, weighted mean */
      vr_cent = (Rmed[i]-Rinf[i])*vrad[l+ns] + (Rsup[i]-Rmed[i])*vrad[l];
      vr_cent /= (Rsup[i]-Rinf[i]);
      /* centered-in-cell azimuthal velocity, arithmetic mean as the azim. velocity
         position with respect to the cell center is symetric */
      if (j < ns-1)
        vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
        vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      /* Total energy is the sum of the kinetic and internal energy. GPE? */
      total += 0.5*Surf[i]*density[l]*(vr_cent*vr_cent + vt_cent*vt_cent) + \
        Surf[i]*energy[l];
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1) MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  } else {
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential) {
     MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
     fulltotal = total;
   }
   return fulltotal;
}
/* <--- */
