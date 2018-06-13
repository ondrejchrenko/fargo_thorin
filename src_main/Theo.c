/** \file Theo.c

A few functions that manipulate the surface density profile.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

extern real      ScalingFactor;

real Sigma(r)
real r;
{
  real cavity=1.0;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; /* This is *not* a steady state */
				/* profile, if a cavity is defined. It first needs */
				/* to relax towards steady state, on a viscous time scale */
  return cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);
}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
   SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/** Initialises the energy array. */
real Energy(r)	/* #THORIN */
real r;
{
  real fac, power, Einit;
  if (ADIABIND == 1.0) {
    mastererr ("The internal energy cannot be initialized with the adiabatic index equal to unity. Terminating now...\n");
    Einit = -1.0;		/* just a randomly chosen irelevant value of energy to prevent a compilation warning  */
    prs_exit (1);
  } else {
    fac = 1.0/(MOLWEIGHT*(ADIABIND-1.0));
    power = -SIGMASLOPE-1.0+2.0*FLARINGINDEX;
    Einit = GASCONST*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,power)*fac;
  }
  return Einit;
}

void FillEnergy()	/* #THORIN */
{
  int i;
  for (i=0; i < NRAD; i++) EnergyMed[i] = Energy(Rmed[i]);
}

void RefillEnergy (Energy)	/* similar to RefillSigma() */
PolarGrid *Energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Energy->Nrad;
  ns = Energy->Nsec;
  field = Energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

real InitCoolingTime(r)
real r;
{
  real coolt;
  coolt = COOLINGTIME*pow(r,2.0+2.0*FLARINGINDEX);
  return coolt;
}

void FillCoolingTime()
{
  int i;
  for (i=0; i<NRAD; i++) CoolingTimeMed[i] = InitCoolingTime(Rmed[i]);
}

real InitQplus(r)
real r;
{
  real viscosity, Qplus;
  viscosity = FViscosity(r);
  Qplus = 2.25*viscosity*SIGMA0*pow(r,-SIGMASLOPE-3.0);
  return Qplus;
}

void FillQplus()
{
  int i;
  for (i=0; i<NRAD; i++) QplusMed[i] = InitQplus(Rmed[i]);
}


void FillVtheta (Vtheta)
PolarGrid *Vtheta;
{
  int i, j, nr, ns, l;
  real *field, vtcorr;
  real moy;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  field = Vtheta->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    vtcorr = Rmed[i]*OmegaFrame;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l] + vtcorr;		/* use values in the inertial frame */
    }
    moy /= (real)ns;
    VthetaMed[i] = moy;
  }
}
/* <--- */

