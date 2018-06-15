/** \file SourceEuler.c 

Contains routines used by the hydrodynamical loop. More specifically,
it contains the main loop itself and all the source term substeps
(with the exception of the evaluation of the viscous force). The
transport substep is treated elsewhere. 

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */
static PolarGrid *TemperInt;
static PolarGrid *VradNew,   *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static real timeCRASH;  
extern boolean Corotating;
extern boolean Restart;		/* #THORIN */

static PolarGrid *EnergyNew, *EnergyInt;	/* #THORIN */

static int AlreadyCrashed = 0, GasTimeStepsCFL;

extern int TimeStep;
extern boolean FastTransport, IsDisk;


boolean DetectCrash (array)
PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) 
	bool = YES;
    }
  }
  return bool;
}
 
void FillPolar1DArrays ()	/* #THORIN */
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
    OmegaInv[i] = pow(Rmed[i],1.5);	/* #THORIN */
    Rmed2[i] = Rmed[i]*Rmed[i];		/* #THORIN */
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Rho, Vr, Vt, En)	/* #THORIN */
PolarGrid *Rho, *Vr, *Vt, *En;
{
  SQRT_ADIABIND_INV = 1.0/sqrt(ADIABIND);
  FillPolar1DArrays (); /* #THORIN: construct the arrays related to the grid (Rinf, Rsup, ...) here */
  FillSigma ();		/* calculate SigmaMed[] and SigmaInf[] based on the surf.density profile */
  FillEnergy ();	/* #THORIN: calculate EnergyMed; see Theo.c */
  if (EnergyEq && IterInitTemper) {
    if (!InitFromFile && !Restart) IterateInitialTemperature ();		/* #THORIN: improve the initial temperature profile by simple iterations */
  }
  InitGasDensityEnergy (Rho, En);	/* #THORIN: we initialize the density and energy
					   here because they are needed to compute
					   the sound speed, pressure and temperature;
					   see Pframeforce.c */
  InitTransport ();	/* allocate fields related to the transport step (transported momenta etc.) */
  InitViscosity ();	/* allocate fields related to the viscous terms (velocity divergence and viscous stress tensor) */
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  GravAccelRad = CreatePolarGrid(NRAD, NSEC, "agr");	/* #THORIN -----> */
  GravAccelTheta = CreatePolarGrid(NRAD, NSEC, "agt");
  if (Pebbles) {					
    PebbleGravAccelRad = CreatePolarGrid(NRAD, NSEC, "pagr");
    PebbleGravAccelTheta = CreatePolarGrid(NRAD, NSEC, "pagt");
  }
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");	/* intermediate energy, calculated in SubStep2, used in SubStep3; defined above as static */
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");	/* calculated in SubStep3, then used to update the energy HD field; defined above as static */
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");	/* see global.h */
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "temper");
  Qplus        = CreatePolarGrid(NRAD, NSEC, "qplus");
  /* note: 'Name' member of the 'PolarGrid' structure is "gas" + "argument" */
  ComputeSoundSpeed (Rho, En);		/* fill the HD pressure field; fill 'globcsvec[]' which is a 1D vector with azimuth-averaged values of sound speed */
  ComputePressureField (Rho, En);
  ComputeTemperatureField (Rho, En);
  InitGasVelocity (Vr, Vt);
  if (DampInit) FillVtheta (Vt);	/* construct the initial field of azim-averaged Vtheta for damping boundary condition */
  if (!ParametricCooling && EnergyEq) InitRadiatDiffusionFields ();		
  if (Damping || Pebbles) SetWaveKillingZones ();		
  if (Pebbles) InitPebbleArrays ();				
  if (TorqueDensity) Torque = CreatePolarGrid(NRAD, NSEC, "torque");	/* #THORIN <----- */
}

real min2 (a,b)
real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}

void AlgoGas (Rho,Vrad,Vtheta,Energy,Label,sys,rsim) 	/* #THORIN non-isothermal gas, pebble disk and rebound simulation added */
PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
PlanetarySystem *sys;
struct reb_simulation *rsim;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl;
  boolean Crashed=NO;
  gastimestepcfl = 1;
  /* ----- */
  if (EnergyEq) {       				/* #THORIN */
    ComputeSoundSpeed (Rho, Energy);
  }
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho,Vrad,Vtheta,Energy,Label); /* see commbound.c*/
    if (Pebbles) SynchronizePebbleDisc ();		  /* #THORIN */
    if (SloppyCFL == YES) 				  /* 2DO not tested */
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI) + 1.e-16;	/* #THORIN: small value added to avoid divisions by zero */
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho,Vrad,Vtheta,Energy,Label); /* #THORIN */
      if (Pebbles) {
	SynchronizePebbleDisc ();	/* #THORIN */  
        PebbleStokesNumbers (Rho);	/* #THORIN: need to update Stokes numbers as the HD background evolves */
      }
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	MinStepForRebound (rsim);			/* #THORIN */
	if (Pebbles) CriticalCharTime (Vrad, Vtheta);	/* #THORIN */
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
      if (PrescribedAccretion) ParametricAccretion (sys, dt);
      if (Pebbles) AccretePebblesOntoPlanets (sys, Rho, Energy, Vtheta, dt);	/* #THORIN */
    }
    dtemp += dt;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      FillForcesArrays (Rho, sys);		/* #THORIN: calculate grav.accelerations (instead of standard potential) */
      AdvanceSystemFromDisk (Rho, sys, dt);	/* #THORIN: use accelerations from the above function + Tanaka&Ward(04) damping */
    }
    AdvanceSystemRebound (sys, rsim, dt);	/* #THORIN */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfoFromRsim (rsim, GET) / dt;	/* #THORIN: FARGO and REBOUND not yet synchronized -> has to get info from 'rsim' */
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES) CorrectVtheta (Vtheta, domega);
      if (Pebbles) CorrectPebblesVtheta (domega);	/* #THORIN */
      OmegaFrame = OmegaNew;
    }
    RotatePsys (rsim, OmegaFrame*dt);		/* #THORIN, see psys.c */
    SynchronizeFargoRebound (sys, rsim);	/* #THORIN */
    if (IsDisk == YES) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);	/* #THORIN */
      Crashed = DetectCrash (Rho);  				/* test for negative density values */
      Crashed = DetectCrash (Energy); 				/* #THORIN: test for negative energy values */
      if (Pebbles) Crashed = DetectCrashPebbles ();		/* #THORIN: test for negative pebble values */
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	  WriteDiskPolar (Temperature, 999); /* #THORIN */
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);
      ComputePressureField (Rho, Energy); 			/* #THORIN */
      if (Pebbles) SourceTermsPebbles (Vrad, Vtheta, dt);	/* #THORIN: this must be before SubStep1 because of backreaction */
      SubStep1 (Vrad, Vtheta, Rho, dt);
      if (Pebbles) {
	SubStep1Pebbles (Vrad, Vtheta, dt);			/* #THORIN !!! has to send non-updated velocities (thus sending Vrad instead of VradInt etc) */
        if (DiffusiveParticles) ParticleDiffusion (Rho);	/* #THORIN */
      }
      SubStep2 (Rho, Energy, dt);				/* #THORIN */
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);	/* #THORIN */
      if (EnergyEq) {
        UpdateDivVelocAndStressTensor (Vrad, Vtheta, Rho);	/* #THORIN: we have to update the viscous stress tensor and div(v) here */
	SubStep3 (Rho, dt);					/* #THORIN */
	ActualiseGas (Energy, EnergyNew); 			/* #THORIN */
      }
      Transport (Rho, Vrad, Vtheta, Energy, Label, dt); 	/* #THORIN */
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);	/* #THORIN */
      ComputeTemperatureField (Rho, Energy); 			/* #THORIN */
      if (Pebbles) EvolvePebbleDisk (dt);			/* #THORIN: advance the disk of pebbles */
    }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void SubStep1 (Vrad, Vtheta, Rho, dt)	/* #THORIN */
PolarGrid *Vrad, *Vtheta, *Rho;
real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho;
  real *vradint, *vthetaint;
  real *accr, *acct, *dragr, *dragt;	/* #THORIN: gas acceleration and backreaction (needed for pebble disk) */
  real gradp, vt2, dxtheta, tmp;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  real *press; 		/* #THORIN */
  real *agr, *agt;	/* #THORIN */
  /* ----- */
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  press = Pressure->Field; 	/* #THORIN: used instead of cs^2*\rho terms */
  agr = GravAccelRad->Field;
  agt = GravAccelTheta->Field;
  if (Pebbles) {	   	/* #THORIN */
    accr  = GasAccelrad->Field;
    acct  = GasAcceltheta->Field;
    if (BackReaction) {
      dragr = DragForceRad->Field;
      dragt = DragForceTheta->Field;
    }
  }
				/* In this substep we take into account     */
				/* the source part of Euler equations       */
				/* (i.e. the R.H.S. in classical notation). */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,invdxtheta,supp_torque,tmp)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	gradp = (press[l]-press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
	vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
	vt2 = vt2*vt2;
	tmp = -gradp+0.5*(agr[l]+agr[lim])+vt2*InvRinf[i];	/* #THORIN */
	vradint[l] = vrad[l]+dt*tmp;
	if (Pebbles) {	/* #THORIN */
          accr[l] = tmp;
	  if (BackReaction) {
	    tmp = dragr[l]*2.0/(rho[l]+rho[lim]);
	    vradint[l] += dt*tmp;
	    accr[l] += tmp;
	  }
	}
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	gradp = (press[l]-press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	tmp = -gradp + 0.5*(agt[l]+agt[ljm]);		/* #THORIN */
	vthetaint[l] = vtheta[l]+dt*tmp;
	vthetaint[l] += dt*supp_torque;
	if (Pebbles) {	/* #THORIN */
          acct[l] = tmp + supp_torque;
	  if (BackReaction) {
            tmp = dragt[l]*2.0/(rho[l]+rho[ljm]);
            vthetaint[l] += dt*tmp;
	    acct[l] += tmp;
	  }
	}
      }
    }
  }
  /* #THORIN: original scheme was wrapped in a single routine ViscousTerms().
   * Now we use separate routines to: 1) update div(v) & stress tensor,
   * 2) modify the velocities using viscous terms, 3) impose Keplerian
   * rotation at the boundaries but only if the Damping condition is not used */
  UpdateDivVelocAndStressTensor (VradInt, VthetaInt, Rho);
  UpdateVelocityWithViscousTerms (VradInt, VthetaInt, Rho, dt);
  /* #THORIN: the following condition leads to poor conservation of ang. momentum;
   * not used anymore
   * if (!Damping) ImposeKeplerianEdges (VthetaInt); */
}

void SubStep2 (Rho, Energy, dt)	/* #THORIN */
PolarGrid *Rho, *Energy;
real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;	/* #THORIN */
  real dxtheta, invdxtheta;
  real dv;
  nr = VradInt->Nrad;
  ns = VradInt->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;	/* note: while intermediate values of vr and vt
				   are already computed at this point, this is
				   the first time we update the energy. We thus
				   calculate the value asociated with the
				   'Intermediate' array here. */
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)	/* #THORIN: dxtheta, ljm and lim deleted below (not needed in this loop) */
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    /* #THORIN: lip, ljm and ljp deleted below (not needed in this loop) */
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    /* #THORIN: lim, lip and ljp deleted below (not needed in this loop) */
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    /* #THORIN: Artifical source term of the energy equation stemming from
     * the artificial viscosity which was designed
     * to suppress shocks on the staggered mesh.
     * (see Stone & Norman 1992). */
    if (EnergyEq) {
#pragma omp for nowait
      for (i=0; i<nr; i++) {
        dxtheta = 2.0*PI/(real)ns*Rmed[i];
        invdxtheta = 1.0/dxtheta;
        for (j = 0; j < ns; j++) {
          l = j+i*ns;
          lip = l+ns;
          ljp = l+1;
          if (j == ns-1) ljp = i*ns;
          energyint[l] = energy[l] -                            \
            dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -       \
            dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
	}
      }
    }
  }
}

/** Numerical step reponsible for the energy
 * update. Calls the energy equation solver. */
void SubStep3 (Rho, dt)		/* #THORIN */
PolarGrid *Rho;
real dt;
{
  real *rho, *energy, *energynew, *qplus, *Trr, *Trp, *Tpp, *divergence;
  int i, j, l, nr, ns;
  real viscosity, r0, r1, r2, qplus1, qplus2, tmp1, tmp2;
  /* --- */
  rho = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  energy = EnergyInt->Field;	/* we use the intermediate value of energy
				   calculated in SubStep2 as the initial value here */
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  /* Stuff needed to extrapolate the heating term for i=0 */
  r0 = Rmed[0];
  r1 = Rmed[1];
  r2 = Rmed[2];
#pragma omp parallel default(none) \
  shared (nr,ns,Rmed,rho,Trr,Trp,Tpp,divergence,qplus,r0,r1,r2,EnergyMed,\
          dt,SigmaMed,CoolingTimeMed,energy,QplusMed,energynew,\
	  ParametricCooling,ADIABIND) \
  private (i,j,l,viscosity,qplus1,qplus2,tmp1,tmp2)
  {
#pragma omp for
  /* As Trp component of the viscosity tensor is defined from i=1,
   * we first calculate the source term caused by viscous heating
   * starting from i=1 */
  for (i=1; i<nr; i++) {
    viscosity = FViscosity (Rmed[i]);
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      if (viscosity != 0.0) {	/* !!! FACTOR 2.0 NEXT TO THE rp TERM BELOW (this was wrong in ADSG version) */
	qplus[l] = 0.5/viscosity/rho[l]*(Trr[l]*Trr[l] +      \
		                         2.0*Trp[l]*Trp[l] +      \
					 Tpp[l]*Tpp[l]);
        qplus[l] += (2.0/9.0)*viscosity*rho[l]*divergence[l]*divergence[l];
      } else {
        qplus[l] = 0.0;
      }
    }
  }
  /* Now we extrapolate to find the heating term for i=0 */
  viscosity = FViscosity (Rmed[0]);
#pragma omp for
  for (j=0; j<ns; j++) {	/* we extrapolate towards the cells in the innermost ring ... */
    qplus1 = qplus[j+ns];	/* ... using two adjacent rings of cells */
    qplus2 = qplus[j+2*ns];
    if (viscosity != 0.0) {
      /* power-law extrapolation */
      qplus[j] = qplus1*exp( log(qplus1/qplus2)*log(r0/r1)/log(r1/r2) );
    } else {
      qplus[j] = 0.0;
    }
  }
  /* If the parametric cooling is used, update energy for all values of 'i' */
  if (ParametricCooling) {	// 2DO deprecated, could be removed
#pragma omp for
    for (i=0; i<nr; i++) {
      for (j=0; j<ns; j++) {
        l = j + i*ns;
        /* Note: QplusMed[] and CoolingTimeMed[] initialized through FillQplus()
         * and FillCoolingTime() in InitGasVelocities(), see Pframeforce.c */
        tmp1 = EnergyMed[i]*dt*rho[l]/SigmaMed[i] + CoolingTimeMed[i]*energy[l] +  \
          dt*CoolingTimeMed[i]*(qplus[l] - QplusMed[i]*rho[l]/SigmaMed[i]);
        tmp2 = dt + CoolingTimeMed[i] + (ADIABIND-1.0)*dt*CoolingTimeMed[i]*divergence[l];
        energynew[l] = tmp1/tmp2;
      }
    }
  }
  }	// end of the omp parallel section
  /* Use an implicit method if radiative terms are accounted for */
  if (!ParametricCooling) ImplicitRadiativeDiffusion (Rho, EnergyInt, EnergyNew, dt);
}
/* <---- */
  		   
int ConditionCFL (Vrad, Vtheta, deltaT)		/* #THORIN */
PolarGrid *Vrad, *Vtheta;
real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, invdt5=1e-30, cs, newdt, dt;
  int ideb=0, jdeb=0;
  real itdbg1=0., itdbg2=0., itdbg3=0., itdbg4=0., itdbg5=0., mdtdbg=0.;
  /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr=0., visct=0.;
  real *soundspeed;	/* #THORIN */
  /* ----- */
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  soundspeed = SoundSpeed->Field;
  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	Vresidual[j] = vt[l];	       /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = soundspeed[l];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      if ( ViscosityAlpha || (VISCOSITY != 0.0) )
	invdt5 = FViscosity(Rmed[i])*4./pow(min2(dxrad,dxtheta),2);
      /*dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+\
			    invdt4*invdt4+invdt5*invdt5); */
      if (!Pebbles) {
	dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+\
                            invdt4*invdt4+invdt5*invdt5+invdtreb_sq);
      } else {
        dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4+invdt5*invdt5+invdtpeb_sq+invdtreb_sq);
      }
      if (dt < newdt) {
	newdt = dt;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2;
	  itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	  itdbg5=1.0/invdt5;
	  mdtdbg = newdt;
	  viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	  visct = dxtheta/dvt/4.0/CVNR/CVNR;
	}
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Artificial viscosity limit     : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Physical viscosity limit       : %g\n", itdbg5);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}

/** Updates the sound speed over the mesh. */
void ComputeSoundSpeed (Rho, Energy)	/* #THORIN */
PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *rho, *energy, *cs;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  energy = Energy->Field;
  cs = SoundSpeed->Field;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      if (!EnergyEq) {
        cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      } else {
        cs[l] = sqrt( ADIABIND*(ADIABIND-1.0)*energy[l]/rho[l] );
      }
    }
  }
  /* this will update the radial field which holds the azimuth-averaged values of sound speed, if needed */
  if (ViscosityAlpha) mpi_make1Dprofile (cs, globcsvec);
}

/** Updates the pressure over the mesh. */
void ComputePressureField (Rho, Energy)	/* #THORIN */
PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *rho, *energy, *pressure, *cs;
  /* ---- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  energy = Energy->Field;
  pressure = Pressure->Field;
  cs = SoundSpeed->Field;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      if (!EnergyEq) {
	pressure[l] = rho[l]*cs[l]*cs[l];
      } else {
        pressure[l] = (ADIABIND-1.0)*energy[l];
      }
    }
  }
}

/** Updates the temperature over the mesh. */
void ComputeTemperatureField (Rho, Energy)	/* #THORIN */
PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real tmp1, tmp2;
  real *rho, *energy, *temperature, *pressure;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  energy = Energy->Field;
  temperature = Temperature->Field;
  pressure = Pressure->Field;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      tmp2 = 1.0/(rho[l]*GASCONST);
      if (EnergyEq) {
        tmp1 = energy[l]*MOLWEIGHT*(ADIABIND-1.0);
        temperature[l] = tmp1*tmp2;
      } else {
        tmp1 = MOLWEIGHT*pressure[l];
	temperature[l] = tmp1*tmp2;
      }
    }
  }
  /* this will update the radial field which holds the azimuth-averaged values of temperature, if needed */
  if (ViscosityAlpha && AlphaFlock) mpi_make1Dprofile (temperature, globtvec);
}
