/**
 * @file Pebbles.c
 *
 * @brief Contains functions reponsible for the pebble disk initialisation,
 * evolution due to source terms and pebble accretion. Also controls several outputs.
 *
 * @author Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>
 *
 * @details The pebble disk equilibrium model is inspired by
 * Lambrechts & Johansen (2014). The evolution follows a standard
 * two-fluid approximation with linear drag coupling and particle diffusion.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Ondřej Chrenko. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

#define TRAPEZMAX 35
#define TRAPEZEPS 1.0e-7

extern boolean Restart, FastTransport;

static PolarGrid *PebbleDensTemp, *PebbleAccelrad, *PebbleAcceltheta;
static PolarGrid *PebbleSize;
static PolarGrid *EtaFaceCentered, *EtaCellCentered;
static PolarGrid *AccretedMass;

static real pebbulkdens;
// static real filtering[50];

// 2DO check everything here if it can work in isothermal case

/** Initialise polar arrays associated with the pebble disk */
void InitPebbleArrays ()
{
  PebbleDens	= CreatePolarGrid (NRAD, NSEC, "pebbledens");
  PebbleVrad	= CreatePolarGrid (NRAD, NSEC, "pebblevrad");
  PebbleVtheta	= CreatePolarGrid (NRAD, NSEC, "pebblevtheta");
  PebbleSize    = CreatePolarGrid (NRAD, NSEC, "pebblesize");
  StokesNumber	= CreatePolarGrid (NRAD, NSEC, "pebblestokes");
  PebbleDensTemp= CreatePolarGrid (NRAD, NSEC, "pebbledensitytemp");
  PebbleAccelrad= CreatePolarGrid (NRAD, NSEC, "pebbleaccelrad");
  PebbleAcceltheta= CreatePolarGrid (NRAD, NSEC, "pebbleacceltheta");
  GasAccelrad   = CreatePolarGrid (NRAD, NSEC, "gasaccelrad");
  GasAcceltheta = CreatePolarGrid (NRAD, NSEC, "gasacceltheta");
  EtaFaceCentered=CreatePolarGrid (NRAD, NSEC, "etafacecentered");
  EtaCellCentered=CreatePolarGrid (NRAD, NSEC, "eta");
  AccretedMass  = CreatePolarGrid (NRAD, NSEC, "accretedmass");
  DragForceRad  = CreatePolarGrid (NRAD, NSEC, "dragforcerad");
  DragForceTheta= CreatePolarGrid (NRAD, NSEC, "createpolargrid");
  pebbulkdens 	= PEBBLEBULKDENS/RHO2CGS;
}

/** Sets up a pebble disk in a coagulation-drift equilibrium */
void EquilPebbleDisk (Rho, Vrad, Vtheta)
PolarGrid *Rho, *Vrad, *Vtheta;
{
  EtaPressureSupport (Vtheta);
  InitPebblesViaFlux (Rho, Vrad);
}

/** Reads arrays necessary to restart the pebble disk */
void RestartPebbleDisk (Rho, index)
PolarGrid *Rho;
int index;
{
  ReadfromFile (PebbleDens, "gaspebbledens", index);
  ReadfromFile (PebbleSize, "gaspebblesize", 0);	// sidenote: stokes number will be calculated in AlgoGas()
  ReadfromFile (PebbleVrad, "gaspebblevrad", index);
  ReadfromFile (PebbleVtheta, "gaspebblevtheta", index);
  BckpFieldsForBC ();
}

/** Imposing the radial mass flux of pebbles in a
 * coagulation-drift equilibrium, initialises their
 * surface density, flow velocity and local Stokes numbers */
void InitPebblesViaFlux (Rho, Vrad)
PolarGrid *Rho, *Vrad;
{
  static real pebbleflux=0.0;
  int i,j,l,nr,ns;
  real *rho, *vr, *etafc, *etacc, *psize, *tau, *pdens, *pvr, *pvt;
  real rho1, vr1, etafc1, etacc1;	// index 1 is for azimuthally averaged vars
  real pdens1, psize1, tau1, pvr1, pvt1;
  real vtcorr, vkepl, fac, etaterm, pvrtemp;
  char command[1024];
  if (pebbleflux==0.0) pebbleflux = PEBBLEFLUX*FLUX2CU;
  pdens = PebbleDens->Field;
  tau = StokesNumber->Field;
  psize = PebbleSize->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  rho = Rho->Field;
  vr = Vrad->Field;
  etafc = EtaFaceCentered->Field;
  etacc = EtaCellCentered->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  for (i=0; i<nr; i++) {
    // prepare azimuthally averaged quantities for constructing 1D radial profiles
    rho1 = 0.0;
    vr1 = 0.0;
    etafc1 = 0.0;
    etacc1 = 0.0;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      rho1 += rho[l];
      vr1 += vr[l];
      etafc1 += etafc[l];
      etacc1 += etacc[l];
    }
    rho1 /= (real)ns;
    vr1 /= (real)ns;
    etafc1 /= (real)ns;
    etacc1 /= (real)ns;
    // 1D steady-state profiles
    tau1 = sqrt(3.0)*PEBBLECOAGULATION*pebbleflux/(32.0*PI*pow(Rmed[i],-1.5)*rho1);	// dominant Stokes number from 2 parameters + hydro quantities
    tau1 = sqrt(tau1)/(Rmed[i]*etacc1);
    psize1 = sqrt(ADIABIND)*tau1*rho1/(2.0*pebbulkdens);     // pebble size corresponding to the dominant Stokes from Epstein law
    vtcorr = Rmed[i]*OmegaFrame;				// steady-state drift velocity part, serves also as the initialization |--------->
    vkepl = pow (Rmed[i], -0.5);
    fac = 1.0/(tau1*tau1 + 1.0);
    etaterm = etafc1*vkepl;
    pvr1 = -2.0*tau1*fac*(etaterm - 0.5/tau1*vr1);	// initial radial velocity (face centered)
    etaterm = etacc1*vkepl;
    pvrtemp = -2.0*tau1*fac*(etaterm - 0.5/tau1*vr1);	// radial velocity in the cell centre (will be used for the flux estimate)
    pvt1 = -fac*(etaterm - 0.5*tau1*vr1);		// initial azim. velocity
    pvt1 += (vkepl - vtcorr);				// transform to corot. <---------|
    pdens1 = pebbleflux/(2.0*PI*Rmed[i]*fabs(pvrtemp));// having the drift velocity, get the density from flux conservation
    // fill HD pebble arrays with 1D profiles
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      psize[l] = psize1;
      tau[l] = tau1;
      pdens[l] = pdens1;
      pvr[l] = pvr1;
      pvt[l] = pvt1;
    }
    // bckp for boundary conditions
    PebDensInit[i] = pdens1;
    PebVradInit[i] = pvr1;
    PebVthetaInit[i] = pvt1 + Rmed[i]*OmegaFrame;	// transform back to intertial value
  }
  WriteDiskPolar (StokesNumber, 0);
  WriteDiskPolar (PebbleSize, 0);
  if (Merge && (CPU_Number > 1)) {
    if (!CPU_Master) return;
    for (i=1; i<CPU_Number; i++) {
      sprintf (command, "cd %s; cat gaspebblestokes0.dat.%05d >> gaspebblestokes0.dat",\
               OUTPUTDIR, i);
      system (command);
      sprintf (command, "cd %s; cat gaspebblesize0.dat.%05d >> gaspebblesize0.dat",\
               OUTPUTDIR, i);
      system (command);
    }
    sprintf (command, "cd %s; rm -f gaspebblestokes0.dat.0*", OUTPUTDIR);
    system (command);
    sprintf (command, "cd %s; rm -f gaspebblesize0.dat.0*", OUTPUTDIR);
    system (command);
  }
}

/** Backs up the initial state of the pebble disk
 * to impose damping boundary conditions later. */
void BckpFieldsForBC ()
{
  real *pdens, *pvr, *pvt;
  int nr, ns, i, j, l;
  pdens = PebbleDens->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  nr = PebbleVrad->Nrad;
  ns = PebbleVtheta->Nsec;
  for (i=0; i<nr; i++) {
    PebDensInit[i] = 0.0;
    PebVradInit[i] = 0.0;
    PebVthetaInit[i] = 0.0;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      PebDensInit[i] += pdens[l];
      PebVradInit[i] += pvr[l];
      PebVthetaInit[i] += pvt[l];
    }
    PebDensInit[i] /= (real)ns;
    PebVradInit[i] /= (real)ns;
    PebVthetaInit[i] /= (real)ns;
    PebVthetaInit[i] += Rmed[i]*OmegaFrame;
  }
}

/** Calculates tha gas rotation parameter eta */
void EtaPressureSupport (Vtheta)
PolarGrid *Vtheta;
{
  int i, j, l, lim, ljp, nr, ns;
  real vtcorr, vkinv, vtcc;
  real *vt, *etafc, *etacc;
  /* ----- */
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  vt = Vtheta->Field;
  etafc = EtaFaceCentered->Field;
  etacc = EtaCellCentered->Field;
#pragma omp parallel default(none) \
  shared (nr,ns,Rmed,OmegaFrame,vt,etacc,etafc) \
  private (i,j,l,vtcorr,vkinv,ljp,lim,vtcc)
  {
#pragma omp for
  for (i=0; i<nr; i++) {
    vtcorr = Rmed[i]*OmegaFrame;
    vkinv = sqrt(Rmed[i]);
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      vtcc = 0.5*(vt[l]+vt[ljp]);
      etacc[l] = 1.0 - (vtcc + vtcorr)*vkinv;
    }
  }
#pragma omp for
  for (i=1; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lim = l - ns;
      etafc[l] = 0.5*(etacc[l]+etacc[lim]);
    }
  }
  } // #end of omp parallel section
  for (j=0; j<ns; j++) {	// adopt nearest values when missing
    etafc[j] = etafc[j+ns];
  }
}

/** Calculates the local Stokes number using the 
 * dominant pebble size, local gas density and parametric
 * pebble material density. */
void PebbleStokesNumbers (Rho)
PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real *rho, *tau, *psize;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  tau = StokesNumber->Field;
  psize = PebbleSize->Field;
#pragma omp parallel for default(none) \
  shared (nr,ns,tau,pebbulkdens,psize,rho,SQRT_ADIABIND_INV) \
  private (i,j,l)
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      tau[l] = 2.0*pebbulkdens*psize[l]/rho[l]*SQRT_ADIABIND_INV;	// given the pebble size, recalculate tau from local density (Epstein)
    }
  }
}

/** Trapezoidal rule to integrate the column mass
 * of pebbles located within the accretion radius.
 * Adopted from Numerical Recipes. */
real Trapzd (n, a, b, H2)
int n;
real a, b, H2;
{
  static real integral;
  real iter, spacing, x, sum;
  int j;
  if (n==1) {
    integral = 0.5*(b-a)*(exp(-0.5*a*a/H2) + exp(-0.5*b*b/H2));
    return integral;
  } else {
    iter = pow(2.0,(real)(n-2));
    spacing = (b-a)/iter;
    x = a + 0.5*spacing;
    sum = 0.0;
    for (j=1; j<=(int)iter; j++) {
      sum += exp(-0.5*x*x/H2);
      x += spacing;
    }
    integral = 0.5*(integral + (b-a)*sum/iter);
    return integral;
  }
}

/** Primitive algorithm to find the mass in a
 * vertical column of pebbles overlapping
 * the accretion radius. */
real IntegrateColumnMass (a, b, Hpeb)
real a, b, Hpeb;
{
  int n;
  real integral, iold, H2;
  H2 = Hpeb*Hpeb;
  for (n=1; n<=TRAPEZMAX; n++) {
    integral = Trapzd (n,a,b,H2);
    if (n > 5) {	// avoid early convergence
      if (fabs(integral-iold) < TRAPEZEPS*fabs(iold) || (integral == 0.0 && iold == 0.0)) return integral;
    }
    iold = integral;
  }
  masterprint ("Warning! Integration of the column mass did not converge. Terminating now...\n");
  prs_exit(1);
  return 0.0;	// will never get here
}

/** Finds the amount of pebbles to be transfered
 * from the pebble disk onto the planets */
void AccretePebblesOntoPlanets (sys, Rho, Energy, Vtheta, dt)
PlanetarySystem *sys;
PolarGrid *Rho, *Energy, *Vtheta;
real dt;
{
  real *pdens, *abs, *ord;
  real *pvr, *pvt, *tau, *cs, *vt, *macc;
  int nr,ns,k,i,j,l,imin,imax,ncell,tempint;
  real dMplanet, dPXplanet, dPYplanet;
  real Xplanet, Yplanet, Zplanet, VXplanet, VYplanet, VZplanet, Mplanet, PXplanet, PYplanet, Rplanet, Rplanet3D, Omegaplanet;
  real vtcorr, xcell, ycell, vrcell, vtcell, vxcell, vycell, flow, flowx, flowy, vrel;
  real ifrac, frac, Mt, Rbondi, tbondi, fac, Reff=0.0, Reff2, Rhill;
  real R2, H, Hpeb, zmin, zmax, dM, mcell, voldenspeb, temp;
  real ang, plradius, plrho, dE, L, tdelay, taper;
  real tau_mean, dM_upper, Hpeb_mean, pdens_mean, veloc, dM_model, limiter;
  /* ----- */
  pdens = PebbleDens->Field;
  nr = PebbleDens->Nrad;
  ns = PebbleDens->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  tau = StokesNumber->Field;
  cs = SoundSpeed->Field;
  macc = AccretedMass->Field;
  vt = Vtheta->Field;
  if (AccretHeating) heatsrc_max = sys->nb;
  mpi_make1Dprofile (vt, vt1D);
  for (k=0; k<sys->nb; k++) {
    dMplanet = dPXplanet = dPYplanet = 0.0;
    Xplanet = sys->x[k];
    Yplanet = sys->y[k];
    Zplanet = sys->z[k];
    VXplanet = sys->vx[k];
    VYplanet = sys->vy[k];
    VZplanet = sys->vz[k];
    Mplanet = sys->mass[k];
    PXplanet = VXplanet*Mplanet;
    PYplanet = VYplanet*Mplanet;
    Rplanet = Xplanet*Xplanet + Yplanet*Yplanet;
    Rplanet3D = sqrt(Rplanet + Zplanet*Zplanet);
    Rplanet = sqrt(Rplanet);
    Rhill = Rplanet3D*pow(1.0/3.0*Mplanet, 1.0/3.0);
    Omegaplanet = pow(Rplanet, -1.5);                                   // transition mass calculation on global grid --->
    ifrac = GetGlobalIFrac (Rplanet);                                   // get radial index w.respect to the global grid
    frac = ifrac-floor(ifrac);
    i = (int)ifrac;
    flow = vt1D[i]*(1.0-frac) + vt1D[i+1]*frac + Rplanet*OmegaFrame;    // flow of gas at planet's orbit (v_rad negleted, azim.average of v_theta, from corot. to inertial)
    flowx = -flow*Yplanet/Rplanet;                                      // conversion from theta to (x,y)
    flowy =  flow*Xplanet/Rplanet;
    vrel = sqrt((VXplanet-flowx)*(VXplanet-flowx) + (VYplanet-flowy)*(VYplanet-flowy) + VZplanet*VZplanet);     // headwind experienced by planet
    Mt = sqrt(1.0/3.0)*pow(vrel,3.0)/Omegaplanet;                       // <---
    imax = Max_or_active-1;						// set range of cells on current CPU that should be checked for accretion --->
    imin = Zero_or_active;
    while ( Rinf[imax] > (Rplanet+Rhill) && imax >= imin) imax--;	// if at the end (imin>imax), it means that current CPU does not contain the critical sphere
    while ( Rsup[imin] < (Rplanet-Rhill) && imin <= imax) imin++;	// and the loop below won't do anything <---
    // STEP 1 - find average tau through the annulus spanned by the Hill sphere
    ncell = 0;
    tau_mean = 0.0;
    for (i=imin; i<=imax; i++) {
      for (j=0; j<ns; j++) {
	l = j + i*ns;
        tau_mean += tau[l];
	ncell++;
      }
    }
    MPI_Allreduce(&tau_mean, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tau_mean = temp;
    MPI_Allreduce(&ncell, &tempint, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ncell = tempint;
    tau_mean /= (real)ncell;
    // STEP 2 - using tau calculated above, find the effective accretion radius
    if (Mplanet < Mt) {			// Bondi regime
      Rbondi = Mplanet/vrel/vrel;	// insert vrel from the headwind experienced by planet
      tbondi = Rbondi/vrel;
      fac = sqrt(tau_mean*pow(Rplanet, 1.5)/tbondi);	// tau/Omega_kepl = t_stopping
      Reff = Rbondi*fac;
    } else {				// Hill regime
      fac = pow(tau_mean/0.1, 1.0/3.0);
      Reff = Rhill*fac;
    }
    if (Reff > Rhill) Reff = Rhill; 	// can't get bigger than Hill (just to be safe)
    Reff2 = Reff*Reff;
    // STEP 3 - find the cells below Reff, find the UPPER LIMIT for ACCRETED MASS
    dM_upper = 0.0;
    Hpeb_mean = 0.0;
    pdens_mean = 0.0;
    ncell = 0;
    for (i=imin; i<=imax; i++) {
      vtcorr = Rmed[i]*OmegaFrame;
      for (j=0; j<ns; j++) {
	l = j + i*ns;
        xcell = abs[l];
        ycell = ord[l];
        R2 = (xcell-Xplanet)*(xcell-Xplanet) + (ycell-Yplanet)*(ycell-Yplanet); // midplane separation cell vs planet
	if (R2 <= Reff2) {
  	  H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
          Hpeb = H*sqrt(PEBBLEALPHA/tau[l]);
	  zmin = Zplanet - sqrt(Reff2 - R2);	// intersections between the vertical column of pebbles and the hill sphere
	  zmax = Zplanet + sqrt(Reff2 - R2);
	  mcell = Surf[i]*pdens[l];
 	  voldenspeb = pdens[l]/Hpeb*SQRT2PI_INV;
          dM = voldenspeb*Surf[i]*IntegrateColumnMass (zmin,zmax,Hpeb);
	  if (dM > mcell) dM = 0.99999999*mcell;
	  macc[l] = dM;				// this is the UPPER LIMIT for the mass transported from a specific cell onto the planet
  	  dM_upper += dM;
	  Hpeb_mean += Hpeb;
	  pdens_mean += pdens[l];
	  ncell++;
        }
      }
    }
    MPI_Allreduce(&dM_upper, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dM_upper = temp;
    MPI_Allreduce(&Hpeb_mean, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Hpeb_mean = temp;
    MPI_Allreduce(&pdens_mean, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    pdens_mean = temp;
    MPI_Allreduce(&ncell, &tempint, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ncell = tempint;
    Hpeb_mean /= (real)ncell;
    pdens_mean /= (real)ncell;
    // STEP 4 - using average values, find the MODEL ACCRETION RATE to reduce the upper
    // limit from STEP 3 (formulae from e.g. Morbidelli et al. 2015)
    if (Mplanet < Mt) {                         // Bondi or Hill regime
      veloc = vrel;				// headwind as the relative velocity
    } else {
      veloc = Reff*Omegaplanet;			// shear velocity
    }
    if (Hpeb_mean < Reff) {        		// 2D or 3D regime
      dM_model = 2.0*Reff*veloc*pdens_mean;
    } else {
      dM_model = PI*Reff*Reff*veloc*pdens_mean/Hpeb_mean*SQRT2PI_INV;
    }
    dM_model *= dt;				// !!! MUST convert the accretion rate to mass change
    if (dM_upper > dM_model) {       		// limit the accretion rate not to exceed the model rate
      limiter = dM_model/dM_upper;
    } else {
      limiter = 1.0;
    }
    // STEP 5 - apply the limiter and finish the accretion
    for (i=imin; i<=imax; i++) {
      vtcorr = Rmed[i]*OmegaFrame;
      for (j=0; j<ns; j++) {
        l = j + i*ns;
        if (macc[l] > 0.0) {			// cells below Rhill have non-zero 'macc'
          xcell = abs[l];
          ycell = ord[l];
          vrcell = pvr[l];
          vtcell = pvt[l] + vtcorr;             // from corot. to inertial
          vxcell = (vrcell*xcell - vtcell*ycell)*InvRmed[i];
          vycell = (vrcell*ycell + vtcell*xcell)*InvRmed[i];
          dM = macc[l]*limiter;                 // limiter applied here
          macc[l] = 0.0;                        // !!! important so that the if-condition works properly next time
          mcell = Surf[i]*pdens[l];
          pdens[l] = (mcell - dM)/Surf[i];      // new surface density
          if (pdens[l] < 0.0) pdens[l] = 1.e-20;  // density floor
          dMplanet += dM;
          dPXplanet += dM*vxcell;
          dPYplanet += dM*vycell;
        }
      }
    }
    MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dMplanet = temp;
    MPI_Allreduce (&dPXplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dPXplanet = temp;
    MPI_Allreduce (&dPYplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dPYplanet = temp;
    // OPTIONAL STEP 6 - account for the accretion heating
    if (AccretHeating) {
      heatsrc_index[k] = -1;	   
      if (Rplanet >= Rinf[0] && Rplanet <= Rsup[nr-1]) {
        tdelay = 4.0*DT*(real)HEATINGDELAY;
        taper = PhysicalTime/tdelay;
        if (taper > 1.0) taper = 1.0;
        plrho = PLANETARYDENSITY/RHO2CGS;
        plradius = pow (3.0*Mplanet/(4.0*PI*plrho), 1.0/3.0);
        L = Mplanet*dMplanet/dt/plradius; // see Benitez et al. (2015); G=1; and we work in terms of LUMINOSITY
        ang = atan2(Yplanet,Xplanet);
        if (ang < 0.0) ang += 2.0*PI;
        i = 0;
        while (Rsup[i] < Rplanet) i++;
        j = floor((real)ns/2.0/PI*ang + 0.5);
        if (j==ns) j=0;
        l = j + i*ns;
	heatsrc_index[k] = l;
	dE = L*taper/Surf[i];	// final dimension: W/m^2 -> average power density in cells closest to the planet
        heatsrc[k] = dE;
      }
    }
    // STEP 7 - Update the planet
    PXplanet += dPXplanet;
    PYplanet += dPYplanet;
    Mplanet  += dMplanet;
    sys->mass[k] = Mplanet;
    sys->vx[k] = PXplanet/Mplanet;
    sys->vy[k] = PYplanet/Mplanet;
  }
}

/** Corrects the azimuthal flow velocity
 * of pebbles to keep up with the frame rotation. */
void CorrectPebblesVtheta (domega)
real domega;
{
  CorrectVtheta (PebbleVtheta, domega);
}

/** Calculates the source terms acting on pebbles */
void SourceTermsPebbles (Vrad,Vtheta,dt)
PolarGrid *Vrad, *Vtheta;
real dt;
{
  real *pvr, *pvt, *tau, *paccr, *pacct, *vr, *vt, *pdens, *dragr, *dragt;
  real *pagr, *pagt;
  real Omega, vt2, drag;
  int nr, ns, i, j, l, lim, ljm, ljp;
  pdens = PebbleDens->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  tau = StokesNumber->Field;
  paccr = PebbleAccelrad->Field;
  pacct = PebbleAcceltheta->Field;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  dragr = DragForceRad->Field;
  dragt = DragForceTheta->Field;
  pagr = PebbleGravAccelRad->Field;
  pagt = PebbleGravAccelTheta->Field;
  DampPebbles (PebbleDens, PebbleVrad, PebbleVtheta, dt);
#pragma omp parallel private(i,j,l,lim,ljm,ljp,vt2,Omega,drag)
  {
#pragma omp for
    for (i = 1; i < nr; i++) {
      Omega = 1.0/OmegaInv[i];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        vt2 = pvt[l]+pvt[ljp]+pvt[lim]+pvt[ljp-ns];
        vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
        vt2 = vt2*vt2;
        paccr[l] = vt2*InvRinf[i]+0.5*(pagr[l]+pagr[lim]);
        if (BackReaction) {
          drag = Omega*(pvr[l] - vr[l])/tau[l];
	  dragr[l] = 0.5*(pdens[l]+pdens[lim])*drag;
	}
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      Omega = 1.0/OmegaInv[i];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
        lim = l-ns;
        ljm = l-1;
        if (j == 0) ljm = i*ns+ns-1;
        ljp = l+1;
        if (j == ns-1) ljp = i*ns;
        pacct[l] = 0.5*(pagt[l]+pagt[ljm]);
	if (BackReaction) {
	  drag = Omega*(pvt[l] - vt[l])/tau[l];
	  dragt[l] = 0.5*(pdens[l]+pdens[ljm])*drag;
	}
      }
    }
  }
}

/** Applies a semi-implicit method to evolve the pebbles
 * dynamically. See Appendix C in Chrenko et al. (2017) */
void SubStep1Pebbles (Vrad,Vtheta,dt)
PolarGrid *Vrad, *Vtheta;
real dt;
{
  real *pvr, *pvt, *tau, *paccr, *pacct, *accr, *acct, *vr, *vt;
  real Omega, efac;
  int nr, ns, i, j, l;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  tau = StokesNumber->Field;
  paccr = PebbleAccelrad->Field;
  pacct = PebbleAcceltheta->Field;
  accr = GasAccelrad->Field;
  acct = GasAcceltheta->Field;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  vr = Vrad->Field;
  vt = Vtheta->Field;
#pragma omp parallel private(j,l,Omega,efac)
  {
#pragma omp for
    for (i = 0; i < nr; i++) {
      Omega = 1.0/OmegaInv[i];
      for (j = 0; j < ns; j++) {
        l = j+i*ns;
	efac = exp(-dt*Omega/tau[l]);
        if (i > 0) {
	  pvr[l] = pvr[l]*efac + accr[l]*dt + \
		 (vr[l] + (paccr[l] - accr[l])*tau[l]/Omega)*(1.0 - efac);
	}
        pvt[l] = pvt[l]*efac + acct[l]*dt + \
		 (vt[l] + (pacct[l] - acct[l])*tau[l]/Omega)*(1.0 - efac);
      }
    }
  }
}

/** Applies the particle diffusion term acting on pebbles.
 * See Eq. (35) in Chrenko et al. (2017) */
void ParticleDiffusion (Rho)
PolarGrid *Rho;
{
  int i,j,l,lim,ljm,nr,ns;
  real corr, dxtheta, invdxtheta;
  real *rho, *pdens, *pvr, *pvt;
  rho = Rho->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  pdens = PebbleDens->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
#pragma omp parallel default(none) \
  shared (nr,ns,Rinf,SCHMIDTNUMBER,rho,pdens,InvDiffRmed,pvr,Rmed,pvt) \
  private (i,j,l,lim,ljm,corr,dxtheta,invdxtheta)
  {
#pragma omp for
  for (i=1; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lim = l - ns;
      corr = -FViscosity(Rinf[i])/SCHMIDTNUMBER;
      corr *= (rho[l] + rho[lim])/(pdens[l] + pdens[lim]);
      corr *= (pdens[l]/rho[l] - pdens[lim]/rho[lim])*InvDiffRmed[i];
      pvr[l] += corr;
    }
  }
#pragma omp for
  for (i=0; i<nr; i++) {
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      corr = -FViscosity(Rmed[i])/SCHMIDTNUMBER;
      corr *= (rho[l] + rho[ljm])/(pdens[l] + pdens[ljm]);
      corr *= (pdens[l]/rho[l] - pdens[ljm]/rho[ljm])*invdxtheta;
      pvt[l] += corr;
    }
  }
  } // #end of omp parallel section
}

/** Calls the transport routines and applies the boundary conditions
 * for pebbles */
void EvolvePebbleDisk (dt)
real dt;
{
  DampPebbles (PebbleDens, PebbleVrad, PebbleVtheta, dt);
  TransportPebbles (PebbleDens, PebbleVrad, PebbleVtheta, dt);
  DampPebbles (PebbleDens, PebbleVrad, PebbleVtheta, dt);
}

/** Synchronises pebble fluid hydrodynamic quantities
 * among the overlapping grid zones. */
void SynchronizePebbleDisc ()
{
  int nr;
  real *pdens, *pvr, *pvt;
  nr = PebbleDens->Nrad;
  pdens = PebbleDens->Field;
  pvr = PebbleVrad->Field;
  pvt = PebbleVtheta->Field;
  if (CPU_Number > 1) {
    SynchronizeOverlapFields (pdens, nr, CPUOVERLAP);	/* we use the function from RadiativeDiffusion.c */
    SynchronizeOverlapFields (pvr, nr, CPUOVERLAP);
    SynchronizeOverlapFields (pvt, nr, CPUOVERLAP);
  }
}

/** Restricts the time step using the CFL
 * condition for the pebble fluid. */
void CriticalCharTime (Vrad, Vtheta)
PolarGrid *Vrad, *Vtheta;
{
  real dxrad, dxtheta, pvtmed=0.0, invdt1, invdt2, invdt3, invdt4, Vresidual, temp;
  real *pvr, *pvt, *vr, *vt;
  int i, j, l, ns;
  pvr = PebbleVrad->Field;
  ns = PebbleVrad->Nsec;
  pvt = PebbleVtheta->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  invdtpeb_sq = -1.e30;
  for (i=One_or_active; i<Max_or_active; i++) {
    dxrad = Rsup[i] - Rinf[i];
    dxtheta = Rmed[i]*2.0*PI/(real)ns;
    if (FastTransport) {
      pvtmed = 0.0;
      for (j=0; j<ns; j++) {
        l = j + i*ns;
	pvtmed += pvt[l];
      }
      pvtmed /= (real)ns;
    }
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      invdt1 = fabs(pvr[l])/dxrad;
      if (FastTransport) {
        Vresidual = pvt[l] - pvtmed;
      } else {
        Vresidual = pvt[l];
      }
      invdt2 = fabs(Vresidual)/dxtheta;
      invdt3 = fabs(pvr[l] - vr[l])/dxrad;	// Rosotti et al. 2011
      invdt4 = fabs(pvt[l] - vt[l])/dxtheta;
      temp = invdt1*invdt1 + invdt2*invdt2 + invdt3*invdt3 + invdt4*invdt4;
      invdtpeb_sq = max2(invdtpeb_sq,temp);
    }
  }
}

/** Outputs the pebble fluid arrays. */
void WritePebbles (index)
int index;
{
  WriteDiskPolar (PebbleDens, index);
  WriteDiskPolar (PebbleVrad, index);
  WriteDiskPolar (PebbleVtheta, index);
  if (Write_Eta) WriteDiskPolar (EtaCellCentered, index);
}

/** Safety check for negative pebble densities. */
boolean DetectCrashPebbles ()
{
  boolean result;
  result = DetectCrash (PebbleDens);
  return result;
}

/** Writes filtering factors if needed. */
/*void UpdateFilteringFile (sys)
PlanetarySystem *sys;
{
  int k;
  char name[256];
  FILE *output;
  if (!CPU_Master) return;
  sprintf (name, "%spebble_filtering.dat", OUTPUTDIR);   
  output = fopenp (name, "a");
  for (k=0; k<sys->nb; k++) {
    fprintf (output, "%d\t%.12g\t%.12g\n", k, PhysicalTime, filtering[k]);
  }
  fclose (output);
}*/

/** Mass accretion onto planets is provided
 * by a parametric prescription using a given mass
 * doubling time. */
void ParametricAccretion (sys, dt)
PlanetarySystem *sys;
real dt;
{
  int k, i, j, l, nr, ns;
  real Mplanet, Xplanet, Yplanet, dMplanet, Rplanet, tdelay, taper;
  real plrho, plradius, L, ang, dE;
  static real doubling=0.0;
  if (doubling == 0.0) doubling = PARAMETRICACCRETION*1000.0*2.0*PI;	// convert the doubling time from kyr to code units
  nr = SoundSpeed->Nrad;
  ns = SoundSpeed->Nsec;
  if (AccretHeating) heatsrc_max = sys->nb;
  for (k=0; k<sys->nb; k++) {
    Mplanet = sys->mass[k];
    Xplanet = sys->x[k];
    Yplanet = sys->y[k];
    dMplanet = Mplanet*dt/doubling;
    Rplanet = sqrt(Xplanet*Xplanet + Yplanet*Yplanet);
    if (AccretHeating) {
      heatsrc_index[k] = -1;
      if (Rplanet >= Rinf[0] && Rplanet <= Rsup[nr-1]) {
        tdelay = 4.0*DT*(real)HEATINGDELAY;
        taper = PhysicalTime/tdelay;
        if (taper > 1.0) taper = 1.0;
        plrho = PLANETARYDENSITY/RHO2CGS;
        plradius = pow (3.0*Mplanet/(4.0*PI*plrho), 1.0/3.0);
        L = Mplanet*Mplanet/doubling/plradius;
        ang = atan2(Yplanet,Xplanet);
        if (ang < 0.0) ang += 2.0*PI;
        i = 0;
        while (Rsup[i] < Rplanet) i++;
        j = floor((real)ns/2.0/PI*ang + 0.5);
        if (j==ns) j=0;
        l = j + i*ns;
        heatsrc_index[k] = l;
        dE = L*taper/Surf[i];   // final dimension: W/m^2 -> average power density in cells closest to the planet
        heatsrc[k] = dE;
      }
    }
    Mplanet += dMplanet;
    sys->mass[k] = Mplanet;
  }
}
