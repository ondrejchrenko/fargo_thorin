/**
 * @file EnergySources.c
 *
 * @brief Subroutines related to the heating/cooling source terms,
 * numerical solver for the energy equation and radiative diffusion.
 *
 * @author Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>
 *
 * @details Calculates the individual energy source terms
 * according to Chrenko et al. (2017). Then solves the energy
 * equation in a linearised implicit form using the successive
 * over-relaxation (SOR) method (see Appendix A in Chrenko et al. 2017).
 *
 * @section     LICENSE
 * Copyright (c) 2017 Ondřej Chrenko. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

/* Control parameters of the SOR-algorithm */
#define SORMAXITERS 1000	/* Max no. of iterations */
#define SOREPS 1.0e-8		/* Precision when comparing the improvement of residues */
/* To avoid convergence issues due to double-precision limitations: */
#define SORMIN 1.0e-15		/* Will be compared with the normalised residue
				   (SOR will finish when close to the DP limit)
				   and the temperature improvements
				   (SOR will finish if the temperature does not improve
				   within the DP limit) */

/* Opacity table according to Bell & Lin (1994), see also Keith & Wardle (2014) */
static real kappa0[7] = {2.0e-4, 2.0e16, 0.1, 2.0e81, 1.0e-8, 1.0e-36, 1.5e20};
static real a[7] = {0.0, 0.0, 0.0, 1.0, 2.0/3.0, 1.0/3.0, 1.0};
static real b[7] = {2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -5.0/2.0};

static PolarGrid *GradTemperRad, *GradTemperTheta, *GradTemperMagnitude;
static PolarGrid *DiffCoefCentered, *DiffCoefIfaceRad, *DiffCoefIfaceTheta;
static PolarGrid *Opacity, *VolumeDensity;
static PolarGrid *Flaring, *Qirradiation;

static PolarGrid *DiscretizationCoefA, *DiscretizationCoefB;
static PolarGrid *MatrixNexttoTemperl, *MatrixNexttoTemperlip, *MatrixNexttoTemperlim;
static PolarGrid *MatrixNexttoTemperljp, *MatrixNexttoTemperljm, *RightHandSide;

static real CV;			// specif. heat capacity
static real omegabest, domega;	// optimum relaxation parameter for SOR and its small increment
static int Niterbest;		// minimum number of iterations reached when minimizing the relaxation parameter
static int jchess1st, jchess2nd;

/** Initialises the polar arrays associated with the heating/cooling
 * processes */
void InitRadiatDiffusionFields ()
{ 
  /* all polar grids are within the scope of this file, apart from Qminus
     which is a global variable (see global.h) */
  Qminus		= CreatePolarGrid (NRAD, NSEC, "qminus");
  GradTemperRad		= CreatePolarGrid (NRAD, NSEC, "gradtemperr");
  GradTemperTheta	= CreatePolarGrid (NRAD, NSEC, "gradtempert");
  GradTemperMagnitude   = CreatePolarGrid (NRAD, NSEC, "gradtemperm");
  DiffCoefCentered	= CreatePolarGrid (NRAD, NSEC, "diffcoefc");
  DiffCoefIfaceRad      = CreatePolarGrid (NRAD, NSEC, "diffcoefr");
  DiffCoefIfaceTheta    = CreatePolarGrid (NRAD, NSEC, "diffcoeft");
  Opacity		= CreatePolarGrid (NRAD, NSEC, "opacity");
  VolumeDensity		= CreatePolarGrid (NRAD, NSEC, "volumedensity");
  Flaring               = CreatePolarGrid (NRAD, NSEC, "flaring");
  Qirradiation		= CreatePolarGrid (NRAD, NSEC, "qirradiation");
  DiscretizationCoefA   = CreatePolarGrid (NRAD, NSEC, "discretcoefa");
  DiscretizationCoefB   = CreatePolarGrid (NRAD, NSEC, "discretcoefb");
  MatrixNexttoTemperl	= CreatePolarGrid (NRAD, NSEC, "matrixl");
  MatrixNexttoTemperlip = CreatePolarGrid (NRAD, NSEC, "matrixlip");
  MatrixNexttoTemperlim = CreatePolarGrid (NRAD, NSEC, "matrixlim");
  MatrixNexttoTemperljp = CreatePolarGrid (NRAD, NSEC, "matrixljp");
  MatrixNexttoTemperljm = CreatePolarGrid (NRAD, NSEC, "matrixljm");
  RightHandSide         = CreatePolarGrid (NRAD, NSEC, "righthandside");
  if (Write_Qbalance) Qbalance = CreatePolarGrid (NRAD, NSEC, "qbalance");
  CV = GASCONST/(MOLWEIGHT*(ADIABIND - 1.0));
}

/** Estimate of the heat loss due to radiation
 * escape in the vertical direction with respect to the midplane.
 * See Eq. (9) in Chrenko et al. (2017). */
void CalculateQminus (Rho)
PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real *rho, *opacity, *qminus, *temper;
  real tau, taueff;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  opacity = Opacity->Field;
  temper = Temperature->Field;
  qminus = Qminus->Field;
#pragma omp parallel for default(none) shared(nr,ns,opacity,rho,temper,qminus,OPACITYDROP) \
  private(i,j,l,tau,taueff)
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      tau = OPACITYDROP*0.5*opacity[l]*rho[l];      /* the surface density is used here (see d'angelo 2003) not the volume density */
      taueff = EffectiveOpticalDepth (tau);
      qminus[l] = 2.0*STEFANBOLTZMANN*pow(temper[l],3.0)/taueff;	/* !!! this is NOT pure qrad (which is propto T^4),
									   but it's a useful form for the implicit inversion of the energy eq. */
    }
  }
}

/** Calculates the stellar irradiation source term.
 * See Eq. (13) in Chrenko et al. (2017). */
void CalculateQirr (Rho)
PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real *rho, *opacity, *qirr, *flaring;
  real Teff4, Rstar2, tau, taueff, Tirr4, r2inv;
  static real Tirr4fac=0.0;
  /* ----- */
  if (Tirr4fac==0.0) {
    Teff4 = EFFECTIVETEMPERATURE/T2SI;
    Teff4 = pow(Teff4, 4.0);
    Rstar2 = 6.957e8/AU_SI;     // 1 solar radius in code units
    Rstar2 = pow(STELLARRADIUS*Rstar2,2.0);
    Tirr4fac = (1.0 - DISCALBEDO)*Teff4*Rstar2;
  }
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  opacity = Opacity->Field;
  qirr = Qirradiation->Field;
  flaring = Flaring->Field;
#pragma omp parallel for default(none) \
  shared(nr,ns,rho,opacity,Rmed2,qirr,Tirr4fac,flaring,OPACITYDROP) \
  private(i,j,l,tau,taueff,Tirr4,r2inv)
  for (i=0; i<nr; i++) {
    r2inv = 1.0/Rmed2[i];
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      if (flaring[l] > 0.0) {
        tau = OPACITYDROP*0.5*opacity[l]*rho[l];
        taueff = EffectiveOpticalDepth (tau);
        Tirr4 = Tirr4fac*r2inv*flaring[l];
        qirr[l] = 2.0*STEFANBOLTZMANN*Tirr4/taueff;
      } else {
        qirr[l] = 0.0;
      }
    }
  }
}

/** Calculates the sine of the grazing angle by reconstructing
 * the surface from the pressure scale height.
 * See Eq. (15) in Chrenko et al. (2017). */
void CalculateFlaring ()
{
  int i,j,l,lim,lip,nr,ns;
  real csin, csout, Hin, Hout, dHdr, H;
  real *flaring, *cs;
  static real Rstar=0.0;
  if (Rstar==0.0) Rstar = STELLARRADIUS*6.957e8/AU_SI;
  flaring = Flaring->Field;
  cs = SoundSpeed->Field;
  nr = SoundSpeed->Nrad;
  ns = SoundSpeed->Nsec;
#pragma omp parallel for default(none) \
  shared(nr,ns,cs,Rmed,Rinf,Rsup,InvDiffRsup,Rstar, \
         flaring,SQRT_ADIABIND_INV,InvRmed,OmegaInv) \
  private(i,j,l,lip,lim,csin,csout,Hin,Hout,dHdr,H)
  for (i=1; i<nr-1; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lip = l + ns;
      lim = l - ns;
      csin = 0.5*(cs[l] + cs[lim]);
      csout = 0.5*(cs[l] + cs[lip]);
      Hin = csin*pow(Rinf[i],1.5)*SQRT_ADIABIND_INV;
      Hout = csout*pow(Rsup[i],1.5)*SQRT_ADIABIND_INV;
      dHdr = (Hout - Hin)*InvDiffRsup[i];
      H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
      flaring[l] = atan(dHdr) - atan((H-0.4*Rstar)*InvRmed[i]);    // Baillié & Charnoz (2014)
      if (flaring[l] >= 0.0) {
        flaring[l] = sin(flaring[l]);
      } else {							// surface parts tilted away from the incoming starlight
        flaring[l] = 0.0;
      }
      if (i==1) {                                               // find values for i=0 and i=nr-1
        H = cs[lim]*OmegaInv[i-1]*SQRT_ADIABIND_INV;
        flaring[lim] = atan(dHdr) - atan((H-0.4*Rstar)/Rmed[i-1]);      // we adopt the closest derivative
        if (flaring[lim] >= 0.0) {
          flaring[lim] = sin(flaring[lim]);
        } else {
          flaring[lim] = 0.0;
        }
      }
      if (i==nr-2) {
        H = cs[lip]*OmegaInv[i+1]*SQRT_ADIABIND_INV;
        flaring[lip] = atan(dHdr) - atan((H-0.4*Rstar)/Rmed[i+1]);      // we adopt the closest derivative
        if (flaring[lip] >= 0.0) {
          flaring[lip] = sin(flaring[lip]);
        } else {
          flaring[lip] = 0.0;
        }
      }
    }
  }
}

/** The main numerical solver of the energy equation. */
void ImplicitRadiativeDiffusion (Rho, EnergyInt, EnergyNew, dt)
PolarGrid *Rho, *EnergyInt, *EnergyNew;
real dt;
{
  real *A, *B, *Ml, *Mlip, *Mlim, *Mljp, *Mljm, *RHS;
  real *temper, *rho, *energynew, *cs, *Dr, *Dt;
  real *qplus, *qminus, *divergence, *qirr;
  int nr, ns, i, j, l, lip, ljp, Niter, k;
  real dxtheta, invdxtheta2, H, afac, bfac;
  static boolean first=YES;
  static int Niterlast;
  static real omega;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  energynew = EnergyNew->Field;
  temper = Temperature->Field;
  cs = SoundSpeed->Field;
  Dr = DiffCoefIfaceRad->Field;
  Dt = DiffCoefIfaceTheta->Field;
  A  = DiscretizationCoefA->Field;
  B  = DiscretizationCoefB->Field;
  Ml = MatrixNexttoTemperl->Field;
  Mlip = MatrixNexttoTemperlip->Field;
  Mlim = MatrixNexttoTemperlim->Field;
  Mljp = MatrixNexttoTemperljp->Field;
  Mljm = MatrixNexttoTemperljm->Field;
  RHS  = RightHandSide->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  qminus = Qminus->Field;
  qirr = Qirradiation->Field;
  ComputeTemperatureField (Rho, EnergyInt);	/* use the intermediate energy (updated in SubStep2) to get T and cs */
  ComputeSoundSpeed (Rho, EnergyInt);
  MidplaneVolumeDensity (Rho);			/* use new cs to convert the surface density to the volume density */
  OpacityProfile ();				/* use new temperature and vol.dens to calc. the opacity */
  CalculateQminus (Rho);			/* estimate the vertical cooling term (z-direction radiation escape) */
  if (StellarIrradiation) {			/* for stellar irradiated discs: */
    CalculateFlaring ();			/* - check which parts are exposed to the incoming radiation */
    CalculateQirr (Rho);			/* - calculate the heating contribution */
  }
  DiffusionCoefs ();				/* calculate the diffusion coefficients */
  /* To update a centered temperature value using SOR, we need the value
   * to be surrounded by 4 average diffusion coefs (one per each interface)
   * and the temperature field in adjacent cells must be known.
   * Update is not possible at the innermost and outermost ring on each CPU.
   * For this reason and for MPI-parallelization of the SOR method, we
   * restrict the computation to the active mesh of each processor
   * (thus everything is fine on middle CPUs, the inner/outer ring problem
   * remains only on the inner/outer CPU).
   * After SOR, the temperature in overlapping zones will be synchronized,
   * leaving old values of T only in the inner ring of inner CPU and in the outer ring
   * of outer CPU. T is updated in these rings by adopting the neighbouring value. */
#pragma omp parallel default(none) \
  shared(nr, ns, One_or_active, MaxMO_or_active, Rmed, Rinf, \
         InvDiffRmed, InvDiffRsup, A, B, Dr, Dt, cs, dt, Ml, Mlip, \
         Mlim, Mljp, Mljm, RHS, CV, divergence, qminus, temper, qplus, ADIABIND, \
         OmegaInv, SQRT_ADIABIND_INV, rho, StellarIrradiation, qirr, AccretHeating, \
	 heatsrc_max, heatsrc_index, heatsrc) \
  private(i, j, l, k, lip, dxtheta, invdxtheta2, ljp, H, afac, bfac)
  {
#pragma omp for
  for (i=1; i<nr; i++) {		// parts of matrices for SOR (no active mesh restriction so far!)
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta2 = 1.0/dxtheta/dxtheta;	// this is 1/(dtheta*Rmed[i])^2
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      A[l] = Dr[l]*Rinf[i]*InvDiffRmed[i];
      H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
      B[l] = Dt[l]*invdxtheta2*2.0*H;
    }
  }
#pragma omp for
  for (i=One_or_active; i<MaxMO_or_active; i++) {	// final matrices for SOR, active mesh restriction
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lip = l + ns;
      ljp = l + 1;
      if (j == ns-1) ljp = i*ns;
      H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
      bfac = dt/(rho[l]*CV);
      afac = bfac*InvDiffRsup[i]*2.0*H/Rmed[i];
      Ml[l] = 1.0 + divergence[l]*(ADIABIND-1.0)*dt + 4.0*bfac*qminus[l] + afac*(A[lip] + A[l]) + bfac*(B[ljp] + B[l]);
      Mlip[l] = - afac*A[lip];		// because of A[lip], the previous cycle is not restricted to active mesh
      Mlim[l] = - afac*A[l];
      Mljp[l] = - bfac*B[ljp];
      Mljm[l] = - bfac*B[l];
      RHS[l]  = temper[l] + bfac*(qplus[l] + 3.0*qminus[l]*temper[l]);
      if (StellarIrradiation) RHS[l] += bfac*qirr[l];
      if (AccretHeating) {
	for (k=0; k<heatsrc_max; k++) {
          if (l==heatsrc_index[k]) RHS[l] += bfac*heatsrc[k];
	}
      }
    }
  }
  } /* #end of omp parallel section */
  if (first) {
    first = NO;
    IterateRelaxationParameter ();
    omega = omegabest;
    Niterlast = Niterbest;
  }
  omega += domega;		// always try to change the relax. parameter a little (domega set in IterateRelaxationParameter())
  if (omega > 1.999999) {	// stay in the [1,2) interval
    omega = 1.999999;
    domega = -domega;
  }
  if (omega < 1.0) {
    omega = 1.0;
    domega = -domega;
  }
  Niter = SuccessiveOverrelaxation (omega, YES);	// solve the implicit eq. with SOR; YES means that the program will crash when exceeding SORMAXITERS
  if (Niter > Niterlast) domega = - domega;		// if the performance is worse than in the previous case, apply the opposite change next time
  Niterlast = Niter;					// (in this way, omega oscillates around the optimum value)
  if (CPU_Number > 1) SynchronizeOverlapFields (temper, nr, CPUOVERLAP);      /* fill the overlapping ghost zones with values from the neighbouring CPU's active mesh */
  if (CPU_Rank == 0) {		// get the temperature in the innermost and outermost ring by adopting the neighbouring value
    for (j=0; j<ns; j++) {
      temper[j] = temper[j+ns];
    }
  }
  if (CPU_Rank == CPU_Number-1) {
    for (j=0; j<ns; j++) {
      temper[j+(nr-1)*ns] = temper[j+(nr-2)*ns];
    }
  }
  for (i=0; i<nr; i++) { // get the energy from the updated and synchronised temperature
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      energynew[l] = rho[l]*temper[l]*CV;
    }
  }
}

/** When solving the energy equation for the first time,
 * the function spans through various values of the SOR
 * parameter in order to find its best value to start with. */
void IterateRelaxationParameter ()	// note: no explicit MPI stuff here -> everything is synchronized in the reduction routines of SOR
{
  real omegamin=1.0, omegamax=1.999999, omega=1.0;	// parameter from interval [1,2)
  int nr, ns, n, i, j, l, nsplit=9, Niter = SORMAXITERS, count=0;
  real *temper, *tempbckp;
  /* ----- */
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
  temper = Temperature->Field;
  tempbckp = (real *) malloc(sizeof(real) * (nr + 3) * (ns + 1) + 5);	// bckp field with the same size as the field of PolarGrid
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      tempbckp[l] = temper[l];			// backup temperature
    }
  }
  domega = (omegamax - omegamin)/nsplit;	// increment to create nsplit itervals bounded by nsplit+1 values
  Niterbest = SORMAXITERS;			// max no. of iterations in SOR is limited by SORMAXITERS (so we use upper limit as the first guess)
  omegabest = 1.0;				// Gauss-Seidel value as the first guess for relax. parameter
  while (YES) {					// infinite loop (till return is called)
    for (n=0; n<=nsplit; n++) {
      Niter = SuccessiveOverrelaxation (omega, NO);	// NO means that the program won't crash when exceeding SORMAXITERS, it will instead try the next omega
      if (Niter == Niterbest) count++;	// count the number of omegas in the sub-interval that reach the best convergence
      if (Niter < Niterbest) {
        Niterbest = Niter;		// keep track of the lowest number of iterations
        omegabest = omega;		// corresponding to the best relax. parameter
      }       
      omega += domega;			// try next value of omega
      for (i=0; i<nr; i++) {		// return the original temperature values so that we solve the same problem again
	for (j=0; j<ns; j++) {
          l = j + i*ns;
	  temper[l] = tempbckp[l];
	}
      }
    }
    if (Niterbest == SORMAXITERS) {	// this means that no call of SOR converged (must increase SORMAXITERS)
      masterprint ("Warning! SOR did not converge when estimating the initial over-relaxation parameter.\n");
      masterprint ("Try using larger SORMAXITERS. Terminating now...\n");
      prs_exit(1);
    }
    if (count==nsplit+1) {		// plateau = all parameters from the tested interval converge equally fast 
      domega = omegamax-omegamin;	// the width of the sub-interval will be used to oscillate omega around the optimum value
      masterprint ("\tInitial over-relaxation parameter: omega = %#.6g, no. of iterations = %d\n", omegabest, Niterbest);
      free (tempbckp);			// deallocate the bckp field
      return;				// leave the routine
    } else {
      count = 0;			// reset counter for the next interval of omegas
    }
    omegamin = omegabest - domega;	// repeat the test on a sub-interval of the previous range
    omegamax = omegabest + domega;
    if (omegamin < 1.0) omegamin = 1.0;
    if (omegamax > 1.999999) omegamax = 1.999999;
    domega = (omegamax - omegamin)/nsplit;
    omega = omegamin;
  }
}

/** The SOR method algorithm inspired by the one from Numerical Recipes */
int SuccessiveOverrelaxation (omega, errcheck)
real omega;
boolean errcheck;
{
  int n, i, j, l, nr, ns, ipass, jchess, temper_changes, tempint;
  int lip, lim, ljp, ljm;
  real resid, normres0 = 0.0, normres, temp, dtemper;
  real *temper, *Ml, *Mlip, *Mlim, *Mljp, *Mljm, *RHS;
  static boolean first=YES;
  /* ----- */
  temper = Temperature->Field;
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
  Ml = MatrixNexttoTemperl->Field;
  Mlip = MatrixNexttoTemperlip->Field;
  Mlim = MatrixNexttoTemperlim->Field;
  Mljp = MatrixNexttoTemperljp->Field;
  Mljm = MatrixNexttoTemperljm->Field;
  RHS  = RightHandSide->Field;
  if (first) {
    first = NO;
    if (CPU_Number > 1) {
      ChessBoardIndexing ();
    } else {
      jchess1st = 0;
      jchess2nd = 1;
    }
  }
  for (i=One_or_active; i<MaxMO_or_active; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lip = l + ns;
      lim = l - ns;
      ljp = l + 1;
      if (j == ns-1) ljp = i*ns;
      ljm = l - 1;
      if (j == 0) ljm = i*ns + ns-1;
      resid = temper[l]*Ml[l] + temper[lip]*Mlip[l] \
        + temper[lim]*Mlim[l] + temper[ljp]*Mljp[l] \
        + temper[ljm]*Mljm[l] - RHS[l];
      normres0 += fabs(resid);
    }
  }
  MPI_Allreduce (&normres0, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  normres0 = temp;
  for (n=1; n<=SORMAXITERS; n++) {
    normres = 0.0;
    temper_changes = 0;
    for (ipass = 1; ipass <= 2; ipass++) {	// array is swept through as a chessboard (blackorwhite cells first, whiteorblack cells second)
      if (ipass == 1) jchess = jchess1st;
      if (ipass == 2) jchess = jchess2nd;
      // jchess alternately attains value 0 or 1, the vallue changes with each increment of index 'i'. The initial value is always different, so that all the fields on the chessboard are covered
#pragma omp parallel for default(none) \
  shared(One_or_active, MaxMO_or_active, ns, temper, Ml, Mlip, \
	 Mlim, Mljp, Mljm, RHS, omega) \
  private(i, j, l, lip, lim, ljp, ljm, resid, dtemper) \
  firstprivate(jchess) reduction(+ : normres, temper_changes)
      for (i=One_or_active; i<MaxMO_or_active; i++) {
        for (j=jchess; j<ns; j+=2) {	// increment of 'j' index is 2
          l = j + i*ns;
          lip = l + ns;
          lim = l - ns;
          ljp = l + 1;
          if (j == ns-1) ljp = i*ns;
          ljm = l - 1;
          if (j == 0) ljm = i*ns + ns-1;
          resid = temper[l]*Ml[l] + temper[lip]*Mlip[l] \
            + temper[lim]*Mlim[l] + temper[ljp]*Mljp[l] \
            + temper[ljm]*Mljm[l] - RHS[l];
          normres += fabs(resid);
	  dtemper = omega*resid/Ml[l];
          temper[l] -= dtemper;
	  if (dtemper < 0.0) dtemper = -dtemper;	// now check if the iterations really improve the temperature within the given precision limit. If yes, sum the number of improvements
	  if (dtemper > temper[l]*SORMIN) temper_changes++;
        }
        jchess = 1 - jchess; // change the starting value of 'j' for the next increment of 'i'
      }
      if (CPU_Number > 1) SynchronizeOverlapFields (temper, nr, 1);	// must update the 1st ring of each ghost zone neighbouring to the active mesh
    }
    MPI_Allreduce (&normres, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    normres = temp;
    if (normres < SOREPS*normres0) return n;	// standard SOR convergence test
    if (normres < SORMIN) return n;	// finish also when the normalised residue is as small as the given precision
    MPI_Allreduce (&temper_changes, &tempint, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    temper_changes = tempint;
    if (temper_changes == 0) return n;	// finish also when temperature in all cells stops being improved
  }
  if (errcheck) {
    masterprint ("Error!\n");
    masterprint ("The SOR solution of the radiative diffusion equation did not converge.\n");
    masterprint ("Terminating now...\n");
    prs_exit (1);
  }
  return n;
}

/** Function ensures the odd-even ordering of the SOR method 
 * when the grid is split on multiple CPUs. */
void ChessBoardIndexing ()
{
  int send[3], recv[3], prevcpu, nextcpu, nractive;
  /* ----- */
  prevcpu = CPU_Rank-1;
  nextcpu = CPU_Rank+1;
  if (CPU_Master) {
    send[0] = 0;	// inner CPU starts with odd cells during the 1st sweep
    send[1] = 1;	// and with even cells during the 2nd sweep
    nractive = NRAD - CPUOVERLAP - 1;	// master has only one ghost zone, but also the innermost ring is skipped in the SOR method
    if (nractive % 2 == 0) {
      send[2] = 0;	// for even number of active rings, next CPU keeps the same indices
    } else {
      send[2] = 1;	// for odd number of active rings, next CPU swaps the indices
    }
    MPI_Send (&send, 3, MPI_INT, nextcpu, 0, MPI_COMM_WORLD);
  }
  if (CPU_Rank > 0) {
    MPI_Recv (&recv, 3, MPI_INT, prevcpu, 0, MPI_COMM_WORLD, &fargostat);
    if (recv[2] == 0) {
      send[0] = recv[0];
      send[1] = recv[1];
    } else {
      send[0] = recv[1];
      send[1] = recv[0];
    }
    if (CPU_Rank < (CPU_Number-1)) {    // outer CPU does not have to send anything
      nractive = NRAD - 2*CPUOVERLAP;	// middle CPUs have two ghost zones
      if (nractive % 2 == 0) {
        send[2] = 0;
      } else {
        send[2] = 1;
      }
      MPI_Send (&send, 3, MPI_INT, nextcpu, 0, MPI_COMM_WORLD);
    }
  }
  jchess1st = send[0];
  jchess2nd = send[1];
}

/** Calculation of the diffusion coefficients */
void DiffusionCoefs ()
{
  int nr, ns, i, j, l, lim, ljm;
  real *voldens, *temper, *opacity, *D, *Dr, *Dt, *gradtmag;
  real gradt, s, fluxlimiter, fac;
  /* ----- */
  voldens = VolumeDensity->Field;
  temper = Temperature->Field;
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
  opacity = Opacity->Field;
  gradtmag = GradTemperMagnitude->Field;
  D = DiffCoefCentered->Field;
  Dr = DiffCoefIfaceRad->Field;
  Dt = DiffCoefIfaceTheta->Field;
  TemperatureGradient ();	/* compute the magnitude of temperature gradient at cell centers */
#pragma omp parallel default(none) \
  shared(nr,ns,voldens,opacity,temper,gradtmag,D,Dr,Dt) \
  private(i,j,l,lim,ljm,gradt,s,fluxlimiter,fac)
  {
#pragma omp for
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      if (i==0) {
	gradt = gradtmag[l+ns];	/* for inner/outer ring, |gradT| is not available. Adopt the radially neighbouring value. */
      } else if (i==nr-1) {	/* this hopefully does not produce a big error cause we use it only to estimate the flux limiter */
	gradt = gradtmag[l-ns];
      } else {
        gradt = gradtmag[l];
      }
      fac = 1.0/(voldens[l]*opacity[l]);
      s = 4.0*gradt*fac/temper[l];
      fluxlimiter = FluxLimiterValue (s);
      D[l] = fluxlimiter*16.0*STEFANBOLTZMANN*pow(temper[l],3.0)*fac;
    }
  }
#pragma omp for
  for (i=1; i<nr; i++) {	/* calc. average diff. coefficient at interfaces where possible */
    for (j=0; j<ns; j++) {	/* (0th ring is skipped for reasons similar to gradT calculation) */
      l = j + i*ns;
      lim = l - ns;
      ljm = l - 1;
      if (j==0) ljm = i*ns+ns-1;
      Dr[l] = 0.5*(D[l] + D[lim]);
      Dt[l] = 0.5*(D[l] + D[ljm]);
    }
  }
  } // #end of the omp parallel section
}

/** Finds the temperature gradients and their magnitude
 * over the mesh. */
void TemperatureGradient ()
{
  int i, j, l, lim, ljm, lip, ljp, ns, nr; 
  real *temper, *gradtr, *gradtt, *gradtmag;
  real dxtheta, invdxtheta, r, t;
  /* ----- */
  temper = Temperature->Field;
  gradtr = GradTemperRad->Field;
  gradtt = GradTemperTheta->Field;
  gradtmag = GradTemperMagnitude->Field;
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
#pragma omp parallel default(none) \
  shared(nr,ns,temper,InvDiffRmed,Rmed,gradtr,gradtt,gradtmag) \
  private(i,j,l,lim,ljm,lip,ljp,dxtheta,invdxtheta,r,t)
  {
#pragma omp for
  for (i=1; i<nr; i++) {		/* on each CPU, radial gradT exists at all interfaces besides the 0th */
    dxtheta = 2.0*PI/(real)ns*Rmed[i];	/* azimuthal gradT exists even for i=0 but we won't need it anyway */
    invdxtheta = 1.0/dxtheta;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lim = l - ns;
      /* dT/dr = ( T(i,j) - T(i-1,j) )/( Rmed[i] - Rmed[i-1] ) */
      gradtr[l] = (temper[l] - temper[lim])*InvDiffRmed[i];
      ljm = l - 1;
      if (j==0) ljm = i*ns+ns-1;
      /* !polar coordinates -->  1/r * dT/dtheta */
      gradtt[l] = (temper[l] - temper[ljm])*invdxtheta; // don't put an extra Rmed[i] here (is already in invdxtheta) !
    }
  }
#pragma omp for
  for (i=1; i<nr-1; i++) {	/* magnitude of gradT exists at all cell centers which are surrounded by 4 gradT values (one per each interface) */
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      lip = l + ns;
      ljp = l + 1;
      if (j == ns-1) ljp = i*ns;
      r = 0.5*(gradtr[lip]+gradtr[l]);
      t = 0.5*(gradtt[ljp]+gradtt[l]);
      gradtmag[l] = sqrt(r*r + t*t);
    }
  }
  } // end of the omp parallel section
}

/** Translates the surface density into
 * the midplane volume density using the local
 * pressure scale height. */
void MidplaneVolumeDensity (Rho)
PolarGrid *Rho;
{
  int nr, ns, i, j, l;
  real *rho, *voldens, *cs, H;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  cs = SoundSpeed->Field;
  voldens = VolumeDensity->Field;
#pragma omp parallel for default(none) shared(nr,ns,cs,Rmed,rho,voldens,SQRT_ADIABIND_INV,OmegaInv) \
  private(i,j,l,H)
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
	/* we estimate the midplane volume density from the vertically integrated
	   surface density (see e.g. Kley & Crida 2008, Pierens 2015) */
      voldens[l] = rho[l]/H*SQRT2PI_INV;
    }
  }
}

/** Fills the opacity polar grid, either with
 * a fixed parametric value or using the Bell & Lin (1994)
 * opacity table. */
void OpacityProfile ()
{
  int nr, ns, i, j, l;
  real *opacity, *temper, *voldens;
  real T, Dens;
  static boolean FillParametricOpacity = YES;
  /* ----- */
  opacity = Opacity->Field;
  temper = Temperature->Field;
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
  voldens = VolumeDensity->Field;
  if (PARAMETRICOPACITY > 0.0) {
    if (FillParametricOpacity) {        /* if the parametric opacity is used ... */
      FillParametricOpacity = NO;       /* fill the array only once ... */
      for (i=0; i<nr; i++) {
        for (j=0; j<ns; j++) {
          l = j + i*ns;
          opacity[l] = PARAMETRICOPACITY*OPA2CU;
        }
      }
    } else {
      return;                           /* ... and exit the function in subsequent calls */
    }
  }
#pragma omp parallel for default(none) \
  shared(nr,ns,temper,voldens,opacity) \
  private(i,j,l,T,Dens)
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      T = temper[l]*T2SI;             /* temperature transformed from code units to Kelvins */
      Dens = voldens[l]*RHO2CGS;      /* transform to cgs */
      opacity[l] = OpacityBellLin (Dens, T);
      opacity[l] *= OPA2CU;
    }
  }
}

/** Calculates the gas opacity in cgs
 * following the prescription from Bell & Lin (1994).
 * The input volume density and temperature of the gas
 * must be also in cgs. */
real OpacityBellLin (Dens, T)
real Dens, T;
{
  real kappa1, kappa2, kappa3, tmp1, tmp2;
  real opacity;
  /* transitions calculated using eq. (12) in Keith & Wardle (2014);
   * comparisons performed in cgs */
  if (T <= 2286.8*pow(Dens,0.0408)) {
    kappa1 = kappa0[0] * pow(Dens,a[0]) * pow(T,b[0]);
    kappa2 = kappa0[1] * pow(Dens,a[1]) * pow(T,b[1]);
    kappa3 = kappa0[2] * pow(Dens,a[2]) * pow(T,b[2]);
    tmp1 = 1.0/(kappa1*kappa1) + 1.0/(kappa2*kappa2);
    tmp1 = 1.0/(tmp1*tmp1);
    tmp2 = pow(kappa3/(1.0 + 1.0e22*pow(T,-10.0)), 4.0);
    opacity = pow(tmp1+tmp2, 0.25);
  }
  else if (T > 2286.8*pow(Dens,0.0408) && T <= 2029.8*pow(Dens,0.0123)) {
    opacity = kappa0[3] * pow(Dens,a[3]) * pow(T,b[3]);
  }
  else if (T > 2029.8*pow(Dens,0.0123) && T <= 1.0e4*pow(Dens,0.0476)) {
    opacity = kappa0[4] * pow(Dens,a[4]) * pow(T,b[4]);
  }
  else if (T > 1.0e4*pow(Dens,0.0476) && T <= (31195.2)*pow(Dens,0.0533)) {
    opacity = kappa0[5] * pow(Dens,a[5]) * pow(T,b[5]);
  }
  else {
    opacity = kappa0[6] * pow(Dens,a[6]) * pow(T,b[6]);
  }
  return opacity;
}

/** Calculates the flux limiter
 * according to Kley (1989) */
real FluxLimiterValue (s)
real s;
{
  real tmp, value;
  if (s <= 2.0) {
    tmp = 3.0 + sqrt(9.0+10.0*s*s);
    value = 2.0/tmp;
  } else {
    tmp = 10.0*s + 9.0 + sqrt(81.0+180.0*s);
    value = 10.0/tmp;
  }
  return value;
}

/** Calculates the effective optical depth
 * in a simple gray model of the disk's vertical
 * structure; see Hubeny (1990). */
real EffectiveOpticalDepth (tau)
real tau;
{
  real value;
  if (StellarIrradiation) {
    value = 0.375*tau + 0.5 + 0.25/tau;
  } else {
    value = 0.375*tau + 0.25*sqrt(3.0) + 0.25/tau;
  }
  return value;
}

/** For a MPI-split grid, synchronizes the values
 * in a requested number 'nsync' of the overlapping
 * radial rings. */
void SynchronizeOverlapFields (field, nr, nsync)
real *field;
int nr, nsync;
{
  MPI_Request req1, req2, req3, req4;
  static real *SendBufferInner, *SendBufferOuter;
  static real *RecvBufferInner, *RecvBufferOuter;
  static boolean allocate=YES;
  int prevcpu, nextcpu, size;
  int sendoffsetin, sendoffsetout, recvoffsetin, recvoffsetout;
  /* ----- */
  if (nsync > CPUOVERLAP) {
    printf ("Error! Requested number of rings to be synchronized by the MPI communication is larger than CPUOVERLAP.\n");
    printf ("Terminating now...\n");
    prs_exit (1);
  }
  if (allocate) {
    allocate=NO;
    SendBufferInner = malloc (CPUOVERLAP*NSEC*sizeof(real));	// buffers allocated with the max size
    SendBufferOuter = malloc (CPUOVERLAP*NSEC*sizeof(real));
    RecvBufferInner = malloc (CPUOVERLAP*NSEC*sizeof(real));
    RecvBufferOuter = malloc (CPUOVERLAP*NSEC*sizeof(real));
  }
  size = nsync*NSEC;
  prevcpu = CPU_Rank-1;
  nextcpu = CPU_Rank+1;
  sendoffsetin = CPUOVERLAP*NSEC;
  sendoffsetout = (nr - CPUOVERLAP - nsync)*NSEC;
  recvoffsetin = (CPUOVERLAP - nsync)*NSEC;
  recvoffsetout = (nr - CPUOVERLAP)*NSEC;
  if (CPU_Rank > 0) memcpy (SendBufferInner, field+sendoffsetin, size*sizeof(real));
  if (CPU_Rank < CPU_Number-1) memcpy (SendBufferOuter, field+sendoffsetout, size*sizeof(real));
  if (CPU_Rank % 2 == 0) {
    if (CPU_Rank > 0) {
      MPI_Isend (SendBufferInner, size, MPI_DOUBLE, prevcpu, 0, MPI_COMM_WORLD, &req1);
      MPI_Irecv (RecvBufferInner, size, MPI_DOUBLE, prevcpu, 0, MPI_COMM_WORLD, &req2);
    }
    if (CPU_Rank < CPU_Number-1) {
      MPI_Isend (SendBufferOuter, size, MPI_DOUBLE, nextcpu, 0, MPI_COMM_WORLD, &req3);
      MPI_Irecv (RecvBufferOuter, size, MPI_DOUBLE, nextcpu, 0, MPI_COMM_WORLD, &req4);
    }
  } else {
    if (CPU_Rank < CPU_Number-1) {
      MPI_Irecv (RecvBufferOuter, size, MPI_DOUBLE, nextcpu, 0, MPI_COMM_WORLD, &req3);
      MPI_Isend (SendBufferOuter, size, MPI_DOUBLE, nextcpu, 0, MPI_COMM_WORLD, &req4);
    }
    if (CPU_Rank > 0) {
      MPI_Irecv (RecvBufferInner, size, MPI_DOUBLE, prevcpu, 0, MPI_COMM_WORLD, &req1);
      MPI_Isend (SendBufferInner, size, MPI_DOUBLE, prevcpu, 0, MPI_COMM_WORLD, &req2);
    }
  }
  if (CPU_Rank > 0) {
    MPI_Wait (&req1, &fargostat);
    MPI_Wait (&req2, &fargostat);
    memcpy (field+recvoffsetin, RecvBufferInner, size*sizeof(real));
  }
  if (CPU_Rank < CPU_Number-1) {
    MPI_Wait (&req3, &fargostat);
    MPI_Wait (&req4, &fargostat);
    memcpy (field+recvoffsetout, RecvBufferOuter, size*sizeof(real));
  }
}

/** Writes an input file for the 'torquemap' code written
 * by Bertram Bitsch. */
void CreateTorqueMapInfile (istep, Surfdens)
int istep;
PolarGrid *Surfdens;
{
  real *cs, *temper, *pressure, *surfdens;
  real CS[MAX1D], T[MAX1D], P[MAX1D], Sigma[MAX1D];	// need to create global azimuthally averaged fields
  real Rho, tmpr, Opacity, Viscalpha, H;
  real Omega, kappa1, kappa2, kappa3, tmp1, tmp2;
  int i;
  char name[256];
  FILE *output;
  // static boolean first=YES;
  /* ---- */
  cs = SoundSpeed->Field;
  temper = Temperature->Field;
  pressure = Pressure->Field;
  surfdens = Surfdens->Field;
  if (CPU_Master) {
    sprintf (name, "%storquemap_infile%d.dat", OUTPUTDIR, istep);
    output = fopenp (name, "w");
  }
  mpi_make1Dprofile (cs, CS);
  mpi_make1Dprofile (temper, T);
  mpi_make1Dprofile (pressure, P);
  mpi_make1Dprofile (surfdens, Sigma);
  if (CPU_Master) {
    masterprint("Writing torque map input file ...");
    for (i=0; i<GLOBALNRAD; i++) {
      Omega = pow(GlobalRmed[i],-1.5);
      H = CS[i]/Omega/sqrt(ADIABIND);
      Rho = Sigma[i]/(sqrt(2.0*PI)*H)*RHO2CGS;
      tmpr = T[i]*T2SI;
      // H[i] = CS[i]/Omega/sqrt(ADIABIND);
      // Rho[i] = Sigma[i]/(2.0*H[i]);
      /* OPACITY PART IN CGS ---> */
      if (tmpr <= 2286.8*pow(Rho,0.0408)) {
        kappa1 = kappa0[0] * pow(Rho,a[0]) * pow(tmpr,b[0]);
        kappa2 = kappa0[1] * pow(Rho,a[1]) * pow(tmpr,b[1]);
        kappa3 = kappa0[2] * pow(Rho,a[2]) * pow(tmpr,b[2]);
        tmp1 = 1.0/(kappa1*kappa1) + 1.0/(kappa2*kappa2);
        tmp1 = 1.0/(tmp1*tmp1);
        tmp2 = pow(kappa3/(1.0 + 1.0e22*pow(tmpr,-10.0)), 4.0);
        Opacity = pow(tmp1+tmp2, 0.25);
      }
      else if (tmpr > 2286.8*pow(Rho,0.0408) && tmpr <= 2029.8*pow(Rho,0.0123)) {
        Opacity = kappa0[3] * pow(Rho,a[3]) * pow(tmpr,b[3]);
      }
      else if (tmpr > 2029.8*pow(Rho,0.0123) && tmpr <= 1.0e5*pow(Rho,0.0476)) {
        Opacity = kappa0[4] * pow(Rho,a[4]) * pow(tmpr,b[4]);
      }
      else if (tmpr > 1.0e5*pow(Rho,0.0476) && tmpr <= (31195.2)*pow(Rho,0.0533)) {
        Opacity = kappa0[5] * pow(Rho,a[5]) * pow(tmpr,b[5]);
      }
      else {
        Opacity = kappa0[6] * pow(Rho,a[6]) * pow(tmpr,b[6]);
      }
      if (PARAMETRICOPACITY > 0.0) Opacity = PARAMETRICOPACITY;
      /* <--- */
      Viscalpha = FViscosity (GlobalRmed[i]);	// this will ALWAYS return the kinematic viscosity
      Viscalpha /= (H*H*Omega);			// this will ALWAYS convert it to alpha
      /* output everything in cgs, only radial coordinate and H should be in AU */
      fprintf (output, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",\
        GlobalRmed[i], Rho, tmpr, P[i]*PRESS2CGS, Opacity, Sigma[i]*SIGMA2CGS, Viscalpha, H);
    }
    fclose (output);
    masterprint(" done\n");
  }
}

/** Initialization routine which improves the temperature
 * profile for a fixed surface density profile. Iteratively
 * equates 'vertical_cooling = simple_viscous_heating'. 
 * Note that the function modifies the initial energy (not the
 * temperature itself) which is used as a starting point for 
 * the initialization. */
void IterateInitialTemperature ()
{
  const real ITEREPS = 1.e-4;
  const int ITERMAX = 10000;
  int i, n;
  real T0, Omega2, H, rho, kappa, tau, taueff, T1, dT;
  real alpha, viscosity, cs, dalpha, tmri, deltat;
  /* ----- */
  masterprint ("\n---------------\n");
  masterprint ("Balancing the vertical cooling and viscous heating\n");
  masterprint ("to improve the initial temperature profile...\n");
  for (i=0; i<NRAD; i++) {
    T0 = EnergyMed[i]*MOLWEIGHT*(ADIABIND-1.0)/(SigmaMed[i]*GASCONST);	// 0th temperature comes from the initial (bad) internal energy profile
    Omega2 = pow(Rmed[i],-3.0);
    n = 0;
    while (YES) {	// iterate until the precision is reached or maximum number of iterations exceeded
      H = sqrt(GASCONST/MOLWEIGHT*T0)*OmegaInv[i];	// get H(T0)
      rho = SigmaMed[i]*SQRT2PI_INV/H;			// get rho(Sigma,H)
      kappa = OpacityBellLin (rho*RHO2CGS, T0*T2SI);	// get kappa from Bell & Lin (1994)
      kappa *= OPA2CU;
      tau = OPACITYDROP*0.5*kappa*SigmaMed[i];
      taueff = EffectiveOpticalDepth (tau);		// get tau_eff from Hubeny (1990)
      if (ViscosityAlpha) {				// get viscosity (from alpha if necessary)
	if (!AlphaFlock) {
          alpha = ALPHAVISCOSITY;
        } else {
          dalpha = ALPHAMRI - ALPHADEAD;
	  tmri = TMRI/T2SI;
	  deltat = TWIDTH/T2SI;
	  alpha = 0.5*dalpha*(1.0 - tanh((tmri - T0)/deltat)) + ALPHADEAD;
        }
	cs = sqrt(ADIABIND*(ADIABIND-1.0)*EnergyMed[i]/SigmaMed[i]);
	viscosity = alpha*cs*cs*OmegaInv[i];
      } else {
        viscosity = VISCOSITY;
      }
      // find new T1 by solving the energy balance for taueff from above and fixed Sigma
      T1 = pow(9.0*SigmaMed[i]*viscosity*Omega2*taueff/(8.0*STEFANBOLTZMANN), 0.25);	// e.g. truncated Eq. (10) from Kretke & Lin (2012) 
      dT = fabs(T1-T0);
      if (fabs(T1-T0)/T0 < ITEREPS) {
	EnergyMed[i] = 0.5*(T0+T1)*SigmaMed[i]*GASCONST/(MOLWEIGHT*(ADIABIND-1.0));
	// printf ("no. of iters: %d\n", n);
        break;
      }	else {
	n++;
	if (n > ITERMAX) {
          printf ("Can't iterate the initial temperature profile within %d iterations.\n", ITERMAX);
	  printf ("Problem for CPU %d, radial distance %f, local radial index %d.\n", \
            CPU_Rank, Rmed[i], i);
	  printf ("Last two temperatures [K] of the iterative process: %#.5g\t%#.5g\n", T0*T2SI, T1*T2SI);
	  printf ("Terminating now...\n");
	  prs_exit (1);
	}
        // to avoid false convergence, new temperature (T0 for the next step) is set randomly between T0,T1 (or T1,T0)
        dT *= (real)rand() / (real)RAND_MAX;
	if (T0 > T1) T0 -= dT;
	if (T0 < T1) T0 += dT;
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  masterprint ("...done on all CPUs\n");
  masterprint ("---------------\n");
}
