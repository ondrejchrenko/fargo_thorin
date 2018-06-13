/** \file global_ex.h

This file is created      
automatically during      
compilation from global.h. Do not edit. 
See perl script           
"varparser.pl" for details

\file global.h

Declares all global variables.
Used to construct automatically
the file global_ex.h. The file
global.h cannot contain any comment,
as it would not be parsed correctly
by varparser.pl
*/                          

extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master;
extern int IMIN;
extern int IMAX;
extern int Zero_or_active;
extern int Max_or_active;
extern int One_or_active;
extern int MaxMO_or_active;		/* MO: Minus One */
extern int GLOBALNRAD;
extern real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
extern real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
extern real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
extern real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper;
extern real OmegaInv[MAX1D], Rmed2[MAX1D];					/* #THORIN */
extern real EnergyMed[MAX1D], globpressvec[MAX1D], globcsvec[MAX1D], WaveKiller[MAX1D], VthetaMed[MAX1D];	/* #THORIN */
extern real QplusMed[MAX1D], CoolingTimeMed[MAX1D];				/* #THORIN */
extern real PebDensInit[MAX1D], PebVradInit[MAX1D], PebVthetaInit[MAX1D];	/* #THORIN */
extern real vt1D[MAX1D], invdtpeb_sq, invdtreb_sq, SQRT_ADIABIND_INV;		/* #THORIN */
extern real OmegaFrame, PhysicalTime, PhysicalTimeInitial;
extern real heatsrc[MAXPLANETS];						/* #THORIN */
extern int heatsrc_max;							/* #THORIN */
extern int TimeStep;
extern boolean EnergyEq, StoreEnergy, ParametricCooling, Damping, DampVrad, DampInit, StellarIrradiation;	/* #THORIN */
extern boolean InitFromFile, Write_Temperature, Write_Energy, Write_Divergence, Write_Qplus, Write_Qbalance;			/*#THORIN*/
extern boolean Collisions, WriteTorque, WriteTorqueMapFile, MonitorNPL, FeelDisk;  				/* #THORIN */
extern boolean Pebbles, Write_Eta, AccretHeating, BackReaction, ActualizeLuminosity, DiffusiveParticles, PrescribedAccretion; 	/* #THORIN */
extern boolean heatsrc_index[MAXPLANETS], TorqueDensity;	/* #THORIN */
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, ViscosityAlpha, RocheSmoothing;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
extern MPI_Status fargostat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt;
extern PolarGrid *Temperature, *Pressure, *SoundSpeed, *Qplus, *Qminus, *Qbalance; 		/* #THORIN */
extern PolarGrid *DivergenceVelocity, *TAURR, *TAURP, *TAUPP;					/* #THORIN */
extern PolarGrid *GasAccelrad, *GasAcceltheta, *DragForceRad, *DragForceTheta;			/* #THORIN */
extern PolarGrid *PebbleDens, *PebbleVrad, *PebbleVtheta, *StokesNumber;
extern PolarGrid *GravAccelRad, *GravAccelTheta, *PebbleGravAccelRad, *PebbleGravAccelTheta;	/* #THORIN */
extern PolarGrid *Torque;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
extern FILE *plout, *discard, *mergers;  /* #THORIN */
