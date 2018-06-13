int CPU_Rank;
int CPU_Number;
boolean CPU_Master;
int IMIN;
int IMAX;
int Zero_or_active;
int Max_or_active;
int One_or_active;
int MaxMO_or_active;		/* MO: Minus One */
int GLOBALNRAD;
real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper;
real OmegaInv[MAX1D], Rmed2[MAX1D];					/* #THORIN */
real EnergyMed[MAX1D], globpressvec[MAX1D], globcsvec[MAX1D], WaveKiller[MAX1D], VthetaMed[MAX1D];	/* #THORIN */
real QplusMed[MAX1D], CoolingTimeMed[MAX1D];				/* #THORIN */
real PebDensInit[MAX1D], PebVradInit[MAX1D], PebVthetaInit[MAX1D];	/* #THORIN */
real vt1D[MAX1D], invdtpeb_sq, invdtreb_sq, SQRT_ADIABIND_INV;		/* #THORIN */
real OmegaFrame, PhysicalTime=0.0, PhysicalTimeInitial;
real heatsrc[MAXPLANETS];						/* #THORIN */
int heatsrc_max;							/* #THORIN */
int TimeStep=0;
boolean EnergyEq, StoreEnergy, ParametricCooling, Damping, DampVrad, DampInit, StellarIrradiation;	/* #THORIN */
boolean InitFromFile, Write_Temperature, Write_Energy, Write_Divergence, Write_Qplus, Write_Qbalance;			/*#THORIN*/
boolean Collisions, WriteTorque, WriteTorqueMapFile, MonitorNPL, FeelDisk;  				/* #THORIN */
boolean Pebbles, Write_Eta, AccretHeating, BackReaction, ActualizeLuminosity, DiffusiveParticles, PrescribedAccretion; 	/* #THORIN */
boolean heatsrc_index[MAXPLANETS], TorqueDensity;	/* #THORIN */
boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
boolean	GotoNextOutput, StoreSigma, ViscosityAlpha, RocheSmoothing;
boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
MPI_Status fargostat;
PolarGrid *CellAbscissa, *CellOrdinate;
PolarGrid *RhoStar, *RhoInt;
PolarGrid *Temperature, *Pressure, *SoundSpeed, *Qplus, *Qminus, *Qbalance; 		/* #THORIN */
PolarGrid *DivergenceVelocity, *TAURR, *TAURP, *TAUPP;					/* #THORIN */
PolarGrid *GasAccelrad, *GasAcceltheta, *DragForceRad, *DragForceTheta;			/* #THORIN */
PolarGrid *PebbleDens, *PebbleVrad, *PebbleVtheta, *StokesNumber;
PolarGrid *GravAccelRad, *GravAccelTheta, *PebbleGravAccelRad, *PebbleGravAccelTheta;	/* #THORIN */
PolarGrid *Torque;
boolean LogGrid;
boolean OverridesOutputdir;
char NewOutputdir[1024];
FILE *plout, *discard, *mergers;  /* #THORIN */
