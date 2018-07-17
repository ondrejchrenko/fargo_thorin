/** \file main.c

Main file of the distribution. Manages the call to initialization
functions, then the main loop.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

boolean         Restart = NO, OpenInner = NO;
int             begin_i = 0, NbRestart = 0;
static int      InnerOutputCounter=0, StillWriteOneOutput;
extern real     LostMass;
extern boolean  Corotating;
real            ScalingFactor = 1.0;

boolean		DumpTorqueNow = NO, DumpTorqueDensNow = NO;	/* #THORIN: output control switches */

int
main(argc, argv)
int argc;
char *argv[];
{
  PolarGrid   *gas_density, *gas_v_rad, *gas_v_theta, *gas_label;
  PolarGrid   *gas_energy;	/* #THORIN */
  int          i;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      TimeToWrite, verbose = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;
  struct reb_simulation *rsim;	/* #THORIN: structure holds a rebound simulation coupled with FARGO */
  int npl;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();			/* Control behavior for floating point
				   exceptions trapping (default is not to do anything) */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-secndovtpfamzibM0123456789") != strlen (argv[i]))	/* #THORIN */
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'c'))
	SloppyCFL = YES;
      if (strchr (argv[i], 'p'))
	Profiling = YES;
      if (strchr (argv[i], 'd'))
	debug = YES;
      if (strchr (argv[i], 'b'))
	CentrifugalBalance = YES;
      if (strchr (argv[i], 'm'))
	Merge = YES;
      if (strchr (argv[i], 'M'))		/* #THORIN */
	MergeDirect = YES;
      if (strchr (argv[i], 'a'))
      if (strchr (argv[i], 'a'))
	MonitorIntegral = YES;
      if (strchr (argv[i], 'z'))
	FakeSequential = YES;
      if (strchr (argv[i], 'i'))
	StoreSigma = YES;
        if (EnergyEq) StoreEnergy = YES;	/* #THORIN */
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
      if ((argv[i][1] >= '1') && (argv[i][1] <= '9')) {
	GotoNextOutput = YES;
	StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'o')) {
	OverridesOutputdir = YES;
	i++;
	sprintf (NewOutputdir, "%s", argv[i]);
      } else {
	if (strchr (argv[i], 'f')) {
	  i++;
	  ScalingFactor = atof(argv[i]);
	  masterprint ("Scaling factor = %g\n", ScalingFactor);
	  if ((ScalingFactor <= 0)) {
	    masterprint ("Incorrect scaling factor\n");
	    PrintUsage (argv[0]);
	  }
	}
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }
  if ((StoreSigma || StoreEnergy) && !(Restart)) {	/* #THORIN */
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("or surface internal energy in a non-restart run.\n");
    mastererr ("Aborted\n");
    prs_exit (1);
  }
  if (Merge && MergeDirect) {
    mastererr ("ERROR - You cannot use -m and -M switches at the same time.\n");
    mastererr ("Restart with only one of them to merge the output files\n");
    mastererr ("from individual CPUs. Terminating now...\n");
    prs_exit (1);
  }
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
  ReadVariables (ParameterFile);	/* #THORIN: InitPlanetarySystem() and ListPlanets() used to be here, replaced by functions of the Rebound interface */
  SplitDomain ();
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  MakeDir (OUTPUTDIR);
  DumpSources (argc, argv);
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_energy         = CreatePolarGrid(NRAD, NSEC, "energy");	/* #THORIN */
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  masterprint ("done.\n");
  /* #THORIN ---> */
  npl = FindNumberOfPlanets (PLANETCONFIG); /* see Psys.c */
  if (CPU_Master) printf ("%d planet(s) found.\n", npl);
  sys = AllocPlanetSystem (npl);        /* see Psys.c */
  sys->nb = npl;
  if (Restart == YES) {
    begin_i = NbRestart * NINTERM;
    rsim = RestartReboundSimulation (sys, NbRestart);
    ListPlanets (sys);
  } else {
    rsim = SetupReboundSimulation (sys, PLANETCONFIG);
    ListPlanets (sys);
  }
  /* <--- #THORIN */
  OmegaFrame = OMEGAFRAME;
  if (Corotating == YES) OmegaFrame = GetPsysInfo (sys, FREQUENCY);
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label); /* #THORIN */
  if (AccretHeating) InitAccretHeatSrc (sys->nb);	/* #THORIN */
  InitComputeAccel ();
  if (Restart) {					/* #THORIN */
    OmegaFrame = GetOmegaFrame (NbRestart);
    GetIterStat (NbRestart);
  } else {
    EmptyIterStat ();
  }
  PhysicalTimeInitial = PhysicalTime;
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    if (InnerOutputCounter == 1) {
      InnerOutputCounter = 0;
      if (WriteTorque) DumpTorqueNow=YES;		/* #THORIN: DumpTorqueNow will allow for the torque calculation in FillForcesArrays ()  */
    }
    if (NINTERM * (TimeStep = (i / NINTERM)) == i)	/* Outputs are done here */ {
      TimeToWrite = YES;
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label);	/* #THORIN: see Output.c */
      OutputNbodySimulation (TimeStep, rsim);		/* #THORIN: output binary file with the N-body part settings */
      if (WriteTorqueMapFile) CreateTorqueMapInfile (TimeStep, gas_density);	/* #THORIN */
      if (TorqueDensity) DumpTorqueDensNow=YES;					/* #THORIN */
      DumpOmegaFrame (TimeStep);			/* #THORIN: print OmegaFrame - it is needed for restart runs */
      DumpIterStat (TimeStep);				/* #THORIN: output info about the SOR setup - it is needed for restart runs */
      if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
	MPI_Finalize();
	return 0;
      }
      StillWriteOneOutput--;
      if (TimeInfo == YES)	/* Time monitoring is done here */
	GiveTimeInfo (TimeStep);
    }
    else {
      TimeToWrite = NO;
    }
				/* Algorithm loop begins here */

				/***********************/
				/* Hydrodynamical Part */
				/***********************/
    InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
    AlgoGas (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys, rsim);	/* #THORIN: see SourceEuler.c */
    GiveSpecificTime (Profiling, t_Hydro);
    if (NOUTELEMENTS * ((i+1) /NOUTELEMENTS) == (i+1)) {	/* #THORIN */
      OutputElements (rsim);
      fflush (plout);
      fflush (discard);
      if (Collisions == YES) fflush (mergers);
    }
    if (MonitorNPL == YES) {	/* #THORIN: If the target number of planets is reached, terminate the run */
      if (sys->nb <= TARGETNPL) {
        masterprint ("Target number of planets %d was reached. Normal termination.\n", TARGETNPL);
	break;
      }
    }
    /* <-- */
    if (MonitorIntegral == YES) {
      masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
      masterprint ("Gas total Energy: %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy)); /* #THORIN: see SideEuler.c */
    }
  }
  reb_free_simulation (rsim);   /* #THORIN */
  FreePlanetary (sys);
  fclose (plout);       	/* #THORIN */
  fclose (discard);
  if (Collisions==YES) fclose (mergers);
  MPI_Finalize ();
  return 0;
}
