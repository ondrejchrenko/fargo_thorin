/**
 * @file ReboundInterface.c
 *
 * @brief Contains the functions interfacing FARGO with
 * the REBOUND package.
 *
 * @author Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>
 *
 * @details The default setup uses the IAS15 integrator to
 * propagate the planets in time and employs
 * the direct collision search with merging, if turned on.
 * This can be easily reprogrammed if needed.
 *
 * Note: This version has no test particles, rsim->N could be
 * used directly as the loop limit. But rsim->N_active is used instead
 * so that test particles could be easily implemented. If rsim->N is used,
 * it indicates that a loop would in principle handle test particles as well.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Ondřej Chrenko. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

extern real OmegaFrame, MassTaper;
extern boolean Corotating;

static real rdenominv;

/** If the REBOUND collision search is successful,
 * this function merges the bodies, outputs information
 * about the merger event and reorganises the particle list. */
int ResolveCollisions (rsim, coll)
struct reb_simulation *rsim;
struct reb_collision coll;
{
  int value=0, pi, pj, pk=0;
  struct reb_particle* const particles = rsim->particles;	/* note: only designed for merging planets at this point !!! */
  pi = coll.p1;		/* array indices of particles participating in the collision have no strict order! both i<j and j<i are possible */
  pj = coll.p2;
  if (CPU_Master) {
    printf ("\n Merger detected. Merging planets %d and %d\n", particles[pi].hash, particles[pj].hash);
    masterprint ("Previous no. of planets: %d\n", rsim->N_active-1);
    fprintf (mergers, "%.12g\t%d\t%d\n", rsim->t, pi, pj);
    fprintf (mergers, "%d\t%f\t%.8g\n", particles[pi].hash, particles[pi].m, particles[pi].r);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pi].x, particles[pi].y, particles[pi].z);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pi].vx, particles[pi].vy, particles[pi].vz);
    fprintf (mergers, "%d\t%f\t%.8g\n", particles[pj].hash, particles[pj].m, particles[pj].r);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pj].x, particles[pj].y, particles[pj].z);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pj].vx, particles[pj].vy, particles[pj].vz);
  }
  value = reb_collision_resolve_merge (rsim, coll);
  rsim->N_active--;		/* change number of massive bodies accordingly */
  if (CPU_Master) {
    switch (value) {
      case 1 :
        ; pk = pj;	/* i-th particle will be discarded, j-th particle replaced by the merger -> output it */
        break;
      case 2 :
        ; pk = pi;	/* vice versa */
        break;
    }
    fprintf (mergers, "%d\t%f\t%.8g\n", particles[pk].hash, particles[pk].m, particles[pk].r);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pk].x, particles[pk].y, particles[pk].z);
    fprintf (mergers, "%#.18g\t%#.18g\t%#.18g\n", particles[pk].vx, particles[pk].vy, particles[pk].vz);
    masterprint ("No. of planets after the merger: %d\n", rsim->N_active-1);
  }
  return value;
}

/** A simple discard routine which looks for
 * planets that were scattered/migrated away from
 * the FARGO grid. Must be called from the heliocentric
 * frame. */
void DiscardParticlesDist (rsim, dt)
struct reb_simulation *rsim;
real dt;
{
  int Nact, i;
  real xr, yr, r2;
  boolean isremoved;
  boolean const keepsorted=YES;
  struct reb_particle* const particles = rsim->particles;
  /* ----- */
  Nact = rsim->N_active;
  for (i=1; i<Nact; i++) {
    xr = particles[i].x;
    yr = particles[i].y;
    r2 = xr*xr + yr*yr;		// check whether the planet's midplane projection left the FARGO grid
    if (r2 < RMIN*RMIN || r2 > RMAX*RMAX) {
      masterprint ("\nPlanet %d discarded. Reason: it is not in the (RMIN, RMAX) domain.\n", particles[i].hash);
      masterprint ("Previous no. of planets: %d\n", rsim->N_active-1);
      if (CPU_Master) {
	fprintf (discard, "%.12g\n", PhysicalTime+dt);
	fprintf (discard, "%d\t%f\t%#.8g\n", particles[i].hash, particles[i].m, particles[i].r);
	fprintf (discard, "%#.18g\t%#.18g\t%#.18g\n", particles[i].x, particles[i].y, particles[i].z);
	fprintf (discard, "%#.18g\t%#.18g\t%#.18g\n", particles[i].vx, particles[i].vy, particles[i].vz);
      }
      isremoved = reb_remove (rsim, i, keepsorted);
      if (isremoved == NO) {
        masterprint ("\nWarning! Failed to remove planet with id %d and loop index %d\n", particles[i].hash, i);
      }
      masterprint ("No. of planets after the discard: %d\n", rsim->N_active-1);		// no need to change N_active manually, reb_remove does this now
    }
  }
}

/** Fills the rebound simulation structure
 * with parameters inherited from FARGO.
 * The integrator type can be easily changed from here. */
void SetupIntegratorParams (rsim)
struct reb_simulation *rsim;
{
  /* rsim->integrator = REB_INTEGRATOR_WHFAST; alternative integrators could be placed here... */
  rsim->integrator = REB_INTEGRATOR_IAS15;
  rsim->ri_ias15.epsilon = IAS15PRECISSION;
  rsim->ri_ias15.min_dt = IAS15MINDT;
  if (Collisions == YES) {
    rsim->collision = REB_COLLISION_DIRECT;
    rsim->collision_resolve = ResolveCollisions;
  } else {
    rsim->collision = REB_COLLISION_NONE;
  }
}

/** Initialises a rebound simulation coupled with FARGO */
struct reb_simulation *SetupReboundSimulation (sys, plfile)
PlanetarySystem *sys;
char *plfile;
{
  FILE *input;
  char s[512], nm[512], filename[256], *s1;
  int npl, i, err;
  float mass, a, e, inc, Omega, omega, f;
  real rho, radius, rad;
  /* ----- */
  struct reb_simulation *rsim = reb_create_simulation ();
  SetupIntegratorParams (rsim);
  npl = sys->nb;
  rsim->N_active = npl + 1;		// N_active = number of massive particles, i.e. planets + primary
  rad = PI/180.0;			// angles to radians conversion
  struct reb_particle primary = {0};	// initialize the primary --->
  primary.hash = 0;			// new version uses 'hash' member instead of 'id' member of reb_particle structure
  primary.m = 1.0;
  reb_add (rsim, primary);		// <---
  input = fopen (plfile, "r");		// initialize the planets; some stuff is similar to InitPlanetarySystem() in psys.c
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'.\n", filename);
    prs_exit (1);
  }
  rho = PLANETARYDENSITY/RHO2CGS;
  rdenominv = 1.0/(4.0*PI*rho);
  i = 1;				// REBOUND has the primary at i=0
  while (fgets(s, 510, input) != NULL) {
    sscanf(s, "%s ", nm);
    if (isalpha(s[0])) {
      s1 = s + strlen(nm);
      // input as the planet mass and orbital elements. Accretion and FEELDISK options are managed by parameters. FEELOTHERS=YES is implicit.
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %f %f %f %f", \
        &mass, &a, &e, &inc, &Omega, &omega, &f);
      err = 0;
      struct reb_particle p = reb_tools_orbit_to_particle_err (G, primary, \
        (real)mass, (real)a, (real)e, (real)inc*rad, (real)Omega*rad, (real)omega*rad, (real)f*rad, &err);
      if (err != 0) {
        masterprint("Error: Planet no. %d cannot be initialized \
	  using the input parameters. Terminating now...\n", i);
	prs_exit (1);
      }
      p.hash = i;
      radius = pow(3.0*(real)mass*rdenominv, 1.0/3.0);
      p.r = radius;			// future 2DO - inflation parameter could be implemented here
      reb_add (rsim, p);	
      sys->mass[i-1] = p.m;		// [i-1] because FARGO has the 1st planet at i=0
      sys->x[i-1] = p.x;
      sys->y[i-1] = p.y;
      sys->z[i-1] = p.z;
      sys->vx[i-1] = p.vx;
      sys->vy[i-1] = p.vy;
      sys->vz[i-1] = p.vz;
      sys->acc[i-1] = ACCRETIONRATE;
      sys->FeelDisk[i-1] = FeelDisk;
      i++;
    }
  }
  fclose (input);
  sprintf (filename, "%snbody.orbits.dat", OUTPUTDIR);
  plout = fopenp (filename, "w");     // "w" should empty the file if it exists 
  sprintf (filename, "%snbody.discard.dat", OUTPUTDIR);
  discard = fopenp (filename, "w");
  if (Collisions == YES) {
    sprintf (filename, "%snbody.mergers.dat", OUTPUTDIR);
    mergers = fopenp (filename, "w");
  }
  return rsim;
}

/** Performs an integration step of the N-body problem */
void AdvanceSystemRebound (sys, rsim, dt)
PlanetarySystem *sys;
struct reb_simulation *rsim;
real dt;
{
  int i;
  real ttarget, mass, radius;
  struct reb_particle com;
  struct reb_particle* const particles = rsim->particles;	/* asterisk means that the pointer is not to be changed, unlike the stuff it points to */
  for (i=1;i<rsim->N_active;i++) {	/* only planets have to be synchronised */
    if (MassTaper < 1.0) {
      mass = MassTaper*(sys->mass[i-1]);			/* masses can change due to accretion etc; mass taper also accounted for */
    } else {
      mass = sys->mass[i-1];
    }
    if (mass != particles[i].m) {				/* if masses change, radii must be updated */
      particles[i].m = mass;
      radius = pow(3.0*mass*rdenominv, 1.0/3.0);
      particles[i].r = radius;
    }
    particles[i].x = sys->x[i-1];
    particles[i].y = sys->y[i-1];
    particles[i].z = sys->z[i-1];
    particles[i].vx = sys->vx[i-1];
    particles[i].vy = sys->vy[i-1];
    particles[i].vz = sys->vz[i-1];
  }
  com = reb_get_com (rsim);		/* a particle representing the centre of mass */
  for (i=0;i<rsim->N;i++) {		/* convert helioc->baryc */
    particles[i].x -= com.x;
    particles[i].y -= com.y;
    particles[i].z -= com.z;
    particles[i].vx -= com.vx;
    particles[i].vy -= com.vy;
    particles[i].vz -= com.vz;
  }
  rsim->t = PhysicalTime;
  rsim->dt = dt;
  ttarget = PhysicalTime + dt;
  if (rsim->integrator==REB_INTEGRATOR_IAS15) reb_integrator_ias15_clear (rsim);	/* dt and planet stuff might have changed since the last step due to e.g. CFL; planet-disk interactions ... */
  reb_integrate (rsim, ttarget);							/* ... one cannot use estimates from the last step and has to reset auxiliary arrays of IAS15 */
  for (i=1;i<rsim->N;i++) {		/* convert back to heliocentric (to apply disk-planet interactions and rotations in FARGO */
    particles[i].x -= particles[0].x;
    particles[i].y -= particles[0].y;
    particles[i].z -= particles[0].z;
    particles[i].vx -= particles[0].vx;
    particles[i].vy -= particles[0].vy;
    particles[i].vz -= particles[0].vz;
  }
  particles[0].x = 0.0;			/* just to be sure that the origin = the primary */
  particles[0].y = 0.0;
  particles[0].z = 0.0;
  particles[0].vx = 0.0;
  particles[0].vy = 0.0;
  particles[0].vz = 0.0;
  DiscardParticlesDist (rsim, dt);	// !has to be called after the helioc. transformation
  /* Here we do not put all planetary stuff back to 'sys', this is done later on
     (after rotations etc) in SynchronizeFargoRebound ().*/
}

/** Synchronises the planetary system between the
 * REBOUND integration and the FARGO simulation. */
void SynchronizeFargoRebound (sys, rsim)
PlanetarySystem *sys;
struct reb_simulation *rsim;
{
  int i;
  boolean lessplanets = NO;
  struct reb_particle* const particles = rsim->particles;
  if (sys->nb > rsim->N_active-1) {	/* if there are more planets in FARGO than in REBOUND */
    lessplanets = YES;			/* something must have happened */
    sys->nb = rsim->N_active-1;		/* so shrink the planetary array in FARGO so all the loops still work */
  }
  if (sys->nb < rsim->N_active-1) {
    printf ("Error! There are more planets in REBOUND than in FARGO! Terminating now...\n");
    prs_exit (1);
  } 
  for (i=1;i<rsim->N_active;i++) {	/* put planetary stuff back */
    sys->x[i-1] = particles[i].x;
    sys->y[i-1] = particles[i].y;
    sys->z[i-1] = particles[i].z;
    sys->vx[i-1] = particles[i].vx;
    sys->vy[i-1] = particles[i].vy;
    sys->vz[i-1] = particles[i].vz;
    if (lessplanets) sys->mass[i-1] = particles[i].m;
    if (MassTaper < 1.0) sys->mass[i-1] = particles[i].m/MassTaper;	/* Mass tapering must be taken back here!!! */
  }
}

/** A simple time step restriction in order
 * not to miss a collision. This should be in principle
 * always be overridden by the IAS15 time step division. */
void MinStepForRebound (rsim)
struct reb_simulation *rsim;
{
  int Nact, i, j;
  real dt, dtnew, r2, rv;
  struct reb_particle* const particles = rsim->particles;
  Nact = rsim->N_active;
  dt = 1.e30;
  for (i=1; i<Nact; i++) {	// loop over planet pairs, make a check as e.g. (5) in Pierens et al. (2013)
    for (j=i+1; j<Nact; j++) {
      r2 = (particles[i].x - particles[j].x)*(particles[i].x - particles[j].x) + \
	   (particles[i].y - particles[j].y)*(particles[i].y - particles[j].y) + \
	   (particles[i].z - particles[j].z)*(particles[i].z - particles[j].z);
      rv = (particles[i].x - particles[j].x)*(particles[i].vx - particles[j].vx) + \
	   (particles[i].y - particles[j].y)*(particles[i].vy - particles[j].vy) + \
	   (particles[i].z - particles[j].z)*(particles[i].vz - particles[j].vz);
      dtnew = fabs(r2/rv);
      if (dtnew < dt) {
        dt = dtnew;
      }
    }
  }
  dt *= 0.25;	// safety factor
  dt *= dt;	// final time needs to be squared before combined in CFL
  invdtreb_sq = 1.0/dt;
}

/** Calculates and outputs the orbital elements. */
void OutputElements (rsim)	// note: plout is a global file name
struct reb_simulation *rsim;
{
  int Nact, i;
  struct reb_orbit orbit;
  struct reb_particle* const particles = rsim->particles;
  if (CPU_Rank != CPU_Number-1) return;
  Nact = rsim->N_active;
  for (i=1; i<Nact; i++) {
    orbit = reb_tools_particle_to_orbit (G, particles[i], particles[0]);
    fprintf (plout, "%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", \
      particles[i].hash, PhysicalTime, orbit.a, orbit.e, orbit.inc, \
      orbit.Omega, orbit.omega, orbit.f, particles[i].m, particles[i].r,
      particles[i].x, particles[i].y, particles[i].z);    
  }
}

/** Stores the entire Rebound simulation
 * in a binary file. Useful for restarts. */
void OutputNbodySimulation (nout, rsim)
int nout;
struct reb_simulation *rsim;
{
  char filename[256];
  if (!CPU_Master) return;
  printf ("Writing the N-body simulation file...\n");
  printf ("\tNumber of planets: %d\n", rsim->N_active-1);
  sprintf (filename, "%snbody.simulation%d.dat", OUTPUTDIR, nout);
  reb_output_binary (rsim, filename);
  printf ("\t\t...done\n");
  fflush (stdout);
}

/** Part of the restart process. */
struct reb_simulation *RestartReboundSimulation (sys, nrestart)
PlanetarySystem *sys;
int nrestart;
{
  int i;
  real rho;
  char filename[256];
  struct reb_particle *particles;
  /* ----- */
  if (CPU_Master) printf ("\n\033[1mRestarting N-body part...\033[0m\n");
  fflush (stdout);
  sprintf (filename, "%snbody.simulation%d.dat", OUTPUTDIR, nrestart);
  struct reb_simulation *rsim = reb_create_simulation_from_binary(filename);
  SetupIntegratorParams (rsim);
  if (CPU_Master) printf ("\033[1mNo problem!\033[0m Function pointers have been reset via the SetupIntegratorParams() function.\n");
  particles = rsim->particles;
  PhysicalTime = rsim->t;
  sys->nb = rsim->N_active - 1;		// get proper 'nb' in case of e.g. previous planetary mergers
  for (i=1; i<rsim->N_active; i++) {
    sys->mass[i-1] = particles[i].m;
    sys->x[i-1] = particles[i].x;
    sys->y[i-1] = particles[i].y;
    sys->z[i-1] = particles[i].z;
    sys->vx[i-1] = particles[i].vx;
    sys->vy[i-1] = particles[i].vy;
    sys->vz[i-1] = particles[i].vz;
    sys->acc[i-1] = ACCRETIONRATE;
    sys->FeelDisk[i-1] = FeelDisk;
  }
  rho = PLANETARYDENSITY/RHO2CGS; 	// get values for updates of planetary radii
  rdenominv = 1.0/(4.0*PI*rho);
  sprintf (filename, "%snbody.orbits.dat", OUTPUTDIR);
  plout = fopenp (filename, "a");     // "a" won't empty the file if it exists
  sprintf (filename, "%snbody.discard.dat", OUTPUTDIR);
  discard = fopenp (filename, "a");
  if (Collisions == YES) {
    sprintf (filename, "%snbody.mergers.dat", OUTPUTDIR);
    mergers = fopenp (filename, "a");
  }
  if (CPU_Master) {
    printf ("%d active planets at restart.\n", rsim->N_active - 1);
    printf ("\n\033[1m...done\n\n\033[0m");
  }
  return rsim;
}
