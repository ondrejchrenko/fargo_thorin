/** \file Psys.c

Contains the functions that set up the planetary system configuration.
In addition, the last two functions allow to track the first planet
(number 0) of the planetary system, in order to perform a calculation
in the frame corotating either with this planet or with its
guiding-center.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

static real Xplanet, Yplanet;

extern boolean GuidingCenter;

int FindNumberOfPlanets (filename)
char *filename;
{
  FILE *input;
  char s[512];
  int Counter=0;
  input = fopen (filename, "r");
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'.\n", filename);
    prs_exit (1);
  }
  while (fgets(s, 510, input) != NULL) {
    if (isalpha(s[0]))
      Counter++;
  }
  fclose (input);
  return Counter;
}

PlanetarySystem *AllocPlanetSystem (nb)	/* #THORIN: 3rd dimension added, acceleration from the disk added */
     int nb;
{
  real *mass, *x, *y, *z, *vx, *vy, *vz, *acc, *ax, *ay, *az;
  boolean *feeldisk, *feelothers;
  int i;
  PlanetarySystem *sys;
  sys  = (PlanetarySystem *)malloc (sizeof(PlanetarySystem));
  if (sys == NULL) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  x    = (real *)malloc (sizeof(real)*(nb+1));
  y    = (real *)malloc (sizeof(real)*(nb+1));
  z    = (real *)malloc (sizeof(real)*(nb+1));
  vx   = (real *)malloc (sizeof(real)*(nb+1));
  vy   = (real *)malloc (sizeof(real)*(nb+1));
  vz   = (real *)malloc (sizeof(real)*(nb+1));
  mass = (real *)malloc (sizeof(real)*(nb+1));
  acc  = (real *)malloc (sizeof(real)*(nb+1));
  if ((x == NULL) || (y == NULL) || (z == NULL) || (vx == NULL) || (vy == NULL) || (vz == NULL) || (acc == NULL) || (mass == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  ax   = (real *)malloc (sizeof(real)*(nb+1));
  ay   = (real *)malloc (sizeof(real)*(nb+1));
  az   = (real *)malloc (sizeof(real)*(nb+1));
  if ((ax == NULL) || (ay == NULL) || (az == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  feeldisk   = (boolean *)malloc (sizeof(real)*(nb+1));
  feelothers = (boolean *)malloc (sizeof(real)*(nb+1));
  if ((feeldisk == NULL) || (feelothers == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  sys->x = x;
  sys->y = y;
  sys->z = z;
  sys->vx= vx;
  sys->vy= vy;
  sys->vz= vz;
  sys->ax= ax;
  sys->ay= ay;
  sys->az= az;
  sys->acc=acc;
  sys->mass = mass;
  sys->FeelDisk = feeldisk;
  sys->FeelOthers = feelothers;
  for (i = 0; i < nb; i++) {
    x[i] = y[i] = z[i] = vx[i] = vy[i] = vz[i] = mass[i] = acc[i] = 0.0;
    ax[i] = ay[i] = az[i] = 0.0;
    feeldisk[i] = feelothers[i] = YES;
  }
  return sys;
}

void FreePlanetary (sys)	/* #THORIN: free also z,vz,ax,ay,az */
PlanetarySystem *sys;
{
  free (sys->x);
  free (sys->vx);
  free (sys->y);
  free (sys->vy);
  free (sys->z);
  free (sys->vz);
  free (sys->ax);
  free (sys->ay);
  free (sys->az);
  free (sys->mass);
  free (sys->acc);
  free (sys->FeelOthers);
  free (sys->FeelDisk);
  free (sys);
}

PlanetarySystem *InitPlanetarySystem (filename)
char *filename;
{
  FILE *input;
  char s[512], nm[512], test1[512], test2[512], *s1;
  PlanetarySystem *sys;
  int i=0, nb;
  float mass, dist, accret;
  boolean feeldis, feelothers;
  nb = FindNumberOfPlanets (filename);
  if (CPU_Master)
    printf ("%d planet(s) found.\n", nb);
  sys = AllocPlanetSystem (nb);
  input = fopen (filename, "r");
  sys->nb = nb;
  while (fgets(s, 510, input) != NULL) {
    sscanf(s, "%s ", nm);
    if (isalpha(s[0])) {
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %s %s", &dist, &mass, &accret, test1, test2);
      sys->mass[i] = (real)mass;
      feeldis = feelothers = YES;
      if (tolower(*test1) == 'n') feeldis = NO;
      if (tolower(*test2) == 'n') feelothers = NO;
      sys->x[i] = (real)dist*(1.0+ECCENTRICITY);
      sys->y[i] = 0.0;
      sys->vy[i] = (real)sqrt((1.0+mass)/dist)*\
	sqrt(1.0-ECCENTRICITY*ECCENTRICITY)/(1.0+ECCENTRICITY);
      sys->vx[i] = -0.0000000001*sys->vy[i];
      sys->acc[i] = accret;
      sys->FeelDisk[i] = feeldis;
      sys->FeelOthers[i] = feelothers;
      i++;
    }
  }
  return sys;
}

void ListPlanets (sys)		/* #THORIN: 3rd dimension added */
PlanetarySystem *sys;
{
  int nb;
  int i;
  nb = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < nb; i++) {
    printf ("Planet number %d\n", i);
    printf ("---------------\n");
    printf ("x = %.10f\ty = %.10f\tz = %.10f\n", sys->x[i],sys->y[i],sys->z[i]);
    printf ("vx = %.10f\tvy = %.10f\tvz = %.10f\n", sys->vx[i],sys->vy[i],sys->vz[i]);
    if (sys->acc[i] == 0.0)
      printf ("Non-accreting.\n");
    else
      printf ("accretion time = %.10f\n", 1.0/(sys->acc[i]));
    if (sys->FeelDisk[i] == YES) {
      printf ("Feels the disk potential\n");
    } else {
      printf ("Doesn't feel the disk potential\n");
    }
    if (sys->FeelOthers[i] == YES) {
      printf ("Feels the other planets potential\n");
    } else {
      printf ("Doesn't feel the other planets potential\n");
    }
    printf ("\n");
  }
}

/** The original function was modified to
 * handle 3D inclined orbits. It is similar
 * to the one used in fargo3D but employs tools
 * from the REBOUND package. */
real GetPsysInfo (sys, action)		/* #THORIN */
PlanetarySystem *sys;
boolean action;
{
  real d1,d2,cross;
  real x,y,z, vx, vy, vz, m, a3;
  real hx,hy,hz,d,Ax,Ay,PerihelionPA;
  real xc, yc;
  struct reb_orbit orbit;
  struct reb_particle p={0};
  static struct reb_particle primary={0};
  /* ----- */
  if (primary.m != 1.0) {
    primary.m = 1.0;
  }
  p.x = xc = x = sys->x[0];
  p.y = yc = y = sys->y[0];
  p.z =      z = sys->z[0];
  p.vx =    vx = sys->vx[0];
  p.vy =    vy = sys->vy[0];
  p.vz =    vz = sys->vz[0];
  p.m  = sys->mass[0];
  m = sys->mass[0]+1.;
  orbit = reb_tools_particle_to_orbit (G, p, primary);
  if (GuidingCenter == YES) {
    if (orbit.e > 1e-8 && orbit.inc > 1e-8) {
 //     if (orbit.e > 1.e-12) {
      hx   = y*vz - z*vy;	/* spec. orbital momentum */
      hy   = z*vx - x*vz;
      hz   = x*vy - y*vx;
      d = sqrt(x*x+y*y+z*z);
      Ax = vy*hz-vz*hy - G*m*x/d;	/* the eccentric vector components */
      Ay = vz*hx-vx*hz - G*m*y/d;
      if (Ax*Ax+Ay*Ay > 0.0) {
        PerihelionPA = atan2(Ay,Ax);	/* angle of peric. w.r.t. the reference direction */
      } else {
        PerihelionPA = atan2(y,x);
      }
      xc = orbit.a*cos(orbit.M+PerihelionPA)*cos(orbit.inc);
      yc = orbit.a*sin(orbit.M+PerihelionPA)*cos(orbit.inc);
    }
  }
  switch (action) {
  case MARK: 
    Xplanet = xc;
    Yplanet = yc;
    return 0.;
    break;
  case GET:
    x = xc;
    y = yc;
    d2 = sqrt(x*x+y*y);
    d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    cross = Xplanet*y-x*Yplanet;
    Xplanet = x;
    Yplanet = y;
    return asin(cross/(d1*d2));
    break;
  case FREQUENCY:
    if (GuidingCenter == YES) {
      a3 = (orbit.a)*(orbit.a)*(orbit.a);
      return sqrt(m/a3);
    } else {
      return (x*vy-y*vx)/(x*x+y*y);
    }
    break;
  }
  return 0.0;
}

/** An analogue of the GetPsysInfo() function;
 * here a rebound simulation structure is used as
 * a call argument. */
real GetPsysInfoFromRsim (rsim, action)		/* #THORIN */
struct reb_simulation *rsim;
boolean action;
{
  struct reb_particle* const particles = rsim->particles;
  real d1,d2,cross;
  real x,y,z, vx, vy, vz, m, a3;
  real hx,hy,hz,d,Ax,Ay,PerihelionPA;
  real xc, yc;
  struct reb_orbit orbit;
  /* ----- */
  xc = x = particles[1].x;	/* first planet has index 1 in REBOUND (but index 0 in FARGO) */
  yc = y = particles[1].y;
  z = particles[1].z;
  vx = particles[1].vx;
  vy = particles[1].vy;
  vz = particles[1].vz;
  m = particles[1].m + particles[0].m;	/* the primary has index 0 in REBOUND (does not exist in FARGO) */
  orbit = reb_tools_particle_to_orbit (G, particles[1], particles[0]);
  if (GuidingCenter == YES) {
    if (orbit.e > 1e-8 && orbit.inc > 1e-8) {
//    if (orbit.e > 1.e-8) {
      hx   = y*vz - z*vy;       /* spec. orbital momentum */
      hy   = z*vx - x*vz;
      hz   = x*vy - y*vx;
      d = sqrt(x*x+y*y+z*z);
      Ax = vy*hz-vz*hy - G*m*x/d;       /* the eccentric vector components */
      Ay = vz*hx-vx*hz - G*m*y/d;
      if (Ax*Ax+Ay*Ay > 0.0) {
        PerihelionPA = atan2(Ay,Ax);    /* angle of peric. w.r.t. the reference direction */
      } else {
        PerihelionPA = atan2(y,x);
      }
      xc = orbit.a*cos(orbit.M+PerihelionPA)*cos(orbit.inc);
      yc = orbit.a*sin(orbit.M+PerihelionPA)*cos(orbit.inc);
    }
  }
  switch (action) {
  case MARK:
    Xplanet = xc;
    Yplanet = yc;
    return 0.;
    break;
  case GET:
    x = xc;
    y = yc;
    d2 = sqrt(x*x+y*y);
    d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    cross = Xplanet*y-x*Yplanet;
    Xplanet = x;
    Yplanet = y;
    return asin(cross/(d1*d2));
    break;
  case FREQUENCY:
    if (GuidingCenter == YES) {
      a3 = (orbit.a)*(orbit.a)*(orbit.a);
      return sqrt(m/a3);
    } else {
      return (x*vy-y*vx)/(x*x+y*y);
    }
    break;
  }
  return 0.0;
}

/** Synchronises the planetary system with the frame.
 * Affects the rebound simulation structure only,
 * 'sys' and 'rsim' are synchronised later on. */
void RotatePsys (rsim, angle)		/* #THORIN */
struct reb_simulation *rsim;
real angle;
{
  int i;
  real sint, cost, xt, yt;
  struct reb_particle* const particles = rsim->particles;
  sint = sin(angle);
  cost = cos(angle);
  for (i=1; i<rsim->N; i++) {
    xt = particles[i].x;
    yt = particles[i].y;
    particles[i].x = xt*cost+yt*sint;
    particles[i].y = -xt*sint+yt*cost;
    xt = particles[i].vx;
    yt = particles[i].vy;
    particles[i].vx = xt*cost+yt*sint;
    particles[i].vy = -xt*sint+yt*cost;
  }
}
