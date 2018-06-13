/**
 * @file Pframeforce.c
 *
 * @brief Calculates the gravitational interactions
 * between the disk and massive bodies in terms of a
 * vertically averaged potential.
 *
 * @author Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>
 *
 * @details The function FillForcesArrays() computes the
 * gravitational acceleration due to planet-disk and star-disk
 * interactions adopting the outline of Muller & Kley (2012).
 * Works both for the gas and pebbles.
 * The indirect terms are included as well as a
 * deep cubic potential with thickness smoothing (see Eq. (37)
 * in Chrenko et al. 2017). There is also a function which
 * provides the inclination damping for 3D planetary orbits.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Ondřej Chrenko. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

extern boolean AllowAccretion, Corotating, Indirect_Term;
static Pair IndirectTerm;
static real vt_int[MAX1D], vt_cent[MAX1D];

extern boolean DumpTorqueNow, DumpTorqueDensNow;	/* #THORIN: new output control */

				/* Below : work in non-rotating frame */
				/* centered on the primary */
/** Using the vertical averaging procedure of Muller &
 * Kley (2012), calculates the acceleration in planet-disk
 * and star-disk interactions for both gas and pebbles. */
void FillForcesArrays (Rho, sys)	/* #THORIN: redesigned using accelerations instead of the potential */
PolarGrid *Rho;
PlanetarySystem *sys;
{
  int i,j,l,nr,ns,k,npl;
  real x, y;
  real xpl, ypl, zpl, mpl, dpl, invdpl3, smoothing;
  int m, Nvert;
  real d2, d, H, H2, zmax, deltaz;
  real rsm[MAXPLANETS], rsm2[MAXPLANETS], sigma, znode, znode2, denom, ax, ay, ar, at, axindir, ayindir;
  real rst, rst2[MAXPLANETS], rst3[MAXPLANETS], rst4[MAXPLANETS], taper, omegainv, r2, r2chk_1, r2chk_2;
  real integstar, sum, tmp, tmpbuf;
  real cs_avr, integstar_avr, sigma_avr, H_avr;
  real tau_avr=0.0, p_H_avr=0.0, p_integstar_avr=0.0, p_integstar=0.0, pax=0.0, pay=0.0, pH=0.0, Hcorr=0.0;
  Pair IndirectTermFromPlanets;
  boolean IntegrateWithPlanet[MAXPLANETS], IntegrateLocally;
  real axstar, aystar, hillcut[MAXPLANETS]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, mcell;
  real rhill[MAXPLANETS], smoothing2[MAXPLANETS], Hpl[MAXPLANETS], s2[MAXPLANETS], integpl[MAXPLANETS], p_integpl[MAXPLANETS], rmorevert2[MAXPLANETS];
  real tmparray[MAXPLANETS]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, axpl[MAXPLANETS]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, aypl[MAXPLANETS]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  real *agr, *agt, *pagr, *pagt, *rho, *tau, *cs, *abs, *ord, *tqwk;
  agr = GravAccelRad->Field;
  agt = GravAccelTheta->Field;
  rho = Rho->Field;
  cs = SoundSpeed->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  if (Pebbles) {
    pagr = PebbleGravAccelRad->Field;
    pagt = PebbleGravAccelTheta->Field;
    tau = StokesNumber->Field;
  }
  if (TorqueDensity) tqwk = Torque->Field;
  npl = sys->nb;
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) {
    agr[i] = 0.0;
    agt[i] = 0.0;
    if (Pebbles) {
      pagr[i] = 0.0;
      pagt[i] = 0.0;
    }
  }
  // pre-calculate planet distances, FULL hill spheres (with no mass tapering) and indirect term caused by planets accelerating the star
  axstar = 0.0;
  aystar = 0.0;
  for (k = 0; k < npl; k++) {
    axpl[k] = 0.0;
    aypl[k] = 0.0;
    xpl = sys->x[k];
    ypl = sys->y[k];
    zpl = sys->z[k];
    mpl = sys->mass[k];
    dpl = sqrt(xpl*xpl+ypl*ypl+zpl*zpl);
    rhill[k] = dpl*pow((1.0/3.0*mpl),1.0/3.0);
    smoothing = ThicknessSmoothing (xpl,ypl); 	// see Force.c
    Hpl[k] = smoothing/THICKNESSSMOOTHING;	// get the thickness along planet's orbit from the smoothing value 
    smoothing2[k] = smoothing*smoothing;
    rst = 0.5*rhill[k];			// specif. lengths for the potential tapered by cubic spline
    rst2[k] = rst*rst;
    rst3[k] = rst2[k]*rst;
    rst4[k] = rst3[k]*rst;
    rsm[k] = 0.05*rhill[k];		// small smoothing length to avoid diverging denom. near planets
    rsm2[k] = rsm[k]*rsm[k];
    rmorevert2[k] = 2.0*rhill[k];
    rmorevert2[k] *= rmorevert2[k];
    if (Indirect_Term == YES) {
      invdpl3 =  1.0/(dpl*dpl*dpl);       // finally, we look for the indirect term if needed
      mpl *= MassTaper;
      axstar += mpl*invdpl3*xpl;
      aystar += mpl*invdpl3*ypl;
    }
  }
  IndirectTermFromPlanets.x = - axstar;
  IndirectTermFromPlanets.y = - aystar;
  axstar = 0.0;	// reset the protostar's acceleration - in the following, its accel. from the disk will be computed
  aystar = 0.0;
#pragma omp parallel for \
  firstprivate (tau_avr,p_H_avr,p_integstar_avr,p_integstar,pax,pay,pH,Hcorr) \
  private (omegainv, r2, Nvert, cs_avr, H_avr, \
H2, zmax, deltaz, sigma_avr, integstar_avr, znode, znode2, sum, \
d2, d, denom, IntegrateLocally, ax, ay, x, y, mcell, H, \
xpl, ypl, zpl, s2, r2chk_1, r2chk_2, IntegrateWithPlanet, integpl, p_integpl, hillcut, \
sigma, integstar, taper, tmp, ar, mpl, at, \
i,j,k,l,m) \
  reduction (+ : axstar,aystar)
  // 2DO axpl[MAXPLANETS], aypl[MAXPLANETS] should also be included in reduction (this should be possible in OpenMP 4.5)
  // 2DO and atomic clauses should be removed
  for (i = 0; i<nr; i++) {
    omegainv = OmegaInv[i];
    r2 = Rmed2[i];
    Nvert = 10;
    // 1|---> in each ring, estimate the acceleration from the star using the Gaussian profile with azim. average H
    cs_avr = 0.0;
    if (Pebbles) tau_avr = 0.0;
    for (j = 0; j<ns; j++) {
      l = j + i*ns;
      cs_avr += cs[l];
      if (Pebbles) tau_avr += tau[l];
    }
    cs_avr /= (real)ns;
    H_avr = cs_avr*omegainv*SQRT_ADIABIND_INV;
    if (Pebbles) {
      tau_avr /= (real)ns;
      p_H_avr = H_avr*sqrt(PEBBLEALPHA/tau_avr);
      p_integstar_avr = 0.0;
    }
    H2 = H_avr*H_avr;
    zmax = 3.0*H_avr;
    deltaz = zmax/(real)Nvert;
    sigma_avr = 0.0;
    integstar_avr = 0.0;
    for (m=1; m<=Nvert; m++) {		// note for pebbles: H_peb = H*sqrt(alpha/tau), thus also znode_peb = znode*sqrt(alpha/tau), thus sum is the same (alpha/tau cancels out)
      znode = ((real)m-0.5)*deltaz;
      znode2 = znode*znode;
      sum = exp(-(0.5*znode2/H2));
      sigma_avr += sum;
      d2 = r2 + znode2;			// only d2 must be recomputed for pebbles as znode_peb != z_node
      d = sqrt(d2);
      denom = d*d2;
      integstar_avr += sum/denom;
      if (Pebbles) {
	znode2 *= PEBBLEALPHA/tau_avr;
        d2 = r2 + znode2;
	d = sqrt(d2);
	denom = d*d2;
	p_integstar_avr += sum/denom;
      }
    }
    // 1|<---
    // now go through individual cells and decide whether to do the integration using
    // the local H. The criteria for this: a) if H is too different compared to H_avr
    // b) if the cell is expected to be strongly influenced by one of the planets
    for (j = 0; j<ns; j++) {
      IntegrateLocally = NO;
      Nvert = 10;
      ax = 0.0;
      ay = 0.0;
      if (Pebbles) {
        pax = 0.0;
        pay = 0.0;
      }
      l = j + i*ns;
      x = abs[l];
      y = ord[l];
      mcell = rho[l]*Surf[i];
      H = cs[l]*omegainv*SQRT_ADIABIND_INV;
      H2 = H*H;
      // 2|---> check for H vs H_avrg
      if (fabs(H-H_avr)/H_avr > 0.05) IntegrateLocally = YES;
      if (Pebbles) {
	Hcorr = sqrt(PEBBLEALPHA/tau[l]);
        pH = H*Hcorr;
	if (fabs(pH-p_H_avr)/p_H_avr > 0.05) IntegrateLocally = YES;
      }
      // 2|<---
      // 3|---> check for the distance w.r.t. planet and thickness related to
      // standard smoothing- Calculate and save cell-planet distances and (if needed) taper
      // to exclude the hill sphere from the torque evaluation
      for (k = 0; k < npl; k++) {
        xpl = sys->x[k];
        ypl = sys->y[k];
        zpl = sys->z[k];
        s2[k] = (x-xpl)*(x-xpl) + (y-ypl)*(y-ypl);
	r2chk_1 = 10.0*Hpl[k];		// 2DO will be moved to parameters
	r2chk_1 *= r2chk_1;
        r2chk_2 = 5.0*rhill[k];		// 2DO will be moved to parameters
        r2chk_2 *= r2chk_2;
        if (MassTaper == 1.0 && (s2[k] < r2chk_1 || s2[k] < r2chk_2)) {
          IntegrateLocally = YES;
          IntegrateWithPlanet[k] = YES;
	  integpl[k] = 0.0;
          if (Pebbles) p_integpl[k] = 0.0;
          d2 = s2[k] + zpl*zpl;
          if (ExcludeHill) {
            hillcut[k] = 1.0/(exp(-(sqrt(d2)/rhill[k]-HILLCUT)/HILLCUT*10.0)+1.0);
          }
        } else {
          IntegrateWithPlanet[k] = NO;
        }
        if (s2[k] < rmorevert2[k]) Nvert = 30;
      }
      // 3|<---
      // 4|---> Perform the local integration if needed
      if (IntegrateLocally) {
        zmax = 3.0*H;
        deltaz = zmax/(real)Nvert;
        sigma = 0.0;
        integstar = 0.0;
        if (Pebbles) p_integstar = 0.0;
        for (m=1; m<=Nvert; m++) {
          znode = ((real)m-0.5)*deltaz;
          znode2 = znode*znode;
          sum = exp(-(0.5*znode2/H2));
          sigma += sum;
          d2 = r2 + znode2;   // part for the star
          d = sqrt(d2);
          denom = d*d2;
          integstar += sum/denom;
          if (Pebbles) {
            d2 = r2 + znode2*Hcorr*Hcorr;
            d = sqrt(d2);
            denom = d*d2;
            p_integstar += sum/denom;
          }
          for (k = 0; k < npl; k++) {	// part for the planets
            if (IntegrateWithPlanet[k]==YES) {
              zpl = sys->z[k];
              d2 = s2[k] + (znode-zpl)*(znode-zpl);
              if (d2 < rst2[k]) {
                d = sqrt(d2);
                taper = -3.0*d2*d2/rst4[k] + 4.0*d*d2/rst3[k];
              } else {
                taper = 1.0;
              }	      
	      d2 += rsm2[k];
              d = sqrt(d2);
              denom = d*d2;
              integpl[k] += sum/denom*taper;		// 2DO condition below should again be parametric
              if (zpl < - 0.05*H || zpl > 0.05*H) {	// For inclined planets, one cannot expect the integral to be symmetric w.r.t. the midplane. 
                d2 = s2[k] + (-znode-zpl)*(-znode-zpl);	// thus integrate also the "mirror" nodes (with -znode vertical coordinate)
                if (d2 < rst2[k]) {
                  d = sqrt(d2);
                  taper = -3.0*d2*d2/rst4[k] + 4.0*d*d2/rst3[k];
                } else {
                  taper = 1.0;
                }
                d2 += rsm2[k];
                d = sqrt(d2);
                denom = d*d2;
                integpl[k] += sum/denom*taper;
              }
              if (Pebbles) {
                d2 = s2[k] + (znode*Hcorr-zpl)*(znode*Hcorr-zpl);
                if (d2 < rst2[k]) {
                  d = sqrt(d2);
                  taper = -3.0*d2*d2/rst4[k] + 4.0*d*d2/rst3[k];
                } else {
                  taper = 1.0;
                }
                d2 += rsm2[k];
                d = sqrt(d2);
                denom = d*d2;
                p_integpl[k] += sum/denom*taper;
                if (zpl < - 0.05*pH || zpl > 0.05*pH) {     // 2DO condition again should be parametric
                  d2 = s2[k] + (-znode*Hcorr-zpl)*(-znode*Hcorr-zpl); // integrate also the "mirror" nodes (with -znode vertical coordinate)
                  if (d2 < rst2[k]) {
                    d = sqrt(d2);
                    taper = -3.0*d2*d2/rst4[k] + 4.0*d*d2/rst3[k];
                  } else {
                    taper = 1.0;
                  }
                  d2 += rsm2[k];
                  d = sqrt(d2);
                  denom = d*d2;
                  p_integpl[k] += sum/denom*taper;
                }
              }
            }
          }
        } // 4|<--- 
      } else {	// 5|---> If local integration is not necessary, use the orbital average
        sigma = sigma_avr;
        integstar = integstar_avr;
        if (Pebbles) p_integstar = p_integstar_avr;
      }	// 5|<---
      // 6|---> cell-star interaction update
      tmp = integstar/sigma;
      ar = - Rmed[i]*tmp;	// acceleration from the star is always radial -> can be updated directly
      agr[l] += ar;
      if (i>=Zero_or_active && i<Max_or_active) {
        axstar += mcell*tmp*x;
        aystar += mcell*tmp*y;
      }
      if (Pebbles) {
        tmp = p_integstar/sigma;
        ar = - Rmed[i]*tmp;
        pagr[l] += ar;
      }
      // 6|<---
      // 7|---> cell-planet interaction update
      for (k = 0; k < npl; k++) {
        xpl = sys->x[k];
        ypl = sys->y[k];
        zpl = sys->z[k];
        mpl = sys->mass[k];
        if (IntegrateWithPlanet[k]==YES) {	// case when local integration was done
          if (zpl < - 0.05*H || zpl > 0.05*H) {
            tmp = (x-xpl)*integpl[k]/(2.0*sigma);	// here integpl is taken from -zmax to zmax, whereas sigma only from 0 to zmax -> must multiply by 2
          } else {
            tmp = (x-xpl)*integpl[k]/sigma;
          }
          ax += - mpl*tmp;
          if (i>=Zero_or_active && i<Max_or_active) {
            if (ExcludeHill) tmp *= hillcut[k];
#pragma omp atomic
            axpl[k] += mcell*tmp;
	    if (TorqueDensity && GETTORQUEFORPLANET==k) tqwk[l] = -ypl*mcell*tmp;
          }
          if (zpl < - 0.05*H || zpl > 0.05*H) {
            tmp = (y-ypl)*integpl[k]/(2.0*sigma);
          } else {
            tmp = (y-ypl)*integpl[k]/sigma;
          } 
          ay += - mpl*tmp;
          if (i>=Zero_or_active && i<Max_or_active) {
            if (ExcludeHill) tmp *= hillcut[k];
#pragma omp atomic
            aypl[k] += mcell*tmp; 
	    if (TorqueDensity && GETTORQUEFORPLANET==k) tqwk[l] += xpl*mcell*tmp;
          }
          if (Pebbles) {
            if (zpl < - 0.05*pH || zpl > 0.05*pH) {
              tmp = (x-xpl)*p_integpl[k]/(2.0*sigma);
            } else {
              tmp = (x-xpl)*p_integpl[k]/sigma;
            }
            pax += -mpl*tmp;
            if (zpl < - 0.05*pH || zpl > 0.05*pH) {
              tmp = (y-ypl)*p_integpl[k]/(2.0*sigma);
            } else {
              tmp = (y-ypl)*p_integpl[k]/sigma;
            }
            pay += -mpl*tmp;
          }          
        } else {  // simple point-mass acceleration otherwise
          mpl *= MassTaper;
          d2 = s2[k] + zpl*zpl + smoothing2[k];
          d = sqrt(d2);
          denom = d*d2;
          tmp = (x-xpl)/denom;
          ax += - mpl*tmp;
          if (i>=Zero_or_active && i<Max_or_active) {
            if (ExcludeHill) tmp *= hillcut[k];
#pragma omp atomic
            axpl[k] += mcell*tmp;
	    if (TorqueDensity && GETTORQUEFORPLANET==k) tqwk[l] = -ypl*mcell*tmp;
          }
          tmp = (y-ypl)/denom;
          ay += - mpl*tmp;
          if (i>=Zero_or_active && i<Max_or_active) {
            if (ExcludeHill) tmp *= hillcut[k];
#pragma omp atomic
            aypl[k] += mcell*tmp;
	    if (TorqueDensity && GETTORQUEFORPLANET==k) tqwk[l] += xpl*mcell*tmp;
          }
          if (Pebbles) {
            d2 += - smoothing2[k] + smoothing2[k]*Hcorr*Hcorr;
            d = sqrt(d2);
            denom = d*d2;
            tmp = (x-xpl)/denom;
            pax += -mpl*tmp;
            tmp = (y-ypl)/denom;
            pay += -mpl*tmp;
          }
        }
      }
      // 7|<---
      ar = (x*ax + y*ay)*InvRmed[i];	// transform ax,ay to ar,at and do the update
      at = (x*ay - y*ax)*InvRmed[i];
      agr[l] += ar;
      agt[l] += at;
      if (Pebbles) {
        ar = (x*pax + y*pay)*InvRmed[i];
        at = (x*pay - y*pax)*InvRmed[i];
        pagr[l] += ar;
        pagt[l] += at;
      }
    }
  }
  // MPI reductions
  MPI_Allreduce (&axstar, &tmpbuf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  IndirectTerm.x = -tmpbuf;
  MPI_Allreduce (&aystar, &tmpbuf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  IndirectTerm.y = -tmpbuf;
  MPI_Allreduce (&axpl, &tmparray, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (k=0; k<npl; k++) {
    sys->ax[k] = tmparray[k];
  }
  MPI_Allreduce (&aypl, &tmparray, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (k=0; k<npl; k++) {
    sys->ay[k] = tmparray[k];
  }
  // correct the disk acceleration using the indirect terms
  axindir = IndirectTerm.x;
  ayindir = IndirectTerm.y;
  if (Indirect_Term == YES) {
    axindir += IndirectTermFromPlanets.x;
    ayindir += IndirectTermFromPlanets.y;
  }
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      x = abs[l];
      y = ord[l];
      ar = (x*axindir + y*ayindir)*InvRmed[i];
      at = (x*ayindir - y*axindir)*InvRmed[i];
      agr[l] += ar;
      agt[l] += at;
      if (Pebbles) {
        pagr[l] += ar;
        pagt[l] += at;
      }
    }
  }
}

/** Updates the planet velocities due to disk forces.
 * Also applies the inclination damping. */
void AdvanceSystemFromDisk (Rho, sys, dt)	/* #THORIN uses directly the accelerations calculated above */
PlanetarySystem *sys;
PolarGrid *Rho;
real dt;		       
{
  int npl, k, i;
  real x,y,z,m,vz,damp;
  char command[1024];
  npl = sys->nb;
  if (DumpTorqueNow) {
    DumpTorqueNow=NO;
    UpdateLog (Rho, sys);
  }
  if (DumpTorqueDensNow) {
    DumpTorqueDensNow=NO;
    WriteDiskPolar (Torque, TimeStep);
    MPI_Barrier (MPI_COMM_WORLD);
    if (Merge && (CPU_Number > 1) && CPU_Master) {
      for (i=1; i<CPU_Number; i++) {
        sprintf (command, "cd %s; cat gastorque%d.dat.%05d >> gastorque%d.dat",\
                 OUTPUTDIR, TimeStep, i, TimeStep);
        system (command);
      }
      sprintf (command, "cd %s; rm -f gastorque%d.dat.0*", OUTPUTDIR, TimeStep);
      system (command);
    }
  }
  for (k = 0; k < npl; k++) {
    if (sys->FeelDisk[k] == YES) {
      m=sys->mass[k];
      x=sys->x[k];
      y=sys->y[k];
      z=sys->z[k];
      vz = sys->vz[k];
      damp = DampingTW04 (Rho, m, x, y, z, vz);
      sys->vx[k] += dt * sys->ax[k];
      sys->vy[k] += dt * sys->ay[k];
      sys->vz[k] += dt * (sys->az[k] + damp);
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

/** Artificial vertical force to damp the orbital
 * inclinations using the Tanaka & Ward (2004)
 * prescription (see also Morbidelli et al. 2007, Pierens & Nelson 2008) */
real DampingTW04 (Rho, m, x, y, z, vz)		/* #THORIN */
PolarGrid *Rho;
real m,x,y,z,vz;
{
  real r, ang, Omega, cs4inv, damp=0.0, rhoavr=0.0, csavr=0.0, tmp;
  real *rho, *cs;
  int i, j, l, ns;
  rho = Rho->Field;
  ns = Rho->Nsec;
  cs = SoundSpeed->Field;
  r = sqrt(x*x + y*y);          /* find position in the midplane */
  if ( (VERTICALDAMPING > 0.0) && (z>0.0) && (r >= Rinf[Zero_or_active] && r <= Rsup[Max_or_active-1]) ) {
     /* condition satisfied only for a CPU which contains the planet within its active zone */
    ang = atan2(y,x);           /* get the damping value on this cpu */
    if (ang < 0.0) ang += 2.0*PI;
    i = Zero_or_active;
    while (Rsup[i] < r) i++;
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      rhoavr += rho[l];
      csavr += cs[l];
    }
    csavr /= (real)ns;
    rhoavr /= (real)ns;
    cs4inv = pow(csavr, -4.0);		/* 2DO isothermal case ? */
    Omega = pow(r, -1.5);               /* Keplerian angular velocity at planets' position */
    damp = VERTICALDAMPING*m*rhoavr*Omega*cs4inv*(-2.176*vz - 0.871*z*Omega);   /* Tanaka&Ward 04 formula */
  }
  if (CPU_Number==1) return damp;       /* we are done in case of a serial build */
  MPI_Allreduce (&damp, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	/* other CPUs have damp=0.0, the sum reduction distributes the result among them */
  damp = tmp;
  return tmp;
}

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}

/** Part of the initialisation */
void InitGasDensityEnergy (Rho, Energy)		/* #THORIN */
PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *rho, *energy;
  /* ----- */
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      rho[l] = SigmaMed[i];	/* SigmaMed[] was initialized through FillSigma()*/
    }
  }
  energy = Energy->Field;
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      energy[l] = EnergyMed[i]; /* EnergyMed[] initialized in FillEnergy(); Theo.c*/
    }
  }
}

void InitGasVelocity (Vr, Vt)		/* #THORIN */
PolarGrid *Vr, *Vt;
{
  int i, j, l, nr, ns;
  real *vr, *vt, *pres, *cs;
  real  r, rg, omega, ri;
  real viscosity, t1, t2, r1, r2;
  /* ----- */
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Vr->Nrad;
  ns  = Vr->Nsec;
  cs  = SoundSpeed->Field;
  pres= Pressure->Field;
  /* loading a file with prescribed soundspeed  profile was possible here,
   * now it is redundant (and thus discarded) */
  mpi_make1Dprofile (pres, globpressvec);/* extract the HD quantity (from all CPUs
					   if necessary) and arrange it into a
					   single 1D array (globpressvec) describing
					   all radial rings by quantities averaged
					   over each ring. */
  /* This part is always calculated, but it is stored in vtheta only if CentrifugalBalance==YES.
     ---->
     !!! 2DO This option has not been tested yet. */
  for (i = 1; i < GLOBALNRAD; i++) {
    /* The equation below is the same as in the original FARGO besides the first
     * term, which is now more general (the pressure difference between two
     * neighbouring rings). */
    vt_int[i]=(globpressvec[i] - globpressvec[i-1])/    \
      (.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1])+    \
      G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    /* The rest of the routine is unchanged, only the density initialization
     * was discarded (it is already done before this routine) */
    vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
  }
  t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
  r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
  t2 = vt_cent[0];
  r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  t1 = t1-r1/(r2-r1)*(t2-t1);
  vt_cent[0] = t1;
  ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[GLOBALNRAD]=vt_cent[GLOBALNRAD-1];
  /* <----- */
  if (ParametricCooling) {
    FillCoolingTime ();	/* for both functions see Theo.c */
    FillQplus ();	/* ! This function calls FViscosity() which might need the azimuth-averaged values of sound speed
			     which are stored in 'globcsvec[]'. It is required that ComputeSoundSpeed() is called before InitGasVelocity(). */
  }
  /* "Usual" velocity initialization */
  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    viscosity = FViscosity (r);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rg = r;
      omega = sqrt(G*1.0/rg/rg/rg);
      /* this one is original: 
      vt[l] = omega*r*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(r,2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX)); */
      /* the initial azimuthal velocity profile is refined using a simple relationship, see e.g. Masset 2002 */
      vt[l] = omega*r*sqrt(1.0-cs[l]*cs[l]*r/ADIABIND);
      /* vt[l] = omega*r*sqrt(1.0-pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX)); */
      vt[l] -= OmegaFrame*r;
      if (CentrifugalBalance)
	vt[l] = vt_cent[i+IMIN];
      if (i == nr) 
	vr[l] = 0.0;
      else {
	vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
	if (ViscosityAlpha) {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
	} else {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
	}
      }
    }
  }
  for (j = 0; j < ns; j++)
    vr[j] = vr[j+ns*nr] = 0.0;
}
/* <---- */
