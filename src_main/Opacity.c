/**
 * @file Opacity.c
 *
 * @brief Contains various opacity laws. 
 *
 * @author Miroslav Brož <mira@sirrah.troja.mff.cuni.cz>
 *
 * @details Opacity laws from various papers merged into
 * a single file. Uses Lin & Papaloizou (1985), Bell & Lin (1994)
 * and Zhu et al. (2012) opacities. The transitions between different
 * regimes are refined to avoid jumps.
 *
 * @section     LICENSE
 * Copyright (c) 2018 Miroslav Brož. See the LICENSE file of the
 * distribution.
 *
 */

#include "fargo.h"

/* kappa_i rho^a_i T^b_i = kappa_j rho^a_j T^b_j */
/* T = (kappa_j/kappa_i rho^(a_j-a_i))^(1/(b_i-b_j)) */
#define TEMP(i,j) pow(kappa[j]/kappa[i] * pow(rho,a[j]-a[i]), 1.0/(b[i]-b[j]))
#define KAPPAFLOOR 1.0e-12


/** Constant opacity function */
real opacity_const (real rho, real T)
{
  return PARAMETRICOPACITY;
}

/** Bell & Lin (1994) opacities. Input values
 * must be in cgs. */
real opacity_BL94 (real rho, real T)
{
  /* 0 ... ice grains */
  /* 1 ... ice evaporation */
  /* 2 ... metal grains */
  /* 3 ... metal evaporation */
  /* 4 ... molecules */
  /* 5 ... H scattering */
  /* 6 ... bound-free, free-free */
  /* 7 ... e- scattering */

  static real kappa[] = {2.0e-4, 2.0e16, 0.1, 2.0e81, 1.0e-8, 1.0e-36, 1.5e20, 0.348};
  static real a[] = {0.0, 0.0, 0.0, 1.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0};
  static real b[] = {2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0};

  real kappa_, kappa1, kappa2, kappa3, tmp1, tmp2;

  kappa_ = KAPPAFLOOR;
  if (T <= TEMP(2,3)) {
    kappa1 = kappa[0] * pow(T,b[0]);  // a[0] = 0.0
    kappa2 = kappa[1] * pow(T,b[1]);  // a[1] = 0.0
    kappa3 = kappa[2] * pow(T,b[2]);  // a[2] = 0.0

    tmp1 = 1.0/(kappa1*kappa1) + 1.0/(kappa2*kappa2);
    tmp1 = 1.0/(tmp1*tmp1);
    tmp2 = pow(kappa3/(1.0 + 1.0e22*pow(T,-10.0)),4.0);
    kappa_ = pow(tmp1+tmp2,0.25);
  }
  else if ((T > TEMP(2,3)) && (T <= TEMP(3,4))) {
    kappa_ = kappa[3] * pow(rho,a[3]) * pow(T,b[3]);
  }
  else if ((T > TEMP(3,4)) && (T <= TEMP(4,5))) {
    kappa_ = kappa[4] * pow(rho,a[4]) * pow(T,b[4]);
  }
  else if ((T > TEMP(4,5)) && (T <= TEMP(5,6))) {
    kappa_ = kappa[5] * pow(rho,a[5]) * pow(T,b[5]);
  }
  else if ((T > TEMP(5,6)) && (T <= TEMP(6,7))) {
    kappa_ = kappa[6] * pow(rho,a[6]) * pow(T,b[6]);
  }
  else {
    kappa_ = kappa[7];
  }
  return kappa_;
}

/** Zhu et al. (1994) opacities. Taken from
 * http://www.physics.unlv.edu/~zhzhu/Opacity.html.
 * Note: There were JUMPS at several transitions,
 * because both cofficients and borders were rounded
 * (corrected for 3->4, 4->5, 5->6, 6->7)!
 * Input values must be in cgs. */
real opacity_ZHU12 (real rho, real T) {

  real xlop, xlp, xlt, rosstabd, pre;

  pre = R_STANDARD*1.0e7/MMW * rho * T;
  xlp = log10(pre);
  xlt = log10(T);

  if(xlt < 2.23567+0.01899*(xlp-5.)){
    xlop = 1.5*(xlt-1.16331)-0.736364;					/* 1 ... water ice grains */
  }
  else if (xlt < 2.30713+0.01899*(xlp-5.)){
    xlop = -3.53154212*xlt+8.767726-(7.24786-8.767726)*(xlp-5.)/16.;	/* 2 ... water ice evaporation */
  }
  else if (xlt < (17.7-0.62+1.5*2.30713)/(1.5+5.832)){
    xlop = 1.5*(xlt-2.30713)+0.62 ;					/* 3 ... non-water grains */
  }
  else if  (xlt < (5.9398+17.7)/(5.832+2.129)){
    xlop = -5.832*xlt+17.7;						/* 4 ... graphite corrosion */
  }
  else if (xlt < (129.88071+5.9398 + (142.996475-129.88071)*0.1*(xlp+4.))/(2.129+42.98075)){
    xlop = 2.129*xlt-5.9398;						/* 5 ... no-graphite grains */
  }
  else if (xlt < (129.88071+15.0125 + (142.996475-129.88071)*0.1*(xlp+4.))/(4.0625+42.98075)){
    xlop = 129.88071-42.98075*xlt+(142.996475-129.88071)*0.1*(xlp+4.0);	/* 6 ... grain evaporation */
  }
  else if (xlt < 3.28+xlp/4.*0.12){
    xlop = -15.0125+4.0625*xlt;						/* 7 ... water vapour */
  }
  else if (xlt < 3.41+0.03328*xlp/4.){
    xlop = 58.9294-18.4808*xlt+(61.6346-58.9294)*xlp/4.;		/* 8 ... water dissociation */
  }
  else if (xlt < 3.76+(xlp-4)/2.*0.03){
    xlop = -12.002+2.90477*xlt+(xlp-4)/4.*(13.9953-12.002);		/* 9 ... molecules */
  }
  else if (xlt < 4.07+(xlp-4)/2.*0.08){
    xlop = -39.4077+10.1935*xlt+(xlp-4)/2.*(40.1719-39.4077);		/* 10 ... H scattering */
  }
  else if (xlt < 5.3715+(xlp-6)/2.*0.5594){
    xlop = 17.5935-3.3647*xlt+(xlp-6)/2.*(17.5935-15.7376);		/* 11 ... bound-free, free-free */
  }
  else{
    xlop = -0.48;							/* 12 ... e- scattering */
  }
  if (xlop < 3.586*xlt-16.85&&xlt<4.){
    xlop = 3.586*xlt-16.85;						/* 13 ... molecules and H scattering */
  }

  rosstabd = pow(10.,xlop);

  return(rosstabd);
}

/** Lin & Papaloizou (1985) opacities.
 * Input values must be in cgs. */
real opacity_LP85 (real rho, real T)
{
  real kappa;

  if (T < pow(2.0e-4/2.0e16,-1.0/9.0)) {
    kappa = 2.0e-4 * T*T;				/* ice grains */
  }
  else if (T < pow(2.0e16/5.0e-3,1.0/8.0)) {
    kappa = 2.0e16 * pow(T,-7.0);			/* ice evaporation */
  }
  else if (T < pow(5.0e-3/2.0e34,-1.0/10.0) * pow(rho,1.0/15.0)) {
    kappa = 5.0e-3 * T;					/* silicate grains */
  }
  else if (T < pow(2.0e34/2.0e-8,1.0/12.0)) {
    kappa = 2.0e34 * pow(rho,2.0/3.0) * pow(T,-9.0);	/* silicate evaporation */
  }
  else if (T < pow(2.0e-8/1.0e-36,1.0/7.0) * pow(rho,1.0/21.0)) {
    kappa = 2.0e-8 * pow(rho,2.0/3.0) * T*T*T;		/* molecules */
  }
  else if (T < pow(1.0e-36/1.5e20,-1.0/12.5) * pow(rho,4.0/75.0)) {
    kappa = 1.0e-36 * pow(rho,1.0/3.0) * pow(T,10.0);	/* H scattering */
  }
  else {
    kappa = 1.5e20 * rho * pow(T,-5.0/2.0);		/* bound-free, free-free */
  }
  return kappa;
}

/** Fills the opacity polar grid */
void OpacityProfile (VolumeDensity, Temperature, Opacity)
PolarGrid *VolumeDensity, *Temperature, *Opacity;
{
  int nr, ns, i, j, l;
  real *opacity, *temper, *voldens;
  real tmp;
  static boolean FillParametricOpacity = YES;

  temper = Temperature->Field;
  nr = Temperature->Nrad;
  ns = Temperature->Nsec;
  voldens = VolumeDensity->Field;
  opacity = Opacity->Field;

  if (PARAMETRICOPACITY > 0.0) {
    if (!FillParametricOpacity) return;
    tmp = PARAMETRICOPACITY*OPA2CU;
    for (i=0; i<nr; i++) {
      for (j=0; j<ns; j++) {
        l = j + i*ns;
        opacity[l] = tmp;
      }
    }
    FillParametricOpacity = NO;
    return;
  }

#pragma omp parallel for default(none) shared(nr,ns,temper,voldens,opacity) private(i,j,l)
  for (i=0; i<nr; i++) {
    for (j=0; j<ns; j++) {
      l = j + i*ns;
      opacity[l] = opacity_func(voldens[l]*RHO2CGS, temper[l]*T2SI) * OPA2CU;
    }
  }
}
