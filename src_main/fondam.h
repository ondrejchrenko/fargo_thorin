/** \file fondam.h

Contains fondamental constants used thorough the code.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#define		G	1.0
#define	       PI	3.14159265358979323844
#define   CPUOVERLAP    5	/* Zeus-like overlap kernel. 2:transport; 2: source, 1:viscous stress  */

/* #THORIN -----> */

#define	GASCONST	1.0
#define MOLWEIGHT	1.0

#define SQRT2PI_INV (1.0/sqrt(2.0*PI))

/* constant values to resolve the opacity table from Bell & Lin (1994), Lin & Papaloizou (1985) respectively,
   and also to properly transform Stefan-Boltzman constant */
#define fac1o15 (1.0/15.0)
#define fac1o21 (1.0/21.0)
#define fac4o75 (4.0/75.0)
#define MSOL_SI 1.98855e30
#define G_SI 6.674e-11
#define GM_SI 1.32712440018e20
#define AU_SI 149597870700.0    /* changing this can modify basic length unit */
#define R_STANDARD 8.3144598
#define MMW 2.4
#define R_SI (R_STANDARD/(MMW*0.001))   /* R_specific = R/mu = 8.3144598 (standard gas constant) / (2.4 (fiducial mol.weight in PPDs) * 1g/mol) = 8.3144598/(2.4*0.001) in SI */
                                        /* changing this means assuming different molar weight */
#define OPA2CU (MSOL_SI*1000.0/pow(AU_SI*100.0,2.0))    /* (kappa in cgs cm^2/g) * OPA2CU = (kappa in code units) */
#define T2SI (GM_SI/R_SI/AU_SI) /* (T in code units) * T2K = (T in Kelvins ~ cgs or SI) */
#define TIME2SI (sqrt(pow(AU_SI,3.0)/GM_SI))    /* (time in code units) * TIME2S = (time in seconds ~ cgs or SI) */
#define RHO2CGS (MSOL_SI*1000.0/pow(AU_SI*100.0,3.0))   /* (VOLUME!!! density in code units) * RHO2CGS = (volume density in g/cm^3 cgs) */
#define STEFANBOLTZMANN (5.670367e-8*(1.0/MSOL_SI)*pow(1.0/TIME2SI,-3.0)*pow(1.0/T2SI,-4.0))
/* below: definition of Stefan constant similar to FARGO3D, using cgs, can be used to check the definition above */
/* old one #define STEFANBOLTZMANNCONTROL (5.670367e-8*pow(1.0/MSOL_SI,-1.5)*pow(1.0/G_SI,-2.5)*pow(1.0/AU_SI,-0.5)*pow(1.0/R_SI,4.0))*/
#define STEFANBOLTZMANNCONTROL (5.670367e-5*pow(1.0/(R_SI*10000.0),4.0)*pow(1.0/(G_SI*1000.0),-2.5)*pow(1.0/(MSOL_SI*1000),-1.5)*pow(1.0/(AU_SI*100.0),-0.5))

/* and some stuff related to the drag coefficient evaluation and to pebble accretion */
#define MOLDIAMETER 2.72        	/* value in angstroms */
#define ANGSTR_CGS 0.00000001   	/* angstrom in centimeters */
#define AMU_CGS 1.660538921e-24     	/* atomic mass unit in grams */
#define MOLCROSSEC_CGS 2.0e-15		/* molecular cross section of H2 in cm^2 */
#define SURFDENS2CGS (MSOL_SI*1000.0/pow(AU_SI*100.0,2.0))	/* (surface density in code units) * SURFDENS2CGS = (surface density in g/cm^2 cgs) */

/* other conversions */
#define SIGMA2CGS (MSOL_SI*1000.0/pow(AU_SI*100.0,2.0))
#define PRESS2CGS (MSOL_SI*1000.0/(AU_SI*100.0*TIME2SI*TIME2SI))
#define FLUX2CU (0.000003/(2.0*PI))	/* (radial flux in earth masses per year) * FLUX2CU = (radial flux in code units) */
/* #THORIN <----- */
