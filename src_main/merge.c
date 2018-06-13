/** \file merge.c

Contains the function that merges the output of different processors.
The resulting merged file is undistinguishable from the file that
would have been produced by a sequential run.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include "fargo.h"

void merge (nb)	/* #THORIN: new outputs */
     int nb;
{
  int i,j;
  char radix[512];
  char command[1024];
  boolean gotonext;
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 3+(AdvecteLabel == YES); j++) {
    switch (j) {
    case 0: strcpy (radix, "dens");
      break;
    case 1: strcpy (radix, "vrad");
      break;
    case 2: strcpy (radix, "vtheta");
      break;
    case 3: strcpy (radix, "label");
      break;
    }
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat",\
	       OUTPUTDIR, radix, nb, i, radix, nb);
      system (command);
    }
    sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",\
	     OUTPUTDIR, radix, nb);
    system (command);
  }
  gotonext = NO;
  for (j = 0; j < 9; j++) {
    switch (j) {
    case 0:
      if (Write_Energy) {
	strcpy (radix, "energy");
      } else {
	gotonext = YES;
      }
      break;
    case 1:
      if (Write_Temperature) {
	strcpy (radix, "temper");
      } else {
	gotonext = YES;
      }
      break;
    case 2:
      if (Write_Divergence) {
	strcpy (radix, "divv");
      } else {
	gotonext = YES;
      }
      break;
    case 3:
      if (Write_Qplus) {
	strcpy (radix, "qplus");
      } else {
	gotonext = YES;
      }
      break;
    case 4:
      if (Write_Qbalance) {
	strcpy (radix, "qbalance");
      } else {
	gotonext = YES;
      }
      break;
    case 5:
      if (Pebbles) {
	strcpy (radix, "pebbledens");
      } else {
	gotonext = YES;
      }
      break;
    case 6:
      if (Pebbles) {
	strcpy (radix, "pebblevrad");
      } else {
	gotonext = YES;
      }
      break;
    case 7:
      if (Pebbles) {
	strcpy (radix, "pebblevtheta");
      } else {
	gotonext = YES;
      }
      break;
    case 8:
      if (Write_Eta) {
	strcpy (radix, "eta");
      } else {
	gotonext = YES;
      }
      break;
    }
    if (gotonext) {
      gotonext = NO;
      continue;
    }
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat",\
	       OUTPUTDIR, radix, nb, i, radix, nb);
      system (command);
    }
    sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",\
	     OUTPUTDIR, radix, nb);
    system (command);
  }
  message ("done\n");
  fflush (stdout);
}
