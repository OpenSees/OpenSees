/*
 *
 *  This file is part of MUMPS 5.4.0, released
 *  on Tue Apr 13 15:26:30 UTC 2021
 *
 *
 *  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#if defined(_WIN32)
#include  "elapse.h"
#include  <time.h>
#include  <sys/timeb.h>
void MUMPS_CALL mumps_elapse(double *val)
{
  time_t	ltime;
  struct    _timeb	tstruct;

  time (&ltime);
  _ftime(&tstruct);
  *val = (double) ltime + (double) tstruct.millitm*(0.001);
}

#else

#include "elapse.h"
#include <sys/time.h>
void mumps_elapse(double *val)
  {
    struct timeval time;
    gettimeofday(&time,(struct timezone *)0);
    *val=time.tv_sec+time.tv_usec*1.e-6;
  }
#endif
