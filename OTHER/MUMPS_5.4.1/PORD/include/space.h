/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: space.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file includes all necessary header files
/
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#ifndef _WIN32
#include <sys/times.h>
#endif
#if defined(__MINGW32__)
#include<sys/time.h>
#endif
#include <math.h>

#ifdef PARIX
#ifdef __EPX
#include <epx/root.h>
#include <epx/link.h>
#include <epx/types.h>
#include <epx/sem.h>
#include <epx/thread.h>
#include <epx/time.h>
#include <epx/logerror.h>
#else
#include <sys/root.h>
#include <sys/link.h>
#include <sys/types.h>
#include <sys/sem.h>
#include <sys/thread.h>
#include <sys/time.h>
#include <sys/logerror.h>
#endif
#include <signal.h>
#endif

#ifdef MPI
#include "mpi.h"
#endif

#include "const.h"
#include "params.h"
#include "macros.h"
#include "types.h"
#include "protos.h"
#include "eval.h"

#define FORTRAN(nu,nl,pl,pc)                     \
void nu ();                                      \
void nl pl                                       \
{ nu pc; }                                       \
void nl##_ pl                                    \
{ nu pc; }                                       \
void nl##__ pl                                   \
{ nu pc; }                                       \
void nu pl
