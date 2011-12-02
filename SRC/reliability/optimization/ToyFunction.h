#ifndef TOYFUNCTION
#define TOYFUNCTION

#include "SnoptProblem.h"

#ifdef __cplusplus
extern "C" {
#endif

  int toyusrf_( integer    *Status, integer *n,    doublereal x[],
		integer    *needF,  integer *neF,  doublereal F[],
		integer    *needG,  integer *neG,  doublereal G[],
		char       *cu,     integer *lencu,
		integer    iu[],    integer *leniu,
		doublereal ru[],    integer *lenru );

  int toyusrfg_( integer    *Status, integer *n,    doublereal x[],
		 integer    *needF,  integer *neF,  doublereal F[],
		 integer    *needG,  integer *neG,  doublereal G[],
		 char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru);

#ifdef __cplusplus
}
#endif


#endif
