#ifndef TOYOPTFUNCTION
#define TOYOPTFUNCTION

#include "snoptanalysis.h"
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

  int toyOptusrf_( integer    *Status, integer *n,    doublereal x[],
		integer    *needF,  integer *neF,  doublereal F[],
		integer    *needG,  integer *neG,  doublereal G[],
		char       *cu,     integer *lencu,
		integer    iu[],    integer *leniu,
		doublereal ru[],    integer *lenru );

  int toyOptusrfg_( integer    *Status, integer *n,    doublereal x[],
		 integer    *needF,  integer *neF,  doublereal F[],
		 integer    *needG,  integer *neG,  doublereal G[],
		 char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru);

#ifdef __cplusplus
}
#endif


#endif
