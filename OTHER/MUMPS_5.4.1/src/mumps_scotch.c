/*
 *
 *  This file is part of MUMPS 5.4.1, released
 *  on Tue Aug  3 09:49:43 UTC 2021
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
/* Interfacing with SCOTCH and pt-SCOTCH */
#include <stdio.h>
#include "mumps_scotch.h"
#if defined(scotch) || defined(ptscotch)
void MUMPS_CALL
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
                    MUMPS_INT * const  petab,
              const MUMPS_INT * const  pfree,
                    MUMPS_INT * const  lentab,
                    MUMPS_INT * const  iwtab,
                    MUMPS_INT * const  nvtab,
                    MUMPS_INT * const  elentab,
                    MUMPS_INT * const  lasttab,
                    MUMPS_INT * const  ncmpa,
                    MUMPS_INT * const  weightused,
                    MUMPS_INT * const  weightrequested )
{
/* weightused(out) = 1 if weight of nodes provided in nvtab are used (esmumpsv is called) 
                   = 0 otherwise
*/
#if ((SCOTCH_VERSION == 6) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 7)
/* esmumpsv prototype with 64-bit integers weights of nodes in the graph are used on entry (nvtab) */
     if ( *weightrequested == 1 )
     {
       *ncmpa = esmumpsv( *n, *iwlen, petab, *pfree,
                          lentab, iwtab, nvtab, elentab, lasttab );
       *weightused=1;
     }
     else
     {
       /* esmumps prototype with standard integers (weights of nodes not used on entry) */
       *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                         lentab, iwtab, nvtab, elentab, lasttab );
       *weightused=0;
     }
#else
     /* esmumps prototype with standard integers (weights of nodes not used on entry) */
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
     *weightused=0;
#endif
}
#endif /* scotch */
#if defined(ptscotch)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr)
{
  MPI_Comm  int_comm;
  int_comm = MPI_Comm_f2c(*comm);
  *ierr = SCOTCH_dgraphInit(graphptr, int_comm);
  return;
}
#endif
