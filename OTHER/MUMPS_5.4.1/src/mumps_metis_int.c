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
#include "mumps_metis_int.h"
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
#  include "metis.h"
#  if defined(parmetis3) || defined(metis4)
     /* IDXTYPEWIDTH not available, use sizeof(idxtype) */
     /* We use metis.h and assume that parmetis datatypes
        will be identical to those of metis.h since it
        does not make senss to combine metis and parmetis
        with different int sizes */
     void MUMPS_CALL
     MUMPS_METIS_IDXSIZE(MUMPS_INT *metis_idx_size)
     {         
       *metis_idx_size=8*sizeof(idxtype);
     }
#  else
     /* Rely on IDXTYPEWIDTH */
     void MUMPS_CALL
     MUMPS_METIS_IDXSIZE(MUMPS_INT *metis_idx_size)
     {
       /* *metis_idx_size=sizeof(idx_t); */
       *metis_idx_size=IDXTYPEWIDTH;
     }
#  endif
#else
  void MUMPS_CALL
  MUMPS_METIS_IDXSIZE(MUMPS_INT *metis_int_size)
  {
    *metis_int_size=-99999;
  }
#endif
