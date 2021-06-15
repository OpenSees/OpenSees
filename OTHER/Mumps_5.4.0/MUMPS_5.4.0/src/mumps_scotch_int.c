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
#include "mumps_scotch_int.h"
#if defined(scotch) || defined(ptscotch)
#  include <stdio.h>
#  include "scotch.h"
#endif
void MUMPS_CALL
MUMPS_SCOTCH_INTSIZE(MUMPS_INT *scotch_intsize)
  {
#    if defined(scotch) || defined(ptscotch)
         *scotch_intsize=8*sizeof(SCOTCH_Num);
#    else
       *scotch_intsize=-99999;
#    endif
}
