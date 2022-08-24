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
#ifndef MUMPS_PORD_H
#define MUMPS_PORD_H
#include "mumps_common.h"
#define MUMPS_PORD_INTSIZE \
  F_SYMBOL(pord_intsize,PORD_INTSIZE)
void MUMPS_CALL MUMPS_PORD_INTSIZE(MUMPS_INT *pord_intsize);
#if defined(pord)
#include <space.h>
MUMPS_INT mumps_pord( PORD_INT, PORD_INT, PORD_INT *, PORD_INT *, PORD_INT * );
#define MUMPS_PORDF \
    F_SYMBOL(pordf,PORDF)
#if defined(INTSIZE64) || defined(PORD_INTSIZE64)
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT8 *nvtx, MUMPS_INT8 *nedges,
             MUMPS_INT8 *xadj, MUMPS_INT8 *adjncy,
             MUMPS_INT8 *nv, MUMPS_INT *ncmpa );
#else
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT *nvtx, MUMPS_INT *nedges,
             MUMPS_INT *xadj, MUMPS_INT *adjncy,
             MUMPS_INT *nv, MUMPS_INT *ncmpa );
#endif
MUMPS_INT mumps_pord_wnd( PORD_INT, PORD_INT, PORD_INT *, PORD_INT *, PORD_INT *, PORD_INT * );
#define MUMPS_PORDF_WND          \
    F_SYMBOL(pordf_wnd,PORDF_WND)
#if defined(INTSIZE64) || defined(PORD_INTSIZE64)
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT8 *nvtx, MUMPS_INT8 *nedges,
                 MUMPS_INT8 *xadj, MUMPS_INT8 *adjncy,
                 MUMPS_INT8 *nv, MUMPS_INT *ncmpa, MUMPS_INT8 *totw );
#else
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT *nvtx, MUMPS_INT *nedges,
                 MUMPS_INT *xadj, MUMPS_INT *adjncy,
                 MUMPS_INT *nv, MUMPS_INT *ncmpa, MUMPS_INT *totw );
#endif
#endif /*PORD*/
#endif /* MUMPS_PORD_H */
