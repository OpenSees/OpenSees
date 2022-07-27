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
#ifndef MUMPS_METIS_H
#define MUMPS_METIS_H
/* Interfacing with 32-bit (par)metis, for METIS 4 or METIS 5 */
#include "mumps_common.h" /* includes mumps_compat.h and mumps_c_types.h */
#if defined(parmetis) || defined(parmetis3)
#include "mpi.h"
#define MUMPS_PARMETIS \
  F_SYMBOL(parmetis,PARMETIS)
void MUMPS_CALL
MUMPS_PARMETIS(MUMPS_INT *first,      MUMPS_INT *vertloctab,
               MUMPS_INT *edgeloctab, MUMPS_INT *numflag,
               MUMPS_INT *options,    MUMPS_INT *order,
               MUMPS_INT *sizes,      MUMPS_INT *comm,
               MUMPS_INT *ierr);
#endif
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
#define MUMPS_METIS_KWAY \
  F_SYMBOL(metis_kway,METIS_KWAY)
void MUMPS_CALL
MUMPS_METIS_KWAY(MUMPS_INT *n,     MUMPS_INT *iptr,
                 MUMPS_INT *jcn,   MUMPS_INT *k,
                 MUMPS_INT *part);
#define MUMPS_METIS_KWAY_AB \
  F_SYMBOL(metis_kway_ab,METIS_KWAY_AB)
void MUMPS_CALL
MUMPS_METIS_KWAY_AB(MUMPS_INT *n,     MUMPS_INT *iptr,
                 MUMPS_INT *jcn,   MUMPS_INT *k,
                 MUMPS_INT *part,  MUMPS_INT *vwgt);
#endif
#endif
