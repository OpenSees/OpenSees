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
#ifndef MUMPS_METIS64_H
#define MUMPS_METIS64_H
/* Interfacing with 64-bit (par)metis, for METIS 4 or METIS 5 */
#include "mumps_common.h" /* includes mumps_compat.h and mumps_c_types.h */
#if defined(parmetis) || defined(parmetis3)
#include "mpi.h"
#define MUMPS_PARMETIS_64 \
  F_SYMBOL(parmetis_64,PARMETIS_64)
void MUMPS_CALL
MUMPS_PARMETIS_64(MUMPS_INT8 *first,      MUMPS_INT8 *vertloctab,
                  MUMPS_INT8 *edgeloctab,
#if defined(parmetis3)
                  MUMPS_INT  *numflag, MUMPS_INT  *options,
#else
                  MUMPS_INT8 *numflag, MUMPS_INT8 *options,
#endif
                  MUMPS_INT8 *order,
                  MUMPS_INT8 *sizes,         MUMPS_INT *comm,
                  MUMPS_INT  *ierr);
#endif
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
#define MUMPS_METIS_KWAY_64 \
  F_SYMBOL(metis_kway_64,METIS_KWAY_64)
void MUMPS_CALL
MUMPS_METIS_KWAY_64(MUMPS_INT8 *n,     MUMPS_INT8 *iptr,
                 MUMPS_INT8 *jcn,   MUMPS_INT8 *k,
                 MUMPS_INT8 *part);
#define MUMPS_METIS_KWAY_AB_64 \
  F_SYMBOL(metis_kway_ab_64,METIS_KWAY_AB_64)
void MUMPS_CALL
MUMPS_METIS_KWAY_AB_64(MUMPS_INT8 *n,     MUMPS_INT8 *iptr,
                 MUMPS_INT8 *jcn,   MUMPS_INT8 *k,
                 MUMPS_INT8 *part,  MUMPS_INT8 *vwgt);
#endif
#endif
