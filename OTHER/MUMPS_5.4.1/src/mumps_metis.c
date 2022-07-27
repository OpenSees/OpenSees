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
#include <stdio.h>  /* For NULL constant (stddef.h) and debug printings */
#include "mumps_metis.h"
#if defined(parmetis) || defined(parmetis3)
/*PARMETIS*/
#if defined(parmetis3)
/* Provide prototype by hand. This is because we are not sure
 * at compilation/preprocessing time whether we use a 32-bit
 * or a 64-bit metis */
  void ParMETIS_V3_NodeND(MUMPS_INT *first, MUMPS_INT *vertloctab, MUMPS_INT *edgeloctab, MUMPS_INT *numflag, MUMPS_INT *options, MUMPS_INT *order, MUMPS_INT *sizes, MPI_Comm *Ccomm);
#else
#include "metis.h"
#include "parmetis.h" /* Prototypes from parmetis.h will be used */
#endif
void MUMPS_CALL
MUMPS_PARMETIS(MUMPS_INT *first,      MUMPS_INT *vertloctab,
               MUMPS_INT *edgeloctab, MUMPS_INT *numflag,
               MUMPS_INT *options,    MUMPS_INT *order,
               MUMPS_INT *sizes,      MUMPS_INT *comm,
               MUMPS_INT *ierr)
{
  MPI_Comm  int_comm;
  int iierr;
  int_comm = MPI_Comm_f2c(*comm);
#if defined(parmetis3)
  ParMETIS_V3_NodeND(first, vertloctab, edgeloctab, numflag, options, order, sizes, &int_comm);
#elif defined(parmetis)
#  if (IDXTYPEWIDTH == 32)
      *ierr=0;
      iierr=ParMETIS_V3_NodeND(first, vertloctab, edgeloctab, numflag, options, order, sizes, &int_comm);
      if(iierr != METIS_OK)
        *ierr=1;
#  else
      /* SHOULD NEVER BE CALLED */
      printf("** Error: ParMETIS version >= 4, IDXTYPE WIDTH !=64, but MUMPS_PARMETIS_64 was called\n");
      *ierr=1;
#  endif
#endif
  return;
}
#endif
#if defined(parmetis) || defined(metis) || defined(parmetis3) || defined(metis4)
#if defined(metis4) || defined(parmetis3) /* parmetis3 comes with metis4 */
/* Provide prototype by hand. This is because we are not sure
 * at compilation/preprocessing time whether we use a 32-bit
 * or a 64-bit metis */
void METIS_PartGraphKway(int *, MUMPS_INT *, MUMPS_INT *, MUMPS_INT *, MUMPS_INT *, int *, int *, int *, int *, int *, MUMPS_INT *);
#else
/* Prototype properly defined in metis.h
 * One can rely on IDXTYPEWIDTH to know at compilation/preprocessing
 * time whether we use a 32-bit or a 64-bit metis */
#include "metis.h"
#endif
/* Interface for metis k-way partitioning with standard ints */
void MUMPS_CALL
MUMPS_METIS_KWAY(MUMPS_INT *n,     MUMPS_INT *iptr,
                 MUMPS_INT *jcn,   MUMPS_INT *k,
                 MUMPS_INT *part)
/* n     -- the size of the graph to be partitioned
   iptr  -- pointer to the beginning of each node's adjacency list
   jcn   -- jcn[iptr[i]:iptr[i+1]-1] contains the list of neighbors of node i
   k     -- the number of parts
   part  -- part[i] is the part node i belongs to */
 {
#if defined(metis4) || defined(parmetis3)
  MUMPS_INT numflag, edgecut, wgtflag, options[8];
  options[0] = 0;
  /* unweighted partitioning */
  wgtflag    = 0;
  /* Use 1-based fortran numbering */
  numflag    = 1;
/* void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); */
  METIS_PartGraphKway(n, iptr, jcn,
                      NULL, NULL, &wgtflag,
                      &numflag, k,
                      options, &edgecut,
                      part);
#else /* METIS >= 5 */
  int ierr;
#  if (IDXTYPEWIDTH == 32)
  MUMPS_INT ncon, edgecut, options[40];
  ierr=METIS_SetDefaultOptions(options);
  options[0]  = 0;
  /* Use 1-based fortran numbering */
  options[17] = 1;
  ncon        = 1;
  ierr = METIS_PartGraphKway(n, &ncon, iptr, jcn,
                             NULL, NULL, NULL,
                             k, NULL, NULL, options,
                             &edgecut, part);
#  else
  /* SHOULD NEVER BE CALLED */
  printf("** Error: METIS version >= 4, IDXTYPE WIDTH !=32, but MUMPS_METIS_KWAY was called\n");
  ierr=1;
#  endif
#endif
  return;
 }
/* Interface for metis k-way partitioning with standard ints and weights on vertices*/
void MUMPS_CALL
MUMPS_METIS_KWAY_AB(MUMPS_INT *n,     MUMPS_INT *iptr,
                 MUMPS_INT *jcn,   MUMPS_INT *k,
                 MUMPS_INT *part, MUMPS_INT *vwgt)
/* n     -- the size of the graph to be partitioned
   iptr  -- pointer to the beginning of each node's adjacency list
   jcn   -- jcn[iptr[i]:iptr[i+1]-1] contains the list of neighbors of node i
   k     -- the number of parts
   part  -- part[i] is the part node i belongs to 
   vwgt  -- weights of the vertices*/
 {
#if defined(metis4) || defined(parmetis3)
  MUMPS_INT numflag, edgecut, wgtflag, options[8];
  options[0] = 0;
  /* unweighted partitioning */
  wgtflag    = 0;
  /* Use 1-based fortran numbering */
  numflag    = 1;
/* void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); */
  METIS_PartGraphKway(n, iptr, jcn,
                      vwgt, NULL, &wgtflag,
                      &numflag, k,
                      options, &edgecut,
                      part);
#else /* METIS >= 5 */
  int ierr;
#  if (IDXTYPEWIDTH == 32)
  MUMPS_INT ncon, edgecut, options[40];
  ierr=METIS_SetDefaultOptions(options);
  options[0]  = 0;
  /* Use 1-based fortran numbering */
  options[17] = 1;
  ncon        = 1;
  ierr = METIS_PartGraphKway(n, &ncon, iptr, jcn,
                             vwgt, NULL, NULL,
                             k, NULL, NULL, options,
                             &edgecut, part);
#  else
  /* SHOULD NEVER BE CALLED */
  printf("** Error: METIS version >= 4, IDXTYPE WIDTH !=32, but MUMPS_METIS_KWAY_AB was called\n");
  ierr=1;
#  endif
#endif
  return;
 }
#endif
