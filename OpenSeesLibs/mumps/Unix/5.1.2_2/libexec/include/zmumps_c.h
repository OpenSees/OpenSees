/*
 *
 *  This file is part of MUMPS 5.1.2, released
 *  on Mon Oct  2 07:37:01 UTC 2017
 *
 *
 *  Copyright 1991-2017 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */

/* Mostly written in march 2002 (JYL) */

#ifndef ZMUMPS_C_H
#define ZMUMPS_C_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mumps_compat.h"
/* Next line defines MUMPS_INT, ZMUMPS_COMPLEX and ZMUMPS_REAL */
#include "mumps_c_types.h"

#ifndef MUMPS_VERSION
/* Protected in case headers of other arithmetics are included */
#define MUMPS_VERSION "5.1.2"
#endif
#ifndef MUMPS_VERSION_MAX_LEN
#define MUMPS_VERSION_MAX_LEN 30
#endif

/*
 * Definition of the (simplified) MUMPS C structure.
 * NB: ZMUMPS_COMPLEX are REAL types in s and d arithmetics.
 */
typedef struct {

    MUMPS_INT      sym, par, job;
    MUMPS_INT      comm_fortran;    /* Fortran communicator */
    MUMPS_INT      icntl[40];
    MUMPS_INT      keep[500];
    ZMUMPS_REAL    cntl[15];
    ZMUMPS_REAL    dkeep[230];
    MUMPS_INT8     keep8[150];
    MUMPS_INT      n;

    MUMPS_INT      nz_alloc; /* used in matlab interface to decide if we
                                free + malloc when we have large variation */

    /* Assembled entry */
    MUMPS_INT      nz;
    MUMPS_INT8     nnz;
    MUMPS_INT      *irn;
    MUMPS_INT      *jcn;
    ZMUMPS_COMPLEX *a;

    /* Distributed entry */
    MUMPS_INT      nz_loc;
    MUMPS_INT8     nnz_loc;
    MUMPS_INT      *irn_loc;
    MUMPS_INT      *jcn_loc;
    ZMUMPS_COMPLEX *a_loc;

    /* Element entry */
    MUMPS_INT      nelt;
    MUMPS_INT      *eltptr;
    MUMPS_INT      *eltvar;
    ZMUMPS_COMPLEX *a_elt;

    /* Ordering, if given by user */
    MUMPS_INT      *perm_in;

    /* Orderings returned to user */
    MUMPS_INT      *sym_perm;    /* symmetric permutation */
    MUMPS_INT      *uns_perm;    /* column permutation */

    /* Scaling (inout but complicated) */
    ZMUMPS_REAL    *colsca;
    ZMUMPS_REAL    *rowsca;
    MUMPS_INT colsca_from_mumps;
    MUMPS_INT rowsca_from_mumps;

    /* RHS, solution, ouptput data and statistics */
    ZMUMPS_COMPLEX *rhs, *redrhs, *rhs_sparse, *sol_loc;
    MUMPS_INT      *irhs_sparse, *irhs_ptr, *isol_loc;
    MUMPS_INT      nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
    MUMPS_INT      schur_mloc, schur_nloc, schur_lld;
    MUMPS_INT      mblock, nblock, nprow, npcol;
    MUMPS_INT      info[40],infog[40];
    ZMUMPS_REAL    rinfo[40], rinfog[40];

    /* Null space */
    MUMPS_INT      deficiency;
    MUMPS_INT      *pivnul_list;
    MUMPS_INT      *mapping;

    /* Schur */
    MUMPS_INT      size_schur;
    MUMPS_INT      *listvar_schur;
    ZMUMPS_COMPLEX *schur;

    /* Internal parameters */
    MUMPS_INT      instance_number;
    ZMUMPS_COMPLEX *wk_user;

    /* Version number: length=14 in FORTRAN + 1 for final \0 + 1 for alignment */
    char version_number[MUMPS_VERSION_MAX_LEN + 1 + 1];
    /* For out-of-core */
    char ooc_tmpdir[256];
    char ooc_prefix[64];
    /* To save the matrix in matrix market format */
    char write_problem[256];
    MUMPS_INT      lwk_user;
    /* For save/restore feature */
    char save_dir[256];
    char save_prefix[256];
} ZMUMPS_STRUC_C;


void MUMPS_CALL
zmumps_c( ZMUMPS_STRUC_C * zmumps_par );

#ifdef __cplusplus
}
#endif

#endif /* ZMUMPS_C_H */

