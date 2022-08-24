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
#include <errno.h>
#include "mumps_common.h"
#include "mumps_c_types.h"
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
# include <pthread.h>
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
extern pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Exported functions */
#define MUMPS_LOW_LEVEL_INIT_ERR_STR \
    F_SYMBOL(low_level_init_err_str,LOW_LEVEL_INIT_ERR_STR)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR( MUMPS_INT *dim, char *err_str, mumps_ftnlen l1 );
/* Export an error to the Fortran layer
   Returns mumps_errno for convenience */
MUMPS_INT mumps_io_error(MUMPS_INT mumps_errno, const char* desc);
/* Export a system error to the Fortran layer (errno must be set)
   Returns mumps_errno for convenience */
MUMPS_INT mumps_io_sys_error(MUMPS_INT mumps_errno, const char* desc);
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
MUMPS_INT mumps_io_init_err_lock();
MUMPS_INT mumps_io_destroy_err_lock();
MUMPS_INT mumps_check_error_th();
MUMPS_INLINE MUMPS_INT mumps_io_protect_err();
MUMPS_INLINE MUMPS_INT mumps_io_unprotect_err();
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
