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
#include "mumps_io_err.h"
#include "mumps_io_basic.h"
#include "mumps_c_types.h"
#if defined( MUMPS_WIN32 )
# include <string.h>
#endif
/* Exported global variables */
char* mumps_err;
MUMPS_INT* dim_mumps_err;
MUMPS_INT mumps_err_max_len;
MUMPS_INT err_flag;
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Functions */
/* Keeps a C pointer to store error description string that will be
   displayed by the Fortran layers.
   * dim contains the size of the Fortran character array to store the
   description.
*/
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR(MUMPS_INT *dim, char* err_str, mumps_ftnlen l1){
  mumps_err = err_str;
  dim_mumps_err = (MUMPS_INT *) dim;
  mumps_err_max_len = (MUMPS_INT) *dim;
  err_flag = 0;
  return;
}
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
MUMPS_INLINE MUMPS_INT
mumps_io_protect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_lock(&err_mutex);
  }
  return 0;
}
MUMPS_INLINE MUMPS_INT
mumps_io_unprotect_err()
{
  if(mumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_unlock(&err_mutex);
  }
  return 0;
}
MUMPS_INT
mumps_io_init_err_lock()
{
  pthread_mutex_init(&err_mutex,NULL);
  return 0;
}
MUMPS_INT
mumps_io_destroy_err_lock()
{
  pthread_mutex_destroy(&err_mutex);
  return 0;
}
MUMPS_INT
mumps_check_error_th()
{
  /* If err_flag != 0, then error_str is set */
  return err_flag;
}
#endif /* MUMPS_WIN32 && WITHOUT_PTHREAD */
MUMPS_INT
mumps_io_error(MUMPS_INT mumps_errno, const char* desc)
{
    MUMPS_INT len;
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag == 0){
    strncpy(mumps_err, desc, mumps_err_max_len);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
    len = (MUMPS_INT) strlen(desc);
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
MUMPS_INT
mumps_io_sys_error(MUMPS_INT mumps_errno, const char* desc)
{
  MUMPS_INT len = 2; /* length of ": " */
  const char* _desc;
  char* _err;
#if defined( MUMPS_WIN32 )
  MUMPS_INT _err_len;
#endif
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_protect_err();
#endif
  if(err_flag==0){
    if(desc == NULL) {
      _desc = "";
    } else {
        len += (MUMPS_INT) strlen(desc);
      _desc = desc;
    }
#if ! defined( MUMPS_WIN32 )
    _err = strerror(errno);
    len += (MUMPS_INT) strlen(_err);
    snprintf(mumps_err, mumps_err_max_len, "%s: %s", _desc, _err);
    /* mumps_err is a FORTRAN string, we do not care about adding a final 0 */
#else
    /* This a VERY UGLY workaround for snprintf: this function has been
     * integrated quite lately into the ANSI stdio: some windows compilers are
     * not up-to-date yet. */
    if( len >= mumps_err_max_len - 1 ) { /* then do not print sys error msg at all */
      len -= 2;
      len = (len >= mumps_err_max_len ) ? mumps_err_max_len - 1 : len;
      _err = strdup( _desc );
      _err[len] = '\0';
      sprintf(mumps_err, "%s", _err);
    } else {
      _err = strdup(strerror(errno));
      _err_len = (MUMPS_INT) strlen(_err);
      /* We will use sprintf, so make space for the final '\0' ! */
      if((len + _err_len) >= mumps_err_max_len) {
        /* truncate _err, not to overtake mumps_err_max_len at the end. */
        _err[mumps_err_max_len - len - 1] = '\0';
        len = mumps_err_max_len - 1;
      } else {
        len += _err_len;
      }
      sprintf(mumps_err, "%s: %s", _desc, _err);
    }
    free(_err);
#endif
    *dim_mumps_err = (len <= mumps_err_max_len ) ? len : mumps_err_max_len;
    err_flag = mumps_errno;
  }
#if ! defined( MUMPS_WIN32 ) && ! defined( WITHOUT_PTHREAD )
  mumps_io_unprotect_err();
#endif
  return mumps_errno;
}
