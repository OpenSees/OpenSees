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
#ifndef MUMPS_IO_H
#define MUMPS_IO_H
#include "mumps_common.h"
#include "mumps_c_types.h"
/*
 *  Two character arrays that are used by low_level_init_prefix
 *  and low_level_init_tmpdir to store intermediate file names.
 *  They are passed to mumps_io_basic.c inside the routine
 *  mumps_low_level_init_ooc_c.
 *  Note that both low_level_init_prefix and low_level_init_tmpdir
 *  MUST be called before low_level_init_ooc_c.
 * 
 */
#define MUMPS_OOC_PREFIX_MAX_LENGTH 63
#define MUMPS_OOC_TMPDIR_MAX_LENGTH 255
static char MUMPS_OOC_STORE_PREFIX[MUMPS_OOC_PREFIX_MAX_LENGTH];
static MUMPS_INT  MUMPS_OOC_STORE_PREFIXLEN=-1;
static char MUMPS_OOC_STORE_TMPDIR[MUMPS_OOC_TMPDIR_MAX_LENGTH];
static MUMPS_INT  MUMPS_OOC_STORE_TMPDIRLEN=-1;
#define MUMPS_DUMPRHSBINARY_C \
    F_SYMBOL(dumprhsbinary_c,DUMP_RHSBINARY_C)
void MUMPS_CALL MUMPS_DUMPRHSBINARY_C ( MUMPS_INT *N, MUMPS_INT *NRHS,
     MUMPS_INT *LRHS, float *RHS, MUMPS_INT *K35, 
     char *filename, mumps_ftnlen l1 );
#define MUMPS_DUMPMATBINARY_C \
    F_SYMBOL(dumpmatbinary_c,DUMP_MATBINARY_C)
void MUMPS_CALL MUMPS_DUMPMATBINARY_C ( MUMPS_INT* N, MUMPS_INT8 *NNZ,
     MUMPS_INT* K35, MUMPS_INT *irn, MUMPS_INT *jcn,
     void *A, MUMPS_INT *is_A_provided,
     char *file_name, mumps_ftnlen l1 );
#define MUMPS_LOW_LEVEL_INIT_PREFIX \
    F_SYMBOL(low_level_init_prefix,LOW_LEVEL_INIT_PREFIX)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_PREFIX(MUMPS_INT * dim, char * str, mumps_ftnlen l1);
#define MUMPS_LOW_LEVEL_INIT_TMPDIR \
    F_SYMBOL(low_level_init_tmpdir,LOW_LEVEL_INIT_TMPDIR)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_TMPDIR(MUMPS_INT * dim, char * str, mumps_ftnlen l1);
MUMPS_INLINE MUMPS_INT
mumps_convert_2fint_to_longlong( MUMPS_INT *, MUMPS_INT *, long long *);
#define MUMPS_LOW_LEVEL_INIT_OOC_C \
    F_SYMBOL(low_level_init_ooc_c,LOW_LEVEL_INIT_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_OOC_C(MUMPS_INT *_myid, MUMPS_INT *total_size_io,MUMPS_INT *size_element,
                           MUMPS_INT *async, MUMPS_INT *k211, MUMPS_INT *nb_file_type,
                           MUMPS_INT *flag_tab , MUMPS_INT* ierr);
#define MUMPS_TEST_REQUEST_C \
    F_SYMBOL(test_request_c,TEST_REQUEST_C)
void MUMPS_CALL
MUMPS_TEST_REQUEST_C(MUMPS_INT *request_id,MUMPS_INT *flag,MUMPS_INT *ierr);
#define MUMPS_WAIT_REQUEST \
    F_SYMBOL(wait_request,WAIT_REQUEST)
void MUMPS_CALL
MUMPS_WAIT_REQUEST(MUMPS_INT *request_id,MUMPS_INT *ierr);
#define MUMPS_LOW_LEVEL_WRITE_OOC_C \
    F_SYMBOL(low_level_write_ooc_c,LOW_LEVEL_WRITE_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_WRITE_OOC_C( const MUMPS_INT * strat_IO, 
                             void * address_block,
                             MUMPS_INT * block_size_int1,
                             MUMPS_INT * block_size_int2,
                             MUMPS_INT * inode,
                             MUMPS_INT * request_arg,
                             MUMPS_INT * type,
                             MUMPS_INT * vaddr_int1,
                             MUMPS_INT * vaddr_int2,
                             MUMPS_INT * ierr);
#define MUMPS_LOW_LEVEL_READ_OOC_C \
    F_SYMBOL(low_level_read_ooc_c,LOW_LEVEL_READ_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_READ_OOC_C( const MUMPS_INT * strat_IO, 
                            void * address_block,
                            MUMPS_INT * block_size_int1,
                            MUMPS_INT * block_size_int2,
                            MUMPS_INT * inode,
                            MUMPS_INT * request_arg,
                            MUMPS_INT * type,
                            MUMPS_INT * vaddr_int1,
                            MUMPS_INT * vaddr_int2,
                            MUMPS_INT * ierr);
#define MUMPS_LOW_LEVEL_DIRECT_READ \
    F_SYMBOL(low_level_direct_read,LOW_LEVEL_DIRECT_READ)
void MUMPS_CALL
MUMPS_LOW_LEVEL_DIRECT_READ(void * address_block,
                            MUMPS_INT * block_size_int1,
                            MUMPS_INT * block_size_int2,
                            MUMPS_INT * type,
                            MUMPS_INT * vaddr_int1,
                            MUMPS_INT * vaddr_int2,
                            MUMPS_INT * ierr);
#define MUMPS_CLEAN_IO_DATA_C \
    F_SYMBOL(clean_io_data_c,CLEAN_IO_DATA_C)
void MUMPS_CALL
MUMPS_CLEAN_IO_DATA_C(MUMPS_INT *myid,MUMPS_INT *step,MUMPS_INT *ierr);
#define MUMPS_GET_MAX_NB_REQ_C \
    F_SYMBOL(get_max_nb_req_c,GET_MAX_NB_REQ_C)
void MUMPS_CALL
MUMPS_GET_MAX_NB_REQ_C(MUMPS_INT *max,MUMPS_INT *ierr);
#define MUMPS_GET_MAX_FILE_SIZE_C \
    F_SYMBOL(get_max_file_size_c,GET_MAX_FILE_SIZE_C)
void MUMPS_CALL
MUMPS_GET_MAX_FILE_SIZE_C(double * max_ooc_file_size);
#define MUMPS_OOC_GET_NB_FILES_C \
    F_SYMBOL(ooc_get_nb_files_c,OOC_GET_NB_FILES_C)
void MUMPS_CALL
MUMPS_OOC_GET_NB_FILES_C(const MUMPS_INT *type, MUMPS_INT *nb_files);
#define MUMPS_OOC_GET_FILE_NAME_C \
    F_SYMBOL(ooc_get_file_name_c,OOC_GET_FILE_NAME_C)
void MUMPS_CALL
MUMPS_OOC_GET_FILE_NAME_C(MUMPS_INT *type, MUMPS_INT *indice, MUMPS_INT *length,
                          char* name, mumps_ftnlen l1);
#define MUMPS_OOC_SET_FILE_NAME_C \
    F_SYMBOL(ooc_set_file_name_c,OOC_SET_FILE_NAME_C)
void MUMPS_CALL
MUMPS_OOC_SET_FILE_NAME_C(MUMPS_INT *type, MUMPS_INT *indice, MUMPS_INT *length, MUMPS_INT *ierr,
                          char* name, mumps_ftnlen l1);
#define MUMPS_OOC_ALLOC_POINTERS_C \
    F_SYMBOL(ooc_alloc_pointers_c,OOC_ALLOC_POINTERS_C)
void MUMPS_CALL
MUMPS_OOC_ALLOC_POINTERS_C(MUMPS_INT *nb_file_type, MUMPS_INT *dim, MUMPS_INT *ierr);
#define MUMPS_OOC_INIT_VARS_C \
    F_SYMBOL(ooc_init_vars_c,OOC_INIT_VARS_C)
void MUMPS_CALL
MUMPS_OOC_INIT_VARS_C(MUMPS_INT *myid_arg, MUMPS_INT *size_element, MUMPS_INT *async,
                      MUMPS_INT *k211, MUMPS_INT *ierr);
#define MUMPS_OOC_START_LOW_LEVEL \
    F_SYMBOL(ooc_start_low_level,OOC_START_LOW_LEVEL)
void MUMPS_CALL
MUMPS_OOC_START_LOW_LEVEL(MUMPS_INT *ierr);
#define MUMPS_OOC_PRINT_STATS \
    F_SYMBOL(ooc_print_stats,OOC_PRINT_STATS)
void MUMPS_CALL
MUMPS_OOC_PRINT_STATS();
#define MUMPS_OOC_REMOVE_FILE_C \
    F_SYMBOL(ooc_remove_file_c,OOC_REMOVE_FILE_C)
void MUMPS_CALL
MUMPS_OOC_REMOVE_FILE_C(MUMPS_INT *ierr, char *name, mumps_ftnlen l1);
#define MUMPS_OOC_END_WRITE_C \
    F_SYMBOL(ooc_end_write_c,OOC_END_WRITE_C)
void MUMPS_CALL
MUMPS_OOC_END_WRITE_C(MUMPS_INT *ierr);
#define MUMPS_OOC_IS_ASYNC_AVAIL \
    F_SYMBOL(ooc_is_async_avail,OOC_IS_ASYNC_AVAIL)
void MUMPS_CALL
MUMPS_OOC_IS_ASYNC_AVAIL(MUMPS_INT *flag);
#endif /* MUMPS_IO_H */
