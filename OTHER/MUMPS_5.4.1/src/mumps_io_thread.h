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
#ifndef MUMPS_IO_THREAD_H
#define MUMPS_IO_THREAD_H
#include "mumps_compat.h"
#include "mumps_c_types.h"
#if ! defined (MUMPS_WIN32) && ! defined (WITHOUT_PTHREAD)
# include <unistd.h>
# include <pthread.h>
# include <sys/types.h>
# include <sys/time.h>
# include <time.h>
# define MAX_IO 20
# define MAX_FINISH_REQ 40
# define IO_FLAG_STOP 1
# define IO_FLAG_RUN 0
# define IO_READ 1
# define IO_WRITE 0
struct request_io{
  MUMPS_INT inode;
  MUMPS_INT req_num; /*request number*/
  void* addr;  /*memory address (either source or dest)*/
  long long size;    /* size of the requested io (unit=size of elementary mumps data)*/
  long long vaddr; /* virtual address for file management */
  MUMPS_INT io_type; /*read or write*/
  MUMPS_INT file_type; /* cb or lu or ... */
  pthread_cond_t local_cond;
  MUMPS_INT int_local_cond;
};
/* Exported global variables */
extern MUMPS_INT io_flag_stop,current_req_num;
extern pthread_t io_thread,main_thread;
extern pthread_mutex_t io_mutex;
extern pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
extern pthread_mutex_t io_mutex_cond;
extern MUMPS_INT int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;
extern MUMPS_INT with_sem;
extern struct request_io *io_queue;
extern MUMPS_INT first_active,last_active,nb_active;
extern MUMPS_INT *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
extern MUMPS_INT mumps_owns_mutex;
extern MUMPS_INT test_request_called_from_mumps;
/* Exported functions */
void* mumps_async_thread_function_with_sem (void* arg);
MUMPS_INT   mumps_is_there_finished_request_th(MUMPS_INT* flag);
MUMPS_INT   mumps_clean_request_th(MUMPS_INT* request_id);
MUMPS_INT   mumps_wait_req_sem_th(MUMPS_INT *request_id);
MUMPS_INT   mumps_test_request_th(MUMPS_INT* request_id,MUMPS_INT *flag);
MUMPS_INT   mumps_wait_request_th(MUMPS_INT *request_id);
MUMPS_INT   mumps_low_level_init_ooc_c_th(MUMPS_INT* async, MUMPS_INT* ierr);
MUMPS_INT   mumps_async_write_th(const MUMPS_INT * strat_IO,void * address_block,long long block_size,
                           MUMPS_INT * inode,MUMPS_INT * request_arg,MUMPS_INT * type,long long vaddr,MUMPS_INT * ierr);
MUMPS_INT   mumps_async_read_th(const MUMPS_INT * strat_IO,void * address_block,long long block_size,MUMPS_INT * inode,MUMPS_INT * request_arg,
                           MUMPS_INT * type,long long vaddr,MUMPS_INT * ierr);
MUMPS_INT mumps_clean_io_data_c_th(MUMPS_INT *myid);
MUMPS_INT mumps_get_sem(void *arg,MUMPS_INT *value);
MUMPS_INT mumps_wait_sem(void *arg,pthread_cond_t *cond);
MUMPS_INT mumps_post_sem(void *arg,pthread_cond_t *cond);
MUMPS_INT mumps_clean_finished_queue_th();
#endif /*_WIN32 && WITHOUT_PTHREAD*/
#endif /* MUMPS_IO_THREAD_H */
