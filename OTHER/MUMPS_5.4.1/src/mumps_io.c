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
#include "mumps_io.h"
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
#include "mumps_c_types.h"
#if ! defined (MUMPS_WIN32) && ! defined (WITHOUT_PTHREAD)
# include "mumps_io_thread.h"
#endif
#if ! defined(MUMPS_WIN32)
double mumps_time_spent_in_sync;
#endif
double read_op_vol,write_op_vol,total_vol;
void MUMPS_CALL MUMPS_DUMPRHSBINARY_C ( MUMPS_INT *N, MUMPS_INT *NRHS,
     MUMPS_INT *LRHS, float *RHS, MUMPS_INT *K35,
     char *filename, mumps_ftnlen l1 )
{
   float *RHSshift; /* float: arbitrary, we use binary content */
   FILE *fd;
   int icol;
   fd=fopen(filename, "w");
   RHSshift=RHS;
   for(icol=0;icol<*NRHS;icol++)
   {
     fwrite(RHSshift, (size_t)(*K35), (size_t)(*N), fd);
     RHSshift=RHSshift+(size_t)(*LRHS)*(size_t)(*K35/sizeof(float));
   }
   fclose(fd);
}
void MUMPS_CALL MUMPS_DUMPMATBINARY_C ( MUMPS_INT *N, MUMPS_INT8 *NNZ,
     MUMPS_INT* K35, MUMPS_INT *irn, MUMPS_INT *jcn,
     void *A, MUMPS_INT *is_A_provided,
     char *filename, mumps_ftnlen l1 )
{
   int64_t i8;
   int32_t myN, tmpi;
   FILE *fd;
   fd=fopen(filename, "w");
   /* cast to int32_t in case MUMPS_INT is 64-bits */
   myN=(int32_t)(*N);
   fwrite( &myN, sizeof(int32_t), 1, fd);
   fwrite( NNZ, sizeof(int64_t), 1, fd);
   if (*NNZ > 0)
   {
     if ( sizeof(MUMPS_INT) == 4 )
     {
       /* write irn and jcn contents directly */
       fwrite( irn, sizeof(int32_t), (size_t)(*NNZ), fd);
       fwrite( jcn, sizeof(int32_t), (size_t)(*NNZ), fd);
     }
     else
     {
       for(i8=0;i8 < *NNZ;i8++)
       {
          tmpi=irn[i8];
          fwrite(&tmpi, sizeof(int32_t), 1, fd);
       }
       for(i8=0;i8 < *NNZ;i8++)
       {
          tmpi=jcn[i8];
          fwrite(&tmpi, sizeof(int32_t), 1, fd);
       }
     }
     if (*is_A_provided)
     {
       fwrite(A, (size_t)(*K35), (size_t)(*NNZ), fd);
     }
   }
   fclose(fd);
}
/* Tests if the request "request_id" has finished. It sets the flag  */
/* argument to 1 if the request has finished (0 otherwise)           */
void MUMPS_CALL
MUMPS_TEST_REQUEST_C(MUMPS_INT *request_id,MUMPS_INT *flag,MUMPS_INT *ierr)
{
  char buf[64]; /* for error message */
  MUMPS_INT request_id_loc;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  MUMPS_INT flag_loc;
#endif
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  request_id_loc=(MUMPS_INT)*request_id;
  switch(mumps_io_flag_async){
  case IO_SYNC:
    /* printf("mumps_test_request_c should not be called with strategy %d\n",mumps_io_flag_async);*/
    /* JY+EA: Allow for this option, since it is similar to wait_request
     * and wait_request is allowed in synchronous mode.
     * We always return TRUE.
     */
    *flag=1;
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *ierr=(MUMPS_INT)mumps_test_request_th(&request_id_loc,&flag_loc);
    *flag=(MUMPS_INT)flag_loc;
    break;
#endif
  default:
    *ierr=-92;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)mumps_io_flag_async);
    mumps_io_error((MUMPS_INT)*ierr,buf);
    return;
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return;
}
/* Waits for the termination of the request "request_id" */
void MUMPS_CALL
MUMPS_WAIT_REQUEST(MUMPS_INT *request_id,MUMPS_INT *ierr)
{
  char buf[64]; /* for error message */
  MUMPS_INT request_id_loc;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  request_id_loc=(MUMPS_INT)*request_id;
  if(*request_id==-1)
    return;
  switch(mumps_io_flag_async){
  case IO_SYNC:
    /* printf("mumps_wait_request should not be called with strategy %d\n",mumps_io_flag_async); */
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *ierr=(MUMPS_INT)mumps_wait_request_th(&request_id_loc);
    break;
#endif
  default:
    *ierr=-92;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)mumps_io_flag_async);
    mumps_io_error((MUMPS_INT)*ierr,buf);
    return;
    /*    printf("Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
          exit (-3);*/
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return;
}
/**
 * Inits the I/O OOC mechanism.
 * Because on some computers, file size is limited, the I/O
 * mechanism must be able to handle a multi-file access to data.
 * Hence, we compute mumps_io_nb_file, which is the the number of files
 * we estimate we need.
 * Because of not exact matching between data packets written and size
 * of files, the recoverment may be imperfect. Consequently, we must
 * be able to reallocate if necessary.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_PREFIX(MUMPS_INT *dim, char *str, mumps_ftnlen l1)
{
  MUMPS_INT i;
  MUMPS_OOC_STORE_PREFIXLEN = *dim;
  if( *dim > MUMPS_OOC_PREFIX_MAX_LENGTH )
      MUMPS_OOC_STORE_PREFIXLEN = MUMPS_OOC_PREFIX_MAX_LENGTH;
  for(i=0;i<MUMPS_OOC_STORE_PREFIXLEN;i++){
    MUMPS_OOC_STORE_PREFIX[i]=str[i];
  }
  return;
}
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_TMPDIR(MUMPS_INT *dim, char *str, mumps_ftnlen l1)
{
  MUMPS_INT i;
  MUMPS_OOC_STORE_TMPDIRLEN=*dim;
  if( *dim > MUMPS_OOC_TMPDIR_MAX_LENGTH )
      MUMPS_OOC_STORE_TMPDIRLEN = MUMPS_OOC_TMPDIR_MAX_LENGTH;
  for(i=0;i<MUMPS_OOC_STORE_TMPDIRLEN;i++){
    MUMPS_OOC_STORE_TMPDIR[i]=str[i];
  }
  return;
}
/* Computes the number of files needed. Uses ceil value. */
/*   mumps_io_nb_file=0; */
/*   mumps_io_last_file_opened=-1; */
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_OOC_C(MUMPS_INT *_myid, MUMPS_INT *total_size_io, MUMPS_INT *size_element,
                           MUMPS_INT *async, MUMPS_INT *k211, MUMPS_INT *nb_file_type,
                           MUMPS_INT *flag_tab, MUMPS_INT *ierr)
{
  char buf[128]; /* for error message */
  MUMPS_INT myid_loc,async_loc,size_element_loc,nb_file_type_loc,*flag_tab_loc;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  MUMPS_INT ierr_loc;
#endif
  long long total_size_io_loc;
  MUMPS_INT i;
  myid_loc=(MUMPS_INT)*_myid;
  async_loc=(MUMPS_INT)*async;
  total_size_io_loc=(long long)*total_size_io;
  size_element_loc=(MUMPS_INT)*size_element;
  nb_file_type_loc=(MUMPS_INT)*nb_file_type;
  flag_tab_loc=(MUMPS_INT *)malloc(nb_file_type_loc*sizeof(MUMPS_INT));
  for (i=0;i<nb_file_type_loc;i++){
    flag_tab_loc[i]=(MUMPS_INT)flag_tab[i];
  }
#if defined(MUMPS_WIN32)
  if(async_loc==IO_ASYNC_AIO||async_loc==IO_ASYNC_TH){
    mumps_io_is_init_called=0;
    *ierr=-92;
    mumps_io_error((MUMPS_INT)*ierr,"Error: Forbidden value of Async flag with _WIN32\n");
    free(flag_tab_loc);
    return;
  }
#endif
#if defined (WITHOUT_PTHREAD)
  if(async_loc==IO_ASYNC_TH){
    mumps_io_is_init_called=0;
    *ierr=-92;
    mumps_io_error((MUMPS_INT)*ierr,"Error: Forbidden value of Async flag with WITHOUT_PTHREAD\n");
    free(flag_tab_loc);
    return;
  }
#endif
  total_vol=0;
  mumps_io_flag_async=async_loc;
  mumps_io_k211=(MUMPS_INT)*k211;
  if (MUMPS_OOC_STORE_PREFIXLEN==-1) {
    *ierr=-92;
    mumps_io_error((MUMPS_INT)*ierr,"Error: prefix not initialized\n");
    free(flag_tab_loc);
    return;
  }
  if (MUMPS_OOC_STORE_TMPDIRLEN==-1) {
    *ierr=-92;
    mumps_io_error((MUMPS_INT)*ierr,"Error: tmpdir not initialized\n");
    free(flag_tab_loc);
    return;
  }
  *ierr=(MUMPS_INT)mumps_init_file_name(MUMPS_OOC_STORE_TMPDIR, MUMPS_OOC_STORE_PREFIX,
                             &MUMPS_OOC_STORE_TMPDIRLEN, &MUMPS_OOC_STORE_PREFIXLEN, &myid_loc);
  if(*ierr<0){
    free(flag_tab_loc);
    return;
  }
  /* Re-initialize lenghts to -1 in order to enable the
   * check on initialization next time this routine is called
   */
  MUMPS_OOC_STORE_PREFIXLEN=-1;
  MUMPS_OOC_STORE_TMPDIRLEN=-1;
  *ierr=(MUMPS_INT)mumps_init_file_structure(&myid_loc,&total_size_io_loc,&size_element_loc,&nb_file_type_loc,flag_tab_loc);
  free(flag_tab_loc);
  if(*ierr<0){
    return;
  }
#if ! defined(MUMPS_WIN32)
  mumps_time_spent_in_sync=0;
#endif
  if(async_loc){
    switch(async_loc){
    case IO_SYNC:
      printf("mumps_low_level_init_ooc_c should not be called with strategy %d\n",(int)mumps_io_flag_async);
      break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      mumps_low_level_init_ooc_c_th(&async_loc,&ierr_loc);
      *ierr=(MUMPS_INT)ierr_loc;
      if(*ierr<0){
        return;
      }
      break;
#endif
    default:
      *ierr=-92;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)*async);
      mumps_io_error((MUMPS_INT)*ierr,buf);
      return;
    }
  }
  mumps_io_is_init_called=1;
  return;
}
/**
 * Writes a contigous block of central memory to the disk.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_WRITE_OOC_C(const MUMPS_INT * strat_IO,
                            void * address_block,
                            MUMPS_INT * block_size_int1,
                            MUMPS_INT * block_size_int2,
                            MUMPS_INT * inode,
                            MUMPS_INT * request_arg,
                            MUMPS_INT * type,
                            MUMPS_INT * vaddr_int1,
                            MUMPS_INT * vaddr_int2,
                            MUMPS_INT * ierr)
{
  MUMPS_INT ret_code=0;
  long long vaddr,block_size;
  char buf[64]; /* for error message */
  MUMPS_INT inode_loc,request_arg_loc,type_loc,ierr_loc,strat_IO_loc;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  inode_loc=(MUMPS_INT)*inode;
  request_arg_loc=(MUMPS_INT)*request_arg;
  type_loc=(MUMPS_INT)*type;
  ierr_loc=(MUMPS_INT)*ierr;
  strat_IO_loc=(MUMPS_INT)*strat_IO;
/* JY 27/2/08: initialize *request_arg to -1 (null request).
 * There were problems of uninitialized requests in the Fortran
 * code. For example when we use the synchronous version, there are
 * still some tests on *request_arg, which is not initialized.*/
  *request_arg=-1;
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
  if(mumps_io_flag_async){
    switch(*strat_IO){
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      ret_code=mumps_async_write_th(&strat_IO_loc, address_block, block_size,&inode_loc,&request_arg_loc,&type_loc,vaddr,&ierr_loc);
      *ierr=(MUMPS_INT)ierr_loc;
      *request_arg=(MUMPS_INT)request_arg_loc;
      if(ret_code<0){
        *ierr=ret_code;
      }
      break;
#endif
    default:
      *ierr=-91;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)*strat_IO);
      mumps_io_error((MUMPS_INT)*ierr,buf);
      return;
    }
  } else {
    ret_code=mumps_io_do_write_block(address_block,block_size,&type_loc,vaddr,&ierr_loc);
    *ierr=(MUMPS_INT)ierr_loc;
    if(ret_code<0){
      *ierr=ret_code;
    }
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  write_op_vol=write_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/**
 * Reads  a contigous block of central memory from the disk.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_READ_OOC_C(const MUMPS_INT * strat_IO,
                           void * address_block,
                           MUMPS_INT * block_size_int1,
                           MUMPS_INT * block_size_int2,
                           MUMPS_INT * inode,
                           MUMPS_INT * request_arg,
                           MUMPS_INT * type,
                           MUMPS_INT * vaddr_int1,
                           MUMPS_INT * vaddr_int2,
                           MUMPS_INT * ierr)
{
  char buf[64]; /* for error message */
  long long vaddr,block_size;
  MUMPS_INT inode_loc,request_arg_loc,type_loc,ierr_loc,strat_IO_loc;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  inode_loc=(MUMPS_INT)*inode;
  request_arg_loc=(MUMPS_INT)*request_arg;
  type_loc=(MUMPS_INT)*type;
  ierr_loc=(MUMPS_INT)*ierr;
  strat_IO_loc=(MUMPS_INT)*strat_IO;  
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
  if(mumps_io_flag_async){
      switch(*strat_IO){
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
      case IO_ASYNC_TH:
        mumps_async_read_th(&strat_IO_loc,address_block,block_size,&inode_loc,&request_arg_loc,&type_loc,vaddr,&ierr_loc);
        *ierr=(MUMPS_INT)ierr_loc;
        *request_arg=(MUMPS_INT)request_arg_loc;
        break;
#endif
      default:
        *ierr=-91;
        sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)*strat_IO);
        mumps_io_error((MUMPS_INT)*ierr,buf);
        return;
      }
  }else{
    mumps_io_do_read_block(address_block,block_size,&type_loc,vaddr,&ierr_loc);
    *ierr=(MUMPS_INT)ierr_loc;
    *request_arg=1;
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  read_op_vol=read_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/* Emergency read from the MUMPS thread during the solve phase.*/
void MUMPS_CALL
MUMPS_LOW_LEVEL_DIRECT_READ(void * address_block,
                            MUMPS_INT * block_size_int1,
                            MUMPS_INT * block_size_int2,
                            MUMPS_INT * type,
                            MUMPS_INT * vaddr_int1,
                            MUMPS_INT * vaddr_int2,
                            MUMPS_INT * ierr)
{
    /*  MUMPS_INT ret_code=0; */
  long long vaddr,block_size;
  MUMPS_INT type_loc,ierr_loc;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  type_loc=(MUMPS_INT)*type;
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    if(mumps_io_flag_async==IO_ASYNC_TH||mumps_io_flag_async==0)
#else
    if(mumps_io_flag_async==0)
#endif
    {
      ierr_loc=mumps_io_do_read_block(address_block,block_size,&type_loc,vaddr,&ierr_loc);
      *ierr=(MUMPS_INT)ierr_loc;
      if(*ierr<0){
         return;
      }
    } else {
    }
#if ! defined(MUMPS_WIN32)
# if ! defined(WITHOUT_PTHREAD)
# endif
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  read_op_vol=read_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/* Cleans the thread/io management data*/
void MUMPS_CALL
MUMPS_CLEAN_IO_DATA_C(MUMPS_INT *myid,MUMPS_INT *step,MUMPS_INT *ierr)
{
  char buf[64]; /* for error message */
  MUMPS_INT step_loc,myid_loc,ierr_loc;
  step_loc=(MUMPS_INT)*step;
  myid_loc=(MUMPS_INT)*myid;
  if(!mumps_io_is_init_called){
    return;
  }
  switch(mumps_io_flag_async){
  case IO_SYNC:
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    ierr_loc=mumps_clean_io_data_c_th(&myid_loc);
    *ierr=(MUMPS_INT)ierr_loc;
    break;
#endif
  default:
    *ierr=-91;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)mumps_io_flag_async);
    mumps_io_error((MUMPS_INT)*ierr,buf);
    return;
  }
  mumps_free_file_pointers(&step_loc);
  mumps_io_is_init_called=0;
  return;
}
void MUMPS_CALL
MUMPS_OOC_PRINT_STATS()
{
#if ! defined(MUMPS_WIN32)
  printf("%d: total time spent in i/o mode = %lf\n",(int)mumps_io_myid,mumps_time_spent_in_sync);
#endif
  printf("%d: Volume of read i/o = %lf\n",(int)mumps_io_myid,read_op_vol);
  printf("%d: Volume of write i/o = %lf\n",(int)mumps_io_myid,write_op_vol);
  total_vol=total_vol+read_op_vol+write_op_vol;
  printf("%d: Total i/o volume = %lf\n",(int)mumps_io_myid,total_vol);
  return;
}
void MUMPS_CALL
MUMPS_GET_MAX_NB_REQ_C(MUMPS_INT *max,MUMPS_INT *ierr)
{
  char buf[64]; /* for error message */
  *ierr=0;
  switch(mumps_io_flag_async){
  case IO_SYNC:
    *max=1;
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *max=MAX_FINISH_REQ+MAX_IO;
    break;
#endif
  default:
    *ierr=-91;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)mumps_io_flag_async);
    mumps_io_error((MUMPS_INT)*ierr,buf);
    return;
  }
  return;
}
void MUMPS_CALL
MUMPS_GET_MAX_FILE_SIZE_C(double * max_ooc_file_size)
{
  *max_ooc_file_size=(double)(MAX_FILE_SIZE);
  return;
}
void MUMPS_CALL
MUMPS_OOC_GET_NB_FILES_C(const MUMPS_INT *type,MUMPS_INT *nb_files)
{
  MUMPS_INT type_loc,nb_files_loc;
  type_loc=(MUMPS_INT)*type;
  mumps_io_get_nb_files(&nb_files_loc,&type_loc);
  *nb_files=(MUMPS_INT)nb_files_loc;
  return;
}
void MUMPS_CALL
MUMPS_OOC_GET_FILE_NAME_C(MUMPS_INT *type,MUMPS_INT *indice,MUMPS_INT *length, char* name, mumps_ftnlen l1)
{
  MUMPS_INT type_loc,indice_loc,length_loc;
  type_loc=(MUMPS_INT)*type;
  indice_loc=(MUMPS_INT)*indice;
  mumps_io_get_file_name(&indice_loc,name,&length_loc,&type_loc);
  *length=(MUMPS_INT)length_loc;
  return;
}
void MUMPS_CALL
MUMPS_OOC_SET_FILE_NAME_C(MUMPS_INT *type, MUMPS_INT *indice, MUMPS_INT *length, MUMPS_INT *ierr,
                          char* name, mumps_ftnlen l1)
{
  MUMPS_INT type_loc,indice_loc,length_loc;
  type_loc=(MUMPS_INT)*type;
  indice_loc=(MUMPS_INT)*indice;
  length_loc=(MUMPS_INT)*length;
  *ierr=(MUMPS_INT)mumps_io_set_file_name(&indice_loc,name,&length_loc,&type_loc);
  return;
}
void MUMPS_CALL
MUMPS_OOC_ALLOC_POINTERS_C(MUMPS_INT *nb_file_type,MUMPS_INT *dim,MUMPS_INT *ierr)
{
  MUMPS_INT i=0;
  MUMPS_INT nb_file_type_loc, *dim_loc;
  nb_file_type_loc=(MUMPS_INT)*nb_file_type;  
  dim_loc=(MUMPS_INT *)malloc(*nb_file_type*sizeof(MUMPS_INT));
  for(i=0;i<nb_file_type_loc;i++){
    dim_loc[i]=(MUMPS_INT)dim[i];
  }
  *ierr=(MUMPS_INT)mumps_io_alloc_pointers(&nb_file_type_loc,dim_loc);
  for(i=0;i<nb_file_type_loc;i++){
    mumps_io_set_last_file(dim_loc+i,&i);
  }
  free(dim_loc);
  return;
}
void MUMPS_CALL
MUMPS_OOC_INIT_VARS_C(MUMPS_INT *myid_arg,
                        MUMPS_INT *size_element,MUMPS_INT *async, MUMPS_INT *k211,
                        MUMPS_INT *ierr)
{
  MUMPS_INT size_element_loc,async_loc,myid_arg_loc;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  mumps_time_spent_in_sync=0;
#endif
  mumps_io_k211=(MUMPS_INT)*k211;
  size_element_loc=(MUMPS_INT)*size_element;
  async_loc=(MUMPS_INT)*async;
  myid_arg_loc=(MUMPS_INT)*myid_arg;
  *ierr=(MUMPS_INT)mumps_io_init_vars(&myid_arg_loc,&size_element_loc,&async_loc);
  return;
}
void MUMPS_CALL
MUMPS_OOC_START_LOW_LEVEL(MUMPS_INT *ierr)
{
  char buf[64]; /* for error message */
  MUMPS_INT ierr_loc;
  read_op_vol=0;
  write_op_vol=0;
  *ierr=(MUMPS_INT)mumps_io_open_files_for_read();
  if(*ierr<0){
    return;
  }
  if(mumps_io_flag_async){
    switch(mumps_io_flag_async){
    case IO_SYNC:
      break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      mumps_low_level_init_ooc_c_th(&mumps_io_flag_async,&ierr_loc);
      *ierr=(MUMPS_INT)ierr_loc;
      if(*ierr<0){
        return;
      }
      break;
#endif
    default:
      *ierr=-91;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",(int)mumps_io_flag_async);
      mumps_io_error((MUMPS_INT)*ierr,buf);
      return;
    }
  }
  mumps_io_is_init_called=1;
  return;
}
void MUMPS_CALL
MUMPS_OOC_REMOVE_FILE_C(MUMPS_INT *ierr, char *name, mumps_ftnlen l1)
{
  char buf[256]; /* for error message, count 256 chars for name */
  *ierr=(MUMPS_INT)remove(name);
  if(*ierr<0){
#if ! defined(MUMPS_WIN32)
    sprintf(buf,"Unable to remove OOC file %s",name);
#else
    sprintf(buf,"Unable to remove OOC file %s with return value %d",name,*ierr);
#endif
    *ierr = -90;
    mumps_io_sys_error((MUMPS_INT)*ierr,buf);
    return;
  }
  return;
}
void MUMPS_CALL
MUMPS_OOC_END_WRITE_C(MUMPS_INT *ierr)
{
  return;
}
void MUMPS_CALL
MUMPS_OOC_IS_ASYNC_AVAIL(MUMPS_INT *flag)
{
#if ( defined (WITHOUT_PTHREAD) || defined(MUMPS_WIN32) ) && ! defined(WITH_AIO)
  *flag=0;
#else
  *flag=1;
#endif
}
/**
 * IMPORTANT:
 * ==========
 *
 *   After every modification of the code of convert_2fint_to_longlong update
 *   the corresponding fortran subroutines MUMPS_OOC_CONVERT_2INTTOVADDR
 *   and MUMPS_OOC_CONVERT_VADDRTO2INT
 */
MUMPS_INLINE MUMPS_INT
mumps_convert_2fint_to_longlong( MUMPS_INT *short_int1, MUMPS_INT *short_int2,
                                 long long * long_int )
{
  *long_int=((long long)(*short_int1)*((long long)1073741824))+(long long)(*short_int2);
  return 0;
}
