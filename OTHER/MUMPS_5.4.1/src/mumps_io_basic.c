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
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
#include "mumps_c_types.h"
/* Exported global variables */
#if ! defined (MUMPS_WIN32)
# if defined(WITH_PFUNC) && ! defined (WITHOUT_PTHREAD)
#  include <pthread.h>
pthread_mutex_t mumps_io_pwrite_mutex;
# endif
/* int* mumps_io_pfile_pointer_array; */
/* int* mumps_io_current_file; */
/* #else /\*MUMPS_WIN32*\/ */
/* FILE** mumps_io_current_file; */
/* FILE** mumps_io_pfile_pointer_array; */
#endif /* MUMPS_WIN32 */
/* mumps_file_struct* mumps_io_pfile_pointer_array;
   mumps_file_struct* mumps_io_current_file; */
mumps_file_type* mumps_files = NULL;
/* int mumps_io_current_file_number; */
char* mumps_ooc_file_prefix = NULL;
/* char** mumps_io_pfile_name; */
/* int mumps_io_current_file_position; */
/* int mumps_io_write_pos; */
/* int mumps_io_last_file_opened; */
MUMPS_INT mumps_elementary_data_size;
MUMPS_INT mumps_io_is_init_called;
MUMPS_INT mumps_io_myid;
MUMPS_INT mumps_io_max_file_size;
/* int mumps_io_nb_file; */
MUMPS_INT mumps_io_flag_async;
MUMPS_INT mumps_io_k211;
/* int mumps_flag_open;*/
MUMPS_INT mumps_directio_flag;
MUMPS_INT mumps_io_nb_file_type;
/* Functions */
MUMPS_INT mumps_set_file(MUMPS_INT type,MUMPS_INT file_number_arg){
  /* Defines the pattern for the file name. The last 6 'X' will be replaced
     so as to name were unique */
  char name[351];
#if ! defined(_WIN32)
  MUMPS_INT fd;
  char buf[64]; /* for error message */
#endif
  mumps_file_struct  *mumps_io_pfile_pointer_array;
  if (file_number_arg > ((mumps_files+type)->mumps_io_nb_file)-1){
    /* mumps_io_nb_file was initialized to the estimated number
    of files inside mumps_io_init_file_struct; this block is
    entered in case of a too small estimation of the required
    number of files. */
    /* We increase the number of files needed and then realloc. */
    ((mumps_files+type)->mumps_io_nb_file)++;
    (mumps_files+type)->mumps_io_pfile_pointer_array=(mumps_file_struct*)realloc((void *)(mumps_files+type)->mumps_io_pfile_pointer_array,((mumps_files+type)->mumps_io_nb_file)*sizeof(mumps_file_struct));
    /* Check for reallocation problem */
    if((mumps_files+type)->mumps_io_pfile_pointer_array==NULL){
      return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
    }
    /* initialize "is_opened", as in mumps_io_init_file_struct */
    ((mumps_files+type)->mumps_io_pfile_pointer_array+((mumps_files+type)->mumps_io_nb_file)-1)->is_opened = 0;
  }
  mumps_io_pfile_pointer_array=(mumps_files+type)->mumps_io_pfile_pointer_array;
  /* 
     Do change the current file 
     Careful: both mumps_io_current_file_number and
              mumps_io_current_file must be changed
   */
  ((mumps_files+type)->mumps_io_current_file_number)=file_number_arg;
  ((mumps_files+type)->mumps_io_current_file)=mumps_io_pfile_pointer_array+file_number_arg;
  if((mumps_io_pfile_pointer_array+file_number_arg)->is_opened!=0){
    /*
       The file already exists and is open.
       The i/o will be performed in the current file (which may not
       be the last one.
     */
    return 0;
  }
/*********************/
/* CREATE A NEW FILE */
/*********************/
/* #if ! defined( MUMPS_WIN32 )*/
/* MinGW does not have a mkstemp function and MinGW defines _WIN32,
 * so we also go in the else branch below with MinGW */
#if ! defined(_WIN32)
  strcpy(name,mumps_ooc_file_prefix);
  fd=mkstemp(name);  
  /* Note that a file name is built by mkstemp and that the file is 
     opened. fd hold the file descriptor to access it.
     We want to close the file that will be opened later
     and might be removed before the end of the processus.
  */
  if(fd < 0) {
    sprintf(buf,"File creation failure");
    return mumps_io_sys_error(-90,buf);
  } else {
    close(fd); 
  }
#else
  sprintf(name,"%s_%d_%d",mumps_ooc_file_prefix,((mumps_files+type)->mumps_io_current_file_number)+1,type);
#endif
/*   *(mumps_io_pfile_pointer_array+mumps_io_current_file_number)=fopen(name,"w+"); */
/*   *(mumps_io_pfile_name+mumps_io_current_file_number)=(char *)malloc((strlen(name)+1)*sizeof(char)); */
/*   if(*(mumps_io_pfile_name+mumps_io_current_file_number)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/*     return -13; */
/*   } */
  strcpy((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->name,name);
  /* See mumps_io_basic.h for comments on the I/O flags passed to open */
#if ! defined( MUMPS_WIN32 )
  (mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file=open(name,(mumps_files+type)->mumps_flag_open,0666); 
  /* 
CPA: for LU factor file: 
(mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file= open(name, O_WRONLY | O_CREAT | O_TRUNC, 0666); */
  if((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file==-1){
    return mumps_io_sys_error(-90,"Unable to open OOC file");
  }
#else
  (mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file=fopen(name,(mumps_files+type)->mumps_flag_open);
  if((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file==NULL){
    return mumps_io_error(-90,"Problem while opening OOC file");
  }
#endif
  (mumps_files+type)->mumps_io_current_file=(mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number);
  ((mumps_files+type)->mumps_io_nb_file_opened)++;
  if((mumps_files+type)->mumps_io_current_file_number>(mumps_files+type)->mumps_io_last_file_opened){
    (mumps_files+type)->mumps_io_last_file_opened=(mumps_files+type)->mumps_io_current_file_number;
  }
  /*  if(*(mumps_io_pfile_pointer_array+mumps_io_current_file_number)==NULL) */
  ((mumps_files+type)->mumps_io_current_file)->write_pos=0;
  ((mumps_files+type)->mumps_io_current_file)->is_opened=1;
  /*  printf("new file created -> num = %d \n", ((mumps_files+type)->mumps_io_last_file_opened));*/
  /*  printf("new file created %d\n",mumps_io_current_file_number);*/
  return 0;
}
void mumps_update_current_file_position(mumps_file_struct* file_arg){
  file_arg->current_pos=file_arg->write_pos;
/*   mumps_io_current_file_position=mumps_io_write_pos; */
}
MUMPS_INT mumps_compute_where_to_write(const double to_be_written,const MUMPS_INT type,long long vaddr,size_t already_written){
  /* Check if the current file has enough memory to receive the whole block*/
  MUMPS_INT ret_code;
  MUMPS_INT file;
  mumps_file_struct *current_file;
  long long vaddr_loc;
  MUMPS_INT pos;
  /* Virtual address based file management scheme */
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size+(long long)already_written;
  mumps_gen_file_info(vaddr_loc,&pos,&file);
  ret_code=mumps_set_file(type,file);
  if(ret_code<0){
    return ret_code;
  }
  current_file=(mumps_files+type)->mumps_io_current_file;
  current_file->write_pos=pos;
  mumps_update_current_file_position(current_file);
  return 0;
}
MUMPS_INT mumps_prepare_pointers_for_write(double to_be_written,MUMPS_INT * pos_in_file, MUMPS_INT * file_number,const MUMPS_INT type,long long vaddr,size_t already_written){
  MUMPS_INT ret_code;
  ret_code=mumps_compute_where_to_write(to_be_written,type,vaddr,already_written);
  if(ret_code<0){
    return ret_code;
  }
  *pos_in_file=((mumps_files+type)->mumps_io_current_file)->current_pos;
  /* should be modified to take into account the file arg */
  *file_number=(mumps_files+type)->mumps_io_current_file_number;
  return 0;
}
MUMPS_INLINE MUMPS_INT mumps_gen_file_info(long long vaddr, MUMPS_INT * pos, MUMPS_INT * file){
  *file=(MUMPS_INT)(vaddr/(long long)mumps_io_max_file_size);
  *pos=(MUMPS_INT)(vaddr%(long long)mumps_io_max_file_size);
  return 0;
}
MUMPS_INT mumps_compute_nb_concerned_files(long long block_size, MUMPS_INT * nb_concerned_files,long long vaddr){
  MUMPS_INT file,pos,available_size;
  long long vaddr_loc;
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size;
  mumps_gen_file_info(vaddr_loc,&pos,&file);
  available_size=mumps_io_max_file_size-pos+1;
  *nb_concerned_files=(MUMPS_INT)my_ceil((double)(my_max(0,((block_size)*(double)(mumps_elementary_data_size))-available_size))/(double)mumps_io_max_file_size)+1;
  return 0;
}
MUMPS_INT mumps_io_do_write_block(void * address_block,
                            long long block_size,
                            MUMPS_INT * type_arg,
                            long long vaddr,
                            MUMPS_INT * ierr){   
  /* Type of fwrite : size_t fwrite(const void *ptr, size_t size, 
                                    *size_t nmemb, FILE *stream); */
  size_t write_size;
  MUMPS_INT i;
  MUMPS_INT nb_concerned_files=0;
  MUMPS_INT ret_code,file_number_loc,pos_in_file_loc;
  double to_be_written;
#if ! defined( MUMPS_WIN32 )
  MUMPS_INT* file;
#else
  FILE** file;
#endif
  MUMPS_INT where;
  void* loc_addr;
  MUMPS_INT type;
  size_t already_written=0;
  char buf[64];
  type=*type_arg;
  loc_addr=address_block;
  mumps_compute_nb_concerned_files(block_size,&nb_concerned_files,vaddr);
  to_be_written=((double)mumps_elementary_data_size)*((double)(block_size));
  /*  printf("nb_concerned -> %d | %lf \n",nb_concerned_files,to_be_written); */
  for(i=0;i<nb_concerned_files;i++){
#if ! defined( MUMPS_WIN32 ) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_protect_pointers();
    }
# endif
#endif
    ret_code=mumps_prepare_pointers_for_write(to_be_written,&pos_in_file_loc,&file_number_loc,type,vaddr,already_written);
    if(ret_code<0){
#if ! defined( MUMPS_WIN32 ) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
# endif
#endif
      return ret_code;
    }
    if((double)(mumps_io_max_file_size-((mumps_files+type)->mumps_io_current_file)->write_pos)>to_be_written){
      write_size=(size_t)to_be_written;
      already_written=(size_t)to_be_written;
    }else{
      write_size=(size_t)((double)(mumps_io_max_file_size-((mumps_files+type)->mumps_io_current_file)->write_pos));
      already_written=already_written+(size_t)write_size;
    }
#if defined( MUMPS_WIN32 )
    write_size=(size_t)(MUMPS_INT)((write_size)/mumps_elementary_data_size);
#endif
    file=&(((mumps_files+type)->mumps_io_current_file)->file);
    where=((mumps_files+type)->mumps_io_current_file)->write_pos;
#if ! defined( MUMPS_WIN32 ) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
# endif
#endif
    ret_code=mumps_io_write__(file,loc_addr,write_size,where,type);
    if(ret_code<0){
      return ret_code;
    }
#if ! defined( MUMPS_WIN32 ) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_protect_pointers();
    }
# endif
#endif
#if ! defined( MUMPS_WIN32 )
    ((mumps_files+type)->mumps_io_current_file)->write_pos=((mumps_files+type)->mumps_io_current_file)->write_pos+((MUMPS_INT)write_size);
    to_be_written=to_be_written-((MUMPS_INT)write_size);
    loc_addr=(void*)((size_t)loc_addr+write_size);
/*     mumps_io_write_pos=mumps_io_write_pos+((MUMPS_INT)write_size); */
/*     to_be_written=to_be_written-((MUMPS_INT)write_size); */
/*     loc_addr=(void*)((size_t)loc_addr+write_size); */
#else
    /* fread and write */
    ((mumps_files+type)->mumps_io_current_file)->write_pos=((mumps_files+type)->mumps_io_current_file)->write_pos+((MUMPS_INT)write_size*mumps_elementary_data_size);
    to_be_written=to_be_written-((MUMPS_INT)write_size*mumps_elementary_data_size);
    loc_addr=(void*)((size_t)loc_addr+(size_t)((MUMPS_INT)write_size*mumps_elementary_data_size));
/*     mumps_io_write_pos=mumps_io_write_pos+((MUMPS_INT)write_size*mumps_elementary_data_size); */
/*     to_be_written=to_be_written-((MUMPS_INT)write_size*mumps_elementary_data_size); */
/*     loc_addr=(void*)((size_t)loc_addr+(size_t)((MUMPS_INT)write_size*mumps_elementary_data_size)); */
#endif
#if ! defined( MUMPS_WIN32 ) && ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
# endif
#endif
  }
  if(to_be_written!=0){
    *ierr = -90;
    sprintf(buf,"Internal (1) error in low-level I/O operation %lf",to_be_written);
    return mumps_io_error(*ierr,buf);
  }
  /* printf("write ok -> %d \n");*/
  return 0;
}
MUMPS_INT mumps_io_do_read_block(void * address_block,
                    long long block_size,
                    MUMPS_INT * type_arg,
                     long long vaddr,
                    MUMPS_INT * ierr){
  size_t size;
#if ! defined( MUMPS_WIN32 )
  MUMPS_INT* file;
#else
  FILE** file;
#endif
  double read_size;
  MUMPS_INT local_fnum,local_offset;
  void *loc_addr;
  long long vaddr_loc;
  MUMPS_INT type;
  type=*type_arg;
  /*  if(((double)(*block_size))*((double)(mumps_elementary_data_size))>(double)mumps_io_max_file_size){
    sprintf(error_str,"Internal error in low-level I/O operation (requested size too big for file system) \n");
    return -90;
    }*/
  if(block_size==0){
    return 0;
  }
  read_size=(double)mumps_elementary_data_size*(double)(block_size);
  /*  if((*file_number<0)&&(read_size<(double)mumps_io_max_file_size)){
    sprintf(error_str,"Internal error (1) in low level read op\n");
    return -90;
    }*/
  loc_addr=address_block;
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size;
  while(read_size>0){
    /* Virtual addressing based management stuff */
    local_fnum=(MUMPS_INT)(vaddr_loc/(long long)mumps_io_max_file_size);
    local_offset=(MUMPS_INT)(vaddr_loc%(long long)mumps_io_max_file_size);
    file=&((((mumps_files+type)->mumps_io_pfile_pointer_array)+local_fnum)->file);
    /* printf("1 read | file -> %d | fnum -> %d | vaddr -> %d \n",*file,local_fnum,(MUMPS_INT)vaddr_loc); */
#if ! defined( MUMPS_WIN32 )
    if(read_size+(double)local_offset>(double)mumps_io_max_file_size){
      size=(size_t)mumps_io_max_file_size-(size_t)local_offset;
    }else{
      size=(size_t)read_size;
    }
#else
    if(read_size+(double)local_offset>(double)mumps_io_max_file_size){
      size=((size_t)mumps_io_max_file_size-(size_t)local_offset)/(size_t)mumps_elementary_data_size;
    }else{
      size=(size_t)(read_size/mumps_elementary_data_size);
    }
#endif
    *ierr=mumps_io_read__(file,loc_addr,size,local_offset,type);
    if(*ierr<0){
      return *ierr;
    }
#if defined( MUMPS_WIN32 )
    size=size*mumps_elementary_data_size;
#endif
    vaddr_loc=vaddr_loc+(long long)size;
    read_size=read_size-(double)size;
    loc_addr=(void*)((size_t)loc_addr+size);
    local_fnum++;
    local_offset=0;
    if(local_fnum>(mumps_files+type)->mumps_io_nb_file){
      *ierr = -90;
      return mumps_io_error(*ierr,"Internal error (2) in low level read op\n");
    }
  }
  return 0;
}
MUMPS_INT mumps_free_file_pointers(MUMPS_INT *step){
  MUMPS_INT i,j,bound,ierr;
/*   Free prefix only for facto  */
  if (*step == 0) free(mumps_ooc_file_prefix);
  if(mumps_files == NULL )
      return 0;
#if ! defined( MUMPS_WIN32 )
#endif
  bound=mumps_io_nb_file_type;
/*   if(*step==0){ */
/*     /\* factorization *\/ */
/*     bound=NB_FILE_TYPE_FACTO; */
/*   }else{ */
/*     /\* solve *\/ */
/*     bound=NB_FILE_TYPE_SOLVE; */
/*   } */
  for(j=0;j<bound;j++){
    if( mumps_files[j].mumps_io_pfile_pointer_array == NULL ) {
      continue;
    }
    for(i=0;i<(mumps_files+j)->mumps_io_nb_file_opened;i++){
#if ! defined( MUMPS_WIN32 )
      ierr=close((((mumps_files+j)->mumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
        return mumps_io_sys_error(-90,"Problem while closing OOC file");
      }
#else
      ierr=fclose((((mumps_files+j)->mumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
        return mumps_io_error(-90,"Problem while closing OOC file\n");
      }    
#endif
      /*     free(*(mumps_io_pfile_name+i)); */
    }
    free((mumps_files+j)->mumps_io_pfile_pointer_array);
  }
/*   free(mumps_io_pfile_name); */
  free(mumps_files);
#if ! defined( MUMPS_WIN32 )
#endif
  return 0;
}
/* Initialize the mumps_file_type structure at <which>th position in
   mumps_files. It only set values with no allocation to avoid any errors. */
void mumps_io_init_file_struct(MUMPS_INT* nb,MUMPS_INT which)
{
  (mumps_files+which)->mumps_io_current_file_number = -1;
  (mumps_files+which)->mumps_io_last_file_opened = -1;
  (mumps_files+which)->mumps_io_nb_file_opened = 0;
  (mumps_files+which)->mumps_io_nb_file=*nb;
  (mumps_files+which)->mumps_io_pfile_pointer_array = NULL;
  (mumps_files+which)->mumps_io_current_file=NULL;
}
/* Allocate the file structures for factor files and initialize the is_opened filed to 0 */
MUMPS_INT mumps_io_alloc_file_struct(MUMPS_INT* nb,MUMPS_INT which)
{
  MUMPS_INT i;
  (mumps_files+which)->mumps_io_pfile_pointer_array=(mumps_file_struct *)malloc((*nb)*sizeof(mumps_file_struct));
  if((mumps_files+which)->mumps_io_pfile_pointer_array==NULL){
    return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
  }
  for(i=0;i<*nb;i++){
    (((mumps_files+which)->mumps_io_pfile_pointer_array)+i)->is_opened=0;
  }
  return 0;
}
MUMPS_INT mumps_init_file_structure(MUMPS_INT* _myid, long long *total_size_io,MUMPS_INT *size_element,MUMPS_INT *nb_file_type,MUMPS_INT *flag_tab)
{
  /* Computes the number of files needed. Uses ceil value. */
  MUMPS_INT ierr;
#if ! defined( MUMPS_WIN32 )
  MUMPS_INT mumps_flag_open;
#endif
  MUMPS_INT i,nb;
  MUMPS_INT mumps_io_nb_file;
  mumps_io_max_file_size=MAX_FILE_SIZE;
  mumps_io_nb_file_type=*nb_file_type;
  mumps_io_nb_file=(MUMPS_INT)((((double)(*total_size_io)*1000000)*((double)(*size_element)))/(double)mumps_io_max_file_size)+1;
  mumps_directio_flag=0;
#if ! defined( MUMPS_WIN32 )
  mumps_flag_open=0;
#endif
  mumps_io_myid=*_myid;
  mumps_elementary_data_size=*size_element;
  /* Allocates the memory necessary to handle the file pointer array.*/
  mumps_files=(mumps_file_type *)malloc(mumps_io_nb_file_type*sizeof(mumps_file_type));
  if(mumps_files==NULL){
    return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
  }
  /* Safe initialization of the mumps_file_type elements */
  for(i=0;i<mumps_io_nb_file_type;i++){
    nb=(flag_tab[i]==0||flag_tab[i]==1) ? mumps_io_nb_file : 1;
    mumps_io_init_file_struct(&nb,i);
  }
  ierr=0;
  for(i=0;i<mumps_io_nb_file_type;i++){
    switch(flag_tab[i]){
    case 0:
#if ! defined( MUMPS_WIN32 )
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_WRONLY|O_CREAT|O_TRUNC;
#  if defined(__MINGW32__)
      /* O_BINARY necessary */
      (mumps_files+i)->mumps_flag_open=(mumps_files+i)->mumps_flag_open|O_BINARY;
#  endif
#else
      strcpy((mumps_files+i)->mumps_flag_open,"wb");
#endif
      break;
    case 1:
#if ! defined( MUMPS_WIN32 )
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDONLY|O_CREAT|O_TRUNC;
#  if defined(__MINGW32__)
      /* O_BINARY necessary */
      (mumps_files+i)->mumps_flag_open=(mumps_files+i)->mumps_flag_open|O_BINARY;
#  endif
#else
      strcpy((mumps_files+i)->mumps_flag_open,"rb");
#endif
      break;
    case 2:
#if ! defined( MUMPS_WIN32 )
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDWR|O_CREAT|O_TRUNC;
#  if defined(__MINGW32__)
      /* O_BINARY necessary */
      (mumps_files+i)->mumps_flag_open=(mumps_files+i)->mumps_flag_open|O_BINARY;
#  endif
#else
      strcpy((mumps_files+i)->mumps_flag_open,"rwb");
#endif
      break;
    default:
      return mumps_io_error(-90,"unknown value of flag_open\n");
    }
    ierr=mumps_io_alloc_file_struct(&nb,i);
    if(ierr<0){
      return ierr;
    }
    ierr=mumps_set_file(i,0);
    if(ierr<0){
      return ierr;
    }
  }
  /* Init the current file.*/
  return 0;
}
MUMPS_INT mumps_init_file_name(char* mumps_dir,char* mumps_file,
                         MUMPS_INT* mumps_dim_dir,MUMPS_INT* mumps_dim_file,MUMPS_INT* _myid){
  MUMPS_INT i;
  char *tmp_dir,*tmp_fname;
  char base_name[20];
  MUMPS_INT dir_flag=0,file_flag=0;
  char mumps_base[10]="mumps_";
  tmp_dir=(char *)malloc(((*mumps_dim_dir)+1)*sizeof(char));
  if(tmp_dir==NULL){
    return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
  }
  tmp_fname=(char *)malloc(((*mumps_dim_file)+1)*sizeof(char));
  if(tmp_fname==NULL){
    return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
  }
  for(i=0;i<*mumps_dim_dir;i++){
    tmp_dir[i]=mumps_dir[i];
  }
  tmp_dir[i]=0;
  for(i=0;i<*mumps_dim_file;i++){
    tmp_fname[i]=mumps_file[i];
  }
  tmp_fname[i]=0;  
  if(strcmp(tmp_dir,UNITIALIZED)==0){
    dir_flag=1;
    free(tmp_dir);
    tmp_dir=getenv("MUMPS_OOC_TMPDIR");
    if(tmp_dir==NULL){
#ifdef _AIX
# ifndef CINES_
      tmp_dir=getenv("TMPDIR");
      if(tmp_dir==NULL){
        tmp_dir=MUMPS_OOC_DEFAULT_DIR;
      }
# else
      tmp_dir=MUMPS_OOC_DEFAULT_DIR;       
# endif       
#else
      tmp_dir=MUMPS_OOC_DEFAULT_DIR;
#endif      
    }
  }
  if(strcmp(tmp_fname,UNITIALIZED)==0){
    free(tmp_fname);
    tmp_fname=getenv("MUMPS_OOC_PREFIX");
    file_flag=1;
  }
  if(tmp_fname!=NULL){
#if ! defined( MUMPS_WIN32 )
      sprintf(base_name,"_%s%d_XXXXXX",mumps_base,(int)*_myid);
#else
      sprintf(base_name,"_%s%d",mumps_base,*_myid);
#endif
      mumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(tmp_fname)+strlen(base_name)+1+1)*sizeof(char));
      if(mumps_ooc_file_prefix==NULL){
        return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
      }
      sprintf(mumps_ooc_file_prefix,"%s%s%s%s",tmp_dir,SEPARATOR,tmp_fname,base_name);
  }else{
#if ! defined( MUMPS_WIN32 )
    sprintf(base_name,"%s%s%d_XXXXXX",SEPARATOR,mumps_base,(int)*_myid);
#else
    sprintf(base_name,"%s%s%d",SEPARATOR,mumps_base,*_myid);
#endif
      mumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(base_name)+1)*sizeof(char));
      if(mumps_ooc_file_prefix==NULL){
        return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
      }
      sprintf(mumps_ooc_file_prefix,"%s%s%s",tmp_dir,SEPARATOR,base_name);
  }  
  if(!dir_flag){
    free(tmp_dir);
  }
  if(!file_flag){
    free(tmp_fname);
  }
  return 0;
}
MUMPS_INT mumps_io_get_nb_files(MUMPS_INT* nb_files, const MUMPS_INT* type){
  *nb_files=((mumps_files+*type)->mumps_io_last_file_opened)+1;
  return 0;
}
MUMPS_INT mumps_io_get_file_name(MUMPS_INT* indice,char* name,MUMPS_INT* length,MUMPS_INT* type){
  MUMPS_INT i;
  i=(*indice)-1;
  strcpy(name,(((mumps_files+*type)->mumps_io_pfile_pointer_array)+i)->name);
  *length=(MUMPS_INT)strlen(name)+1;
  return 0;  
}
MUMPS_INT mumps_io_alloc_pointers(MUMPS_INT* nb_file_type,MUMPS_INT * dim){
  MUMPS_INT ierr;
  MUMPS_INT i;
  /* This is called by solve step, we have only one type of files */
  mumps_io_nb_file_type=*nb_file_type;
  mumps_files=(mumps_file_type *)malloc(mumps_io_nb_file_type*sizeof(mumps_file_type));
  if(mumps_files==NULL){
    return mumps_io_error(-13,"Allocation problem in low-level OOC layer\n");
  }
  for(i=0;i<mumps_io_nb_file_type;i++){
    mumps_io_init_file_struct(dim+i,i);
    ierr=mumps_io_alloc_file_struct(dim+i,i);
    if(ierr<0){
      return ierr;
    }
  }
  return 0;
}
MUMPS_INT mumps_io_init_vars(MUMPS_INT* myid_arg,MUMPS_INT* size_element,MUMPS_INT* async_arg){
#if ! defined( MUMPS_WIN32 )
  MUMPS_INT mumps_flag_open;
#endif
  MUMPS_INT i;
  mumps_io_max_file_size=MAX_FILE_SIZE;
  mumps_directio_flag=0;
#if ! defined( MUMPS_WIN32 )
  mumps_flag_open=0;
#endif
  /* must be changed when we will have more than one file type during solve step */
  for(i=0;i<mumps_io_nb_file_type;i++){
#if ! defined( MUMPS_WIN32 )
    (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDONLY;
#else
    strcpy((mumps_files+i)->mumps_flag_open,"rb");
#endif
  }
  mumps_io_myid=*myid_arg;
  mumps_elementary_data_size=*size_element;
  mumps_io_flag_async=*async_arg;
  return 0;
}
MUMPS_INT mumps_io_set_file_name(MUMPS_INT* indice,char* name,MUMPS_INT* length,MUMPS_INT* type){
  MUMPS_INT i;
  i=(*indice)-1;
/*   *(mumps_io_pfile_name+i)=(char *) malloc((*length)*strlen(name)); */
/*   if(*(mumps_io_pfile_name+i)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer"); */
/*     return -13; */
/*   } */
  strcpy((((mumps_files+*type)->mumps_io_pfile_pointer_array)+i)->name,name);
  return 0;  
}
MUMPS_INT mumps_io_open_files_for_read(){
  MUMPS_INT i,j;
  mumps_file_struct  *mumps_io_pfile_pointer_array;
#if defined (sgi) || defined (__sgi)
  struct dioattr dio;
#endif
  for(j=0;j<mumps_io_nb_file_type;j++){
    mumps_io_pfile_pointer_array=(mumps_files+j)->mumps_io_pfile_pointer_array;
    for(i=0;i<(mumps_files+j)->mumps_io_nb_file;i++){
#if ! defined( MUMPS_WIN32 )
      (mumps_io_pfile_pointer_array+i)->file=open((mumps_io_pfile_pointer_array+i)->name,(mumps_files+j)->mumps_flag_open);
      if((mumps_io_pfile_pointer_array+i)->file==-1){
        return mumps_io_sys_error(-90,"Problem while opening OOC file");
      }
#else
      (mumps_io_pfile_pointer_array+i)->file=fopen((mumps_io_pfile_pointer_array+i)->name,(mumps_files+j)->mumps_flag_open);      
      if((mumps_io_pfile_pointer_array+i)->file==NULL){
        return mumps_io_error(-90,"Problem while opening OOC file");
      }
      (mumps_io_pfile_pointer_array+i)->is_opened=1;
#endif
    }
  }
  return 0;
}
MUMPS_INT mumps_io_set_last_file(MUMPS_INT* dim,MUMPS_INT* type){
  (mumps_files+*type)->mumps_io_last_file_opened=*dim-1;
  (mumps_files+*type)->mumps_io_nb_file_opened=*dim;
  return 0;
}
#if ! defined( MUMPS_WIN32 ) &&  ! defined (WITHOUT_PTHREAD)
# ifdef WITH_PFUNC
MUMPS_INT mumps_io_protect_pointers(){
  pthread_mutex_lock(&mumps_io_pwrite_mutex);
  return 0;
}
MUMPS_INT mumps_io_unprotect_pointers(){
  pthread_mutex_unlock(&mumps_io_pwrite_mutex);
  return 0;
}
MUMPS_INT mumps_io_init_pointers_lock(){
  pthread_mutex_init(&mumps_io_pwrite_mutex,NULL);
  return 0;
}
MUMPS_INT mumps_io_destroy_pointers_lock(){
  pthread_mutex_destroy(&mumps_io_pwrite_mutex);
  return 0;
}
# endif /*WITH_PFUNC*/
#endif /* _WIN32 && WITHOUT_PTHREAD */
 MUMPS_INT mumps_io_read__(void * file,void * loc_addr,size_t size,MUMPS_INT local_offset,MUMPS_INT type){
  MUMPS_INT ret_code;
#if ! defined( MUMPS_WIN32 )
  if(!mumps_directio_flag){
    ret_code=mumps_io_read_os_buff__(file,loc_addr, size,local_offset);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=mumps_io_read_win32__(file,loc_addr, size,local_offset);
  if(ret_code<0){
    return ret_code;
  }
#endif  
  return 0;
}
#if ! defined( MUMPS_WIN32 )
MUMPS_INT mumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,MUMPS_INT local_offset){
  size_t ret_code;
  /* printf("Read with buff %d %d %d\n",(MUMPS_INT) size, local_offset,*((MUMPS_INT *)file)); */
# ifdef WITH_PFUNC
  ret_code=pread(*(MUMPS_INT *)file,loc_addr,size,local_offset);
# else
  lseek(*(MUMPS_INT *)file,(long) local_offset,SEEK_SET);
  ret_code=read(*(MUMPS_INT *)file,loc_addr,size);
# endif
  if((MUMPS_INT) ret_code==-1){
    return mumps_io_sys_error(-90,"Problem with low level read");
  }
  return 0;
}
#endif
#if defined( MUMPS_WIN32 )
MUMPS_INT mumps_io_read_win32__(void * file,void * loc_addr,size_t size,MUMPS_INT local_offset){
  size_t ret_code;
  fseek(*(FILE **)file,(long) local_offset,SEEK_SET);
  ret_code=fread(loc_addr,mumps_elementary_data_size,size,*(FILE **)file);
  if((ret_code!=size)||(ferror(*(FILE**)file))){
    return mumps_io_error(-90,"Problem with I/O operation\n");
  }
  return 0;
}
#endif
MUMPS_INT mumps_io_write__(void *file, void *loc_addr, size_t write_size, MUMPS_INT where,MUMPS_INT type){
  MUMPS_INT ret_code;
#if ! defined( MUMPS_WIN32 )
  if(!mumps_directio_flag){
    ret_code=mumps_io_write_os_buff__(file,loc_addr, write_size,where);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=mumps_io_write_win32__(file,loc_addr, write_size,where);
  if(ret_code<0){
    return ret_code;
  }
#endif
  return 0;
}
#if ! defined( MUMPS_WIN32 )
MUMPS_INT mumps_io_write_os_buff__(void *file, void *loc_addr, size_t write_size, MUMPS_INT where){
  size_t ret_code;
  /* printf("write with buff %d %d %d\n",(MUMPS_INT) write_size, where,*((MUMPS_INT *)file)); */
# ifdef WITH_PFUNC
  ret_code=pwrite(*(MUMPS_INT *)file,loc_addr,write_size,where);
# else
  /*in this case all the I/O's are made by the I/O thread => we don't
    need to protect the file pointer.*/
  lseek(*(MUMPS_INT *)file,(long)where,SEEK_SET); 
  ret_code=write(*(MUMPS_INT *)file,loc_addr,write_size);
# endif
  if((MUMPS_INT)ret_code==-1){
    return mumps_io_sys_error(-90,"Problem with low level write");
  } else if(ret_code!=write_size){
    return mumps_io_error(-90,"Error not enough space on disk \n");
  }
  return 0;
}
#endif
#if defined( MUMPS_WIN32 )
MUMPS_INT mumps_io_write_win32__(void *file, void *loc_addr, size_t write_size, MUMPS_INT where){
  size_t ret_code;
  fseek(*(FILE **)file,(long)where,SEEK_SET);  
  ret_code=fwrite(loc_addr,mumps_elementary_data_size, write_size,*(FILE**)file);
  if((ret_code!=write_size)||(ferror(*(FILE**)file))){
    return mumps_io_error(-90,"Problem with I/O operation\n");
  }
  return 0;
}
#endif
MUMPS_INT mumps_compute_file_size(void *file,size_t *size){
  /* Compute the size of the file pointed by file and return it in
     size */
#if defined(MUMPS_WIN32)
  /* This works well as soon as we don't use threads under WIN32 */
  MUMPS_INT ret_code;
  long pos=0;
  /* Get the current position */
  pos=ftell(*(FILE **)file);
  /* Move the file pointer to the end of the file */
  fseek(*(FILE **)file,0,SEEK_END);
  /* Get the current position which is in fact the size */
  *size=(size_t)ftell(*(FILE **)file);
  /* Restore the old position */
  fseek(*(FILE **)file,pos,SEEK_SET);
#else
  struct stat file_info;
  /* fstat does everything :-) */
  fstat(*(MUMPS_INT *)file, &file_info);
  *size = (size_t)file_info.st_size;
#endif
  return 0;
}
