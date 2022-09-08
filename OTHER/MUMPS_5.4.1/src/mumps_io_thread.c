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
#include "mumps_io_thread.h"
#include "mumps_c_types.h"
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
/* Exported global variables */
MUMPS_INT io_flag_stop,current_req_num;
pthread_t io_thread,main_thread;
pthread_mutex_t io_mutex;
pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
pthread_mutex_t io_mutex_cond;
MUMPS_INT int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;
MUMPS_INT with_sem;
struct request_io *io_queue;
MUMPS_INT first_active,last_active,nb_active;
MUMPS_INT *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
MUMPS_INT mumps_owns_mutex;
MUMPS_INT test_request_called_from_mumps;
/* Other global variables */
double inactive_time_io_thread;
MUMPS_INT time_flag_io_thread;
struct timeval origin_time_io_thread;
/**
 * Main function of the io thread when semaphores are used.
 */
void*  mumps_async_thread_function_with_sem (void* arg){
   struct request_io *current_io_request;
   MUMPS_INT ierr,_sem_stop;
   struct timeval start_time,end_time;
   MUMPS_INT ret_code;
   for (;;){      
     gettimeofday(&start_time,NULL);
       if(with_sem==2){
         mumps_wait_sem(&int_sem_io,&cond_io);
       }
     /*     sem_wait(&sem_io);  */
     gettimeofday(&end_time,NULL);
     if(time_flag_io_thread){
       inactive_time_io_thread=inactive_time_io_thread+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
     }else{
       inactive_time_io_thread=((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)origin_time_io_thread.tv_sec+((double)origin_time_io_thread.tv_usec/1000000));
     }
     if(!time_flag_io_thread){
       time_flag_io_thread=1;
     }
     /* Check if the main thread ordered to stop this slave thread */
     /*     sem_getvalue(&sem_stop,&_sem_stop); */
       if(with_sem==2){
         mumps_get_sem(&int_sem_stop,&_sem_stop);
       }
     if(_sem_stop==IO_FLAG_STOP){
       /* The thread must stop */
       break; /* Breaks the while loop. */
     }
      current_io_request=&io_queue[first_active];
      switch(current_io_request->io_type)
        {
         case IO_WRITE:
           ret_code=mumps_io_do_write_block(current_io_request->addr,
                                     current_io_request->size,
                                     &(current_io_request->file_type),
                                     current_io_request->vaddr,
                                     &ierr);
           if(ret_code<0){
             goto end;
           }
           break;
         case IO_READ:
           ret_code=mumps_io_do_read_block(current_io_request->addr,
           current_io_request->size,
           &(current_io_request->file_type),
           current_io_request->vaddr,
           &ierr);
           if(ret_code<0){
             goto end;
           }
           break;
         default:
           printf("Error : Mumps_IO : Operation %d is neither READ nor WRITE\n",current_io_request->io_type);
           exit (-3);
         }
      /* Notify that the IO was performed */
      /* Wait that finished_requests queue could register 
       the notification */
       if(with_sem==2){
         mumps_wait_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
       }
      pthread_mutex_lock(&io_mutex);
      /* Updates active queue bounds */
      /* Register the notification in finished_requests queue
       and updates its bounds. */
      finished_requests_id[last_finished_requests]=current_io_request->req_num;
      finished_requests_inode[last_finished_requests]=current_io_request->inode;
      last_finished_requests=(last_finished_requests+1)%(MAX_FINISH_REQ); /* ??? */
      nb_finished_requests++;      
      /* Realeases the lock : ***UNLOCK*** */
      nb_active--;      
      if(first_active<MAX_IO-1){
        first_active++; 
      }
      else{
        first_active=0;
      }
      if(with_sem==2){
        mumps_post_sem(&(current_io_request->int_local_cond),&(current_io_request->local_cond));
      }
      pthread_mutex_unlock(&io_mutex);  
      /* Finally increases the number of free active requests.*/
      /*      sem_post(&sem_nb_free_active_requests); */
          mumps_post_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
   }
 end:
   /* The main thread ordered the end of the IO thread (it changed sem_stop).
    We exit. */
   pthread_exit(NULL);
/* FIXME Not reached */
/*   return NULL; */
}
MUMPS_INT mumps_test_request_th(MUMPS_INT* request_id,MUMPS_INT *flag){
  /* Tests if the request "request_id" has finished. It sets the flag  */
  /* argument to 1 if the request has finished (0 otherwise)           */
  MUMPS_INT request_pos;
  MUMPS_INT i;
  i=mumps_check_error_th();
  if(i!=0){
    return i;
  }
  pthread_mutex_lock(&io_mutex);
  /*  printf("entering test !!! \n"); */
  if(*request_id < smallest_request_id){
    *flag=1;
    /*    exit (-2); */
  }else{
    if(nb_finished_requests==0){
      *flag=0;
    }else{
      request_pos=(first_finished_requests+nb_finished_requests-1)%(MAX_IO*2);    
      if(*request_id > finished_requests_id[request_pos]){
        /*the request has not been treated yet since it is not in 
        the list of treated requests*/
        i=0;
        /*this loop is only for checking (no special treatment is done*/
        while(i<nb_active){
          request_pos=(first_active+i)%(MAX_IO);
          if(io_queue[request_pos].req_num==*request_id){
            break;
          }
          i++;
        }
        if(i==nb_active){
          return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_test_request_th (1))\n");
        }
        *flag=0;
      }else{
        /*the request has been treated yet since it is the list of
          treated requests*/
        i=0;
        while(i<nb_finished_requests){
          request_pos=(first_finished_requests+i)%(MAX_IO*2);
          if(finished_requests_id[request_pos]==*request_id){
            break;
          }
          i++;
        }
        if(i==nb_finished_requests){
          return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_test_request_th (2))\n");
        }      
        *flag=1;
      }
    }
  }
  /*  printf("je vais essayer de nettoyer %d\n",nb_finished_requests);
    pthread_mutex_unlock(&io_mutex);
    sleep (10); */
  mumps_owns_mutex=1;
  /* printf("ici nb_finished_requests=%d\n",nb_finished_requests);*/
  mumps_clean_finished_queue_th();
  mumps_owns_mutex=0;
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
MUMPS_INT mumps_wait_req_sem_th(MUMPS_INT *request_id){
  MUMPS_INT i,j;
  j=first_active;
  for(i=0;i<nb_active;i++){
    if(io_queue[j].req_num==*request_id) break;
    j=(j+1)%MAX_IO;
  }
  if(i<nb_active){
    mumps_wait_sem(&(io_queue[j].int_local_cond),&(io_queue[j].local_cond));
  }
  return 0;
}
MUMPS_INT mumps_wait_request_th(MUMPS_INT *request_id){
  /* waits for the termination of the request "request_id" */
  MUMPS_INT flag=0,ierr;
  if(with_sem!=2){
    while(!flag){
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }else{
    ierr=mumps_test_request_th(request_id,&flag);
    if(ierr!=0)return ierr;
    if(!flag){
      mumps_wait_req_sem_th(request_id);
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }
  return 0;
}
MUMPS_INT mumps_is_there_finished_request_th(MUMPS_INT* flag){
  if(!mumps_owns_mutex) pthread_mutex_lock(&io_mutex);
  if(nb_finished_requests==0){
    *flag=0;
  }else{
    *flag=1;
    /*    printf("finished : %d\n",nb_finished_requests);     */
  }
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  return 0;
}
MUMPS_INT mumps_clean_finished_queue_th(){
   /* Cleans the finished request queue. On exit, the queue is empty.*/
   MUMPS_INT local_flag;
   MUMPS_INT cur_req;
   MUMPS_INT loc_owned_mutex=0,ierr;
   if(!mumps_owns_mutex){
      pthread_mutex_lock(&io_mutex);
      mumps_owns_mutex=1;
      loc_owned_mutex=1;
  }
  /* this block of code is designed for avoiding deadlocks between
     the two threads*/
   mumps_is_there_finished_request_th(&local_flag);
   while(local_flag){
     ierr=mumps_clean_request_th(&cur_req);
     if(ierr!=0){
       return ierr;
     }
     mumps_is_there_finished_request_th(&local_flag);
   }
   if((!mumps_owns_mutex)||(loc_owned_mutex)){
      pthread_mutex_unlock(&io_mutex);
      mumps_owns_mutex=0;
   }
   return 0;
}
MUMPS_INT mumps_clean_request_th(MUMPS_INT* request_id){
  MUMPS_INT ierr;
  ierr=mumps_check_error_th();
  if(ierr!=0){
    return ierr;
  }
  if(!mumps_owns_mutex)pthread_mutex_lock(&io_mutex);
  *request_id=finished_requests_id[first_finished_requests];
  if(smallest_request_id!=finished_requests_id[first_finished_requests]){
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_clean_request_th)\n");
  }
  finished_requests_id[first_finished_requests]=-9999;
  first_finished_requests=(first_finished_requests+1)%(MAX_FINISH_REQ);
  nb_finished_requests--;
  /*we treat the io requests in their arrival order => we just have to
    increase smallest_request_id*/
  smallest_request_id++;
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  if(with_sem) {
      if(with_sem==2){
        mumps_post_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
      }
  }
  return 0;
}
MUMPS_INT mumps_low_level_init_ooc_c_th(MUMPS_INT* async, MUMPS_INT* ierr){
  MUMPS_INT i, ret_code;    
  char buf[128];
  /* Computes the number of files needed. Uses ceil value. */
  *ierr=0;
  current_req_num=0;
  with_sem=2;
  first_active=0;
  last_active=0;
  nb_active=0;
  first_finished_requests=0;
  last_finished_requests=0;
  nb_finished_requests=0;
  smallest_request_id=0;
  mumps_owns_mutex=0;
  inactive_time_io_thread=0;
  time_flag_io_thread=0;
  gettimeofday(&origin_time_io_thread,NULL);
  /*  mumps_io_flag_async=*async; */
  if(*async!=IO_ASYNC_TH){
    *ierr = -91;
    sprintf(buf,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
    return mumps_io_error(*ierr,buf);
  }
  if(*async){
    pthread_mutex_init(&io_mutex,NULL);
    mumps_io_init_err_lock();
#ifdef WITH_PFUNC
    mumps_io_init_pointers_lock();
#endif
    io_queue=(struct request_io *)malloc(MAX_IO*sizeof(struct request_io));
    if(with_sem==2){
      for(i=0;i<MAX_IO;i++){
        pthread_cond_init(&(io_queue[i].local_cond),NULL);
        io_queue[i].int_local_cond=0;
      }
    }
    finished_requests_id=(MUMPS_INT *)malloc(MAX_IO*2*sizeof(MUMPS_INT));
    finished_requests_inode=(MUMPS_INT *)malloc(MAX_IO*2*sizeof(MUMPS_INT));
    for(i=0;i<MAX_IO*2;i++){
      finished_requests_id[i]=-9999;
      finished_requests_inode[i]=-9999;
    }
    if(with_sem){
      switch(with_sem){
      case 2:
        int_sem_io=0;
        int_sem_stop=0;
        int_sem_nb_free_finished_requests=MAX_FINISH_REQ;
        int_sem_nb_free_active_requests=MAX_IO;
        pthread_cond_init(&cond_stop,NULL);
        pthread_cond_init(&cond_io,NULL);
        pthread_cond_init(&cond_nb_free_active_requests,NULL);
        pthread_cond_init(&cond_nb_free_finished_requests,NULL);
        pthread_mutex_init(&io_mutex_cond,NULL);
        break;
      default:
        *ierr = -92;
        sprintf(buf,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
        return mumps_io_error(*ierr,buf);
      }
      ret_code=pthread_create(&io_thread,NULL,mumps_async_thread_function_with_sem,NULL);
    }
    if(ret_code!=0){
      errno = ret_code;
      return mumps_io_sys_error(-92,"Unable to create I/O thread");
    }
    main_thread=pthread_self();
  }
  return 0;
}
MUMPS_INT mumps_async_write_th(const MUMPS_INT * strat_IO, 
                        void * address_block,
                        long long block_size,
                        MUMPS_INT * inode,
                        MUMPS_INT * request_arg,
                        MUMPS_INT * type,
                        long long vaddr,
                        MUMPS_INT * ierr){
  MUMPS_INT cur_req;
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){ 
    mumps_clean_finished_queue_th();
      if(with_sem==2){
        mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }
  if(nb_active<=MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }
    else{
      last_active=(last_active+1)%MAX_IO;
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode; 
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=0;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    *ierr = -91;
    return mumps_io_error(*ierr,"Internal error in OOC Management layer (mumps_async_write_th)\n");
    /*    exit(-3);*/
  }
  pthread_mutex_unlock(&io_mutex);
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
        mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  return 0;
}
MUMPS_INT mumps_async_read_th(const MUMPS_INT * strat_IO, 
                       void * address_block,
                       long long  block_size,
                       MUMPS_INT * inode,
                       MUMPS_INT * request_arg,
                        MUMPS_INT * type,
                       long long vaddr,
                       MUMPS_INT * ierr){
  MUMPS_INT cur_req;  
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){
    mumps_clean_finished_queue_th();
    /* end of the block*/
      if(with_sem==2){
        mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }
  if(nb_active<MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }else{
      last_active=(last_active+1)%MAX_IO;
      /*        if(last_active<MAX_IO-1){
                cur_req=last_active+1;
                last_active++;
                }else{
                cur_req=0;
                last_active=0;
                }*/
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode;
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=1;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    *ierr = -91;
    return mumps_io_error(*ierr,"Internal error in OOC Management layer (mumps_async_read_th)\n");
  }
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
        mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
MUMPS_INT mumps_clean_io_data_c_th(MUMPS_INT *myid){
  MUMPS_INT i;
  /* cleans the thread/io management data*/
  if(mumps_io_flag_async){
    /*we can use either signals or mutexes for this step */
    if(with_sem){
        if(with_sem==2){
          mumps_post_sem(&int_sem_stop,&cond_stop);
          mumps_post_sem(&int_sem_io,&cond_io);
        }
    }else{
      pthread_mutex_lock(&io_mutex);
      io_flag_stop=1;
      pthread_mutex_unlock(&io_mutex);
    }
    pthread_join(io_thread,NULL);
    pthread_mutex_destroy(&io_mutex);
    mumps_io_destroy_err_lock();
#ifdef WITH_PFUNC
    mumps_io_destroy_pointers_lock();
#endif
    if(with_sem){
        if(with_sem==2){
          pthread_cond_destroy(&cond_stop);
          pthread_cond_destroy(&cond_io);
          pthread_cond_destroy(&cond_nb_free_active_requests);
          pthread_cond_destroy(&cond_nb_free_finished_requests);
          pthread_mutex_destroy(&io_mutex_cond);
        }
    }
  }
  if(with_sem==2){
    for(i=0;i<MAX_IO;i++){
      pthread_cond_destroy(&(io_queue[i].local_cond));
    }
  }
  free(io_queue);
  free(finished_requests_id);
  free(finished_requests_inode);
  return 0;
}
MUMPS_INT mumps_get_sem(void *arg,MUMPS_INT *value){
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    *value=*((MUMPS_INT *)arg);
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps__get_sem)\n");
  }
  return 0;
}
MUMPS_INT mumps_wait_sem(void *arg,pthread_cond_t *cond){
  MUMPS_INT *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(MUMPS_INT *)arg;
    while(*tmp_pointer==0){
      pthread_cond_wait(cond, &io_mutex_cond);
    }
    (*tmp_pointer)--;
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_wait_sem)\n");
  }
  return 0;  
}
MUMPS_INT mumps_post_sem(void *arg,pthread_cond_t *cond){
  MUMPS_INT *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(MUMPS_INT *)arg;
    (*tmp_pointer)++;
    if(*tmp_pointer==1){
      pthread_cond_broadcast(cond);
    }
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    return mumps_io_error(-91,"Internal error in OOC Management layer (mumps_post_sem)\n");
  }
  return 0;  
}
#endif /* MUMPS_WIN32 && WITHOUT_PTHREAD */
