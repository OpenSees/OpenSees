/* 
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */


#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#elif _WIN32

double SuperLU_timer_() {
    return 0.0;
}

#else

double SuperLU_timer_()
{
    double dclock();
    return (dclock());
}

#endif


#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double usertimer_()
{
    struct tms use;
    double tmp;
    times(&use);
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / CLK_TCK;
}

#if 0
/*
 * It uses the system call gethrtime(3C), which is accurate to 
 * nanoseconds.  This routine is MT-safe.
 *
 */
#include <sys/time.h>
 
double timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#endif


#if 0

/* Get user, syste & wall-clock time */

#include <sys/time.h>
#include <sys/resource.h>
#include </usr/ucbinclude/sys/rusage.h>

double last_user_time = -1;
double last_system_time = -1;
double last_wall_time = -1;

void inittime() 
{
  struct timeval tp;
  struct timezone tzp;
  struct rusage use ;

  gettimeofday(&tp,&tzp);
  getrusage(RUSAGE_SELF, &use) ;
  last_user_time = ((double) use.ru_utime.tv_sec
		    + (double) use.ru_utime.tv_usec/1000000.) ;
  last_system_time = ((double) use.ru_stime.tv_sec
		      + (double) use.ru_stime.tv_usec/1000000. );  
  last_wall_time = tp.tv_sec + (double) tp.tv_usec / 1000000. ;
}

double gettimes(double *utime, double *stime, double *wtime)
{
  struct rusage use ;
  double user_time;
  double system_time;
  double wall_time;
  struct timeval tp;
  struct timezone tzp;

  getrusage(RUSAGE_SELF, &use) ;
  gettimeofday(&tp,&tzp);
  user_time = ((double) use.ru_utime.tv_sec  + 
	       (double) use.ru_utime.tv_usec / 1000000.)
    - last_user_time;
  system_time = ((double) use.ru_stime.tv_sec + 
		 (double) use.ru_stime.tv_usec / 1000000.)
    - last_system_time;
  wall_time = (double) tp.tv_sec + (double) tp.tv_usec / 1000000.
    - last_wall_time;
  last_user_time = last_system_time = last_wall_time = -1;


  *utime = user_time;
  *stime = system_time;
  *wtime = wall_time;
  return ;
}

#endif

