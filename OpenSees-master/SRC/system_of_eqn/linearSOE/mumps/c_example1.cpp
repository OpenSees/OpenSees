/* This file is part of MUMPS VERSION 4.7.3
This Version was built on Fri May 4 15:54:01 2007

*/
/* $Id: c_example1.cpp,v 1.1 2007-06-01 17:44:35 fmk Exp $ */
/* Example program using the C interface to the 
 * double precision version of MUMPS, dmumps_c.
 * We solve the system A x = RHS with
 *   A = diag(1 2) and RHS = [1 4]^T
 * Solution is [1 2]^T */
#include <stdio.h>
#include <mpi.h>
extern "C" {
#include <dmumps_c.h>
}

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654


#include <OPS_Globals.h>
#include <StandardStream.h>
#include <mpi.h>
#include <SimulationInformation.h>
 
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
SimulationInformation simulationInfo;
  
double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;



#if defined(MAIN_COMP)
/*
 * Some Fortran compilers (COMPAQ fort) define main inside
 * their runtime library while a Fortran program translates
 * to MAIN_ or MAIN__ which is then called from "main". This
 * is annoying because MAIN__ has no arguments and we must
 * define argc/argv arbitrarily !!
 */
int MAIN__();
int MAIN_()
  {
    return MAIN__();
  }

int MAIN__()
{
  int argc=1;
  char * name = "c_example";
  char ** argv ;
#else
int main(int argc, char ** argv)
{
#endif
  DMUMPS_STRUC_C id;

  int n = 6;
  int nz = 21;
  int irnL[] = {1,2,3,4,5,6,2,3,4,5,6,3,4,5,6,4,5,6,5,6,6};
  int jcnL[] = {1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5,6};
  int irnS[] = {1,2,3,2,3,3};
  int jcnS[] = {1,1,1,2,2,3};

  int *irn = irnL;
  int *jcn = jcnL;

  double a[21];
  double rhs[6];

  int myid, ierr, numP;
#if defined(MAIN_COMP)
  argv = &name;
#endif
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numP);
  
  double nump = numP * 1.0;

  /* Define A and rhs */
  rhs[0]=2.0;rhs[1]=2.0;rhs[2]=2.0;rhs[3]=2.0;rhs[4]=2.0;rhs[5]=2.0;

  a[0]=4169.95/nump;
  a[1]=0.0;
  a[2]=10075.0/nump;
  a[3]=-4030;
  a[4]=0.0;
  a[5]=0.0;
  a[6]=10084.0/nump;
  a[7]=-1612.0/nump;
  a[8]=0.0;
  a[9]=-8.9556;
  a[10]=-1612;
  a[11]=1354080.0/nump;
  a[12]=0.0;
  a[13]=1612.0;
  a[14]=193440.0;
  a[15]=4169.93;
  a[16]=0.0;
  a[17]=10075;
  a[18]=10084;
  a[19]=1612;
  a[20]=1354000;

  if (myid != 0) {
      nz = 6;

      a[0]=4169.95/nump;
      a[1]=0.0;
      a[2]=10075.0/nump;

      a[3]=10084.0/nump;
      a[4]=-1612.0/nump;

      a[5]=1354080.0/nump;

      irn = irnS;
      jcn = jcnS;
  }

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  id.job=JOB_INIT; id.par=1; id.sym=2; id.comm_fortran=USE_COMM_WORLD;

  /* parallel solver; distributed i/p matrix A */
  id.ICNTL(5)=0; id.ICNTL(18)=3;  
  dmumps_c(&id);

  /* Define the problem on the host */
  /*
  id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
  id.a = a; id.rhs = rhs;
  */

  /* parallel solver; distributed i/p matrix A */

  id.ICNTL(5)=0; id.ICNTL(18)=3;  


  id.n = n; id.nz_loc =nz; id.irn_loc=irn; id.jcn_loc=jcn;
  id.a_loc = a; id.rhs = rhs;

/* No outputs */
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0; 

        id.job=1;
       dmumps_c(&id);

  id.job=5;
  dmumps_c(&id);

  if (myid == 0) {
    printf("Solution is : (%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e \n", rhs[0],rhs[1], rhs[2], rhs[3], rhs[4], rhs[5]);
  }

  rhs[0]=2.0;rhs[1]=2.0;rhs[2]=2.0;rhs[3]=1.0;rhs[4]=1.0;rhs[5]=1.0;

  id.job=3;
  dmumps_c(&id);

  /* Terminate instance */
  id.job=JOB_END; 
  dmumps_c(&id); 

  if (myid == 0) {
    printf("Solution is : (%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e \n", rhs[0],rhs[1], rhs[2], rhs[3], rhs[4], rhs[5]);
  }

  ierr = MPI_Finalize();
  return 0;
}














