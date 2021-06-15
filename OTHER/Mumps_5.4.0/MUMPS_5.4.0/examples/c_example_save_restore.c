/*
 *
 *  This file is part of MUMPS 5.4.0, released
 *  on Tue Apr 13 15:26:30 UTC 2021
 *
 */
/* Example program using the C interface to the 
 * double real arithmetic version of MUMPS, dmumps_c.
 * We solve the system A x = RHS with
 *   A = diag(1 2) and RHS = [1 4]^T
 * Solution is [1 2]^T */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#if defined(MAIN_COMP)
/*
 * Some Fortran compilers (COMPAQ fort) define "main" in
 * their runtime library while a Fortran program translates
 * to MAIN_ or MAIN__ which is then called from "main".
 * We defined argc/argv arbitrarily in that case.
 */
int MAIN__();
int MAIN_()
  {
    return MAIN__();
  }

int MAIN__()
{
  int argc=1;
  char * name = "c_example_save_restore";
  char ** argv ;
#else
int main(int argc, char ** argv)
{
#endif
  DMUMPS_STRUC_C id_save,id_restore;
  MUMPS_INT n = 2;
  MUMPS_INT8 nnz = 2;
  MUMPS_INT irn[] = {1,2};
  MUMPS_INT jcn[] = {1,2};
  double a[2];
  double rhs[2];

  int error = 0;
/* When compiling with -DINTSIZE64, MUMPS_INT is 64-bit but MPI
   ilp64 versions may still require standard int for C interface. */
/* MUMPS_INT myid, ierr; */
  int myid, ierr;
#if defined(MAIN_COMP)
  argv = &name;
#endif
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  /* Define A and rhs */
  rhs[0]=1.0;rhs[1]=4.0;
  a[0]=1.0;a[1]=2.0;

  /* Initialize MUMPS save instance. Use MPI_COMM_WORLD */
  id_save.comm_fortran=USE_COMM_WORLD;
  id_save.par=1; id_save.sym=0;
  id_save.job=JOB_INIT;
  dmumps_c(&id_save);
  /* Define the problem on the host */
  if (myid == 0) {
    id_save.n = n; id_save.nnz =nnz; id_save.irn=irn; id_save.jcn=jcn;
    id_save.a = a;
  }
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
  id_save.ICNTL(1)=-1; id_save.ICNTL(2)=-1;
  id_save.ICNTL(3)=-1; id_save.ICNTL(4)=0;
  /* Call the MUMPS package on the save instance (analyse and factorization). */
  id_save.job=4;
  dmumps_c(&id_save);

  /* MUMPS save feature on the save instance. */
  strcpy(id_save.save_prefix,"csave_restore");
  strcpy(id_save.save_dir,"/tmp");
  if (myid == 0) {
    printf("Saving MUMPS instance in %s with prefix %s.\n",
        id_save.save_dir, id_save.save_prefix);
  }
  id_save.job=7;
  dmumps_c(&id_save);
  if (id_save.infog[0]<0) {
    printf("\n (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
        myid, id_save.infog[0], id_save.infog[1]);
    error = 1;
  } else if (myid == 0) {
    printf("                                                    DONE\n\n");
  }

  /* Terminate the save instance. */
  id_save.job=JOB_END;
  dmumps_c(&id_save);



  if (!error) {
    /* Initialize MUMPS restore instance. Use MPI_COMM_WORLD */
    id_restore.comm_fortran=USE_COMM_WORLD;
    id_restore.par=1; id_restore.sym=0;
    id_restore.job=JOB_INIT;
    dmumps_c(&id_restore);
    /* Define the rhs on the host */
    if (myid == 0) {
      id_restore.rhs = rhs;
    }

    /* No outputs */
    id_save.ICNTL(1)=-1; id_save.ICNTL(2)=-1;
    id_save.ICNTL(3)=-1; id_save.ICNTL(4)=0;

    /* MUMPS restore feature on restore instance. */
    if (myid == 0) {
      printf("Restoring MUMPS instance in %s with prefix %s.\n",
          id_save.save_dir, id_save.save_prefix);
    }
    strcpy(id_restore.save_prefix,"csave_restore");
    strcpy(id_restore.save_dir,"/tmp");
    id_restore.job=8;
    dmumps_c(&id_restore);
    if (id_save.infog[0]<0) {
      printf("\n (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
          myid, id_save.infog[0], id_save.infog[1]);
      error = 1;
    } else if (myid == 0) {
      printf("                                                    DONE\n\n");
    }
  }

  if (!error) {
    /* Call the MUMPS package on restore instance (solve). */
    if (myid == 0) {
      printf("Calling MUMPS package (solve).\n");
    }
    id_restore.job=3;
    dmumps_c(&id_restore);
    if (id_save.infog[0]<0) {
      printf("=> (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
          myid, id_save.infog[0], id_save.infog[1]);
      error = 1;
    } else if (myid == 0) {
      printf("                                                    DONE\n\n");
    }

    /* Deletes the saved and the OOC files. */
    if (myid == 0) {
      printf("Removing save files.\n");
    }
    id_restore.job=-3;
    dmumps_c(&id_restore);
    if (id_save.infog[0]<0) {
      printf("=> (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
          myid, id_save.infog[0], id_save.infog[1]);
      error = 1;
    } else if (myid == 0) {
      printf("                                                    DONE\n\n");
    }

    /* Terminate the restore instance. */
    id_restore.job=JOB_END;
    dmumps_c(&id_restore);
  }

  if (myid == 0) {
    if (!error) {
      printf("Solution is : (%8.2f  %8.2f)\n", rhs[0],rhs[1]);
    } else {
      printf("An error has occured, please check error code returned by MUMPS.\n");
    }
  }
  ierr = MPI_Finalize();
  return 0;
}
