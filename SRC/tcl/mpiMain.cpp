/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2006-01-13 19:35:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/mpiMain.cpp,v $

/* 
 * tclAppInit.c --
 *
 *	Provides a default version of the main program and Tcl_AppInit
 *	procedure for Tcl applications (without Tk).
 *
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994-1997 Sun Microsystems, Inc.
 * Copyright (c) 1998-1999 by Scriptics Corporation.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: mpiMain.cpp,v 1.1 2006-01-13 19:35:27 fmk Exp $
 */

extern "C" {
#include <tcl.h>
#include <tk.h>
}


// #include <mpi.h>
#include "commands.h"

/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */

#ifdef _KAI
extern "C" int matherr();
#endif

#ifdef _UNIX
#include <math.h>

//int *tclDummyMathPtr = (int *)matherr;
#endif


#ifdef TCL_TEST

extern "C" {
#include "tclInt.h"
}



extern int		Procbodytest_Init _ANSI_ARGS_((Tcl_Interp *interp));
extern int		Procbodytest_SafeInit _ANSI_ARGS_((Tcl_Interp *interp));
extern int		TclObjTest_Init _ANSI_ARGS_((Tcl_Interp *interp));
extern int		Tcltest_Init _ANSI_ARGS_((Tcl_Interp *interp));
#ifdef TCL_THREADS
extern int		TclThread_Init _ANSI_ARGS_((Tcl_Interp *interp));
#endif

#endif /* TCL_TEST */

#ifdef TCL_XT_TEST
extern void		XtToolkitInitialize _ANSI_ARGS_((void));
extern int		Tclxttest_Init _ANSI_ARGS_((Tcl_Interp *interp));
#endif

/*
 *----------------------------------------------------------------------
 *
 * main --
 *
 *	This is the main program for the application.
 *
 * Results:
 *	None: Tcl_Main never returns here, so this procedure never
 *	returns either.
 *
 * Side effects:
 *	Whatever the application does.
 *
 *----------------------------------------------------------------------
 */

extern void g3TclMain(int argc, char **argv, Tcl_AppInitProc *appInitProc);
#include <stdio.h>
#include <string.h>

#include <PartitionedDomain.h>
#include <MPI_MachineBroker.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <FEM_ObjectBroker.h>

extern PartitionedDomain theDomain;

extern int OPS_PARALLEL_PROCESSING;
extern int OPS_NUM_SUBDOMAINS;
extern bool OPS_PARTITIONED;
extern FEM_ObjectBroker *OPS_OBJECT_BROKER;
extern MachineBroker    *OPS_MACHINE;
extern bool OPS_USING_MAIN_DOMAIN;
extern int OPS_MAIN_DOMAIN_PARTITION_ID;

static MPI_MachineBroker *theMachineBroker = 0;

int
main(int argc, char **argv)
{
  FEM_ObjectBroker theBroker;
  MPI_MachineBroker theMachine(&theBroker, argc, argv);
  theMachineBroker = &theMachine;

  int rank = theMachine.getPID();
  int np = theMachine.getNP();

  //
  // depending on rank we do something
  //
  if (rank != 0) {

    //
    // on slave processes we spin waiting to create & run actors
    //
    fprintf(stderr, "Slave Process Running\n");
    theMachine.runActors();
    fprintf(stderr, "Slave Process DONE %d\n", rank);

  } else {

    //
    // on process 0 we create some ShadowSubdomains & then start the OpenSees interpreter
    //
    fprintf(stderr, "Master Process Running OpenSees Interpreter\n");   

    //
    // set some global parameters
    //
    OPS_OBJECT_BROKER = &theBroker;
    OPS_MACHINE = &theMachine;
    OPS_PARALLEL_PROCESSING = np;

    if (np%2 == 0) {
      OPS_NUM_SUBDOMAINS = np;
      OPS_USING_MAIN_DOMAIN = true;
      OPS_MAIN_DOMAIN_PARTITION_ID = 1;
    } else
      OPS_NUM_SUBDOMAINS = np - 1;

    OPS_PARTITIONED    = false;
    
    //
    // the rest straightr out of regular tclMain to start our interpreter
    //

    /*
     * The following #if block allows you to change the AppInit
     * function by using a #define of TCL_LOCAL_APPINIT instead
     * of rewriting this entire file.  The #if checks for that
     * #define and uses Tcl_AppInit if it doesn't exist.
     */

#ifndef TCL_LOCAL_APPINIT
#define TCL_LOCAL_APPINIT Tcl_AppInit    
#endif
    
    /* fmk - comment out the following block to get to compile 
       extern "C" int TCL_LOCAL_APPINIT _ANSI_ARGS_((Tcl_Interp *interp));
       fmk - end commented block */

    /*
     * The following #if block allows you to change how Tcl finds the startup
     * script, prime the library or encoding paths, fiddle with the argv,
     * etc., without needing to rewrite Tcl_Main()
     */
    
#ifdef TCL_LOCAL_MAIN_HOOK
    extern int TCL_LOCAL_MAIN_HOOK _ANSI_ARGS_((int *argc, char ***argv));
#endif
    
#ifdef TCL_XT_TEST
    XtToolkitInitialize();
#endif
    
#ifdef TCL_LOCAL_MAIN_HOOK
    TCL_LOCAL_MAIN_HOOK(&argc, &argv);
#endif

    g3TclMain(argc, argv, TCL_LOCAL_APPINIT);

    // some clean up to shut the remotes down if still running
    theDomain.clearAll();
    
    // shutdown the remote machines
    theMachine.shutdown();
  }
  
  //
  // mpi clean up
  //
  

  fprintf(stderr, "Process Terminating %d\n", rank);
  
  return 0;
}

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in the interp's result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */


int Tcl_AppInit(Tcl_Interp *interp)
{
    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

#ifdef TCL_TEST
#ifdef TCL_XT_TEST
     if (Tclxttest_Init(interp) == TCL_ERROR) {
	 return TCL_ERROR;
     }
#endif
    if (Tcltest_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    Tcl_StaticPackage(interp, "Tcltest", Tcltest_Init,
            (Tcl_PackageInitProc *) NULL);
    if (TclObjTest_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
#ifdef TCL_THREADS
    if (TclThread_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
#endif
    if (Procbodytest_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    Tcl_StaticPackage(interp, "procbodytest", Procbodytest_Init,
            Procbodytest_SafeInit);
#endif /* TCL_TEST */

    /*
     * Call the init procedures for included packages.  Each call should
     * look like this:
     *
     * if (Mod_Init(interp) == TCL_ERROR) {
     *     return TCL_ERROR;
     * }
     *
     * where "Mod" is the name of the module.
     */

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    if (g3AppInit(interp) < 0)
	return TCL_ERROR;

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}


int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theDomain.clearAll();

  //
  // mpi clean up
  //

  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
  }

  Tcl_Exit(0);

  return 0;
}

