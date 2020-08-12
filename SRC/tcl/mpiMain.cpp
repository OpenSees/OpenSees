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
                                                                        
// $Revision: 1.13 $
// $Date: 2010-04-23 23:01:14 $
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
 * See Tcl/Tk License Terms for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: mpiMain.cpp,v 1.13 2010-04-23 23:01:14 fmk Exp $

Tcl/Tk License Terms:
This software is copyrighted by the Regents of the University of
California, Sun Microsystems, Inc., Scriptics Corporation, ActiveState
Corporation and other parties.  The following terms apply to all files
associated with the software unless explicitly disclaimed in
individual files.

The authors hereby grant permission to use, copy, modify, distribute,
and license this software and its documentation for any purpose, provided
that existing copyright notices are retained in all copies and that this
notice is included verbatim in any distributions. No written agreement,
license, or royalty fee is required for any of the authorized uses.
Modifications to this software may be copyrighted by their authors
and need not follow the licensing terms described here, provided that
the new terms are clearly indicated on the first page of each file where
they apply.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS.

GOVERNMENT USE: If you are acquiring this software on behalf of the
U.S. government, the Government shall have only "Restricted Rights"
in the software and related documentation as defined in the Federal 
Acquisition Regulations (FARs) in Clause 52.227.19 (c) (2).  If you
are acquiring the software on behalf of the Department of Defense, the
software shall be classified as "Commercial Computer Software" and the
Government shall have only "Restricted Rights" as defined in Clause
252.227-7013 (c) (1) of DFARs.  Notwithstanding the foregoing, the
authors grant the U.S. Government and others acting in its behalf
permission to use and distribute the software in accordance with the
terms specified in this license. 

 */

extern "C" {
#include <tcl.h>
}


// #include <mpi.h>
#include "commands.h"

/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */

#ifdef _KAI
//extern "C" int matherr();
#endif

#ifdef _UNIX
//#include <math.h>

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

extern void g3TclMain(int argc, char **argv, Tcl_AppInitProc *appInitProc, int rank, int np);
#include <stdio.h>
#include <string.h>

#include <PartitionedDomain.h>
#include <MPI_MachineBroker.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <FEM_ObjectBrokerAllClasses.h>

#include <DomainPartitioner.h>

extern PartitionedDomain theDomain;
extern int OPS_PARALLEL_PROCESSING;
extern int OPS_NUM_SUBDOMAINS;
extern bool OPS_PARTITIONED;
extern FEM_ObjectBroker *OPS_OBJECT_BROKER;
extern MachineBroker    *OPS_MACHINE;
extern bool OPS_USING_MAIN_DOMAIN;
extern int OPS_MAIN_DOMAIN_PARTITION_ID;

extern DomainPartitioner *OPS_DOMAIN_PARTITIONER;
extern GraphPartitioner  *OPS_GRAPH_PARTITIONER;
extern LoadBalancer      *OPS_BALANCER;
extern FEM_ObjectBroker  *OPS_OBJECT_BROKER;
extern MachineBroker     *OPS_MACHINE;
extern Channel          **OPS_theChannels;

extern MachineBroker *theMachineBroker;
extern Channel **theChannels;
extern int numChannels;
extern int OPS_rank;
extern int OPS_np;

#include <mpi.h>

int
main(int argc, char **argv)
{
  theMachineBroker = new MPI_MachineBroker(0, argc, argv);
  FEM_ObjectBrokerAllClasses theBroker;
  theMachineBroker->setObjectBroker(&theBroker);

  OPS_rank = theMachineBroker->getPID();
  OPS_np = theMachineBroker->getNP();

  //
  // depending on rank we do something
  //
  if (OPS_rank != 0) {

    //
    // on secondary processes we spin waiting to create & run actors
    //
    fprintf(stderr, "Secondary Process Running %d\n", OPS_rank);
    theMachineBroker->runActors();

  } else {

    //
    // on process 0 we create some ShadowSubdomains & then start the OpenSees interpreter
    //
    fprintf(stderr, "Primary Process Running OpenSees Interpreter %d\n", OPS_rank);   

    //
    // set some global parameters
    //
    OPS_OBJECT_BROKER = &theBroker;
    //    OPS_MACHINE = &theMachine;
    OPS_MACHINE = theMachineBroker;
    OPS_PARALLEL_PROCESSING = OPS_np;
    
    /* only use p0 if even number of partitions
       if (OPS_np%2 == 0) {
       OPS_NUM_SUBDOMAINS = OPS_np;
       OPS_USING_MAIN_DOMAIN = true;
       OPS_MAIN_DOMAIN_PARTITION_ID = 1;
       } else
       OPS_NUM_SUBDOMAINS = OPS_np - 1;
    */
    // always use p0 even if ODD number of partitions
    OPS_NUM_SUBDOMAINS = OPS_np;
    OPS_USING_MAIN_DOMAIN = true;
    OPS_MAIN_DOMAIN_PARTITION_ID = 1;
    
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

    g3TclMain(argc, argv, TCL_LOCAL_APPINIT, 1, 0);

    // some clean up to shut the remotes down if still running
    theDomain.clearAll();
    
    // shutdown the remote machines
    theMachineBroker->shutdown();
  }
  
  //
  // mpi clean up
  //

  fprintf(stderr, "Process Terminating %d\n", OPS_rank);

  if (theMachineBroker != 0)
    delete theMachineBroker;
  
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

    if (OpenSeesAppInit(interp) < 0)
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


