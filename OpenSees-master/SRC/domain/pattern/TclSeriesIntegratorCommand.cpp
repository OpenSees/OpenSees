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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: fmk 
// Created: 11/00
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the groundMotion command in the interpreter.

#include <tcl.h>
#include <string.h>

#include <TrapezoidalTimeSeriesIntegrator.h>
#include <SimpsonTimeSeriesIntegrator.h>


// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in 
//   one array, which is why Tcl_Free needs only be called on the array.
static void cleanup(TCL_Char **argv)  {
    Tcl_Free((char *) argv);
}


TimeSeriesIntegrator *
TclSeriesIntegratorCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char *arg)
{
    int argc;
    TCL_Char **argv;
    
    // split the list
    if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
        opserr << "WARNING could not split series integrator list " << arg << endln;
        return 0;
    }
    
    TimeSeriesIntegrator *theSeriesIntegrator = 0;
    
    if (strcmp(argv[0],"Trapezoidal") == 0) {
        theSeriesIntegrator = new TrapezoidalTimeSeriesIntegrator();
    }
    
    else if (strcmp(argv[0],"Simpson") == 0) {
        theSeriesIntegrator = new SimpsonTimeSeriesIntegrator();
    }
    
    else {
        // type of load pattern type unknown
        opserr << "WARNING unknown TimeSeriesIntegrator type " << argv[0] << " - ";
        opserr << " SeriesIntegratorType <type args>\n\tvalid types: Trapezoidal or Simpson\n";
        cleanup(argv);
        return 0;
    }
    
    cleanup(argv);
    return theSeriesIntegrator;
}
