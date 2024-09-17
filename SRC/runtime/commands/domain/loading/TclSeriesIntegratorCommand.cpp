//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the groundMotion command in the interpreter.
//
// Written: fmk
// Created: 11/00
//
#include <tcl.h>
#include <string.h>

#include <TrapezoidalTimeSeriesIntegrator.h>
#include <SimpsonTimeSeriesIntegrator.h>

// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in
//   one array, which is why Tcl_Free needs only be called on the array.
static void
cleanup(TCL_Char ** argv)
{
  Tcl_Free((char *)argv);
}

TimeSeriesIntegrator *
TclDispatch_newSeriesIntegrator(ClientData clientData, Tcl_Interp* interp, TCL_Char * const arg)
{
  int argc;
  TCL_Char ** argv;

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    opserr << "WARNING could not split series integrator list " << arg << "\n";
    return 0;
  }

  TimeSeriesIntegrator *theSeriesIntegrator = nullptr;

  if (strcmp(argv[0], "Trapezoidal") == 0) {
    theSeriesIntegrator = new TrapezoidalTimeSeriesIntegrator();
  }

  else if (strcmp(argv[0], "Simpson") == 0) {
    theSeriesIntegrator = new SimpsonTimeSeriesIntegrator();
  }

  else {
    // type of load pattern type unknown
    opserr << "WARNING unknown TimeSeriesIntegrator type " << argv[0] << " - ";
    opserr << " SeriesIntegratorType <type args>\n\tvalid types: Trapezoidal "
              "or Simpson\n";
    cleanup(argv);
    return 0;
  }

  cleanup(argv);
  return theSeriesIntegrator;
}
