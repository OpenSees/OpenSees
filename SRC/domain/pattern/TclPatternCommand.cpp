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

// $Revision: 1.12 $
// $Date: 2006-02-08 19:29:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclPatternCommand.cpp,v $

// File: ~/domain/pattern/TclPatternComand.C
//
// Written: fmk
// Created: 07/99
//
// Modified: fmk 11/00 - removed TimeSeries stuff from file, now an external procedure
// Revision: B
//
//  Modified:  Nov.    2002,  Zhaohui  Yang and Boris Jeremic added Plastic Bowl
//  loading (aka Domain Reduction Method) commands
//
// Description: This file contains the function invoked when the user invokes
// the Pattern command in the interpreter. It is invoked by the
// TclModelBuilder_addPattern function in the TclModelBuilder.C file. Current
// valid Pattern types are:

// What: "@(#) TclPatternCommand.C, revA"

#include <TclModelBuilder.h>

#include <tcl.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <RectangularSeries.h>
#include <TrigSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <UniformExcitation.h>
#include <MultiSupportPattern.h>
#include <GroundMotion.h>
#include <GroundMotionRecord.h>
#include <TimeSeriesIntegrator.h>
#include <PBowlLoading.h>

#include <string.h>



LoadPattern *theTclLoadPattern =0;
MultiSupportPattern *theTclMultiSupportPattern =0;


extern TimeSeriesIntegrator *
TclSeriesIntegratorCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char *arg);

extern TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char *arg);


int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv, Domain *theDomain)
{
  LoadPattern *thePattern = 0;

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
      opserr << "WARNING invalid command - want: pattern type ";
      opserr << " <type args> {list of load and sp constraints commands}\n";
      opserr << "           valid types: Plain, UniformExcitation, MultiSupport\n";
      return TCL_ERROR;
  }

  TimeSeries *theSeries = 0;
  int patternID =0;

  if (Tcl_GetInt(interp, argv[2], &patternID) != TCL_OK) {
    opserr << "WARNING invalid patternID: pattern type " << argv[2]
	<< "<type args>\n";
    return TCL_ERROR;
  }

  int  commandEndMarker = 0;

  if (strcmp(argv[1],"Plain") == 0) {

      thePattern = new LoadPattern(patternID);
      theSeries = TclSeriesCommand(clientData, interp, argv[3]);


      if (thePattern == 0 || theSeries == 0) {

    if (thePattern == 0) {
        opserr << "WARNING - out of memory creating LoadPattern ";
        opserr << patternID << endln;
    } else {
        opserr << "WARNING - problem creating TimeSeries for LoadPattern ";
        opserr << patternID << endln;
    }

    // clean up the memory and return an error
    if (thePattern != 0)
        delete thePattern;
    if (theSeries != 0)
        delete theSeries;
    return TCL_ERROR;
      }

      thePattern->setTimeSeries(theSeries);

  }

  else if (strcmp(argv[1],"UniformExcitation") == 0) {

  int dir;

  if (Tcl_GetInt(interp, argv[3], &dir) != TCL_OK) {
    opserr << "WARNING invalid patternID: pattern type " << argv[2]
	<< "<type args>\n";
    return TCL_ERROR;
  }

  dir--; // subtract 1 for c indexing

  TimeSeries *accelSeries = 0;
  TimeSeries *velSeries = 0;
  TimeSeries *dispSeries = 0;
  TimeSeriesIntegrator *seriesIntegrator = 0;
  double vel0 = 0.0;
  
  int currentArg = 4;
  bool doneSeries = false;
  while (currentArg < argc-1 && doneSeries == false) {

  if ((strcmp(argv[currentArg],"-vel0") == 0) ||
      (strcmp(argv[currentArg],"-initialVel") == 0)) {

    currentArg++;
    if ((currentArg < argc) &&
        (Tcl_GetDouble(interp, argv[currentArg], &vel0) != TCL_OK)) {
      opserr << "WARNING invalid vel0: pattern type UniformExciation\n";
      return TCL_ERROR;
    }

    currentArg++;
  }


  else if ((strcmp(argv[currentArg],"-accel") == 0) ||
     (strcmp(argv[currentArg],"-acceleration") == 0)) {

    currentArg++;
    accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
    if (accelSeries == 0) {
      opserr << "WARNING invalid accel series: pattern UniformExcitation -accel <args>\n";
      return TCL_ERROR;
    }
    
    currentArg++;

  } else if ((strcmp(argv[currentArg],"-vel") == 0) ||
       (strcmp(argv[currentArg],"-velocity") == 0)) {

    currentArg++;
    velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

    if (velSeries == 0) {
      opserr << "WARNING invalid vel series: " << argv[currentArg];
      opserr << " pattern UniformExcitation -vel {series}\n";
      return TCL_ERROR;
    }

    currentArg++;

  } else if ((strcmp(argv[currentArg],"-disp") == 0) ||
       (strcmp(argv[currentArg],"-displacement") == 0)) {

    currentArg++;
    dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

    if (dispSeries == 0) {
      opserr << "WARNING invalid vel series: " << argv[currentArg];
      opserr << " pattern UniformExcitation -vel {series}\n";
      return TCL_ERROR;
    }

    currentArg++;
  } else if ((strcmp(argv[currentArg],"-int") == 0) ||
       (strcmp(argv[currentArg],"-integrator") == 0)) {

    currentArg++;
    seriesIntegrator = TclSeriesIntegratorCommand(clientData, interp,
              argv[currentArg]);
    if (seriesIntegrator == 0) {
      opserr << "WARNING invalid series integrator: " << argv[currentArg];
      opserr << " - pattern UniformExcitation -int {Series Integrator}\n";
      return TCL_ERROR;
    }
    currentArg++;
  }

  else
    doneSeries = true;
      }

      if (dispSeries == 0 && velSeries == 0 && accelSeries == 0) {
    opserr << "WARNING invalid series, want - pattern UniformExcitation";
    opserr << "-disp {dispSeries} -vel {velSeries} -accel {accelSeries} ";
    opserr << "-int {Series Integrator}\n";
    return TCL_ERROR;
      }
      
      GroundMotion *theMotion = new GroundMotion(dispSeries, velSeries,
             accelSeries, seriesIntegrator);

      if (theMotion == 0) {
    opserr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
    opserr << patternID << endln;

    return TCL_ERROR;
      }

  // create the UniformExcitation Pattern
  thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0);

  if (thePattern == 0) {
      opserr << "WARNING ran out of memory creating load pattern - pattern UniformExcitation ";
      opserr << patternID << endln;

      // clean up memory allocated up to this point and return an error
      if (theMotion != 0)
    delete theMotion;

      return TCL_ERROR;
  }

  // Added by MHS to prevent call to Tcl_Eval at end of this function
  commandEndMarker = currentArg;
  }

  else if (strcmp(argv[1],"Uniform") == 0) {

  // First search for file name and time step
  int numInputs = argc;
  TCL_Char *accelFileName = 0;
  double dt = 0.0;
  for (int i = 5; i < argc; i++)
  {
    if (strcmp (argv[i],"-accel") == 0 && i+2 < argc)
      {
      // Read the input file name
            accelFileName = argv[i+1];

      // Read the time interval
            if (Tcl_GetDouble(interp, argv[i+2], &dt) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "time interval - pattern UniformExcitation: "
               << patternID << endln;
          return TCL_ERROR;
      }
      numInputs -= 3;

    }
    }

    if (numInputs < 5) {
      opserr << "WARNING insufficient number of arguments - want: pattern ";
      opserr << "UniformExcitation " << patternID << " dir factor\n";
      return TCL_ERROR;
    }

  int dir;
  if (Tcl_GetInt(interp, argv[3], &dir) != TCL_OK) {

    char inputDir;
    inputDir = argv[3][0];
    switch (inputDir)
      {
      case 'X': case 'x': case '1': // Global X
	dir = 0;
	break;
      case 'Y': case 'y': case '2': // Global Y
	dir = 1;
	break;
      case 'Z': case 'z': case '3': // Global Z
	dir = 2;
      break;
      default:
	opserr << "WARNING cannot read direction for excitation \n";
	opserr << "UniformExcitation " << patternID << " dir factor" << endln;
	return TCL_ERROR;
	break;
      }
  } else
    dir--;    // change to c++ indexing

  double factor;
  if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) {
    opserr << "WARNING insufficient number of arguments - want: pattern ";
    opserr << "UniformExcitation " << patternID << " dir factor\n";
    return TCL_ERROR;
  }

    GroundMotionRecord *theMotion;


    // read in the ground motion
  if (accelFileName == 0)
  {
    opserr << "WARNING -- No ground motion data provided\n";
    opserr << "UniformExcitation tag: " << patternID << endln;
    return TCL_ERROR;
  }

    theMotion = new GroundMotionRecord(accelFileName, dt, factor);

  if (theMotion == 0) {
    opserr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
    opserr << patternID << endln;

    return TCL_ERROR;
    }

    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID);

    if (thePattern == 0) {
    opserr << "WARNING ran out of memory creating load pattern - pattern UniformExcitation ";
    opserr << patternID << endln;

    // clean up memory allocated up to this point and return an error
    if (theMotion != 0)
      delete theMotion;

    return TCL_ERROR;
    }
  }

  else if ((strcmp(argv[1],"MultipleSupportExcitation") == 0) ||
     (strcmp(argv[1],"MultipleSupport") == 0) ||
     (strcmp(argv[1],"MultiSupport") == 0)) {

      theTclMultiSupportPattern = new MultiSupportPattern(patternID);
      thePattern = theTclMultiSupportPattern;

      if (thePattern == 0) {
    opserr << "WARNING ran out of memory creating load pattern - pattern MultipleSupportExcitation ";
    opserr << patternID << endln;

      }
      commandEndMarker = 2;
  }

//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
//Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002
  else if (strcmp(argv[1],"PBowlLoading") == 0) {

      // First search for file name and time step
      TCL_Char * PBEleFileName = 0;
      TCL_Char * accelFileName = 0;
      TCL_Char * displFileName = 0;

      double dt = 0.02;
      double cf = 1.00;
      double xp = 0.0;
      double xm = 0.0;
      double yp = 0.0;
      double ym = 0.0;
      double zp = 0.0;
      double zm = 0.0;

      int currentArg = 3;

      while (currentArg < argc-1)
      {
    if ((strcmp(argv[currentArg],"-pbele") == 0) ||
        (strcmp(argv[currentArg],"-PBEle") == 0)) {

      currentArg++;
      PBEleFileName = argv[currentArg];
      currentArg++;
    }
    else if ((strcmp(argv[currentArg],"-acce") == 0) ||
        (strcmp(argv[currentArg],"-acceleration") == 0)) {

      currentArg++;
      accelFileName = argv[currentArg];
      currentArg++;
    }
    else if ((strcmp(argv[currentArg],"-disp") == 0) ||
         (strcmp(argv[currentArg],"-displacement") == 0)) {

      currentArg++;
      displFileName = argv[currentArg];
      currentArg++;
    }
    else if (strcmp(argv[currentArg],"-dt") == 0)  {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &dt) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "time interval - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }
    else if ((strcmp(argv[currentArg],"-factor") == 0)||
       (strcmp(argv[currentArg],"-Factor") == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &cf) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "load factor - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-xminus") == 0)||
       (strcmp(argv[currentArg],"-xm"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &xm) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Left x  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-xplus") == 0)||
       (strcmp(argv[currentArg],"-xp"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &xp) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Right x  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-yminus") == 0)||
       (strcmp(argv[currentArg],"-ym"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &ym) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Left y  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-yplus") == 0)||
       (strcmp(argv[currentArg],"-yp"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &yp) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Right y  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-zminus") == 0)||
       (strcmp(argv[currentArg],"-zm"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &zm) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Left z  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

    else if ((strcmp(argv[currentArg],"-zplus") == 0)||
       (strcmp(argv[currentArg],"-zp"   ) == 0 )) {

      currentArg++;
          if (Tcl_GetDouble(interp, argv[currentArg], &zp) != TCL_OK)
      {
        opserr << "WARNING problem reading ground motion "
           << "Right y  - pattern PBowlLoading: "
               << patternID << endln;
        return TCL_ERROR;
      }
      currentArg++;
    }

      } //End of while loop

      //test cout << "Tcl parameters dt " << dt << " xp " << xp << " xm " << xm << " yp " << yp << " ym " << ym << " zm " << zm << " done...\n";
      thePattern = new PBowlLoading(patternID, PBEleFileName, displFileName, accelFileName, dt, cf, xp, xm, yp, ym, zp, zm);
      if (thePattern == 0) {
      opserr << "WARNING ran out of memory creating load pattern - pattern PBowlLoading ";
      opserr << patternID << endln;

      }
// Added by Joey Yang to prevent call to Tcl_Eval at end of this function
      commandEndMarker = currentArg;
  }

  else {
      opserr << "WARNING unknown pattern type " << argv[1];
      opserr << " - want: pattern patternType " << patternID ;
      opserr << " \t valid types: Plain, UniformExcitation, MultiSupportExciatation \n";
      return TCL_ERROR;
  }

  // now add the load pattern to the modelBuilder
  if (theDomain->addLoadPattern(thePattern) == false) {
    opserr << "WARNING could not add load pattern to the domain " << *thePattern;
    delete thePattern; // free up the memory, pattern destroys the time series
    return TCL_ERROR;
  }

  theTclLoadPattern = thePattern;

  // use TCL_Eval to evaluate the list of load and single point constraint commands
  if (commandEndMarker < (argc-1)) {
    if (Tcl_Eval(interp, argv[argc-1]) != TCL_OK) {
  opserr << "WARNING - error reading load pattern information in { } ";
  return TCL_ERROR;
    }
  }

  theTclMultiSupportPattern = 0;

  return TCL_OK;
}




