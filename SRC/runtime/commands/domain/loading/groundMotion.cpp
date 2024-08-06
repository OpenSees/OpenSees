//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the GroundMotion command in the interpreter. 
//
//
// Written: fmk
// Created: 11/00
//
#include <tcl.h>
#include <GroundMotionRecord.h>

#include <string.h>
#include <MultiSupportPattern.h>
#include <InterpolatedGroundMotion.h>
#include <TimeSeriesIntegrator.h>

extern TimeSeries *TclSeriesCommand(ClientData clientData, Tcl_Interp *interp,
                                    TCL_Char * const arg);

extern TimeSeriesIntegrator *TclDispatch_newSeriesIntegrator(ClientData clientData,
                                                        Tcl_Interp *interp,
                                                        TCL_Char * const arg);

static int
TclCommand_newGroundMotion(G3_Runtime* rt,
                       int argc,
                       TCL_Char ** const argv,
                       MultiSupportPattern *thePattern);

int
TclCommand_addGroundMotion(ClientData clientData, Tcl_Interp *interp,
                           int argc, TCL_Char ** const argv)

{
  G3_Runtime *rt = G3_getRuntime(interp);
  MultiSupportPattern* pattern = 
    (MultiSupportPattern *)Tcl_GetAssocData(interp,"theTclMultiSupportPattern", NULL);

  if (pattern == nullptr) {
    opserr << "ERROR no multi-support pattern\n";
    return TCL_ERROR;
  }
  return TclCommand_newGroundMotion(rt, argc, argv, pattern);
}


static int
TclCommand_newGroundMotion(G3_Runtime* rt, int argc,
                       TCL_Char ** const argv, MultiSupportPattern *thePattern)
{

  int gMotionTag;
  GroundMotion *theMotion = nullptr;
  Tcl_Interp *interp = G3_getInterpreter(rt);

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
    opserr << "WARNING invalid command - want: groundMotion tag type <args>\n";
    opserr << "           valid types: AccelRecord and Interpolated \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[1], &gMotionTag) != TCL_OK) {
    opserr << "WARNING invalid tag: groundMotion tag  type <args>\n";
    return TCL_ERROR;
  }

  int startArg = 2;

  if ((strcmp(argv[startArg], "Series") == 0) ||
      (strcmp(argv[startArg], "Plain") == 0)) {

    TimeSeries *accelSeries = nullptr;
    TimeSeries *velSeries   = nullptr;
    TimeSeries *dispSeries  = nullptr;
    TimeSeriesIntegrator *seriesIntegrator = nullptr;

    int currentArg = startArg + 1;
    double dtInt = 0.01;
    double fact = 1.0;

    while (currentArg < argc - 1) {
      if ((strcmp(argv[currentArg], "-accel") == 0) ||
          (strcmp(argv[currentArg], "-acceleration") == 0)) {

        currentArg++;
        accelSeries = TclSeriesCommand((ClientData)0, interp, argv[currentArg]);

        if (accelSeries == 0) {
          opserr << "WARNING invalid accel series: " << argv[currentArg];
          opserr << " groundMotion tag Series -accel {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-vel") == 0) ||
                 (strcmp(argv[currentArg], "-velocity") == 0)) {

        currentArg++;
        velSeries = TclSeriesCommand((ClientData)0, interp, argv[currentArg]);

        if (velSeries == 0) {
          opserr << "WARNING invalid vel series: " << argv[currentArg];
          opserr << " groundMotion tag Series -vel {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-disp") == 0) ||
                 (strcmp(argv[currentArg], "-displacement") == 0)) {

        currentArg++;
        dispSeries = TclSeriesCommand((ClientData)0, interp, argv[currentArg]);

        if (dispSeries == 0) {
          opserr << "WARNING invalid disp series: " << argv[currentArg];
          opserr << " groundMotion tag Series -disp {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-int") == 0) ||
                 (strcmp(argv[currentArg], "-integrator") == 0)) {

        currentArg++;
        seriesIntegrator =
            TclDispatch_newSeriesIntegrator((ClientData)0, interp, argv[currentArg]);
        if (seriesIntegrator == 0) {
          opserr << "WARNING invalid series integrator: " << argv[currentArg];
          opserr << " - groundMotion tag Series -int {Series Integrator}\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-dtInt") == 0) ||
                 (strcmp(argv[currentArg], "-dtIntegrator") == 0) ||
                 (strcmp(argv[currentArg], "-deltaT") == 0)) {

        currentArg++;
        if (Tcl_GetDouble(interp, argv[currentArg], &dtInt) != TCL_OK) {
          opserr << "WARNING invalid dtInt: " << argv[currentArg];
          opserr << " - groundMotion tag Series -dtInt dt\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-fact") == 0) ||
                 (strcmp(argv[currentArg], "-scale") == 0) ||
                 (strcmp(argv[currentArg], "-factor") == 0)) {

        currentArg++;
        if (Tcl_GetDouble(interp, argv[currentArg], &fact) != TCL_OK) {
          opserr << "WARNING invalid factor: " << argv[currentArg];
          opserr << " - groundMotion tag Series -fact factor\n";
          return TCL_ERROR;
        }

        currentArg++;
      }
    }

    theMotion = new GroundMotion(dispSeries, velSeries, accelSeries,
                                 seriesIntegrator, dtInt, fact);

    if (theMotion == 0) {
      opserr << "WARNING ran out of memory creating ground motion - pattern "
                "UniformExcitation ";
      opserr << gMotionTag << endln;

      return TCL_ERROR;
    }
  }

  else if (strcmp(argv[startArg], "Interpolated") == 0) {

    int endMotionIDs = startArg + 1;
    int motionID;
    while (Tcl_GetInt(interp, argv[endMotionIDs], &motionID) == TCL_OK)
      endMotionIDs++;

    int numMotions = endMotionIDs - startArg - 1;
    GroundMotion **theMotions;

    if (numMotions != 0) {
      // TODO: LEAKED?
      theMotions = new GroundMotion *[numMotions];
      ID motionIDs(numMotions);
      for (int i = 3; i < endMotionIDs; ++i) {

        if (Tcl_GetInt(interp, argv[i], &motionID) != TCL_OK)
          return TCL_ERROR;

        motionIDs[i - 3] = motionID;

        GroundMotion *theMotion1 = thePattern->getMotion(motionID);
        if (theMotion1 == nullptr) {
          opserr << "WARNING no groundMotion with tag " << motionID << " :";
          opserr << " pattern MultiSupport gMotion1? gMotion? .. ";
          opserr << "-fact fact1? fact2? .. \n";
          return TCL_ERROR;

        } else
          theMotions[i - 3] = theMotion1;
      }
    } else {
      opserr << "WARNING no gMotionTags want :";
      opserr << " pattern MultiSupport gMotion1? gMotion? .. ";
      opserr << "-fact fact1? fact2? .. \n";
      return TCL_ERROR;
    }
    Vector facts(numMotions);
    for (int i = 0; i < numMotions; ++i) {
      double fact;
      if (Tcl_GetDouble(interp, argv[i + endMotionIDs + 1], &fact) != TCL_OK)
        return TCL_ERROR;
      facts[i] = fact;
    }

    theMotion = new InterpolatedGroundMotion(theMotions, facts, false);

  }

  else {
    opserr << "WARNING unknown pattern type " << argv[1];
    opserr << " - want: pattern patternType " << gMotionTag;
    opserr << " \t valid types: Plain, UniformExcitation \n";
    return TCL_ERROR;
  }

  // now add the load pattern to the modelBuilder
  if (theMotion != nullptr) {
    if (thePattern->addMotion(*theMotion, gMotionTag) < 0) {
      opserr << "WARNING could not add ground motion with tag " << gMotionTag;
      opserr << " to pattern\n ";
      delete theMotion; // free the memory, pattern destroys the time series
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}
