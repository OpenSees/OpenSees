//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
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
#include <BasicModelBuilder.h>

extern TimeSeries *TclSeriesCommand(ClientData clientData, Tcl_Interp *interp,
                                    TCL_Char * const arg);

extern TimeSeriesIntegrator *TclDispatch_newSeriesIntegrator(ClientData clientData,
                                                        Tcl_Interp *interp,
                                                        TCL_Char * const arg);

static int
TclCommand_newGroundMotion(ClientData, Tcl_Interp*,
                       int argc,
                       TCL_Char ** const argv,
                       MultiSupportPattern *thePattern);

int
TclCommand_addGroundMotion(ClientData clientData, Tcl_Interp *interp,
                           int argc, TCL_Char ** const argv)

{
  MultiSupportPattern* pattern = 
    (MultiSupportPattern *)Tcl_GetAssocData(interp,"theTclMultiSupportPattern", NULL);

  if (pattern == nullptr) {
    opserr << "ERROR no multi-support pattern\n";
    return TCL_ERROR;
  }
  return TclCommand_newGroundMotion(clientData, interp, argc, argv, pattern);
}


static int
TclCommand_newGroundMotion(ClientData clientData, Tcl_Interp* interp, int argc,
                       TCL_Char ** const argv, MultiSupportPattern *thePattern)
{

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "invalid command - want: groundMotion tag type <args>\n";
    opserr << "           valid types: AccelRecord and Interpolated \n";
    return TCL_ERROR;
  }

  int gMotionTag;
  if (Tcl_GetInt(interp, argv[1], &gMotionTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid tag: groundMotion tag  type <args>\n";
    return TCL_ERROR;
  }

  int startArg = 2;

  GroundMotion *theMotion = nullptr;
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
        accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (accelSeries == nullptr) {
          opserr << OpenSees::PromptValueError << "invalid accel series: " << argv[currentArg];
          opserr << " groundMotion tag Series -accel {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-vel") == 0) ||
                 (strcmp(argv[currentArg], "-velocity") == 0)) {

        currentArg++;
        velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (velSeries == nullptr) {
          opserr << OpenSees::PromptValueError << "invalid vel series: " << argv[currentArg];
          opserr << " groundMotion tag Series -vel {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-disp") == 0) ||
                 (strcmp(argv[currentArg], "-displacement") == 0)) {

        currentArg++;
        dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (dispSeries == nullptr) {
          opserr << OpenSees::PromptValueError << "invalid disp series: " << argv[currentArg];
          opserr << " groundMotion tag Series -disp {series}\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-int") == 0) ||
                 (strcmp(argv[currentArg], "-integrator") == 0)) {

        currentArg++;
        seriesIntegrator =
            TclDispatch_newSeriesIntegrator((ClientData)0, interp, argv[currentArg]);
        if (seriesIntegrator == nullptr) {
          opserr << OpenSees::PromptValueError << "invalid series integrator: " << argv[currentArg];
          opserr << " - groundMotion tag Series -int {Series Integrator}\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-dtInt") == 0) ||
                 (strcmp(argv[currentArg], "-dtIntegrator") == 0) ||
                 (strcmp(argv[currentArg], "-deltaT") == 0)) {

        currentArg++;
        if (Tcl_GetDouble(interp, argv[currentArg], &dtInt) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid dtInt: " << argv[currentArg];
          opserr << " - groundMotion tag Series -dtInt dt\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-fact") == 0) ||
                 (strcmp(argv[currentArg], "-scale") == 0) ||
                 (strcmp(argv[currentArg], "-factor") == 0)) {

        currentArg++;
        if (Tcl_GetDouble(interp, argv[currentArg], &fact) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid factor: " << argv[currentArg];
          opserr << " - groundMotion tag Series -fact factor\n";
          return TCL_ERROR;
        }

        currentArg++;
      }
    }

    theMotion = new GroundMotion(dispSeries, velSeries, accelSeries,
                                 seriesIntegrator, dtInt, fact);
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
          opserr << OpenSees::PromptValueError << "no groundMotion with tag " << motionID << " :";
          opserr << " pattern MultiSupport gMotion1? gMotion? .. ";
          opserr << "-fact fact1? fact2? .. \n";
          return TCL_ERROR;

        } else
          theMotions[i - 3] = theMotion1;
      }
    } else {
      opserr << OpenSees::PromptValueError << "no gMotionTags want :";
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
    opserr << OpenSees::PromptValueError << "unknown pattern type " << argv[1];
    opserr << " - want: pattern patternType " << gMotionTag;
    opserr << " \t valid types: Plain, UniformExcitation \n";
    return TCL_ERROR;
  }

  // now add the load pattern to the modelBuilder
  if (theMotion != nullptr) {
    if (thePattern->addMotion(*theMotion, gMotionTag) < 0) {
      opserr << OpenSees::PromptValueError << "could not add ground motion with tag " << gMotionTag;
      opserr << " to pattern\n ";
      delete theMotion; // free the memory, pattern destroys the time series
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}
