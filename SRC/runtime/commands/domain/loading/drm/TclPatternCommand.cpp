/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the function invoked when the user invokes
// the Pattern command in the interpreter. It is invoked by the
// TclBasicBuilder_addPattern function in the TclBasicBuilder.C file. Current
// valid Pattern types are:
//
// Written: fmk
// Created: 07/99
//
// Modified: fmk 11/00 - removed TimeSeries stuff from file, now an external
// procedure
//
// Modified:  Nov.    2002,  Zhaohui  Yang and Boris Jeremic added Plastic Bowl
// loading (aka Domain Reduction Method) commands
//
#include <BasicModelBuilder.h>

#include <runtimeAPI.h>
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

#include <DRMLoadPattern.h>
#include <DRMInputHandler.h>
#include <PlaneDRMInputHandler.h>

#include <string.h>

#include <SimulationInformation.h>
extern SimulationInformation simulationInfo;
extern const char *getInterpPWD(Tcl_Interp *interp); // interpreter.cpp

LoadPattern *theTclLoadPattern = 0;
MultiSupportPattern *theTclMultiSupportPattern = 0;

extern TimeSeriesIntegrator *TclDispatch_newSeriesIntegrator(ClientData clientData,
                                                        Tcl_Interp *interp,
                                                        TCL_Char *arg);

extern TimeSeries *TclSeriesCommand(ClientData clientData, Tcl_Interp *interp,
                                    TCL_Char *arg);

int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char **argv, Domain *theDomain)
{
  LoadPattern *thePattern = nullptr;

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "invalid command - want: pattern type ";
    opserr << " <type args> {list of load and sp constraints commands}\n";
    opserr
        << "           valid types: Plain, UniformExcitation, MultiSupport\n";
    return TCL_ERROR;
  }

  TimeSeries *theSeries = nullptr;
  int patternID = 0;

  if (Tcl_GetInt(interp, argv[2], &patternID) != TCL_OK) {
    opserr << "WARNING invalid patternID: pattern type " << argv[2]
           << "<type args>\n";
    return TCL_ERROR;
  }

  int commandEndMarker = 0;

  if (strcmp(argv[1], "Plain") == 0) {

    thePattern = new LoadPattern(patternID);
    theSeries = TclSeriesCommand(clientData, interp, argv[3]);

    if (theSeries == nullptr) {

      opserr << "WARNING - problem creating TimeSeries for LoadPattern ";
      opserr << patternID << "\n";

      // clean up the memory and return an error
      if (thePattern != nullptr)
        delete thePattern;

      if (theSeries != nullptr)
        delete theSeries;

      return TCL_ERROR;
    }

    thePattern->setTimeSeries(theSeries);

  }

  else if (strcmp(argv[1], "UniformExcitation") == 0) {

    int dir;

    if (Tcl_GetInt(interp, argv[3], &dir) != TCL_OK) {
      opserr << "WARNING invalid patternID: pattern type " << argv[2]
             << "<type args>\n";
      return TCL_ERROR;
    }

    dir--; // subtract 1 for c indexing

    TimeSeries *accelSeries = nullptr;
    TimeSeries *velSeries = nullptr;
    TimeSeries *dispSeries = nullptr;
    TimeSeriesIntegrator *seriesIntegrator = nullptr;
    double vel0 = 0.0;

    int currentArg = 4;
    bool doneSeries = false;
    while (currentArg < argc - 1 && doneSeries == false) {

      if ((strcmp(argv[currentArg], "-vel0") == 0) ||
          (strcmp(argv[currentArg], "-initialVel") == 0)) {

        currentArg++;
        if ((currentArg < argc) &&
            (Tcl_GetDouble(interp, argv[currentArg], &vel0) != TCL_OK)) {
          opserr << "WARNING invalid vel0: pattern type UniformExcitation\n";
          return TCL_ERROR;
        }

        currentArg++;
      }

      else if ((strcmp(argv[currentArg], "-accel") == 0) ||
               (strcmp(argv[currentArg], "-acceleration") == 0)) {

        currentArg++;
        accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
        if (accelSeries == 0) {
          opserr << "WARNING invalid accel series: pattern UniformExcitation "
                    "-accel <args>\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-vel") == 0) ||
                 (strcmp(argv[currentArg], "-velocity") == 0)) {

        currentArg++;
        velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (velSeries == 0) {
          opserr << "WARNING invalid vel series: " << argv[currentArg];
          opserr << " pattern UniformExcitation -vel {series}\n";
          return TCL_ERROR;
        }

        currentArg++;

      } else if ((strcmp(argv[currentArg], "-disp") == 0) ||
                 (strcmp(argv[currentArg], "-displacement") == 0)) {

        currentArg++;
        dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (dispSeries == 0) {
          opserr << "WARNING invalid vel series: " << argv[currentArg];
          opserr << " pattern UniformExcitation -vel {series}\n";
          return TCL_ERROR;
        }

        currentArg++;
      } else if ((strcmp(argv[currentArg], "-int") == 0) ||
                 (strcmp(argv[currentArg], "-integrator") == 0)) {

        currentArg++;
        seriesIntegrator =
            TclDispatch_newSeriesIntegrator(clientData, interp, argv[currentArg]);
        if (seriesIntegrator == nullptr) {
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

    GroundMotion *theMotion =
        new GroundMotion(dispSeries, velSeries, accelSeries, seriesIntegrator);

    if (theMotion == 0) {
      opserr << "WARNING ran out of memory creating ground motion - pattern "
                "UniformExcitation ";
      opserr << patternID << "\n";

      return TCL_ERROR;
    }

    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0);

    if (thePattern == 0) {
      opserr << "WARNING ran out of memory creating load pattern - pattern "
                "UniformExcitation ";
      opserr << patternID << "\n";

      // clean up memory allocated up to this point and return an error
      if (theMotion != 0)
        delete theMotion;

      return TCL_ERROR;
    }

    // Added by MHS to prevent call to Tcl_Eval at end of this function
    commandEndMarker = currentArg;
  }

  else if (strcmp(argv[1], "Uniform") == 0) {

    // First search for file name and time step
    int numInputs = argc;
    TCL_Char *accelFileName = 0;
    double dt = 0.0;
    for (int i = 5; i < argc; ++i) {
      if (strcmp(argv[i], "-accel") == 0 && i + 2 < argc) {
        // Read the input file name
        accelFileName = argv[i + 1];

        // Read the time interval
        if (Tcl_GetDouble(interp, argv[i + 2], &dt) != TCL_OK) {
          opserr << "WARNING problem reading ground motion "
                 << "time interval - pattern UniformExcitation: " << patternID
                 << "\n";
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
      switch (inputDir) {
      case 'X':
      case 'x':
      case '1': // Global X
        dir = 0;
        break;
      case 'Y':
      case 'y':
      case '2': // Global Y
        dir = 1;
        break;
      case 'Z':
      case 'z':
      case '3': // Global Z
        dir = 2;
        break;
      default:
        opserr << "WARNING cannot read direction for excitation \n";
        opserr << "UniformExcitation " << patternID << " dir factor" << "\n";
        return TCL_ERROR;
        break;
      }
    } else
      dir--; // change to c++ indexing

    double factor;
    if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) {
      opserr << "WARNING insufficient number of arguments - want: pattern ";
      opserr << "UniformExcitation " << patternID << " dir factor\n";
      return TCL_ERROR;
    }

    GroundMotionRecord *theMotion;

    // read in the ground motion
    if (accelFileName == 0) {
      opserr << "WARNING -- No ground motion data provided\n";
      opserr << "UniformExcitation tag: " << patternID << "\n";
      return TCL_ERROR;
    }

    theMotion = new GroundMotionRecord(accelFileName, dt, factor);

    if (theMotion == 0) {
      opserr << "WARNING ran out of memory creating ground motion - pattern "
                "UniformExcitation ";
      opserr << patternID << "\n";
      return TCL_ERROR;
    }

    const char *pwd = getInterpPWD(interp);
    simulationInfo.addInputFile(accelFileName, pwd);

    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID);

    if (thePattern == 0) {
      opserr << "WARNING ran out of memory creating load pattern - pattern "
                "UniformExcitation ";
      opserr << patternID << "\n";

      // clean up memory allocated up to this point and return an error
      if (theMotion != 0)
        delete theMotion;

      return TCL_ERROR;
    }
  }

  else if ((strcmp(argv[1], "MultipleSupportExcitation") == 0) ||
           (strcmp(argv[1], "MultipleSupport") == 0) ||
           (strcmp(argv[1], "MultiSupport") == 0)) {

    theTclMultiSupportPattern = new MultiSupportPattern(patternID);
    thePattern = theTclMultiSupportPattern;

    if (thePattern == 0) {
      opserr << "WARNING ran out of memory creating load pattern - pattern "
                "MultipleSupportExcitation ";
      opserr << patternID << "\n";
    }
    commandEndMarker = 2;
  }

  // Added Joey Yang and Boris Jeremic at UC Davis 10/31/2002

  //////// //////// ///////// ////////// /////  // DRMLoadPattern add BEGIN
  else if (strcmp(argv[1], "DRMLoadPattern") == 0) {
    TCL_Char *InputDataFileName = 0;

    if ((strcmp(argv[2], "-inputdata") == 0) ||
        (strcmp(argv[2], "-InputData") == 0)) {
      InputDataFileName = argv[3];
    }

    // now parse the input file name to extract the pattern input data
    std::ifstream ifile(InputDataFileName);
    int num_steps;
    ifile >> num_steps;
    double dt;
    ifile >> dt;
    int steps_cached;
    ifile >> steps_cached;
    double *ele_d = new double[3];
    ifile >> ele_d[0];
    ifile >> ele_d[1];
    ifile >> ele_d[2];
    double *drm_box_crds = new double[6];
    for (int i = 0; i < 6; ++i)
      ifile >> drm_box_crds[i];
    int n1;
    ifile >> n1;
    int n2;
    ifile >> n2;

    std::string inps;
    int nf = 6;
    char **files = new char *[nf];
    int *f_d = new int[3 * (nf - 1)];
    int ne1, ne2;
    for (int i = 0; i < nf; ++i) {
      ifile >> inps;
      files[i] = (char *)inps.c_str();
      if (i < (nf - 1)) {
        ifile >> ne1;
        ifile >> ne2;
        f_d[3 * i] = (ne1 + 1) * (ne2 + 1);
        f_d[3 * i + 1] = ne1;
        f_d[3 * i + 2] = ne2;
      }
    }

    Mesh3DSubdomain *myMesher = new Mesh3DSubdomain(theDomain);
    PlaneDRMInputHandler *patternhandler = new PlaneDRMInputHandler(
        1.0, files, nf, dt, 0, num_steps, f_d, 15, n1, n2, drm_box_crds,
        drm_box_crds, ele_d, myMesher, steps_cached, theDomain);
    DRMLoadPattern *ptr = new DRMLoadPattern(1, 1.0, patternhandler, theDomain);
    ptr->setMaps();
    thePattern = ptr;

  } // end else if DRMLoadPattern

  ///////// ///////// ///////// //////// //////  // DRMLoadPattern add END

  else {
    opserr << "WARNING unknown pattern type " << argv[1];
    opserr << " - want: pattern patternType " << patternID;
    opserr << " \t valid types: Plain, UniformExcitation, "
              "MultiSupportExciatation \n";
    return TCL_ERROR;
  }

  // now add the load pattern to the modelBuilder
  if (theDomain->addLoadPattern(thePattern) == false) {
    opserr << "WARNING could not add load pattern to the domain "
           << *thePattern;
    delete thePattern; // free up the memory, pattern destroys the time series
    return TCL_ERROR;
  }

  theTclLoadPattern = thePattern;

  // use TCL_Eval to evaluate the list of load and single point constraint
  // commands
  if (commandEndMarker < (argc - 1)) {
    if (Tcl_Eval(interp, argv[argc - 1]) != TCL_OK) {
      opserr << "WARNING - error reading load pattern information in { } ";
      return TCL_ERROR;
    }
  }

  theTclMultiSupportPattern = 0;

  return TCL_OK;
}
