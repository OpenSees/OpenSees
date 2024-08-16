//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user invokes
// the "pattern" command in the interpreter. 
//
// The command takes the form:
//
//     pattern <type> <tag> {
//      <load-like>...
//     };
// 
// SIDE EFFECTS:
// - Create an instance of class <type> with tag <tag>.
// - Add this instance to the domain
//
// Written: fmk
// Created: 07/99
//
// Modified: fmk 11/00 - removed TimeSeries stuff from file, now an external
// procedure
//
#include <tcl.h>
#include <assert.h>
#include <BasicModelBuilder.h>
#include <runtimeAPI.h>

#include <runtimeAPI.h>
#include <G3_Logging.h>

#include <Domain.h>
#include <LoadPattern.h>

#include <UniformExcitation.h>
#include <MultiSupportPattern.h>
#include <GroundMotion.h>
#include <GroundMotionRecord.h>
#include <TimeSeriesIntegrator.h>

#ifdef OPSDEF_DRM
#  include <DRMLoadPattern.h>
#  include <DRMInputHandler.h>
#  include <PlaneDRMInputHandler.h>
#  include <DRMLoadPatternWrapper.h>
#endif

#ifdef _H5DRM
#  include <H5DRM.h>
#endif

#include <NodalThermalAction.h>   //L.Jiang [SIF]
#include <NodalLoad.h>

#include <string.h>



Tcl_CmdProc TclCommand_addSP;
Tcl_CmdProc TclCommand_addNodalLoad;

extern TimeSeriesIntegrator *TclDispatch_newSeriesIntegrator(ClientData clientData,
                                                        Tcl_Interp *interp,
                                                        TCL_Char * const arg);

extern TimeSeries *TclSeriesCommand(ClientData clientData, Tcl_Interp *interp,
                                    TCL_Char * const arg);

//
// This command creates a scope where the following commands
// behave differently:
// - load
// - sp
int
TclCommand_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,
                      TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain* domain = builder->getDomain();
  LoadPattern *thePattern = nullptr;

  // make sure at least one other argument to contain integrator
  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "invalid command - want: pattern type ";
    opserr << " <type args> {list of load and sp constraints commands}\n";
    opserr << "           valid types: Plain, UniformExcitation, MultiSupport\n";
    return TCL_ERROR;
  }

  int commandEndMarker = 3;

  TimeSeries *theSeries = nullptr;

  int patternID;
  if (Tcl_GetInt(interp, argv[2], &patternID) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid patternID: " << argv[2] << "\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "Plain") == 0) {
    if (argc < 4) {
        opserr << G3_ERROR_PROMPT << "Invalid command for Plain pattern.\n";
        return TCL_ERROR;
    }

    double fact = 1.0;
    const char* series_arg = nullptr;

    while (commandEndMarker < argc) {
      if ((strcmp(argv[commandEndMarker], "-fact") == 0) ||
          (strcmp(argv[commandEndMarker], "-factor") == 0)) {

        if (Tcl_GetDouble(interp, argv[++commandEndMarker], &fact) != TCL_OK) {
          opserr << G3_ERROR_PROMPT
                 << "invalid factor: " << argv[commandEndMarker] << "\n";
          return TCL_ERROR;
        }
        commandEndMarker++;

      } else if (series_arg == nullptr) {
        series_arg = argv[commandEndMarker];
        commandEndMarker++;

      } else {
        break;
      }
    }

    thePattern = new LoadPattern(patternID, fact);
    theSeries = TclSeriesCommand(clientData, interp, series_arg);

    if (theSeries == nullptr) {
      opserr << G3_ERROR_PROMPT << "problem creating TimeSeries for LoadPattern "
             << patternID << endln;

      // clean up the memory and return an error
      if (thePattern != nullptr)
        delete thePattern;
      return TCL_ERROR;
    }

    thePattern->setTimeSeries(theSeries);
  }

  else if (strcmp(argv[1], "UniformExcitation") == 0) {

    int dir;

    if (Tcl_GetInt(interp, argv[3], &dir) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid patternID: pattern type " << argv[2]
             << "<type args>\n";
      return TCL_ERROR;
    }

    dir--; // subtract 1 for C indexing

    TimeSeries *accelSeries = nullptr;
    TimeSeries *velSeries   = nullptr;
    TimeSeries *dispSeries  = nullptr;
    TimeSeriesIntegrator *seriesIntegrator = nullptr;
    double vel0 = 0.0;
    double fact = 1.0;

    int currentArg = 4;
    bool doneSeries = false;
    while (currentArg < argc - 1 && doneSeries == false) {

      if ((strcmp(argv[currentArg], "-vel0") == 0) ||
          (strcmp(argv[currentArg], "-initialVel") == 0)) {

        currentArg++;
        if ((currentArg < argc) &&
            (Tcl_GetDouble(interp, argv[currentArg], &vel0) != TCL_OK)) {
          opserr << OpenSees::PromptValueError << "invalid vel0: pattern type UniformExcitation\n";
          return TCL_ERROR;
        }

        currentArg++;
      }

      else if ((strcmp(argv[currentArg], "-fact") == 0) ||
               (strcmp(argv[currentArg], "-factor") == 0)) {

        currentArg++;
        if ((currentArg < argc) &&
            (Tcl_GetDouble(interp, argv[currentArg], &fact) != TCL_OK)) {
          opserr << OpenSees::PromptValueError << "invalid fact: pattern type UniformExcitation\n";
          return TCL_ERROR;
        }

        currentArg++;
      }

      else if ((strcmp(argv[currentArg], "-accel") == 0) ||
               (strcmp(argv[currentArg], "-acceleration") == 0)) {

        currentArg++;
        accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (accelSeries == nullptr) {
          // Assume TclSeriesCommand printed error info
          opserr << "      in <series> for -accel flag of: '" << argv[currentArg];
          opserr << "'\n      pattern UniformExcitation ... -accel <series>\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-vel") == 0) ||
                 (strcmp(argv[currentArg], "-velocity") == 0)) {

        currentArg++;
        velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (velSeries == nullptr) {
          // Assume TclSeriesCommand printed error info
          opserr << "      in <series> for -vel[ocity] flag of: '" << argv[currentArg];
          opserr << "'\n      pattern UniformExcitation ... -velocity <series>\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-disp") == 0) ||
                 (strcmp(argv[currentArg], "-displacement") == 0)) {

        currentArg++;
        dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

        if (dispSeries == nullptr) {
          // Assume TclSeriesCommand printed error info
          opserr << "      in <series> for -disp[lacement] flag '" << argv[currentArg] << "' of ";
          opserr << "'\n      pattern UniformExcitation ... -displacement <series>\n";
          return TCL_ERROR;
        }
        currentArg++;

      } else if ((strcmp(argv[currentArg], "-int") == 0) ||
                 (strcmp(argv[currentArg], "-integrator") == 0)) {

        currentArg++;
        seriesIntegrator =
            TclDispatch_newSeriesIntegrator(clientData, interp, argv[currentArg]);
        if (seriesIntegrator == nullptr) {
          opserr << OpenSees::PromptValueError << "invalid series integrator: " << argv[currentArg];
          opserr << " - pattern UniformExcitation -int {Series Integrator}\n";
          return TCL_ERROR;
        }
        currentArg++;
      }
      else
        doneSeries = true;
    }

    if (dispSeries == nullptr && velSeries == nullptr && accelSeries == nullptr) {
      opserr << G3_ERROR_PROMPT << "invalid series, expected:\n    pattern UniformExcitation";
      opserr << "-disp {dispSeries} -vel {velSeries} -accel {accelSeries} ";
      opserr << "-int {Series Integrator}" << "\n";
      return TCL_ERROR;
    }

    GroundMotion *theMotion =
        new GroundMotion(dispSeries, velSeries, accelSeries, seriesIntegrator);

    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0, fact);

    // Added by MHS to prevent call to Tcl_Eval at end of this function
    commandEndMarker = currentArg+1;
  }

  else if (strcmp(argv[1], "Uniform") == 0) {

    // First search for file name and time step
    int numInputs = argc;
    TCL_Char *accelFileName = nullptr;
    double dt = 0.0;
    for (int i = 5; i < argc; ++i) {
      if (strcmp(argv[i], "-accel") == 0 && i + 2 < argc) {
        // Read the input file name
        accelFileName = argv[i + 1];

        // Read the time interval
        if (Tcl_GetDouble(interp, argv[i + 2], &dt) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "problem reading ground motion "
                 << "time interval - pattern UniformExcitation: " << patternID
                 << endln;
          return TCL_ERROR;
        }
        numInputs -= 3;
      }
    }

    if (numInputs < 5) {
      opserr << OpenSees::PromptValueError << "insufficient number of arguments - want: pattern ";
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
        opserr << OpenSees::PromptValueError << "cannot read direction for excitation \n";
        opserr << "UniformExcitation " << patternID << " dir factor" << endln;
        return TCL_ERROR;
        break;
      }
    } else
      dir--; // change to c++ indexing

    double factor;
    if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "insufficient number of arguments - want: pattern ";
      opserr << "UniformExcitation " << patternID << " dir factor\n";
      return TCL_ERROR;
    }

    GroundMotionRecord *theMotion = nullptr;

    // Read in the ground motion
    if (accelFileName == 0) {
      opserr << G3_ERROR_PROMPT << "No ground motion data provided\n";
      opserr << "UniformExcitation tag: " << patternID << endln;
      return TCL_ERROR;
    }

    theMotion = new GroundMotionRecord(accelFileName, dt, factor);

    // Create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID);
    Tcl_SetAssocData(interp,"theTclMultiSupportPattern", NULL, (ClientData)0);

  }

  else if ((strcmp(argv[1], "MultipleSupportExcitation") == 0) ||
           (strcmp(argv[1], "MultipleSupport") == 0) ||
           (strcmp(argv[1], "MultiSupport") == 0)) {

    MultiSupportPattern *theTclMultiSupportPattern = new MultiSupportPattern(patternID);
    Tcl_SetAssocData(interp,"theTclMultiSupportPattern", NULL, (ClientData)theTclMultiSupportPattern);
    thePattern = theTclMultiSupportPattern;

//  commandEndMarker = 2;
  }


#ifdef OPSDEF_DRM
  //////// //////// ///////// ////////// /////  // DRMLoadPattern add BEGIN
  else if (strcmp(argv[1], "DRMLoadPattern") == 0) {
    TCL_Char *InputDataFileName = nullptr;

    if ((strcmp(argv[2], "-inputdata") == 0) ||
        (strcmp(argv[2], "-InputData") == 0)) {
      InputDataFileName = argv[3];

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

      Mesh3DSubdomain *myMesher = new Mesh3DSubdomain(domain);
      PlaneDRMInputHandler *patternhandler = new PlaneDRMInputHandler(
          1.0, files, nf, dt, 0, num_steps, f_d, 15, n1, n2, drm_box_crds,
          drm_box_crds, ele_d, myMesher, steps_cached, domain);
      DRMLoadPattern *ptr =
          new DRMLoadPattern(1, 1.0, patternhandler, domain);
      ptr->setMaps();
      thePattern = ptr;
      theTclMultiSupportPattern = 0;
      Tcl_SetAssocData(interp,"theTclMultiSupportPattern", NULL, (ClientData)0);
    } else {

      //     TCL_Char * ifp = 0;
      double INVALID = 0.7111722273337;
      int c_arg = 3;
      int end = argc - 1;
      double dt = INVALID;
      double *ele_d = new double[3];
      double *drm_box_crds = new double[6];
      for (int i = 0; i < 3; ++i) {
        ele_d[i] = INVALID;
        drm_box_crds[2 * i] = INVALID;
        drm_box_crds[2 * i + 1] = INVALID;
      }

      int nf = 6;
      char **files = new char *[nf];
      files[5] = "./NONE";
      int *f_d = new int[15];
      int num_steps = 1;
      int steps_cached = 10;
      int n1, n2;
      n2 = 0;
      double factor = 1.0;

      while (c_arg < end) {

        if ((strcmp(argv[c_arg], "-dt") == 0) ||
            (strcmp(argv[c_arg], "-deltaT") == 0)) {
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &dt) != TCL_OK) {
            opserr << " Error reading deltaT for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        if ((strcmp(argv[c_arg], "-numSteps") == 0) ||
            (strcmp(argv[c_arg], "-numberOfSteps") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &num_steps) != TCL_OK) {
            opserr << " Error reading number of steps for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-stepsCached") == 0) ||
                 (strcmp(argv[c_arg], "-cache") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &steps_cached) != TCL_OK) {
            opserr << " Error reading number of steps for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-gridSize") == 0) ||
                 (strcmp(argv[c_arg], "-eleSize") == 0)) {
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &ele_d[0]) != TCL_OK) {
            opserr << " Error reading deltaT for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &ele_d[1]) != TCL_OK) {
            opserr << " Error reading deltaT for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &ele_d[2]) != TCL_OK) {
            opserr << " Error reading deltaT for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-gridDataFace1") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[1]) != TCL_OK) {
            opserr << " Error reading grid data f1 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[2]) != TCL_OK) {
            opserr << " Error reading grid data f1 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          f_d[0] = (f_d[1] + 1) * (f_d[2] + 1);
        }

        else if ((strcmp(argv[c_arg], "-gridDataFace2") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[4]) != TCL_OK) {
            opserr << " Error reading grid data f2 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[5]) != TCL_OK) {
            opserr << " Error reading grid data f2 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          f_d[3] = (f_d[4] + 1) * (f_d[5] + 1);
        }

        else if ((strcmp(argv[c_arg], "-gridDataFace3") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[7]) != TCL_OK) {
            opserr << " Error reading grid data f3 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[8]) != TCL_OK) {
            opserr << " Error reading grid data f3 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          f_d[6] = (f_d[7] + 1) * (f_d[8] + 1);
        }

        else if ((strcmp(argv[c_arg], "-gridDataFace4") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[10]) != TCL_OK) {
            opserr << " Error reading grid data f4 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[11]) != TCL_OK) {
            opserr << " Error reading grid data f4 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          f_d[9] = (f_d[10] + 1) * (f_d[11] + 1);
        }

        else if ((strcmp(argv[c_arg], "-gridDataFace5") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[13]) != TCL_OK) {
            opserr << " Error reading grid data f5 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &f_d[14]) != TCL_OK) {
            opserr << " Error reading grid data f5 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          f_d[12] = (f_d[13] + 1) * (f_d[14] + 1);
        }

        else if ((strcmp(argv[c_arg], "-filePathFace1") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[0] = new char[tmp.size() + 1];
          strcpy(files[0], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-filePathFace2") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[1] = new char[tmp.size() + 1];
          strcpy(files[1], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-filePathFace3") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[2] = new char[tmp.size() + 1];
          strcpy(files[2], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-filePathFace4") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[3] = new char[tmp.size() + 1];
          strcpy(files[3], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-filePathFace5a") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[4] = new char[tmp.size() + 1];
          strcpy(files[4], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-filePathFace5b") == 0)) {
          c_arg++;
          std::string tmp(argv[c_arg]);
          files[5] = new char[tmp.size() + 1];
          strcpy(files[5], tmp.c_str());
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-fileFace5aGridPoints") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &n1) != TCL_OK) {
            opserr << " Error reading grid data f5 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-fileFace5bGridPoints") == 0)) {
          c_arg++;
          if (Tcl_GetInt(interp, argv[c_arg], &n2) != TCL_OK) {
            opserr << " Error reading grid data f5 for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-factor") == 0) ||
                 (strcmp(argv[c_arg], "-Factor") == 0)) {
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &factor) != TCL_OK) {
            opserr << " Error reading number of steps for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }

        else if ((strcmp(argv[c_arg], "-DRMBoxCrds") == 0)) {
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[0]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, xmin for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[1]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, xmax for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[2]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, ymin for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[3]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, ymax for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[4]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, zmin for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
          if (Tcl_GetDouble(interp, argv[c_arg], &drm_box_crds[5]) != TCL_OK) {
            opserr << " Error reading DRM box Crds, zmax for DRMLoadPattern \n";
            return TCL_ERROR;
          }
          c_arg++;
        }
      }

      thePattern = new DRMLoadPatternWrapper(patternID, factor, files, nf, dt,
                                             num_steps, f_d, 15, n1, n2,
                                             drm_box_crds, ele_d, steps_cached);
      // theTclMultiSupportPattern = 0;
      Tcl_SetAssocData(interp,"theTclMultiSupportPattern", NULL, (ClientData)0);
      commandEndMarker = c_arg+1;
    }

  }    // end else if DRMLoadPattern
#endif // OPSDEF_DRM

#ifdef _H5DRM
  else if ((strcmp(argv[1], "H5DRM") == 0) || 
           (strcmp(argv[1], "h5drm") == 0)) {

    int tag = 0;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "insufficient number of arguments - want: pattern ";
      opserr << "H5DRM tag filename factor\n";
      return TCL_ERROR;
    }

    std::string filename = argv[3];

    double factor = 1.0;
    if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "insufficient number of arguments - want: pattern ";
      opserr << "H5DRM " << patternID << " filename factor\n";
      return TCL_ERROR;
    }

    opserr << "Creating H5DRM tag = " << tag
           << " filename = " << filename.c_str() << " factor = " << factor
           << endln;

    thePattern = new H5DRM(tag, filename, factor);

    opserr << "Done! Creating H5DRM tag = " << tag
           << " filename = " << filename.c_str() << " factor = " << factor
           << endln;

    domain->addLoadPattern(thePattern);
    return TCL_OK;
  }
#endif

  else {
    opserr << OpenSees::PromptValueError << "unknown pattern type " << argv[1];
    opserr << " \t valid types: Plain, UniformExcitation, "
              "MultiSupportExciatation \n";
    return TCL_ERROR;
  }

  // now add the load pattern to the modelBuilder
  if (domain->addLoadPattern(thePattern) == false) {
    opserr << OpenSees::PromptValueError << "could not add load pattern to the domain "
           << *thePattern;
    delete thePattern;
    return TCL_ERROR;
  }


  builder->setEnclosingPattern(thePattern);

  // use TCL_Eval to evaluate the list of load and single point constraint
  // commands

  if (commandEndMarker < argc) {
    // Set the Pattern for "sp" command
//  Tcl_CmdInfo info;
//  assert(Tcl_GetCommandInfo(interp, "sp", &info)==1);
//  info.clientData = (ClientData)thePattern;
//  Tcl_SetCommandInfo(interp, "sp", &info);

    // Tcl_CreateCommand(interp, "nodalLoad", TclCommand_addNodalLoad, (ClientData)thePattern, NULL);
    Tcl_Eval(interp, "rename load opensees::import;");
    Tcl_Eval(interp, "rename nodalLoad load;");
    if (Tcl_Eval(interp, argv[commandEndMarker]) != TCL_OK) {
      // opserr << OpenSees::PromptValueError << "- error reading load pattern information in { }";
      opserr << G3_ERROR_PROMPT << Tcl_GetStringResult(interp);
//    Tcl_Eval(interp, "puts $errorInfo; flush stdout;");
//    Tcl_Exit(TCL_ERROR);
      return TCL_ERROR;
    }

    Tcl_SetAssocData(interp,"theTclMultiSupportPattern", NULL, (ClientData)0);
//  info.clientData = (ClientData)builder;
//  Tcl_SetCommandInfo(interp, "sp", &info);

    Tcl_Eval(interp, "rename load nodalLoad;");
    Tcl_Eval(interp, "rename opensees::import load;");
  }


  return TCL_OK;
}

int
TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  
  // TODO
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  LoadPattern *theTclLoadPattern = builder->getEnclosingPattern();
  int nodeLoadTag = builder->getNodalLoadTag();

  int ndf = builder->getNDF(); // argc - 2;
  NodalLoad *theLoad = nullptr;

  bool isLoadConst = false;
  bool explicitPatternPassed = false;
  int  loadPatternTag = 0;

  if (true) {
    // make sure at least one other argument to contain type
    if (argc < (2 + ndf)) {
      opserr << OpenSees::PromptValueError << "bad command - want: load nodeId " << ndf << " forces\n";
      return TCL_ERROR;
    }

    // get the id of the node
    int nodeId;
    if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nodeId: " << argv[1];
      opserr << " - load nodeId " << ndf << " forces\n";
      return TCL_ERROR;
    }

    // get the load vector
    Vector forces(ndf);
    for (int i = 0; i < ndf; ++i) {
      double theForce;
      if (Tcl_GetDouble(interp, argv[2 + i], &theForce) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid force " << i + 1 << " in load " << nodeId;
        opserr << ", got " << ndf << " forces\n";
        return TCL_ERROR;
      } else
        forces(i) = theForce;
    }

    // allow some additional options at end of command
    int endMarker = 2 + ndf;
    while (endMarker != argc) {
      if (strcmp(argv[endMarker], "-const") == 0) {
        // allow user to specify const load
        isLoadConst = true;
      } else if (strcmp(argv[endMarker], "-pattern") == 0) {
        // allow user to specify load pattern other than current
        endMarker++;
        explicitPatternPassed = true;
        if (endMarker == argc ||
            Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

          opserr << OpenSees::PromptValueError << "invalid patternTag - load " << nodeId << " ";
          opserr << ndf << " forces pattern patterntag\n";
          return TCL_ERROR;
        }
      }
      endMarker++;
    }

    // get the current pattern tag if no tag given in i/p
    if (explicitPatternPassed == false) {
      if (theTclLoadPattern == nullptr) {
        opserr << OpenSees::PromptParseError << "no current load pattern - load " << nodeId;
        opserr << " " << ndf << " forces\n";
        return TCL_ERROR;
      } else
        loadPatternTag = theTclLoadPattern->getTag();
    }

    // create the load
    theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
  }

  // add the load to the domain
  if (builder->getDomain()->addNodalLoad(theLoad, loadPatternTag) == false) {
    opserr << OpenSees::PromptValueError << "BasicModelBuilder - could not add load to domain\n";
    delete theLoad;
    return TCL_ERROR;
  }
  builder->incrNodalLoadTag();

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}

