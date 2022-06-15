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
// $Date: 2007/09/29 01:54:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclPatternCommand.cpp,v $

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

#include <FireLoadPattern.h>  //Added by UoE openSees Group

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

#include <FireLoadPattern.h>  //Added by UoE openSees Group

#include <DRMLoadPattern.h>
#include <DRMInputHandler.h>
#include <PlaneDRMInputHandler.h>
#include <DRMLoadPatternWrapper.h>

#ifdef _H5DRM
#include <H5DRM.h>
#endif

#include <string.h>

#include <SimulationInformation.h>
extern SimulationInformation simulationInfo;
//extern const char * getInterpPWD(Tcl_Interp *interp);  // commands.cpp

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

    double fact = 1.0;
    if (argc==7 && ((strcmp(argv[4],"-fact") == 0) ||
        (strcmp(argv[4],"-factor") == 0))) {
            if (Tcl_GetDouble(interp, argv[5], &fact) != TCL_OK) {
                opserr << "WARNING invalid fact: pattern type Plain\n";
                return TCL_ERROR;
            }
    }
    thePattern = new LoadPattern(patternID,fact);
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

  //--Adding FireLoadPattern:[BEGIN] by UoE OpenSees Group--//
   else if (strcmp(argv[1],"Fire") == 0) {
     FireLoadPattern *theFirePattern = new FireLoadPattern(patternID);
     thePattern = theFirePattern;
     TimeSeries *theSeries1 = TclSeriesCommand(clientData, interp, argv[3]);
     TimeSeries *theSeries2 = TclSeriesCommand(clientData, interp, argv[4]);
     TimeSeries *theSeries3 = TclSeriesCommand(clientData, interp, argv[5]);
     TimeSeries *theSeries4 = TclSeriesCommand(clientData, interp, argv[6]);
     TimeSeries *theSeries5 = TclSeriesCommand(clientData, interp, argv[7]);
     TimeSeries *theSeries6 = TclSeriesCommand(clientData, interp, argv[8]);
     TimeSeries *theSeries7 = TclSeriesCommand(clientData, interp, argv[9]);
     TimeSeries *theSeries8 = TclSeriesCommand(clientData, interp, argv[10]);
     TimeSeries *theSeries9 = TclSeriesCommand(clientData, interp, argv[11]);

     //opserr << "series1 ";
     //*theSeries1->Print;
     // if (thePattern == 0 || theSeries == 0) {
     
     if (thePattern == 0) {
       opserr << "WARNING - out of memory creating LoadPattern ";
       opserr << patternID << endln;
     } 
     else if (theSeries1 == 0) {
       opserr << "WARNING - problem creating TimeSeries1 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries2 == 0) {
       opserr << "WARNING - problem creating TimeSeries2 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries3 == 0) {
       opserr << "WARNING - problem creating TimeSeries3 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries4 == 0) {
       opserr << "WARNING - problem creating TimeSeries4 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries5 == 0) {
       opserr << "WARNING - problem creating TimeSeries5 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries6 == 0) {
       opserr << "WARNING - problem creating TimeSeries6 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries7 == 0) {
       opserr << "WARNING - problem creating TimeSeries7 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries8 == 0) {
       opserr << "WARNING - problem creating TimeSeries8 for LoadPattern ";
       opserr << patternID << endln;
     }
     else if (theSeries9 == 0){
       opserr << "WARNING - problem creating TimeSeries9 for LoadPattern ";
       opserr << patternID << endln;
     }
     
     // clean up the memory and return an error
     //if (thePattern != 0)
     //   delete thePattern;
     //if (theSeries != 0)
     //   delete theSeries;
     //return TCL_ERROR;
     // }
     
     theFirePattern->setFireTimeSeries(theSeries1, theSeries2, 
                       theSeries3, theSeries4, theSeries5, 
                       theSeries6, theSeries7, theSeries8, theSeries9);
     
   }
  //--Adding FireLoadPattern:[END] by UoE OpenSees Group--//

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
    double fact = 1.0;
    
    int currentArg = 4;
    bool doneSeries = false;
    while (currentArg < argc-1 && doneSeries == false) {
      
      if ((strcmp(argv[currentArg],"-vel0") == 0) ||
      (strcmp(argv[currentArg],"-initialVel") == 0)) {
    
    currentArg++;
    if ((currentArg < argc) &&
        (Tcl_GetDouble(interp, argv[currentArg], &vel0) != TCL_OK)) {
      opserr << "WARNING invalid vel0: pattern type UniformExcitation\n";
      return TCL_ERROR;
    }
    
    currentArg++;
      }
      
      else if ((strcmp(argv[currentArg],"-fact") == 0) || (strcmp(argv[currentArg],"-factor") == 0)) {
    
    currentArg++;
    if ((currentArg < argc) &&
        (Tcl_GetDouble(interp, argv[currentArg], &fact) != TCL_OK)) {
      opserr << "WARNING invalid fact: pattern type UniformExcitation\n";
      return TCL_ERROR;
    }
    
    currentArg++;
      }
      
      
      else if ((strcmp(argv[currentArg],"-accel") == 0) ||
           (strcmp(argv[currentArg],"-acceleration") == 0)) {
    
    currentArg++;
    accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
    
    if (accelSeries == 0) {
      opserr << "WARNING invalid accel series: " << argv[currentArg];
      opserr << " pattern UniformExcitation -accel {series}\n";
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
      opserr << "WARNING invalid disp series: " << argv[currentArg];
      opserr << " pattern UniformExcitation -disp {series}\n";
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
    thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0, fact);
    
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
    if (accelFileName == 0) {
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

    //    const char *pwd = getInterpPWD(interp);
    //    simulationInfo.addInputFile(accelFileName, pwd);  

    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID);
    theTclMultiSupportPattern = 0;

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
  
#ifdef _H5DRM
  else if ((strcmp(argv[1],"H5DRM") == 0) ||
       (strcmp(argv[1],"h5drm") == 0) ) 
  {
      int tag = 0;
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
          opserr << "WARNING insufficient number of arguments - want: pattern ";
          opserr << "H5DRM tag filename\n";
          return TCL_ERROR;
      }

      std::string filename = argv[3];

      double factor=1.0;
      if(argc > 4)
      {
        if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) 
        {
              opserr << "WARNING insufficient number of arguments - want: pattern ";
              opserr << "H5DRM " << patternID << " filename factor\n";
              return TCL_ERROR;
        }
      }

      double crd_scale=1.0;
      if(argc > 5)
      {
        if (Tcl_GetDouble(interp, argv[5], &crd_scale) != TCL_OK) 
        {
              opserr << "WARNING insufficient number of arguments - want: pattern ";
              opserr << "H5DRM " << patternID << " filename factor crd_scale\n";
              return TCL_ERROR;
        }
      }
      double distance_tolerance=1e-3;
      if(argc > 6)
      {
        if (Tcl_GetDouble(interp, argv[6], &distance_tolerance) != TCL_OK) 
        {
              opserr << "WARNING insufficient number of arguments - want: pattern ";
              opserr << "H5DRM " << patternID << " filename factor crd_scale distance_tolerance\n";
              return TCL_ERROR;
        }
      }

      int do_coordinate_transformation=true;
      if(argc > 7)
      {
        if (Tcl_GetInt(interp, argv[7], &do_coordinate_transformation) != TCL_OK) 
        {
              opserr << "WARNING insufficient number of arguments - want: pattern ";
              opserr << "H5DRM " << patternID << " filename factor crd_scale distance_tolerance do_coordinate_transformation\n";
              return TCL_ERROR;
        }
      }
      double stuff[12];
      if(argc > 7+12)
      {
        for (int i = 0; i < 12; ++i)
        {
          if (Tcl_GetDouble(interp, argv[7+i+1], &stuff[i]) != TCL_OK) 
          {
                opserr << "WARNING insufficient number of arguments - want: pattern ";
                opserr << "H5DRM " << patternID << " filename factor crd_scale distance_tolerance do_coordinate_transformation stuff\n";
                return TCL_ERROR;
          }
        }
      }

      // opserr << "!!Creating H5DRM tag = " << tag 
      // << "\n   filename = " << filename.c_str() 
      // << "\n   factor = " << factor 
      // << "\n   crd_scale = " << crd_scale 
      // << "\n   distance_tolerance = " << distance_tolerance 
      // << "\n   do_coordinate_transformation = " << do_coordinate_transformation 
      // << endln;
      double T00=stuff[0]; double T01=stuff[1]; double T02=stuff[2];
      double T10=stuff[3]; double T11=stuff[4]; double T12=stuff[5];
      double T20=stuff[6]; double T21=stuff[7]; double T22=stuff[8];
      double x00=stuff[9]; double x01=stuff[10]; double x02=stuff[11];

      // opserr << "T = " << endln;
      // opserr << T00 << " " << T01 << " " << T02 << endln;
      // opserr << T10 << " " << T11 << " " << T12 << endln;
      // opserr << T20 << " " << T21 << " " << T22 << endln;
      // opserr << "x0 = " << endln;
      // opserr << x00 << " " << x01 << " " << x02 << endln;
      
      
      
      

      
      
      
      
      
      

      thePattern = new H5DRM(tag, filename, factor,crd_scale, distance_tolerance, do_coordinate_transformation, T00, T01, T02, T10, T11, T12, T20, T21, T22, x00, x01, x02);

      // opserr << "Done! Creating H5DRM tag = " << tag << " filename = " << filename.c_str() << " factor = " << factor << endln;

      theDomain->addLoadPattern(thePattern);
      return TCL_OK;
    }
#endif

  //////// //////// ///////// ////////// /////  // DRMLoadPattern add BEGIN
  else if (strcmp(argv[1],"DRMLoadPattern") == 0) {
    TCL_Char * InputDataFileName = 0;
    
    if ((strcmp(argv[2],"-inputdata") == 0) || (strcmp(argv[2],"-InputData") == 0)) {
      InputDataFileName = argv[3];
    
      // now parse the input file name to extract the pattern input data
      std::ifstream ifile(InputDataFileName);
      int num_steps; ifile >> num_steps;
      double dt; ifile >> dt;
      int steps_cached; ifile >> steps_cached;
      double* ele_d = new double[3];
      ifile >> ele_d[0];
      ifile >> ele_d[1];
      ifile >> ele_d[2];
      double* drm_box_crds = new double[6];
      for (int i=0; i<6; i++)
    ifile >> drm_box_crds[i];
      int n1; ifile >> n1;
      int n2; ifile >> n2;
      
      std::string inps;
      int nf = 6;
      char** files = new char*[nf];
      int* f_d = new int[3*(nf-1)];
      int ne1,ne2;
      for (int i=0; i<nf; i++) {
    ifile >> inps;
    files[i] = (char*) inps.c_str();
    if (i <(nf-1)) {
      ifile >> ne1;
      ifile >> ne2; 
      f_d[3*i] = (ne1+1)*(ne2+1);
      f_d[3*i+1] = ne1;
      f_d[3*i+2] = ne2;
    }
      }
      
      Mesh3DSubdomain * myMesher = new Mesh3DSubdomain(theDomain);
      PlaneDRMInputHandler* patternhandler = new PlaneDRMInputHandler(1.0,files,nf,dt,0,num_steps,f_d,15,n1,n2,
                                      drm_box_crds,drm_box_crds,ele_d,
                                      myMesher, steps_cached,theDomain);
      DRMLoadPattern* ptr = new DRMLoadPattern(1,1.0,patternhandler,theDomain);
      ptr->setMaps();
      thePattern = ptr;
      theTclMultiSupportPattern = 0;
    } else {

      //     TCL_Char * ifp = 0;
      double INVALID = 0.7111722273337;
      int c_arg = 3;
      int end = argc-1;
      double dt=INVALID; 
      double* ele_d = new double[3];
      double* drm_box_crds = new double[6];
      for (int i=0; i<3; i++) {
    ele_d[i] = INVALID;
    drm_box_crds[2*i] = INVALID;
    drm_box_crds[2*i+1] = INVALID;
      }
      
      int nf =6;
      char** files = new char*[nf];
      files[5] = "./NONE";
      int* f_d = new int[15];
      int num_steps=1;
      int steps_cached=10;
      int n1,n2;
      n2=0;
      double factor =1.0;
      
      
      while ( c_arg < end ) {
    
    if ((strcmp(argv[c_arg],"-dt") == 0) || (strcmp(argv[c_arg],"-deltaT") == 0) ) {
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &dt) != TCL_OK) {
        opserr << " Error reading deltaT for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    if ((strcmp(argv[c_arg],"-numSteps") == 0) || (strcmp(argv[c_arg],"-numberOfSteps") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &num_steps) != TCL_OK) {
        opserr << " Error reading number of steps for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-stepsCached") == 0) || (strcmp(argv[c_arg],"-cache") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &steps_cached) != TCL_OK) {
        opserr << " Error reading number of steps for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-gridSize") == 0) || (strcmp(argv[c_arg],"-eleSize") == 0) ) {
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &ele_d[0]) != TCL_OK) {
        opserr << " Error reading deltaT for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &ele_d[1]) != TCL_OK) {
        opserr << " Error reading deltaT for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &ele_d[2]) != TCL_OK) {
        opserr << " Error reading deltaT for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-gridDataFace1") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[1]) != TCL_OK) {
        opserr << " Error reading grid data f1 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[2]) != TCL_OK) {
        opserr << " Error reading grid data f1 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      f_d[0] = (f_d[1]+1)*(f_d[2]+1);
    }
    
    else if ((strcmp(argv[c_arg],"-gridDataFace2") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[4]) != TCL_OK) {
        opserr << " Error reading grid data f2 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[5]) != TCL_OK) {
        opserr << " Error reading grid data f2 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      f_d[3] = (f_d[4]+1)*(f_d[5]+1);
    }
    
    else if ((strcmp(argv[c_arg],"-gridDataFace3") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[7]) != TCL_OK) {
        opserr << " Error reading grid data f3 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[8]) != TCL_OK) {
        opserr << " Error reading grid data f3 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      f_d[6] = (f_d[7]+1)*(f_d[8]+1);
    }
    
    else if ((strcmp(argv[c_arg],"-gridDataFace4") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[10]) != TCL_OK) {
        opserr << " Error reading grid data f4 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[11]) != TCL_OK) {
        opserr << " Error reading grid data f4 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      f_d[9] = (f_d[10]+1)*(f_d[11]+1);
    }
    
    else if ((strcmp(argv[c_arg],"-gridDataFace5") == 0) ) {
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[13]) != TCL_OK) {
        opserr << " Error reading grid data f5 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetInt(interp,argv[c_arg], &f_d[14]) != TCL_OK) {
        opserr << " Error reading grid data f5 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      f_d[12] = (f_d[13]+1)*(f_d[14]+1);
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace1") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[0] = new char[tmp.size()+1];
      strcpy(files[0],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace2") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[1] = new char[tmp.size()+1];
      strcpy(files[1],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace3") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[2] = new char[tmp.size()+1];
      strcpy(files[2],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace4") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[3] = new char[tmp.size()+1];
      strcpy(files[3],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace5a") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[4] = new char[tmp.size()+1];
      strcpy(files[4],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-filePathFace5b") == 0) ) {
      c_arg++;
      std::string tmp(argv[c_arg]);
      files[5] = new char[tmp.size()+1];
      strcpy(files[5],tmp.c_str());
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-fileFace5aGridPoints") == 0) ) {
      c_arg++;  
      if (Tcl_GetInt(interp,argv[c_arg], &n1) != TCL_OK) {
        opserr << " Error reading grid data f5 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-fileFace5bGridPoints") == 0) ) {
      c_arg++;  
      if (Tcl_GetInt(interp,argv[c_arg], &n2) != TCL_OK) {
        opserr << " Error reading grid data f5 for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-factor") == 0) || (strcmp(argv[c_arg],"-Factor") == 0) ) {
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &factor) != TCL_OK) {
        opserr << " Error reading number of steps for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
    else if ((strcmp(argv[c_arg],"-DRMBoxCrds") == 0) ) {
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[0]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, xmin for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[1]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, xmax for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[2]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, ymin for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[3]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, ymax for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[4]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, zmin for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
      if (Tcl_GetDouble(interp,argv[c_arg], &drm_box_crds[5]) != TCL_OK) {
        opserr << " Error reading DRM box Crds, zmax for DRMLoadPattern \n";
        exit(-1);
      }
      c_arg++;
    }
    
      }
      
      thePattern = new DRMLoadPatternWrapper(patternID,factor,files,nf,dt,num_steps,f_d,15,n1,n2,
                         drm_box_crds,ele_d,
                         steps_cached);
      theTclMultiSupportPattern = 0;      
      commandEndMarker = c_arg;

    }
    
  }// end else if DRMLoadPattern


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
  
  //  theTclMultiSupportPattern = 0;

  return TCL_OK;
}



