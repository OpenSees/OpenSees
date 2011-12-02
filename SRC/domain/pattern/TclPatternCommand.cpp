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
                                                                        
// $Revision: 1.3 $
// $Date: 2000-12-14 08:41:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclPatternCommand.cpp,v $

// File: ~/domain/pattern/TclPatternComand.C
// 
// Written: fmk 
// Created: 07/99
//
// Modified: fmk 11/00 - removed TimeSeries stuff from file, now an external procedure
// Revision: B
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

#include <string.h>



LoadPattern *theTclLoadPattern =0;
MultiSupportPattern *theTclMultiSupportPattern =0;


extern TimeSeriesIntegrator *
TclSeriesIntegratorCommand(ClientData clientData, Tcl_Interp *interp, char *arg);

extern TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, char *arg);
		 

int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			      int argc, char **argv, Domain *theDomain)
{
  LoadPattern *thePattern = 0;

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
      cerr << "WARNING invalid command - want: pattern type ";
      cerr << " <type args> {list of load and sp constraints commands}\n";
      cerr << "           valid types: Plain, UniformExcitation, MultiSupport\n";
      return TCL_ERROR;
  }    

  TimeSeries *theSeries = 0;
  int patternID =0;

  if (Tcl_GetInt(interp, argv[2], &patternID) != TCL_OK) {
    cerr << "WARNING invalid patternID: pattern type " << argv[2]
	<< "<type args>\n";
    return TCL_ERROR;		
  }

  int  commandEndMarker = 0; 	

  if (strcmp(argv[1],"Plain") == 0) {

      thePattern = new LoadPattern(patternID);       
      theSeries = TclSeriesCommand(clientData, interp, argv[3]);

	
      if (thePattern == 0 || theSeries == 0) {

	  if (thePattern == 0) {
	      cerr << "WARNING - out of memory creating LoadPattern ";
	      cerr << patternID << endl;
	  } else {
	      cerr << "WARNING - problem creating TimeSeries for LoadPattern ";
	      cerr << patternID << endl;
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
	char inputDir;
	inputDir = argv[3][0];
	switch (inputDir)  {

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
			cerr << "WARNING cannot read direction for excitation \n";
			cerr << "UniformExcitation " << patternID << " dir factor" << endl;
			return TCL_ERROR;
			break;
	}



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
	    cerr << "WARNING invalid vel0: pattern type UniformExciation\n";
	    return TCL_ERROR;		
	  }

	  currentArg++;
	}


	else if ((strcmp(argv[currentArg],"-accel") == 0) ||
		 (strcmp(argv[currentArg],"-acceleration") == 0)) {
	  
	  currentArg++;
	  accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
	  if (accelSeries == 0) {
	    cerr << "WARNING invalid accel series: pattern UniformExcitation -accel <args>\n";
	    return TCL_ERROR;		
	  }

	  currentArg++;

	} else if ((strcmp(argv[currentArg],"-vel") == 0) ||
		   (strcmp(argv[currentArg],"-velocity") == 0)) {
	  
	  currentArg++;
	  velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

	  if (velSeries == 0) {
	    cerr << "WARNING invalid vel series: " << argv[currentArg];
	    cerr << " pattern UniformExcitation -vel {series}\n";
	    return TCL_ERROR;		
	  }
	  
	  currentArg++;

	} else if ((strcmp(argv[currentArg],"-disp") == 0) ||
		   (strcmp(argv[currentArg],"-displacement") == 0)) {
	  
	  currentArg++;
	  dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

	  if (dispSeries == 0) {
	    cerr << "WARNING invalid vel series: " << argv[currentArg];
	    cerr << " pattern UniformExcitation -vel {series}\n";
	    return TCL_ERROR;		
	  }
	  
	  currentArg++;
	} else if ((strcmp(argv[currentArg],"-int") == 0) ||
		   (strcmp(argv[currentArg],"-integrator") == 0)) {
	  
	  currentArg++;
	  seriesIntegrator = TclSeriesIntegratorCommand(clientData, interp, 
							argv[currentArg]);
	  if (seriesIntegrator == 0) {
	    cerr << "WARNING invalid series integrator: " << argv[currentArg];
	    cerr << " - pattern UniformExcitation -int {Series Integrator}\n";
	    return TCL_ERROR;		
	  }	  
	  currentArg++;
	}

	else
	  doneSeries = true;
      }
      
      if (dispSeries == 0 && velSeries == 0 && accelSeries == 0) {
	  cerr << "WARNING invalid series, want - pattern UniformExcitation";
	  cerr << "-disp {dispSeries} -vel {velSeries} -accel {accelSeries} ";
	  cerr << "-int {Series Integrator}\n";
	  return TCL_ERROR;			  
      }

      GroundMotion *theMotion = new GroundMotion(dispSeries, velSeries, 
						 accelSeries, seriesIntegrator);

      if (theMotion == 0) {
	  cerr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
	  cerr << patternID << endl;
	
	  return TCL_ERROR;	      
      }
    
	// create the UniformExcitation Pattern
	thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0);

	if (thePattern == 0) {
	    cerr << "WARNING ran out of memory creating load pattern - pattern UniformExcitation ";
	    cerr << patternID << endl;
	
	    // clean up memory allocated up to this point and return an error
	    if (theMotion != 0)
		delete theMotion;

	    return TCL_ERROR;	      
	}
  }

  else if (strcmp(argv[1],"Uniform") == 0) {

	// First search for file name and time step
	int numInputs = argc;
	char *accelFileName = 0;
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
			  cerr << "WARNING problem reading ground motion "
				   << "time interval - pattern UniformExcitation: "
		           << patternID << endl;
		      return TCL_ERROR;
			}
			numInputs -= 3;

		}
    }

    if (numInputs < 5) {
      cerr << "WARNING insufficient number of arguments - want: pattern ";
      cerr << "UniformExcitation " << patternID << " dir factor\n";
      return TCL_ERROR;
    }    

	int dir;
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
			cerr << "WARNING cannot read direction for excitation \n";
			cerr << "UniformExcitation " << patternID << " dir factor" << endl;
			return TCL_ERROR;
			break;
	}

	double factor;
	if (Tcl_GetDouble(interp, argv[4], &factor) != TCL_OK) {
	  cerr << "WARNING insufficient number of arguments - want: pattern ";
	  cerr << "UniformExcitation " << patternID << " dir factor\n";
	  return TCL_ERROR;
	} 
    
    GroundMotionRecord *theMotion;


    // read in the ground motion
	if (accelFileName == 0)
	{
		cerr << "WARNING -- No ground motion data provided\n";
		cerr << "UniformExcitation tag: " << patternID << endl;
		return TCL_ERROR;
	}
	
    theMotion = new GroundMotionRecord(accelFileName, dt, factor);

	if (theMotion == 0) {
	  cerr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
	  cerr << patternID << endl;
	
	  return TCL_ERROR;	      
    }
    
    // create the UniformExcitation Pattern
    thePattern = new UniformExcitation(*theMotion, dir, patternID);

    if (thePattern == 0) {
	  cerr << "WARNING ran out of memory creating load pattern - pattern UniformExcitation ";
	  cerr << patternID << endl;
	
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
	  cerr << "WARNING ran out of memory creating load pattern - pattern MultipleSupportExcitation ";
	  cerr << patternID << endl;
      
      }      
      commandEndMarker = 2;
  } 

  else { 
      cerr << "WARNING unknown pattern type " << argv[1];
      cerr << " - want: pattern patternType " << patternID ;
      cerr << " \t valid types: Plain, UniformExcitation, MultiSupportExciatation \n";
      return TCL_ERROR;      
  }

  // now add the load pattern to the modelBuilder
  if (theDomain->addLoadPattern(thePattern) == false) {
    cerr << "WARNING could not add load pattern to the domain " << *thePattern;
    delete thePattern; // free up the memory, pattern destroys the time series
    return TCL_ERROR;
  }

  theTclLoadPattern = thePattern;

  // use TCL_Eval to evaluate the list of load and single point constraint commands
  if (commandEndMarker < (argc-1)) {
    Tcl_Eval(interp, argv[argc-1]);
  }

  theTclMultiSupportPattern = 0;

  return TCL_OK;
}




