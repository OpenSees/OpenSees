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
                                                                        
// $Revision: 1.1 $
// $Date: 2000-12-12 07:33:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/TclGroundMotionCommand.cpp,v $

// File: ~/domain/groundMotion/TclGroundMotionComand.C
// 
// Written: fmk 
// Created: 11/00
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the GroundMotion command in the interpreter. It is invoked by the 
// TclModelBuilder_addGroundMotion function in the TclModelBuilder.C file. Current 
// valid GroundMotion types are:

// What: "@(#) TclGroundMotionCommand.C, revA"

#include <TclModelBuilder.h>

#include <tcl.h>
#include <GroundMotionRecord.h>

#include <string.h>
#include <MultiSupportPattern.h>
#include <InterpolatedGroundMotion.h>
#include <TimeSeriesIntegrator.h>


extern TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, char *arg);

extern TimeSeriesIntegrator *
TclSeriesIntegratorCommand(ClientData clientData, Tcl_Interp *interp, char *arg);

int
TclGroundMotionCommand(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char **argv, MultiSupportPattern *thePattern)
{
  GroundMotion *theMotion = 0;
  int gMotionTag;
    
  // make sure at least one other argument to contain integrator
  if (argc < 4) {
    cerr << "WARNING invalid command - want: groundMotion tag type <args>\n";
    cerr << "           valid types: AccelRecord and Interpolated \n";
    return TCL_ERROR;
  }    


  if (Tcl_GetInt(interp, argv[1], &gMotionTag) != TCL_OK) {
    cerr << "WARNING invalid tag: groundMotion tag  type <args>\n";
    return TCL_ERROR;		
  }
  int startArg = 2;

  if ((strcmp(argv[startArg],"Series") == 0) ||
      (strcmp(argv[startArg],"Plain") == 0)) {

    TimeSeries *accelSeries = 0;
    TimeSeries *velSeries = 0;
    TimeSeries *dispSeries = 0;
    TimeSeriesIntegrator *seriesIntegrator = 0;
    
    int currentArg = startArg+1;	

    while (currentArg < argc-1) {
      if ((strcmp(argv[currentArg],"-accel") == 0) ||
	  (strcmp(argv[currentArg],"-acceleration") == 0)) {
	  
	currentArg++;
	accelSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
	if (accelSeries == 0) {
	  cerr << "WARNING invalid accel series: groundMotion tag Series -accel <args>\n";
	  return TCL_ERROR;		
	}

	currentArg++;

      } else if ((strcmp(argv[currentArg],"-vel") == 0) ||
		 (strcmp(argv[currentArg],"-velocity") == 0)) {
	  
	currentArg++;
	velSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);
	
	if (velSeries == 0) {
	  cerr << "WARNING invalid vel series: " << argv[currentArg];
	  cerr << " groundMotion tag Series -vel {series}\n";
	  return TCL_ERROR;		
	}
	
	currentArg++;

      } else if ((strcmp(argv[currentArg],"-disp") == 0) ||
		 (strcmp(argv[currentArg],"-displacement") == 0)) {
	  
	currentArg++;
	dispSeries = TclSeriesCommand(clientData, interp, argv[currentArg]);

	if (dispSeries == 0) {
	  cerr << "WARNING invalid vel series: " << argv[currentArg];
	  cerr << " groundMotion tag Series -vel {series}\n";
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
	  cerr << " - groundMotion tag Series -int {Series Integrator}\n";
	  return TCL_ERROR;		
	}	  
	currentArg++;
      }

    }

    theMotion = new GroundMotion(dispSeries, velSeries, 
				 accelSeries, seriesIntegrator);

    if (theMotion == 0) {
      cerr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
      cerr << gMotionTag << endl;
      
      return TCL_ERROR;	      
    }      
  }      

  else if (strcmp(argv[startArg],"Interpolated") == 0) {

    int endMotionIDs = startArg+1;
    int motionID;
    while (Tcl_GetInt(interp, argv[endMotionIDs], &motionID) == TCL_OK) {
      endMotionIDs++;
    }
    int numMotions = endMotionIDs - startArg - 1;
    GroundMotion **theMotions;
    if (numMotions != 0) {
      theMotions = new GroundMotion *[numMotions];	
      ID motionIDs(numMotions); 
      for (int i=3; i<endMotionIDs; i++) {
	if (Tcl_GetInt(interp, argv[i], &motionID) != TCL_OK)	
	  return TCL_ERROR;	
	motionIDs[i-3] = motionID;	  
	GroundMotion *theMotion1 = 
	  thePattern->getMotion(motionID);
	
	if (theMotion1 == 0) {
	  cerr << "WARNING no groundMotion with tag " << motionID <<" :";
	  cerr << " pattern MultiSupport gMotion1? gMotion? .. ";
	  cerr << "-fact fact1? fact2? .. \n";
	  return TCL_ERROR;		    
	} else
	  theMotions[i-3] = theMotion1;
      }
    } else {
      cerr << "WARNING no gMotionTags want :";
      cerr << " pattern MultiSupport gMotion1? gMotion? .. ";
      cerr << "-fact fact1? fact2? .. \n";
      return TCL_ERROR;
    }
    Vector facts(numMotions);
    for (int i=0; i<numMotions; i++) {
      double fact;
      if (Tcl_GetDouble(interp, argv[i+startArg+1+numMotions+2], &fact) != TCL_OK)	
	return TCL_ERROR;		    
      facts[i] = fact;
    }
    
    theMotion = new InterpolatedGroundMotion(theMotions, facts, false);
	
  }
    
  else { 
    cerr << "WARNING unknown pattern type " << argv[1];
    cerr << " - want: pattern patternType " << gMotionTag ;
    cerr << " \t valid types: Plain, UniformExcitation \n";
    return TCL_ERROR;      
  }


  // now add the load pattern to the modelBuilder
  cerr << *thePattern;
  if (theMotion != 0) {
    if (thePattern->addMotion(*theMotion, gMotionTag) < 0) {
      cerr << "WARNING could not add ground motion with tag " << gMotionTag;
      cerr << " to pattern\n ";
      delete theMotion; // free up the memory, pattern destroys the time series
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}




