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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclPatternCommand.cpp,v $

// File: ~/domain/pattern/TclPatternComand.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
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
#include <GroundMotionRecord.h>

#include <string.h>

LoadPattern *theTclLoadPattern =0;

int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			      int argc, char **argv, Domain *theDomain)
{
  LoadPattern *thePattern = 0;

  // make sure at least one other argument to contain integrator
  if (argc < 4) {
      cerr << "WARNING invalid command - want: pattern type ";
      cerr << " <type args> {list of load and sp constraints commands}\n";
      cerr << "           valid types: Plain, UniformExcitation \n";
      return TCL_ERROR;
  }    

  TimeSeries *theSeries = 0;
  int patternID =0;

  if (Tcl_GetInt(interp, argv[2], &patternID) != TCL_OK) {
    cerr << "WARNING invalid patternID: pattern type " << argv[2]
	<< "<type args>\n";
    return TCL_ERROR;		
  }

  bool commandEnd = true;
  int  commandEndMarker = 0; 	

  if (strcmp(argv[1],"Plain") == 0) {
    //
    // check argv[3] for type of Pattarn (pattern & time series)  
    // and create the object 
    //
	if (strcmp(argv[3],"Constant") == 0) {
		// LoadPattern and ConstantSeries - read args & create LinearSeries object
		double cFactor = 1.0;
	
		int endMarker = 4;
		while (endMarker != argc) {
			if (strcmp(argv[endMarker],"-factor") == 0) {
			// allow user to specify the factor
				endMarker++;
				if (endMarker == argc || 
					Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
  
					cerr << "WARNING invalid cFactor - pattern Plain " << patternID
						<< " Constant -factor cFactor\n";
					return TCL_ERROR;
				}
			}
			else { 
				// it must be a string of loads and sp constraint commands
				commandEnd = false;
				commandEndMarker = endMarker;
			}
			endMarker++;
		}
		
		thePattern = new LoadPattern(patternID);       
		theSeries = new ConstantSeries(cFactor);       	
	
		if (thePattern == 0 || theSeries == 0) {
			cerr << "WARNING ran out of memory - pattern Plain Constant "
				<< patternID << endl;

			// clean up the memory and return an error
			if (thePattern != 0)
				delete thePattern;
			if (theSeries != 0)
				delete theSeries;
			return TCL_ERROR;	      
		}

		thePattern->setTimeSeries(theSeries);
	}

	else if (strcmp(argv[3],"Trig") == 0 || strcmp(argv[3],"Sine") == 0) {
		// LoadPattern and TrigSeries - read args & create TrigSeries object
		double cFactor = 1.0;
		double tStart, tFinish, period;
		double shift = 0.0;
		
		if (argc < 7) {
			cerr << "WARNING invalid command - pattern Plain " << patternID << endl;
			cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
			return TCL_ERROR;	
		}	
		if (Tcl_GetDouble(interp, argv[4], &tStart) != TCL_OK) {
			cerr << "WARNING invalid tStart - pattern Plain " << patternID << endl;
			cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
			return TCL_ERROR;				
		}
		if (Tcl_GetDouble(interp, argv[5], &tFinish) != TCL_OK) {
			cerr << "WARNING invalid tFinish - pattern Plain " << patternID << endl;
			cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
			return TCL_ERROR;	
		}     
		if (Tcl_GetDouble(interp, argv[6], &period) != TCL_OK) {
			cerr << "WARNING invalid period - pattern Plain " << patternID << endl;
			cerr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
			return TCL_ERROR;	
		}     


		int endMarker = 7;
		while (endMarker != argc) {
			if (strcmp(argv[endMarker],"-factor") == 0) {
			// allow user to specify the factor
				endMarker++;
				if (endMarker == argc || 
					Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
  
					cerr << "WARNING invalid cFactor - pattern Plain " << patternID
						<< " Linear  -factor cFactor\n";
					return TCL_ERROR;
				}
			}

			else if (strcmp(argv[endMarker],"-shift") == 0) {
			// allow user to specify phase shift
				endMarker++;
				if (endMarker == argc || 
					Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
  
					cerr << "WARNING invalid phase shift - pattern Plain " << patternID
						<< " Trig -shift shift\n";
					return TCL_ERROR;
				}
			}

			else { 
				// it must be a string of loads and sp constraint commands
				commandEnd = false;
				commandEndMarker = endMarker;
			}
			endMarker++;
		}
		
		thePattern = new LoadPattern(patternID);       
		theSeries = new TrigSeries(tStart, tFinish, period, shift, cFactor);
	
		if (thePattern == 0 || theSeries == 0) {
			cerr << "WARNING ran out of memory - pattern Plain Trig "
				<< patternID << endl;

			// clean up the memory and return an error
			if (thePattern != 0)
				delete thePattern;
			if (theSeries != 0)
				delete theSeries;
			return TCL_ERROR;	      
		}

		thePattern->setTimeSeries(theSeries);
	}

    else if (strcmp(argv[3],"Linear") == 0) {
      // LoadPattern and LinearSeries - read args & create LinearSeries object
      double cFactor = 1.0;
      
      int endMarker = 4;
      while (endMarker != argc) {

	if (strcmp(argv[endMarker],"-factor") == 0) {
	  // allow user to specify the factor
	  endMarker++;
	  if (endMarker == argc || 
	      Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	    cerr << "WARNING invalid cFactor - pattern Plain " << patternID
		<< " Linear  -factor cFactor\n";
	    return TCL_ERROR;
	  }
	}  
	else { 
	  // it must be a string of loads and sp constraint commands
	  commandEnd = false;
	  commandEndMarker = endMarker;
	}
	endMarker++;
      }

      thePattern = new LoadPattern(patternID);       
      theSeries = new LinearSeries(cFactor);       	
      if (thePattern == 0 || theSeries == 0) {
	cerr << "WARNING ran out of memory - pattern Plain Linear "
	    << patternID << endl;
	
	// clean up the memory and return an error
	if (thePattern != 0)
	  delete thePattern;
	if (theSeries != 0)
	  delete theSeries;
	return TCL_ERROR;	      
      }
      thePattern->setTimeSeries(theSeries);    
    }
    
    else if (strcmp(argv[3],"Rectangular") == 0) {
      // LoadPattern and RectangularSeries - read args and create RectangularSeries object
      double tStart, tFinish;
      double cFactor = 1.0;
      if (argc < 6) {
	cerr << "WARNING invalid command - pattern Plain Rectangular "
	    << patternID << " Rectangular tStart tFinish <cFactor>\n";
	return TCL_ERROR;	
      }	
      if (Tcl_GetDouble(interp, argv[4], &tStart) != TCL_OK) {
	cerr << "WARNING invalid tStart - pattern Plain Rectangular ";
	cerr << patternID << " Rectangular tStart tFinish cFactor\n";
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &tFinish) != TCL_OK) {
	cerr << "WARNING invalid tStart - pattern Plain Rectangular ";
	cerr << patternID << " Rectangular tStart tFinish cFactor\n";
	return TCL_ERROR;	
      }     

      int endMarker = 6;
      while (endMarker != argc) {
	
	if (strcmp(argv[endMarker],"-factor") == 0) {
	  // allow user to specify the factor
	  endMarker++;
	  if (endMarker == argc || 
	      Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	    cerr << "WARNING invalid cFactor - pattern Plain Rectangular";
	    cerr << patternID << " Rectangular -factor cFactor\n";
	    return TCL_ERROR;
	  }
	}  
	else { 
	  // it must be a string of loads and sp constraint commands
	  commandEnd = false;
	  commandEndMarker = endMarker;
	}
	endMarker++;
      }

      thePattern = new LoadPattern(patternID);       
      theSeries = new RectangularSeries(tStart, tFinish, cFactor); 
      
      if (thePattern == 0 || theSeries == 0) {
	cerr << "WARNING ran out of memory - pattern Plain " << patternID;
	cerr << " Rectangular tStart tFinish cFactor\n";

	// clean up the memory and return an error
	if (thePattern != 0)
	  delete thePattern;
	if (theSeries != 0)
	  delete theSeries;
	return TCL_ERROR;	      
      }
      
      thePattern->setTimeSeries(theSeries);
    }
    
    else if (strcmp(argv[3],"Series") == 0) {
      // LoadPattern and PathSeries - read args and create RectangularSeries object
      double cFactor = 1.0;
      if (argc < 6) {
	cerr << "WARNING invalid command - pattern Plain  " << patternID;
	cerr << " Series -dt timeIncr -values {list of points } {list of ";
	cerr << " load and sp_constraint commands}\n";
	return TCL_ERROR;	
      }

      double timeIncr = 1.0;
      int endMarker = 4;
      bool done = false;
      char *fileTimeName = 0;
      char *filePathName = 0;
      Vector *dataPath = 0;
      Vector *dataTime = 0;

      commandEndMarker = argc-1;

      while (endMarker != argc && done == false) {
	
	if (strcmp(argv[endMarker],"-dt") == 0) {
	  // allow user to specify the time increment
	  endMarker++;
	  if (endMarker == argc || 
	      Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {
	    
	    cerr << "WARNING invalid dt - pattern Plain " << patternID;
	    cerr << " Series -dt dt ... \n";
	    return TCL_ERROR;
	  }
	} 

	else if (strcmp(argv[endMarker],"-factor") == 0) {
	  // allow user to specify the factor
	  endMarker++;
	  if (endMarker == argc || 
	      Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	    
	    cerr << "WARNING invalid cFactor - pattern Plain " << patternID;
	    cerr << " Series -factor ... \n";
	    return TCL_ERROR;
	  }
	} 

	else if (strcmp(argv[endMarker],"-filePath") == 0) {
	  // allow user to specify the file name containg the data points
	  endMarker++;
	  if (endMarker != argc) {
	    filePathName = argv[endMarker];
	  }
	}

	else if (strcmp(argv[endMarker],"-fileTime") == 0) {
	  // allow user to specify the file name containg the data points
	  endMarker++;
	  if (endMarker != argc) {
	    fileTimeName = argv[endMarker];
	  }
	}

	else if (strcmp(argv[endMarker],"-values") == 0) {
	  // allow user to specify the data points in tcl list
	  endMarker++;
	  if (endMarker != argc) {
	    int pathSize;
	    char **pathStrings;
	  
	    if (Tcl_SplitList(interp, argv[endMarker], 
			      &pathSize, &pathStrings) != TCL_OK) {

	      cerr << "WARNING problem reading path - pattern Plain " << patternID;
	      cerr << " Series -values {path} ... \n";
	      return TCL_ERROR;
	    }
	  
	    dataPath = new Vector(pathSize);
	    for (int i=0; i<pathSize; i++) {
	      double value;
	      if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
		cerr << "WARNING problem reading path - pattern Plain ";
		cerr << patternID << " Series -values {path} ... \n";
		return TCL_ERROR;
	      }
	      (*dataPath)(i) = value;
	    }
	    // free up the array of pathsStrings .. see tcl man pages as to why
#ifdef TCL_Free
	    Tcl_Free((char *) pathStrings);
#endif
	  }
	}

	else if (strcmp(argv[endMarker],"-time") == 0) {
	  // allow user to specify the data points in tcl list
	  endMarker++;
	  if (endMarker != argc) {
	    int pathSize;
	    char **pathStrings;
	  
	    if (Tcl_SplitList(interp, argv[endMarker], 
			      &pathSize, &pathStrings) != TCL_OK) {

	      cerr << "WARNING problem reading path - pattern Plain " << patternID;
	      cerr << " Series -time {times} ... \n";
	      return TCL_ERROR;
	    }
	  
	    dataTime = new Vector(pathSize);
	    for (int i=0; i<pathSize; i++) {
	      double value;
	      if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
		cerr << "WARNING problem reading path - pattern Plain ";
		cerr << patternID << " Series -values {path} ... \n";
		return TCL_ERROR;
	      }
	      (*dataTime)(i) = value;
	    }
	    // free up the array of pathsStrings .. see tcl man pages as to why
#ifdef TCL_Free	    
	    Tcl_Free((char *) pathStrings);
#endif	    
	  }
	}
      
	else 
	  done = true;
      
	endMarker++;
      }
    

      if (argc == endMarker)
	commandEnd = false;
    
      thePattern = new LoadPattern(patternID);       
      if (filePathName != 0 && fileTimeName == 0)
	theSeries = new PathSeries(filePathName, timeIncr, cFactor);
      else if (filePathName != 0 && fileTimeName != 0)
	theSeries = new PathTimeSeries(filePathName, fileTimeName, cFactor); 
      else if (dataPath != 0 && dataTime == 0) {
	  theSeries = new PathSeries(*dataPath, timeIncr, cFactor); 
	  delete dataPath;
      } else if (dataPath != 0 && dataTime != 0) {
	theSeries = new PathTimeSeries(*dataPath, *dataTime, cFactor);  
	delete dataPath;      
	delete dataTime;      
      } else {
	cerr << "WARNING choice of options invalid - valid options for ";
	cerr << "pattern Plain patternID Path:\n";
	cerr << " \t -fileT fileTimeName -fileP filePathName \n";
	cerr << " \t -dt constTimeIncr -file filePathName\n";
	cerr << " \t -dt constTimeIncr -values {list of points on path}\n";
	cerr << " \t -time {list of time points} -values {list of points on path}\n";
	if (thePattern == 0)
	  delete thePattern;
	return TCL_ERROR;
      }      

      if (thePattern == 0 || theSeries == 0) {
	cerr << "WARNING ran out of memory - pattern Plain " << patternID;
	cerr << " Path ... \n";
	
	// clean up the memory and return an error
	if (thePattern != 0)
	  delete thePattern;
	if (theSeries != 0)
	  delete theSeries;
	return TCL_ERROR;	      
      }
      
      thePattern->setTimeSeries(theSeries);
  }

    else {
	// type of load pattern type unknown
	cerr << "WARNING unknown Series type " << argv[3] << " - want: pattern Plain ";
	cerr << patternID;
	cerr << " seriesType <type args>\n\tvalid types: Linear, Rectangular, Path\n";
	return TCL_ERROR;
    }
  }

  else if (strcmp(argv[1],"UniformExcitation") == 0) {

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
  
  else { 
      cerr << "WARNING unknown pattern type " << argv[1];
      cerr << " - want: pattern patternType " << patternID ;
      cerr << " \t valid types: Plain, UniformExcitation \n";
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
  if (commandEnd == false) 
    Tcl_Eval(interp, argv[commandEndMarker]);

  return TCL_OK;
}




