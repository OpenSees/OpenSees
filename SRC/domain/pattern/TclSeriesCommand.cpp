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
                                                                        
// $Revision: 1.11 $
// $Date: 2003-10-30 22:46:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TclSeriesCommand.cpp,v $

// Written: fmk 
// Created: 11/00
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the Pattern command in the interpreter. It is invoked by the 
// TclModelBuilder_addPattern function in the TclModelBuilder.C file. Current 
// valid Pattern types are:

// What: "@(#) TclPatternCommand.C, revA"

#include <tcl.h>
#include <Domain.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <RectangularSeries.h>
#include <TrigSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <string.h>

/////////////// terje //////////////////////////
#include <DiscretizedRandomProcessSeries.h>
#include <SimulatedRandomProcessSeries.h>
#include <Spectrum.h>
#include <RandomNumberGenerator.h>
#include <ReliabilityDomain.h>
extern ReliabilityDomain *theReliabilityDomain;
extern RandomNumberGenerator *theRandomNumberGenerator;
////////////////////////////////////////////////

// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in 
//   one array, which is why Tcl_Free needs only be called on the array.
static void cleanup(TCL_Char **argv) {
	  Tcl_Free((char *) argv);
}


TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char *arg)
{

  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    opserr << "WARNING could not split series list " << arg << endln;
    return 0;
  }
			    
  TimeSeries *theSeries = 0;

  if (strcmp(argv[0],"Constant") == 0) {
    // LoadPattern and ConstantSeries - read args & create LinearSeries object
    double cFactor = 1.0;
	
    int endMarker = 1;
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	    
	opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	opserr << " Constant -factor cFactor\n";
	cleanup(argv);
	return 0;
      }
      endMarker++;
    }

    theSeries = new ConstantSeries(cFactor);       	
    
  } else if (strcmp(argv[0],"Trig") == 0 || 
	     strcmp(argv[0],"Sine") == 0) {
    // LoadPattern and TrigSeries - read args & create TrigSeries object
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double shift = 0.0;
      
    if (argc < 4) {
      opserr << "WARNING not enough TimeSeries args - ";
      opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }	
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << " - ";
      opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;				
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tFinish " << argv[2] << " - ";
      opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }     
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      opserr << "WARNING invalid period " << argv[3] << " - ";
      opserr << " Trig tStart tFinish period <-shift shift> <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }     
    
    int endMarker = 4;
    
    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker],"-factor") == 0) {
	// allow user to specify the factor
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	  opserr << "WARNING invalid cFactor " << argv[endMarker] << " -";
	  opserr << " Trig  tStart tFinish period -factor cFactor\n";
	  cleanup(argv);
	  return 0;
	}
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
	// allow user to specify phase shift
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
	    
	  opserr << "WARNING invalid phase shift " << argv[endMarker] << " - ";
	  opserr << " Trig tStart tFinish period -shift shift\n";
	  cleanup(argv);
	  return 0;
	}
      }
      endMarker++;
    }

    theSeries = new TrigSeries(tStart, tFinish, period, shift, cFactor);
	
  }	

  else if (strcmp(argv[0],"Linear") == 0) {
    // LoadPattern and LinearSeries - read args & create LinearSeries object
    double cFactor = 1.0;
      
    int endMarker = 1;
    
    if ((endMarker < argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	
	opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	opserr << " Linear  -factor cFactor\n";
	cleanup(argv);
	return 0;
      }	
      endMarker++;
    }

    theSeries = new LinearSeries(cFactor);       	
  }

// AddingSensitivity:BEGIN ///////////////////////////////////

#ifdef _RELIABILITY

	else if (strcmp(argv[0],"DiscretizedRandomProcess") == 0) {


		double mean, maxStdv;
		ModulatingFunction *theModFunc;

		if (Tcl_GetDouble(interp, argv[1], &mean) != TCL_OK) {
			opserr << "WARNING invalid input: random process mean \n";
			cleanup(argv);
			return 0;
		}

		if (Tcl_GetDouble(interp, argv[2], &maxStdv) != TCL_OK) {
			opserr << "WARNING invalid input: random process max stdv \n";
			cleanup(argv);
			return 0;
		}


		// Number of modulating functions
		int argsBeforeModList = 3;
		int numModFuncs = argc-argsBeforeModList;

		// Create an array to hold pointers to modulating functions
		ModulatingFunction **theModFUNCS = new ModulatingFunction *[numModFuncs];

		// For each modulating function, get the tag and ensure it exists
		int tagI;
		for (int i=0; i<numModFuncs; i++) {
			if (Tcl_GetInt(interp, argv[i+argsBeforeModList], &tagI) != TCL_OK) {
				opserr << "WARNING invalid modulating function tag. " << endln;
				cleanup(argv);
				return 0;
			}

			theModFunc = 0; 
			theModFunc = theReliabilityDomain->getModulatingFunction(tagI);

			if (theModFunc == 0) {
				opserr << "WARNING modulating function number "<< argv[i+argsBeforeModList] << "does not exist...\n";
				delete [] theModFUNCS;
				cleanup(argv);
				return 0;
			}
			else {
				theModFUNCS[i] = theModFunc;
			}

		}	

		// Parsing was successful, create the random process series object
		theSeries = new DiscretizedRandomProcessSeries(numModFuncs,theModFUNCS,mean,maxStdv);       	
	}

 



	else if (strcmp(argv[0],"SimulatedRandomProcess") == 0) {

		int spectrumTag, numFreqIntervals;
		double mean;

		if (Tcl_GetInt(interp, argv[1], &spectrumTag) != TCL_OK) {
			opserr << "WARNING invalid input to SimulatedRandomProcess: spectrumTag" << endln;
			cleanup(argv);
			return 0;
		}

		if (Tcl_GetDouble(interp, argv[2], &mean) != TCL_OK) {
			opserr << "WARNING invalid input to SimulatedRandomProcess: mean" << endln;
			cleanup(argv);
			return 0;
		}

		if (Tcl_GetInt(interp, argv[3], &numFreqIntervals) != TCL_OK) {
			opserr << "WARNING invalid input to SimulatedRandomProcess: numFreqIntervals" << endln;
			cleanup(argv);
			return 0;
		}

		// Check that the random number generator exists
		if (theRandomNumberGenerator == 0) {
			opserr << "WARNING: A random number generator must be instantiated before SimulatedRandomProcess." << endln;
			cleanup(argv);
			return 0;
		}

		// Check that the spectrum exists
		Spectrum *theSpectrum = 0;
		theSpectrum = theReliabilityDomain->getSpectrum(spectrumTag);
		if (theSpectrum == 0) {
			opserr << "WARNING: Could not find the spectrum for the SimulatedRandomProcess." << endln;
			cleanup(argv);
			return 0;
		}


		// Parsing was successful, create the random process series object
		theSeries = new SimulatedRandomProcessSeries(theRandomNumberGenerator,theSpectrum,numFreqIntervals,mean);
	}

#endif
  
// AddingSensitivity:END /////////////////////////////////////




  else if (strcmp(argv[0],"Rectangular") == 0) {
    // LoadPattern and RectangularSeries - read args and create RectangularSeries object
    double tStart, tFinish;
    double cFactor = 1.0;
    if (argc < 3) {
      opserr << "WARNING not enogh args - ";
      opserr << " Rectangular tStart tFinish <-factor cFactor>\n";
      cleanup(argv);
      return 0;	
    }	
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << " - ";
      opserr << " Rectangular tStart tFinish <-factor factor>\n";
      cleanup(argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[2] << " - ";
      opserr << " Rectangular tStart tFinish <-factor fcator>\n";
      cleanup(argv);
      return 0;	
    }     

    int endMarker =  3;
    
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	
	opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	opserr << " Rectangular tStart tFinish -factor cFactor\n";
	cleanup(argv);
	return 0;
      }
      endMarker++;
    }  

    theSeries = new RectangularSeries(tStart, tFinish, cFactor); 
      
  }
    
  else if ((strcmp(argv[0],"Series") == 0) ||
	   (strcmp(argv[0],"Path") == 0)) {

    // LoadPattern and PathSeries - read args and create RectangularSeries object
    double cFactor = 1.0;
    if (argc < 3) {

      opserr << "WARNING not enough args - ";
      opserr << " Series -dt timeIncr -values {list of points }\n"; 
      cleanup(argv);
      return 0;	
    }

    double timeIncr = 0.0;
    int endMarker =  1;
    bool done = false;
    int fileName = 0;
    int fileTimeName = 0;
    int filePathName = 0;
    Vector *dataPath = 0;
    Vector *dataTime = 0;

    while (endMarker < argc && done == false) {
	
      if (strcmp(argv[endMarker],"-dt") == 0) {
	// allow user to specify the time increment
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {
	  
	  opserr << "WARNING invalid dt " << argv[endMarker] << " - ";
	  opserr << " Series -dt dt ... \n";
	  cleanup(argv);
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-factor") == 0) {
	// allow user to specify the factor
	endMarker++;
	if (endMarker == argc || 
	    Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
	  
	  opserr << "WARNING invalid cFactor " << argv[endMarker] << " - ";
	  opserr << " Series -factor ... \n";
	  cleanup(argv);
	  return 0;
	}
      } 

      else if (strcmp(argv[endMarker],"-file") == 0) {
	// allow user to specify the file name containg time and data points
	endMarker++;
	if (endMarker != argc) {
	  fileName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-filePath") == 0) {
	// allow user to specify the file name containg the data points
	endMarker++;
	if (endMarker != argc) {
	  filePathName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-fileTime") == 0) {
	// allow user to specify the file name containg the data points
	endMarker++;
	if (endMarker != argc) {
	  fileTimeName = endMarker; // argv[endMarker];
	}
      }

      else if (strcmp(argv[endMarker],"-values") == 0) {
	// allow user to specify the data points in tcl list
	endMarker++;
	if (endMarker != argc) {
	  int pathSize;
	  TCL_Char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
		      
	    opserr << "WARNING problem splitting path list " << argv[endMarker] << " - ";
	    opserr << " Series -values {path} ... \n";
	    cleanup(argv);
	    return 0;
	  }
	  
	  dataPath = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      opserr << "WARNING problem reading path data value " << pathStrings[i] << " - ";
	      opserr << " Series -values {path} ... \n";
	      cleanup(argv);
	      cleanup(pathStrings);
	      return 0;
	    }
	    (*dataPath)(i) = value;
	  }
	  // free up the array of pathsStrings .. see tcl man pages as to why
	  cleanup(pathStrings);
	}
      }

      else if (strcmp(argv[endMarker],"-time") == 0) {
	// allow user to specify the data points in tcl list
	endMarker++;
	if (endMarker != argc) {
	  int pathSize;
	  TCL_Char **pathStrings;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &pathSize, &pathStrings) != TCL_OK) {
			
	    opserr << "WARNING problem spltting time path " << argv[endMarker] << " - ";
	    opserr << " Series -time {times} ... \n";
	    cleanup(argv);
	    return 0;
	  }
	  
	  dataTime = new Vector(pathSize);
	  for (int i=0; i<pathSize; i++) {
	    double value;
	    if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
	      opserr << "WARNING problem reading time path value " << pathStrings[i] << " - ";
	      opserr << " Series -values {path} ... \n";

	      cleanup(argv);
	      cleanup(pathStrings);
	      return 0;
	    }
	    (*dataTime)(i) = value;
	  }
	  // free up the array of pathsStrings .. see tcl man pages as to why
	  cleanup(pathStrings);
	}
      }
      
      endMarker++;
    }
    

    if (filePathName != 0 && fileTimeName == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(argv[filePathName], timeIncr, cFactor);
    }

    else if (fileName != 0) {
      theSeries = new PathTimeSeries(argv[fileName], cFactor);

    } else if (filePathName != 0 && fileTimeName != 0)
      theSeries = new PathTimeSeries(argv[filePathName], argv[fileTimeName], cFactor); 

    else if (dataPath != 0 && dataTime == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(*dataPath, timeIncr, cFactor); 
      delete dataPath;

    } else if (dataPath != 0 && dataTime != 0) {
      theSeries = new PathTimeSeries(*dataPath, *dataTime, cFactor);  
      delete dataPath;      
      delete dataTime;      

    } else {
      opserr << "WARNING choice of options for Path Series invalid - valid options for ";
      opserr << " Path are\n";
      opserr << " \t -fileT fileTimeName -fileP filePathName \n";
      opserr << " \t -dt constTimeIncr -file filePathName\n";
      opserr << " \t -dt constTimeIncr -values {list of points on path}\n";
      opserr << " \t -time {list of time points} -values {list of points on path}\n";
      cleanup(argv);
      return 0;
    }      
	
  } 
	
  else {
    // type of load pattern type unknown
    opserr << "WARNING unknown Series type " << argv[0] << " - ";
    opserr << " valid types: Linear, Rectangular, Path, Constant, Trig, Sine\n";
    cleanup(argv);
    return 0;
  }

  // clean up after ourselves and return the series
  cleanup(argv);
  return theSeries;
}




